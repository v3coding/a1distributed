#include "core/graph.h"
#include "core/utils.h"
#include <future>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

static int numberOfThreads;
static Graph g;

typedef struct threadStats{
    pthread_t threadID;
    double runtime;
    int threadTriangles;
    double threadRuntime;
    uintV *array1;
    uintE len1;
    uintV *array2;
    uintE len2;
    uintV u;
    uintV v;
    uintV n;
} threadStats;

typedef struct threadObject{
    pthread_mutex_t* writeMutex;
    pthread_t* threads;
    int totalTriangles;
    int uniqueTriangles;
    double totalRuntime;
    threadStats* threadStatistics;
    int totalThreads;
    int threadsRunning;
    int trianglesRemaining;
    Graph* g;
}   threadObject;

long countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2,
                     uintV u, uintV v) {

  uintE i = 0, j = 0; // indexes for array1 and array2
  long count = 0;

  if (u == v)
    return count;

  while ((i < len1) && (j < len2)) {
    if (array1[i] == array2[j]) {
      if ((array1[i] != u) && (array1[i] != v)) {
        count++;
      }
      i++;
      j++;
    } else if (array1[i] < array2[j]) {
      i++;
    } else {
      j++;
    }
  }
  return count;
}

void* countTrianglesP(void *_arg){
  timer t1;
  t1.start();
  threadObject* threadData = (threadObject*) _arg;
  int triangle_count;
  int index;


  for(int i = 0; i < numberOfThreads; i++){
    if(pthread_self() == threadData->threadStatistics[i].threadID){
      index = i;
    }
  }

  uintV u = threadData->threadStatistics[index].u;
  uintV n = threadData->threadStatistics[index].n;
  int threadTriangles = 0;
    for (; u < n; u++) {
    uintE out_degree = g.vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g.vertices_[u].getOutNeighbor(i);
        threadTriangles += countTriangles(g.vertices_[u].getInNeighbors(),
                                         g.vertices_[u].getInDegree(),
                                         g.vertices_[v].getOutNeighbors(),
                                         g.vertices_[v].getOutDegree(), u, v);
      }
    }
    double endtime = t1.stop();
    pthread_mutex_lock(threadData->writeMutex);
    threadData->threadStatistics[index].threadTriangles = threadTriangles;
    threadData->totalTriangles += threadTriangles;
    threadData->totalRuntime += endtime;
    threadData->threadStatistics[index].threadRuntime = endtime;
    pthread_mutex_unlock(threadData->writeMutex);
  return 0;
}

void triangleCountSerial(Graph &g) {
  uintV n = g.n_;
  long triangle_count = 0;
  double time_taken = 0.0;
  timer t1;

  pthread_t threads[numberOfThreads];
  threadObject threadHolder;
  threadStats* stats = new threadStats[numberOfThreads];
  threadHolder.threadStatistics = stats;
  threadHolder.writeMutex = new pthread_mutex_t;
  threadHolder.threadsRunning = 0;
  threadHolder.totalThreads = numberOfThreads;
  threadHolder.totalTriangles = 0;
  threadHolder.totalRuntime = 0;
  threadHolder.uniqueTriangles = 0;
  threadHolder.threadStatistics->threadTriangles = 0;
  pthread_mutex_init(threadHolder.writeMutex,NULL);

  uintV startIndex = 0;
  uintV endIndex = 0;
  uintV numberRemaining = g.n_;
  uintV numberPerThread = numberRemaining / numberOfThreads;
  int threadsRemaining = numberOfThreads;

  //std::cout << "Before Thread Spawn" << std::endl;

// u is the start, n is the finish
  for(int i = 0; i < numberOfThreads; i++){
    pthread_create(&threads[i],NULL,countTrianglesP,&threadHolder);
    threadHolder.threadStatistics[i].threadID = threads[i];
    threadHolder.threadStatistics[i].u = startIndex;
    endIndex = startIndex + numberPerThread;
    if(threadsRemaining > 1){
      //std::cout << "EndIndex = " << endIndex << std::endl;
      threadHolder.threadStatistics[i].n = endIndex;
      //std::cout << "thread " << i << " End Index in Thread = " << threadHolder.threadStatistics[i].n << std::endl;
    }
    else{
      //std::cout << "EndIndex = " << endIndex << std::endl;
      threadHolder.threadStatistics[i].n = g.n_;
      //std::cout << "thread " << i << " End Index in Thread = " << threadHolder.threadStatistics[i].n << std::endl;
    }
    startIndex += numberPerThread;
    threadsRemaining--;
  }

  //std::cout << "HERE2" << std::endl;
  //std::cout << "Number of total nodes = " << g.n_ << std::endl;

    for(int i = 0; i < numberOfThreads; i++){
      
        pthread_join(threads[i],NULL);
        //std::cout << "JOINED : " << i << std::endl;
    }

  //std::cout << "HERE3" << std::endl;


  // The outNghs and inNghs for a given vertex are already sorted

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  // Process each edge <u,v>
  //for (uintV u = 0; u < n; u++) {
    // For each outNeighbor v, find the intersection of inNeighbor(u) and
    // outNeighbor(v)
   // uintE out_degree = g.vertices_[u].getOutDegree();
   // for (uintE i = 0; i < out_degree; i++) {
    //  uintV v = g.vertices_[u].getOutNeighbor(i);
    //  triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
    //                                   g.vertices_[u].getInDegree(),
    //                                   g.vertices_[v].getOutNeighbors(),
    //                                   g.vertices_[v].getOutDegree(), u, v);
   // }
  //}
  time_taken = t1.stop();
  // -------------------------------------------------------------------
  // Here, you can just print the number of non-unique triangles counted by each
  // thread std::cout << "thread_id, triangle_count, time_taken\n"; Print the
  // above statistics for each thread Example output for 2 threads: thread_id,
  // triangle_count, time_taken 1, 102, 0.12 0, 100, 0.12

  for(int i = 0; i < numberOfThreads; i++){
    std::cout << i << ", " << threadHolder.threadStatistics[i].threadTriangles << ", " << threadHolder.threadStatistics[i].threadRuntime << std::endl;
  }

  // Print the overall statistics
  std::cout << "Number of triangles : " << threadHolder.totalTriangles << "\n";
  std::cout << "Number of unique triangles : " << threadHolder.totalTriangles / 3 << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << threadHolder.totalRuntime << "\n";
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "triangle_counting_serial",
      "Count the number of triangles using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  numberOfThreads = n_workers;
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";
  //Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";

  triangleCountSerial(g);

  return 0;
}
