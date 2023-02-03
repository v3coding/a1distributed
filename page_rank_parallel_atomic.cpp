#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif

static int numberOfThreads;
static int maxIterations;
static Graph g;

typedef struct threadStats{
    pthread_t threadID;
    PageRankType threadPageRank;
    double threadRuntime;
    uintV u;
    uintV v;
    uintV n;
} threadStats;

typedef struct threadObject{
    pthread_mutex_t* writeMutex;
    PageRankType totalPageRank;
    CustomBarrier* barrier;
    double totalRuntime;
    threadStats* threadStatistics;
    int resetOnce;
    PageRankType *pr_curr;
    PageRankType *pr_next;
}   threadObject;

void* pageRankParallel(void *_arg){
  timer t1;
  t1.start();

 // std::cout << pthread_self() << std::endl;

  double time_taken = 0.0;
  int index;
  threadObject* threadData = (threadObject*) _arg;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------

    for(int i = 0; i < numberOfThreads; i++){
        if(pthread_self() == threadData->threadStatistics[i].threadID){
            index = i;
        }
    }
   // std::cout << "Thread " << index << " alive.. " << std::endl;

    //threadData->resetOnce = 0;

    uintV u = threadData->threadStatistics[index].u;
    uintV n = threadData->threadStatistics[index].n;

  PageRankType *pr_curr = threadData->pr_curr;
  PageRankType *pr_next = threadData->pr_next;

 // std::cout << "HERE" << std::endl;

  
  //PageRankType *pr_currLocal = new PageRankType[n];
  PageRankType *pr_nextLocal = new PageRankType[n];

 // std::cout << "HERE2" << std::endl;

 for (uintV i = 0; i < n; i++) {
   // pr_currLocal[i] = INIT_PAGE_RANK;
    pr_nextLocal[i] = 0.0;
  }

  //std::cout << "HERE3" << std::endl;
  for (int iter = 0; iter < maxIterations; iter++) {
   // std::cout << "Iteration " << iter << std::endl;
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = index; u <= n - numberOfThreads; u += numberOfThreads) {
      uintE out_degree = g.vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g.vertices_[u].getOutNeighbor(i);
        pr_nextLocal[v] += (pr_curr[u] / out_degree);
      }
      //threadData->barrier->wait();
    }
    for (uintV v = 0; v < n; v++) {
      pr_nextLocal[v] = PAGE_RANK(pr_nextLocal[v]);
      // reset pr_curr for the next iteration
      pr_curr[v] = pr_nextLocal[v];
      pr_nextLocal[v] = 0.0;
    }
  }
 // std::cout << "Breaking Main Loop " << std::endl;

  threadData->barrier->wait();

  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
    time_taken = t1.stop();
    pthread_mutex_lock(threadData->writeMutex);
    threadData->threadStatistics[index].threadPageRank = sum_of_page_ranks;
    threadData->totalPageRank += sum_of_page_ranks;
    //threadData->totalRuntime += time_taken;
    threadData->threadStatistics[index].threadRuntime = time_taken;
    pthread_mutex_unlock(threadData->writeMutex);
    return 0;
}

void pageRankSerial(Graph &g, int max_iters) {
  uintV n = g.n_;

  timer t;
  t.start();

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

threadObject threadHolder;
  threadHolder.writeMutex = new pthread_mutex_t;
  threadHolder.barrier = new CustomBarrier(numberOfThreads);
  threadStats* stats = new threadStats[numberOfThreads];
  threadHolder.threadStatistics = stats;
  threadHolder.totalPageRank = 0;
  threadHolder.totalRuntime = 0;
  pthread_mutex_init(threadHolder.writeMutex,NULL);
  pthread_t threads[numberOfThreads];

  threadHolder.pr_curr = pr_curr;
  threadHolder.pr_next = pr_next;

// u is the start, n is the finish
  for(int i = 0; i < numberOfThreads; i++){
    pthread_create(&threads[i],NULL,pageRankParallel,&threadHolder);
    threadHolder.threadStatistics[i].threadID = threads[i];
    threadHolder.threadStatistics[i].n = n;
  }

  for(int i = 0; i < numberOfThreads; i++){
    pthread_join(threads[i],NULL);
  //  std::cout << "Thread " << i << " Joined " << std::endl;
  }

  // -------------------------------------------------------------------
  // std::cout << "thread_id, time_taken\n";
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12
  for(int i = 0; i < numberOfThreads; i++){
    std::cout << i << ", " << threadHolder.threadStatistics[i].threadRuntime << ", " << threadHolder.threadStatistics[i].threadPageRank << std::endl;
  }

  threadHolder.totalRuntime = t.stop();

  std::cout << "Sum of page rank : " << threadHolder.totalPageRank << "\n";
  std::cout << "Time taken (in seconds) : " << threadHolder.totalRuntime << "\n";
  delete[] pr_curr;
  delete[] pr_next;

}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_push",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  numberOfThreads = n_workers;
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using FLOAT\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";

  //Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  maxIterations = max_iterations;
  pageRankSerial(g, max_iterations);

  return 0;
}
