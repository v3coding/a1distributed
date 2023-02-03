#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_POINTS "12345678"

static int numberOfThreads;
static int points;

typedef struct threadStats{
    pthread_t threadID;
    int threadProduced;
    int threadHits;
    double threadTime;
    double runtime;
} threadStats;

typedef struct threadObject{
    pthread_mutex_t* writeMutex;
    pthread_t* threads;
    int retval;
    int n;
    int T;
    int pointsCalculated;
    threadStats* threadStatistics;
}   threadObject;

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;
}

uint get_points_in_circle(uint n, uint random_seed) {
  uint circle_count = 0;
  double x_coord, y_coord;
  for (uint i = 0; i < n; i++) {
    x_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    y_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    if ((sqr(x_coord) + sqr(y_coord)) <= 1.0)
      circle_count++;
  }
  return circle_count;
}

void* getPoints(void *_arg) {
    timer t1;
    t1.start();
    threadObject* threadData = (threadObject*) _arg;
    pthread_mutex_lock(threadData->writeMutex);
    int n = threadData->n;
    int T = threadData->T;
    pthread_mutex_unlock(threadData->writeMutex);
    
    int pointsCalculated = 0;
    uint random_seed = rand();

  uint circle_count = 0;
  double x_coord, y_coord;
  for (uint i = 0; i < n/T; i++) {
    x_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    y_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    if ((sqr(x_coord) + sqr(y_coord)) <= 1.0){
      circle_count++;
    }
    pointsCalculated++;
  }
    pthread_mutex_lock(threadData->writeMutex);
    threadData->retval += circle_count;
    threadData->pointsCalculated += pointsCalculated;
    pthread_mutex_unlock(threadData->writeMutex);

    for(int i = 0; i < numberOfThreads; i++){
        if(pthread_self() == threadData->threadStatistics[i].threadID){
        threadData->threadStatistics[i].runtime = t1.stop();
        threadData->threadStatistics[i].threadHits = circle_count;
        threadData->threadStatistics[i].threadProduced = pointsCalculated;
    }
  }

    return 0;
}

void piCalculation(uint n) {
  timer serial_timer;
  double time_taken = 0.0;
  uint random_seed;

  // Create threads and distribute the work across T threads
  threadObject threadHolder;
  threadHolder.writeMutex = new pthread_mutex_t;
  threadStats* stats = new threadStats[numberOfThreads];
  threadHolder.threadStatistics = stats;
  pthread_mutex_init(threadHolder.writeMutex,NULL);
  pthread_t threads[numberOfThreads];
  threadHolder.n = points;
  //std::cout << "HERE : " << threadHolder.n << std::endl;
  threadHolder.T = numberOfThreads;
  threadHolder.retval = 0;

//
  serial_timer.start();
  // ------------------------------------------------------------------
  for(int i = 0; i < numberOfThreads; i++){
    pthread_create(&threads[i],NULL,getPoints,&threadHolder);
    threadHolder.threadStatistics[i].threadID = threads[i];
  }
  for(int i = 0; i < numberOfThreads; i++){
    pthread_join(threads[i],NULL);
  }
  //uint circle_points = get_points_in_circle(n, random_seed);
  double pi_value = 4.0 * (double)threadHolder.retval / (double)threadHolder.n;
  // -------------------------------------------------------------------
  time_taken = serial_timer.stop();

  std::cout << "thread_id, points_generated, circle_points, time_taken\n";
  for(int i = 0; i < numberOfThreads; i++){
    std::cout << i << ", " << threadHolder.threadStatistics[i].threadProduced << ", " << threadHolder.threadStatistics[i].threadHits << ", " << std::setprecision(TIME_PRECISION) << threadHolder.threadStatistics[i].runtime << std::endl;
  }
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, points_generated, circle_points, time_taken
  // 1, 100, 90, 0.12
  // 0, 100, 89, 0.12

  // Print the overall statistics
  std::cout << "Total points generated : " << threadHolder.pointsCalculated << "\n";
  std::cout << "Total points in circle : " << threadHolder.retval << "\n";
  std::cout << "Result : " << std::setprecision(VAL_PRECISION) << pi_value
            << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";

  delete threadHolder.writeMutex;
  delete stats;
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("pi_calculation",
                           "Calculate pi using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nPoints", "Number of points",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_points = cl_options["nPoints"].as<uint>();
  uint n_workers = cl_options["nWorkers"].as<uint>();
  std::cout << std::fixed;
  std::cout << "Number of points : " << n_points << "\n";
  std::cout << "Number of workers : " << n_workers << "\n";
  numberOfThreads = n_workers;
  points = n_points;

  piCalculation(n_points);

  return 0;
}
