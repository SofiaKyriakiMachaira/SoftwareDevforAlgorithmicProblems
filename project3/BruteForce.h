// Basic header file for BruteForce.cpp

#include "common.h" // which has struct POINTDISTANCE
#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <climits>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include <cstring>

class BRUTEFORCE {
  private:
    double **all_p; // All points
    int dimensions; // Dimensions (28*28)
    int N; // Amount of neighbors to find
    int r; // Radius
    int how_many; // Number of images

  double EuclidianDistance(double*, double*, int);
  public:
    BRUTEFORCE(double**, int, int, int, int); // Constructor
    ~BRUTEFORCE(); // Destructor
    vector<POINTDISTANCE> NNearestNeighbors(double*);
    vector<POINTDISTANCE> RangeSearch(double*);
};