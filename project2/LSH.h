// Header file for LSH.cpp

#include "common.h"
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

class LSH {
  private:
  double **all_p; // All points
  int table_size;
  int w;
  int *g_temp; // Temporary g
  vector<int> *gvector; // Vector of buckets of Gs
  vector<int> *idvector; // IDs
  int *id_temp; // Temporary id
  int dimensions;
  int N; // For NearestNeighbors
  int k;
  int r; // Radius
  int l;
  int how_many; // Number of images

  double EuclidianDistance(double*, double*, int);
  void HashTable(int , int , int , double*); // Creating hashtables

  public:
  LSH(double**, int, int, int , int, int , int);
  ~LSH(); // Destructor
  vector<POINTDISTANCE> NNearestNeighbors(double*,int );
  vector<POINTDISTANCE> RangeSearch(double*, int);
};