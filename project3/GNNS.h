// Header file for Graph Nearest Neighbor.cpp
#include "LSH.h"
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
#include <cstdio>
#include <fstream>
#include <bits/stdc++.h>

class GNNS{
  private:
  double **all_p; // All points
  int how_many; 
  int dimensions;
  int k; // Parameters
  int R;
  int E;
  int N;
  vector<POINTDISTANCE> *edges; // The graph is stored here

  double EuclidianDistance(double*, double*, int);
  void GraphCreate(); // Construction of graph
  public:
  GNNS(double**, int, int, int, int, int, int); // Constructor
  ~GNNS(); // Destructor
  vector<POINTDISTANCE> GraphNNSearch(double *); // GNNS algorithm
};
