// Header file for MRNG.cpp
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
using namespace std;

// Make new POINTDISTANCE that also includes whether or not this point has been checked

struct POINTDISTANCECHECK { // Saves position of point and its euclidian distance
  int point;
  double distance;
  int checked;
};

struct PDCompareFunctionDistanceCheck { // Compares distance between 2 points
    inline bool operator() (const POINTDISTANCECHECK& pd1, const POINTDISTANCECHECK& pd2){
        return (pd1.distance < pd2.distance);
    }
};

struct PDCompareFunctionPointCheck { // Compares position between 2 points
    inline bool operator() (const POINTDISTANCECHECK& pd1, const POINTDISTANCECHECK& pd2){
        return (pd1.point < pd2.point);
    }
};

struct PDCompareFunctionEqualPointCheck { // Checks equality between two points of POINTDISTANCE
  explicit PDCompareFunctionEqualPointCheck(const int& s) : point(s) {}
  inline bool operator() (const POINTDISTANCECHECK& pd1){
    return (pd1.point == point);
  }
  int point;
};


class MRNG{
  private:
  double **all_p; // All points
  int how_many;
  int dimensions;
  int k; // Parameters
  int R;
  int E;
  int N;
  int L;
  vector<POINTDISTANCE> *edges; // Graph is stored here

  double EuclidianDistance(double*, double*, int);
  void GraphCreate(); // Graph construction
  public:
  MRNG(double**, int, int, int, int, int, int, int); // Constructor
  ~MRNG(); // Destructor
  vector<POINTDISTANCE> SearchGraph(double *); // Search On Graph algorithm
};