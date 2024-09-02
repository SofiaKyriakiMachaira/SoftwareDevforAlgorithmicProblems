// Header file for CentralClustering.cpp

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
#include <iterator>
#include <cfloat>

class TEAMSCENTRAL{ // To store results (points/cluster leaders + euclidian distances of points within its cluster)
  public:
    vector<POINTDISTANCE> *teams_euclidian;
    vector<double*> leaders;
};

class CENTRALCLUSTERING{
  private:
    double **all_p; // All points
    int dimensions;
    int k;
    int how_many; // Number of images
    TEAMSCENTRAL centrals;
    double EuclidianDistance(double *, double *, int);
    int FindSum(int);
    bool CheckIfDifferent(int *,int);
    double* AveragePoint(vector <POINTDISTANCE>);
    double AverageDistance(vector <POINTDISTANCE>, int);
    double AverageSilhouette(double*);
    TEAMSCENTRAL Initialization(); // Initializing the first random clusters
  public:
    CENTRALCLUSTERING(double **, int, int, int);
    ~CENTRALCLUSTERING();
    TEAMSCENTRAL Reverse(int, int, int, int, int, int);
    TEAMSCENTRAL Classic();
    double Silhouette(TEAMSCENTRAL, FILE*); // Can also write in output file, returns stotal, average silhouette
};