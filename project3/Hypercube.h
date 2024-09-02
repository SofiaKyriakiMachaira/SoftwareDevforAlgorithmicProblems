// Header file for Hypercube.cpp
#include "common.h"
#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <climits>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <list>
#include <cstring>
#include <cstdlib>

class HYPERCUBE{
  private :
    class H10{
    public:
      vector<int> h;
      vector<int> to10; // To 1 or 0
    };

    H10 hto10;
    double **all_p; // All points
    vector <int> * fvector; // F function
    int dimensions; // Parameters
    int N;
    int k;
    int M;
    int r; // Radius
    int probes;
    int w;
    int how_many; // Number of images
    int *fq; // F function for query
    int fqdecimal;
    int check_m;
    vector <int> previus_probe;
    vector<int> FindProbe(int, int, vector<int>, int*, int*, int); // Searches the next neighbor of original probe
    vector<int> ProbeSearch(int, vector<int>, int); // Searches within probe its neighbors
    double EuclidianDistance(double*, double*, int);
    int CheckH10(int);
    void HashTable(int, int, int, double*); // Create hashtables
  public :
    HYPERCUBE(double**, int, int, int, int, int, int, int);
    ~HYPERCUBE();
    vector<POINTDISTANCE> NNearestNeighbors(double*,int); // RandomProjections
    vector<POINTDISTANCE> RangeSearch(double*, int);
};