#include "BruteForce.h"
using namespace std;

BRUTEFORCE::BRUTEFORCE(double **all_p1, int how_many1, int N1, int dimensions1, int r1){ // Set parameters, create a copy of input points
  this->all_p = new double*[how_many1];
  for(int i = 0; i < how_many1; i++)
    this->all_p[i] = new double[dimensions1];

  for (int i = 0; i < how_many1; i++)
    memcpy(this->all_p[i], all_p1[i], dimensions1*sizeof(double));
  
  this->dimensions = dimensions1;
  this->N = N1;
  this->r = r1;
  this->how_many = how_many1;
}

BRUTEFORCE::~BRUTEFORCE(){ // Free points
  for(int i = 0; i < how_many; i++)
    delete[] all_p[i];

  delete[] all_p;
}

double BRUTEFORCE::EuclidianDistance(double *p, double *q, int dimensions){ // Calculate distance between 2 points and return it
  double dist = 0;
  for (int i = 1; i < dimensions; i++)
    dist += pow(p[i] - q[i], 2);

  return sqrt(dist);
}

vector<POINTDISTANCE> BRUTEFORCE::NNearestNeighbors(double *q){ // Exhaustive search
  vector <POINTDISTANCE> point_distance;
  for(int i = 0; i < this->how_many; i++){
    POINTDISTANCE pd;
    pd.distance = EuclidianDistance(this->all_p[i], q, this->dimensions);
    pd.point = i;
    point_distance.push_back(pd);
  }
  sort(point_distance.begin(), point_distance.end(), PDCompareFunctionDistance()); // Sort based on distance and return the closest
  vector <POINTDISTANCE> results;
  for(int i = 0; i < this->N; i++)
    results.push_back(point_distance[i]);
  return results;
}

vector<POINTDISTANCE> BRUTEFORCE::RangeSearch(double *q){ // Exhaustive search
  vector <POINTDISTANCE> point_distance;
  for(int i = 0; i < this->how_many; i++){
    POINTDISTANCE pd;
    pd.distance = EuclidianDistance(this->all_p[i], q, this->dimensions);
    pd.point = i;
    if (pd.distance < (double)(this->r)) // Within radius
      point_distance.push_back(pd);
  }

  sort(point_distance.begin(), point_distance.end(), PDCompareFunctionDistance()); // Sort based on distance
  return point_distance;
}