#include "MRNG.h"

MRNG::MRNG(double **all_p1, int how_many1, int dimensions1 = 784, int k1 = 50, int N1 = 1, int R1 = 1, int E1 = 30, int L1 = 20){
  this->all_p = new double*[how_many1];
  for(int i = 0; i < how_many1; i++)
    this->all_p[i] = new double[dimensions1];

  for (int i = 0; i < how_many1; i++)
    memcpy(this->all_p[i], all_p1[i], dimensions1*sizeof(double));

  this->dimensions = dimensions1;
  this->how_many = how_many1;
  this->k = k1;
  this->N = N1;
  this->R = R1;
  this->E = E1;
  this->L = L1;
  this->edges = new vector<POINTDISTANCE>[how_many1]; // To store the edges every point has
  GraphCreate();
}

MRNG::~MRNG(){ // For leaks
  for(int i = 0; i < how_many; i++)
    delete[] all_p[i];
  delete[] all_p;
  delete[] edges;
}

double MRNG::EuclidianDistance(double *p, double *q, int dimensions){ // Calculate distance between input point p and query point q
  double dist = 0;
  for (int i = 1; i < dimensions; i++)
    dist += pow(p[i] - q[i], 2);
  return sqrt(dist);
}

void MRNG::GraphCreate(){ // MRNG Construction
  // FILE NAME PURPOSEFULLY ALTERED IN ORDER TO RUN main3.cpp
  FILE * graph_file = fopen("mrnggrapfedsh.txt", "r"); // Check if file already exists
  // CHANGE BACK TO "gnnsgraph.txt" for normal results
  if (!graph_file){ // If not, create it now
    graph_file = fopen("mrnggraph.txt", "w");
    // Pseudoalgorithm's S set is our this->all_p array
    vector <POINTDISTANCE> Lp[this->how_many];
    for (int i = 0; i < this->how_many; i++){

      vector <POINTDISTANCE> Rp;
      for (int j = 0; j < this->how_many; j++){
        if (j != i){
          POINTDISTANCE p;
          p.point = j;
          p.distance = EuclidianDistance(this->all_p[j], this->all_p[i], this->dimensions);
          Rp.push_back(p);
        }
      }
      sort(Rp.begin(), Rp.end(), PDCompareFunctionDistance());
      Lp[i].push_back(Rp[0]);
      int min = Rp[0].distance;
      for (int q = 1; q < Rp.size(); q++){
        if (min == Rp[q].distance){ // Store in Lp all Rp elements with the minimum distance
          Lp[i].push_back(Rp[q]);
        } else {
          break;
        }
      }
      int j = 0;
      for (j = 0; j < Rp.size(); j++){

        bool condition = true;
        if (count_if(Lp[i].begin(), Lp[i].end(), PDCompareFunctionEqualPoint(Rp[j].point)) == 0){ // If Rp[j].point doesn't exist in Lp[i]
          int Lpsize = Lp[i].size();
          for (int q = 0; q < Lpsize; q++){
            double pr = Rp[j].distance;
            double rt = EuclidianDistance(this->all_p[Lp[i][q].point], this->all_p[Rp[j].point], this->dimensions);
            double pt = Lp[i][q].distance;
            if ((pr >= rt) && (pr >= pt)){
              condition = false;
              break;
            }
          }
        } else {
          continue;
        }
        if (condition == true)
          Lp[i].push_back(Rp[j]);
      }
    }
    for (int i = 0; i < this->how_many; i++){
      for (int j = 0; j < Lp[i].size(); j++){
        POINTDISTANCE p;
        p.point = Lp[i][j].point;
        p.distance = Lp[i][j].distance;
        this->edges[i].push_back(p);
      }
    }
    for (int i = 0; i < this->how_many; i++){
      fprintf(graph_file, "%ld\n", this->edges[i].size());
      for (int j = 0; j < this->edges[i].size(); j++)
        fprintf(graph_file, "%d %f\n", (this->edges[i][j]).point, (this->edges[i][j]).distance);
    }
    fclose(graph_file);
  } else { // If graph text file pre-exists, load it in array of vectors this->edges
    ifstream graph;
    graph.open("mrnggraph.txt");
    
    for (int i = 0; i < this->how_many; i++){
      string firstline;
      getline(graph, firstline);
      string neighbors;
      stringstream getwords(firstline);
      getwords >> neighbors;
      int count_neighbors = 0;
      int neighborsi = atoi(neighbors.c_str());
      while (count_neighbors < neighborsi){
        POINTDISTANCE p;
        count_neighbors++;
        string line;
        getline(graph, line);
        string current;
        stringstream getwords(line);
        getwords >> current;
        p.point = atoi(current.c_str());
        getwords >> current;
        p.distance = atof(current.c_str());
        this->edges[i].push_back(p);
      }
    }
    graph.close();
  }
}

vector <POINTDISTANCE> MRNG::SearchGraph(double * q){ // Search On Graph algorithm from slides
  vector <POINTDISTANCECHECK> R; // Use POINTDISTANCECHECK which also has an int variable "checked", equal 0 if not-checked and equal 1 if checked.
  const int range_from = 0, range_to = this->how_many - 1;
  random_device rand_dev;
  mt19937 generator(rand_dev());
  uniform_int_distribution<int> distr(range_from, range_to); // Choose point uniformly
  int start_node = 0;
  start_node = distr(generator);
  POINTDISTANCECHECK p1;
  p1.point = start_node;
  p1.distance = EuclidianDistance(this->all_p[start_node], q, this->dimensions);
  p1.checked = 0;
  R.push_back(p1);
  int i = 1;
  while (i < this->L){
    POINTDISTANCECHECK p;
    for (int j = 0; j < R.size(); j++){
      if (R[j].checked == 0){ // If not checked
        R[j].checked = 1;
        p = R[j];
        break;
      }
    }
    for (int j = 0; j < edges[p.point].size(); j++){
      if (count_if(R.begin(), R.end(), PDCompareFunctionEqualPointCheck(edges[p.point][j].point)) == 0){ // If no edge exists
        POINTDISTANCECHECK p2;
        p2.point = edges[p.point][j].point;
        p2.distance = EuclidianDistance(this->all_p[edges[p.point][j].point], q, this->dimensions);
        p2.checked = 0;
        R.push_back(p2);
        i++;
      }
      sort(R.begin(), R.end(), PDCompareFunctionDistanceCheck());
    }
  }
  vector <POINTDISTANCE> results;
  for(int i1 = 0; i1 < this->k; i1++){ // To return only k neighbors
    POINTDISTANCE point;
    point.point = R[i1].point;
    point.distance = R[i1].distance;
    results.push_back(point);
  }
  return results;
}