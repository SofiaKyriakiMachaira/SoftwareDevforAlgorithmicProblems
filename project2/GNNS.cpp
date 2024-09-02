#include "GNNS.h"
//#include "Hypercube.h"

using namespace std;

GNNS::GNNS(double **all_p1, int how_many1, int dimensions1 = 784, int k1 = 50, int N1 = 1, int R1 = 1, int E1 = 30){
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
  this->edges = new vector<POINTDISTANCE>[how_many1]; // Create vector to store the edges each point has
  GraphCreate();
}

GNNS::~GNNS(){ // For leaks
  for(int i = 0; i < how_many; i++)
    delete[] all_p[i];
  delete[] all_p;
  delete[] edges;
}

double GNNS::EuclidianDistance(double *p, double *q, int dimensions){ // Calculate distance between input point p and query point q
  double dist = 0;
  for (int i = 1; i < dimensions; i++)
    dist += pow(p[i] - q[i], 2);
  return sqrt(dist);
}

void GNNS::GraphCreate(){ // Opens a file and stores the graph inside for easier access upon re-running it
  FILE * graph_file = fopen("gnnsgraph.txt", "r");
  if (!graph_file){
    LSH lsh(this->all_p, this->how_many, this->k, this->dimensions, 4, 10000, 5); // Default values
    // HYPERCUBE lsh(this->all_p, this->how_many, this->k, 100, this->dimensions, 14, 10000, 10); // To run with Hypercube
    graph_file = fopen("gnnsgraph.txt", "w");
    for (int i = 0; i < this->how_many; i++){
      lsh.NNearestNeighbors(this->all_p[i], -1).swap(this->edges[i]);
      fprintf(graph_file, "%ld\n", this->edges[i].size());
      for (int j = 0; j < this->edges[i].size(); j++)
        fprintf(graph_file, "%d %f\n", (this->edges[i][j]).point, (this->edges[i][j]).distance);
    }
    fclose(graph_file);
  }
  else {
    ifstream graph;
    graph.open("gnnsgraph.txt");
    
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

vector<POINTDISTANCE> GNNS::GraphNNSearch(double * q){
  int T = 50;
  const int range_from = 0, range_to = how_many - 1;
  random_device rand_dev;
  mt19937 generator(rand_dev());
  uniform_int_distribution<int> distr(range_from, range_to);
  int Ynow = 0;
  vector <POINTDISTANCE> S;
  for(int ra = 0; ra < this->R; ra++){
    POINTDISTANCE p;
    Ynow = distr(generator);
    p.point = Ynow;
    p.distance = EuclidianDistance(this->all_p[Ynow], q, this->dimensions);
    double Ynow_Distance = 0.0;
    double YPrevious = p.distance;
    if (count_if(S.begin(), S.end(), PDCompareFunctionEqualPoint(p.point)) == 0)
      S.push_back(p);
    for (int times = 1; ((times <= T) && (abs(Ynow_Distance - YPrevious) >= (0.01*YPrevious))); times++){ 
      // T has not been exceeded and difference between new neighbor is more than 1/100
      vector<POINTDISTANCE> min;
      for(int i = 0; i < this->E && i < this->edges[Ynow].size(); i++){
        auto it = find_if(S.begin(), S.end(), PDCompareFunctionEqualPoint(edges[Ynow][i].point)); // Find if point exists in S
        if (it != S.end()){ // if it exists, add it to vector min
          min.push_back((*it)); 
        } else { // If not, add it
          POINTDISTANCE po;
          po.point = edges[Ynow][i].point;
          po.distance = EuclidianDistance(this->all_p[edges[Ynow][i].point], q, this->dimensions);
          S.push_back(po);
          min.push_back(po);
        }
      }
      sort(min.begin(), min.end(), PDCompareFunctionDistance()); // Sort points based on distance
      YPrevious = Ynow_Distance; 
      Ynow = min[0].point;
      Ynow_Distance = min[0].distance;
    }
  }
  sort(S.begin(), S.end(), PDCompareFunctionDistance()); // Sort points based on distance
  vector <POINTDISTANCE> results;
  for(int i = 0; i < this->N; i++) // To return only N neighbors
    results.push_back(S[i]);
  return results;
}