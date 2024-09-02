#include "LSH.h"
using namespace std;

// Constructor sets default values too (just in case) or takes given parameters and creates a copy of all input points
LSH::LSH(double **all_p1, int how_many1, int N1 = 1, int dimensions1 = 784, int k1 = 4, int r1 = 10000, int l1 = 5){
  this->all_p = new double*[how_many1];
  for(int i = 0; i < how_many1; i++)
    this->all_p[i] = new double[dimensions1];

  for (int i = 0; i < how_many1; i++)
    memcpy(this->all_p[i], all_p1[i], dimensions1*sizeof(double));

  this->dimensions = dimensions1;
  this->k = k1;
  this->r = r1;
  this->l = l1;
  this->N = N1;
  this->how_many = how_many1;
  this->w = 100;
  this->table_size = this->how_many / 1000;
  this->gvector = new vector<int>[this->table_size];
  this->idvector = new vector<int>[this->table_size];
  this->id_temp = new int[this->how_many];
  this->g_temp = new int[this->how_many];
  for (int i = 0; i < this->l; i++){
    HashTable(this->w, this->how_many, this->k, NULL); // Call hashtable to create every h, g, id function
  }
}

LSH::~LSH(){ // For leaks
  for(int i = 0; i < how_many; i++)
    delete[] all_p[i];

  delete[] id_temp;
  delete[] g_temp;
  delete[] all_p;
  delete[] gvector;
  delete[] idvector;
}

double LSH::EuclidianDistance(double *p, double *q, int dimensions){ // Calculate distance between input point p and query point q
  double dist = 0;
  for (int i = 1; i < dimensions; i++)
    dist += pow(p[i] - q[i], 2);
  return sqrt(dist);
}

void LSH::HashTable(int w, int n, int k, double *q) {
  int h[n][k];
  unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
  default_random_engine e(seed);
  random_device generator;
  normal_distribution<double> distribution (0.0, 1.0);
  double v[this->dimensions]; // To store normal distribution
  for(int i = 0; i < this->dimensions; i++)
    v[i] = distribution(generator);
  int w1 = w - 1;
  uniform_real_distribution<> U(0.0, w1); // For t
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < k; i++){

      double t = U(e);
      double innerproduct;
      if(q == NULL){ // If no query point has been given
        innerproduct = inner_product(this->all_p[j], this->all_p[j] + this->dimensions, v, 0.0);
      } else { // If query point has been given
        innerproduct = inner_product(q, q + this->dimensions, v, 0.0);
      }
      h[j][i] = floor((abs(innerproduct) + t) / w); // [((p*v)+t)/w]
    }
  }
  // Classic hash function
  int M = UINT_MAX - 5;
  int r[k];
  int id[n]; // For querying trick
  for (int i = 0; i < k; i++)
    r[i] = rand();
  int g[n];
  for (int j = 0; j < n; j++){
    int current = 0;
    for (int i = 0; i < k; i++)
      current += (r[i] * h[j][i]);
    id[j] = abs(current) % M;
    g[j] = id[j] % table_size;
    if(q == NULL){ // If no query point has been given
      if (count(gvector[g[j]].begin(), gvector[g[j]].end(), j) == 0) // If j isn't in gvector 
        gvector[g[j]].push_back(j); // Add it in both gvector and idvector
        idvector[g[j]].push_back(id[j]);
    }
  }
  memcpy(this->id_temp, id, n*sizeof(int));
  memcpy(this->g_temp, g, n*sizeof(int));
}
// custom N every time to add if needed
vector<POINTDISTANCE> LSH::NNearestNeighbors(double* q, int N1){
  vector <int> number_of_index_equal_id; // Index with equal ID as q in input file points
  vector <int> number_of_index_equal_g; // Index with equal G as q in input file points
  vector <POINTDISTANCE> point_distance;
  if (N1 > 0)
    this->N = N1;
  
  HashTable(this->w, 1, this->k, q); // Call hashtable to find G and ID of q
  int g_q = this->g_temp[0];
  int id_q = this->id_temp[0];

    for(int j = 0; j < gvector[g_q].size(); j++){
      if (idvector[g_q][j] == id_q){
          number_of_index_equal_id.push_back(gvector[g_q][j]);
      } else {
          number_of_index_equal_g.push_back(gvector[g_q][j]);
      }
  }
  // Only calculate distance between q and points with equal ID or G as q
  for(int i = 0; i < number_of_index_equal_id.size(); i++){
    POINTDISTANCE pd;
    pd.distance = EuclidianDistance(this->all_p[number_of_index_equal_id[i]], q, this->dimensions);
    pd.point = number_of_index_equal_id[i];
    point_distance.push_back(pd);
  }

  if(point_distance.size() < this->N){ // If we haven't exceeded the amount of neighbors required
    for(int i = 0; i < number_of_index_equal_g.size(); i++){
      POINTDISTANCE pd;
      pd.distance = EuclidianDistance(all_p[number_of_index_equal_g[i]], q, this->dimensions);
      pd.point = number_of_index_equal_g[i];
      point_distance.push_back(pd);
    }
  }

  sort(point_distance.begin(), point_distance.end(), PDCompareFunctionDistance()); // Sort based on distance
  vector <POINTDISTANCE> results;
  int N2 = this->N;
  for(int i = 0; ((i < N2) && (i < point_distance.size()) && (results.size() < this->N)); i++){
    if ((N1 == -1) && (point_distance[i].distance == (double)0.0)){
      N2++;
    } else {
      results.push_back(point_distance[i]);
    }
  }
  return results;
}

vector<POINTDISTANCE> LSH::RangeSearch(double * q, int R){
  vector <int> number_of_index_equal_g;
  vector <POINTDISTANCE> point_distance;
  if (R > 0) // To change radius (used for clustering)
    this->r = R;

  HashTable(this->w, 1, this->k, q); // Same technique as NNearestNeighbors
  int g_q = this->g_temp[0];
  int id_q = this->id_temp[0];

  for(int j = 0; j < gvector[g_q].size(); j++)
    number_of_index_equal_g.push_back(gvector[g_q][j]);

  for(int i = 0; i < number_of_index_equal_g.size(); i++){
    POINTDISTANCE pd;
    pd.distance = EuclidianDistance(all_p[number_of_index_equal_g[i]], q, this->dimensions);
    if(pd.distance < (double)(this->r)){ // Only store points within radius
      pd.point = number_of_index_equal_g[i];
      point_distance.push_back(pd);
    }
  }

  sort(point_distance.begin(), point_distance.end(), PDCompareFunctionDistance()); // Sort points based on distance
  return point_distance;
}
