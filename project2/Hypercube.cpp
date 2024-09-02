#include "Hypercube.h"
using namespace std;

HYPERCUBE::HYPERCUBE(double **all_p1, int how_many1, int N1, int M1, int dimensions1, int k1, int r1, int probes1){
  this->all_p = new double*[how_many1];
  for(int i = 0; i < how_many1; i++)
    this->all_p[i] = new double[dimensions1];

  for (int i = 0; i < how_many1; i++)
    memcpy(this->all_p[i], all_p1[i], dimensions1*sizeof(double));
  // After copying all input points, set parameters and create hashtable
  int k2 = pow(2, k1);
  this->fvector = new vector<int>[k2];
  this->M = M1;
  this->fq = new int[k1];
  this->dimensions = dimensions1;
  this->k = k1;
  this->r = r1;
  this->probes = probes1;
  this->N = N1;
  this->how_many = how_many1;
  this->w = rand() % 6 + 2;
  HashTable(this->w, this->how_many, this->k, NULL);
}

HYPERCUBE::~HYPERCUBE(){ // For leaks
  for(int i = 0; i < how_many; i++)
    delete[] all_p[i];

  delete[] all_p;
  delete[] fq;
  delete[] fvector;
}

double HYPERCUBE::EuclidianDistance(double *p, double *q, int dimensions){
  double dist = 0;
  for (int i = 1; i < dimensions; i++)
    dist += pow(p[i] - q[i], 2);

  return sqrt(dist);
}

int HYPERCUBE::CheckH10(int h){ // Check if H has already a value (0 or 1 or set it)
  int ret = -1;
  auto it = find(hto10.h.begin(), hto10.h.end(), h);

  // If element was found, return it
  if (it != hto10.h.end())  {
    int index = it - hto10.h.begin();
    return hto10.to10[index];
  } else { // If not, do "coin flip" to set it to 0 or 1
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
    default_random_engine e(seed);
    uniform_int_distribution<> U(0,1);
    ret = U(e);
    hto10.h.push_back(h);
    hto10.to10.push_back(ret);
    return ret;
  }
}

void HYPERCUBE::HashTable(int w, int n, int k, double *q) {
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
    int binary = 0;
    int power = k - 1;
    for (int i = 0; i < k; i++){

      double t = U(e);
      double innerproduct;
      if(q == NULL){ // If no query point has been given
        innerproduct = inner_product(this->all_p[j], this->all_p[j] + this->dimensions, v, 0.0);
      } else { // If query point has been given
        innerproduct = inner_product(q, q + this->dimensions, v, 0.0);
      }

      h[j][i] = floor((abs(innerproduct) + t) / w); // [((p*v)+t)/w]
      if(q == NULL){ // If no query point has been given
        binary = CheckH10(h[j][i])*pow(2, power) + binary;
      } else { // If query point has been given
        this->fq[i] = CheckH10(h[j][i]);
        binary = CheckH10(h[j][i])*pow(2, power) + binary;
      }
      power--;
    }

    if(q == NULL){ // If no query point has been given
      fvector[binary].push_back(j);
    } else { // If query point has been given
      fqdecimal = binary;
    }
  }
}

vector<int> HYPERCUBE::FindProbe(int count_check, int count_max, vector<int> number_of_index_equal, int *probesfinder, int *i_have_check, int probe){
  // Searches the next neighbor of original probe
  for(int i = 0; i < this->k && this->check_m < this->M; i++){
    int newprobe = probe;
    if(count_check == count_max){ // Final change required
      if(i_have_check[i] == -1){ // If previously checked or not
          newprobe = newprobe + probesfinder[i]; // Makes new probe
          number_of_index_equal = ProbeSearch(count_max, number_of_index_equal, newprobe);
      }
    } else {
      if(i_have_check[i] == -1){
        newprobe = newprobe + probesfinder[i];
        i_have_check[i] = 1;
        count_check++;
        number_of_index_equal = FindProbe(count_check, count_max, number_of_index_equal, probesfinder, i_have_check, newprobe);
        count_check--;
        i_have_check[i] = -1;
      }
    }
  }
  return number_of_index_equal;
}

vector<int> HYPERCUBE::ProbeSearch(int count_max, vector<int> number_of_index_equal, int probe){ // Searches within probe its neighbors
  // M and probes have not been exceeded
  if(count(this->previus_probe.begin(), this->previus_probe.end(),probe) == 0){
    this->previus_probe.push_back(probe);
    for(int j = 0; j < this->fvector[probe].size() && this->check_m < this->M; j++){
      number_of_index_equal.push_back(this->fvector[probe][j]);
      this->check_m++;
    }
  }
  return number_of_index_equal;
}

// custom N every time to add if needed
vector<POINTDISTANCE> HYPERCUBE::NNearestNeighbors(double *q,int N1){
  HashTable(this->w, 1, this->k, q); // Call hashtable for query q
  vector <int> number_of_index_equal;
  vector <POINTDISTANCE> point_distance;
  int self = 0;
  if (N1 > 0){
    this->N = N1;
  } else if (N1 < 0){
    self = 1;
  }

  int probefinder[this->k];
  int i_have_check[this->k];
  int power = this->k - 1;
  for(int i = 0; i < this->k; i++){
    if(fq[i] == 0) {
      probefinder[i] = pow(2, power);
    } else {
      probefinder[i] = -pow(2, power);
    }
    i_have_check[i] = -1;
    power--;
  }
  this->check_m = 0;
  int count_max = 0;
  this->previus_probe.clear();
  number_of_index_equal = ProbeSearch(count_max, number_of_index_equal, this->fqdecimal);
  count_max = 1;

  while (count_max <= this->probes &&  this->check_m < this->M){
    number_of_index_equal = FindProbe(1, count_max, number_of_index_equal, probefinder, i_have_check, this->fqdecimal);
    count_max++;
  }

  if(number_of_index_equal.size() == 0){
    cout << "No close point with q." << endl;
    return point_distance;
  }
  for(int i = 0; i < number_of_index_equal.size(); i++){ // Similar method with LSH
    POINTDISTANCE pd;
    pd.distance = EuclidianDistance(all_p[number_of_index_equal[i]], q, this->dimensions);
    pd.point = number_of_index_equal[i];
    point_distance.push_back(pd);
  }
  vector <POINTDISTANCE> results;

  sort(point_distance.begin(), point_distance.end(), PDCompareFunctionDistance()); // Sort based on distance

  int N2 = this->N;
  int M2 = this->M;
  for(int i = 0; ((i < N2) && (i < M2) && (i < point_distance.size()) && (results.size() < this->N)); i++){
    if ((self == 1) && (point_distance[i].distance == (double)0.0)){
      N2++;
      M2++;
    } else {
      results.push_back(point_distance[i]);
    }
  }

  return results;
}

vector<POINTDISTANCE> HYPERCUBE::RangeSearch(double *q, int R){
  HashTable(this->w, 1, this->k, q);
  vector <int> number_of_index_equal;
  vector <POINTDISTANCE> point_distance;
  if (R > 0) // To set different radius
    this->r = R;

  int probefinder[this->k];
  int i_have_check[this->k];
  int power = this->k-1;
  for(int i = 0; i < this->k; i++){
    if(fq[i] == 0){
      probefinder[i] = pow(2, power);
    }else {
      probefinder[i] = -pow(2, power);
    }
    i_have_check[i] = -1;
    power--;
  }

  this->check_m = 0;
  int count_max = 0;
  this->previus_probe.clear();
  number_of_index_equal = ProbeSearch(count_max, number_of_index_equal, this->fqdecimal);
  count_max = 1;
  while (count_max <= this->probes && this->check_m < this->M){
    number_of_index_equal = FindProbe(1, count_max, number_of_index_equal, probefinder, i_have_check, this->fqdecimal);
    count_max++;
  }

  if(number_of_index_equal.size() == 0){
    cout << "no close point with q" << endl;
    return point_distance;
  }

  for(int i = 0; i < number_of_index_equal.size(); i++){ // Find euclidian distance for points with equal index in f as query q
    POINTDISTANCE pd;
    pd.distance = EuclidianDistance(all_p[number_of_index_equal[i]], q, this->dimensions);
    if(pd.distance < (double)(this->r)){
      pd.point = number_of_index_equal[i];
      point_distance.push_back(pd);
    }
  }

  sort(point_distance.begin(), point_distance.end(), PDCompareFunctionDistance()); // Sort results and return them
  return point_distance;
}