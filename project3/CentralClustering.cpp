#include "CentralClustering.h"
#include "LSH.h"
#include "Hypercube.h"
#define LLOYDSLIMIT 20
using namespace std;

CENTRALCLUSTERING::CENTRALCLUSTERING(double **all_p1, int how_many1, int dimensions1, int k1){ // Initialization with KMeans++ method
  this->all_p = new double*[how_many1]; // Initialization of points dimension etc.
  for(int i = 0; i < how_many1; i++)
    this->all_p[i] = new double[dimensions1];

  for (int i = 0; i < how_many1; i++)
    memcpy(this->all_p[i], all_p1[i], dimensions1*sizeof(double));

  this->dimensions = dimensions1;
  this->k = k1;
  this->how_many = how_many1;
  this->centrals = Initialization(); // K-means++;
}

CENTRALCLUSTERING::~CENTRALCLUSTERING(){ // For leaks
  for(int i = 0; i < how_many; i++)
    delete[] all_p[i];

  delete[] all_p;
}

TEAMSCENTRAL CENTRALCLUSTERING::Initialization(){
  TEAMSCENTRAL centrals1;
  centrals1.teams_euclidian = new vector<POINTDISTANCE>[this->k]; // The euclidian Distance of all points to center i for k clusters
  int possibilities[this->how_many];
  int sum = FindSum(this->how_many);
  int random_factor = 0;

  for(int i = 0; i < this->how_many; i++){ // Create possibilities of i/sum that is the same of i^2/S(i^2) because S i=0->i=n(i/sum) =1
    random_factor = random_factor + i; // We don't use power to not exceed MAX_INT in case of large amount input points
    possibilities[i] = random_factor;
  }

  int leaders_central[k]; // The first centrals of k clusters (we save the points)
  leaders_central[0] = rand() % this->how_many; // Take a random point for first central
  for(int i = 1; i < k; i++)
    leaders_central[i] = -1;

  int lead = 0;
  while(lead < this->k){ // Create 2 to k centers

    vector<POINTDISTANCE> random_cent; // The euclidian average distance of 1 to lead centrals of all points
    for(int i = 0; i < this->how_many; i++){
      POINTDISTANCE pd;
      pd.distance = EuclidianDistance(this->all_p[leaders_central[lead]], this->all_p[i], this->dimensions);
      pd.point = i;
      centrals1.teams_euclidian[lead].push_back(pd);
      int new_how_many = this->how_many;
      if(CheckIfDifferent(leaders_central, i) == true){ // If i is not a leader in a cluster then find the average distance
        double distance_all = pd.distance;
        for(int j = 0; j < lead; j++)
          distance_all = distance_all + centrals1.teams_euclidian[j][i].distance;
        POINTDISTANCE pd_all;
        pd_all.distance = distance_all;
        pd_all.point = i;
        random_cent.push_back(pd_all);
      }
    }
    lead++;
    sort(random_cent.begin(), random_cent.end(), PDCompareFunctionDistance()); // Take from smaller distance to bigger distance to match the possibilities
    sum = FindSum((this->how_many - lead)); // n-t lead=t page 43/53
    random_factor = (rand() % sum); // Find the i/sum possibility
    for(int i = 0; i < this->how_many; i++){
      if(random_factor < possibilities[i]){
        leaders_central[lead] = random_cent[i-1].point;
        break;
      }
    }
  }

  for(int i = 0; i < this->k; i++){
    centrals1.leaders.push_back(all_p[leaders_central[i]]); // Add centers of cluster to return class object
    sort(centrals1.teams_euclidian[i].begin(), centrals1.teams_euclidian[i].end(), PDCompareFunctionPoint());
  }

  return centrals1;
}

int CENTRALCLUSTERING::FindSum(int n){ // Si=1->i=n i the sum of 1+2+3...+n
  int sum = 0;
  for (int i = 1; i <= n; i++)
    sum = sum + i;
  return sum;
}

bool CENTRALCLUSTERING::CheckIfDifferent(int *n, int check){ // Check if check doesn't exist in array n
  for (int i = 0; i <= this->k; i++){
    if(n[i] == check)
      return false;
  }
  return true;
}

double* CENTRALCLUSTERING::AveragePoint(vector<POINTDISTANCE> central){ // Find average point in vector
  double *all_points = new double[this->dimensions];
  for(int j = 0; j < this->dimensions; j++){
    all_points[j] = 0;
    for(int i = 0; i < central.size(); i++)
      all_points[j] += all_p[central[i].point][j];

    all_points[j] = (all_points[j]/central.size());
  }
  return all_points;
}

double CENTRALCLUSTERING::AverageDistance(vector<POINTDISTANCE> central, int itself){ // Find average distance in vector
  double all_points = 0;
  for(int i = 0; i < central.size(); i++){
    if (i == itself)
      continue;
    all_points += central[i].distance;
  }
  all_points = (all_points / central.size());
  return all_points;
}

double CENTRALCLUSTERING::AverageSilhouette(double* central){ // Find average silhouette
  double all_points = 0;
  for(int i = 0; i < this->how_many; i++)
    all_points += central[i];

  all_points = (all_points / this->how_many);
  return all_points;
}

double CENTRALCLUSTERING::EuclidianDistance(double *p, double *q, int dimensions){
  double dist = 0;
  for (int i = 1; i < dimensions; i++)
    dist += pow(p[i] - q[i], 2);
  return sqrt(dist);
}

TEAMSCENTRAL CENTRALCLUSTERING::Classic(){ // Lloyds method

  int check_flag = 0;
  int check_time = 0;
  vector <POINTDISTANCE> team_euclidian_new[this->k];
  int old_leaders_size[this->k];
  do {
    for (int i = 0; i< this->k; i++)
      team_euclidian_new[i].clear();

    for(int i = 0; i < centrals.teams_euclidian[0].size(); i++){ // We use team_euclidian_new so centrals.teams_euclidian sizes remain the same

      double min_distance = centrals.teams_euclidian[0][i].distance; // Find minimum distance
      int min_point = 0;
      for(int j = 1; j < this->k; j++){

        if(centrals.teams_euclidian[j][i].distance < min_distance){
          min_point = j;
          min_distance = centrals.teams_euclidian[j][i].distance;
        }
      }
      team_euclidian_new[min_point].push_back(centrals.teams_euclidian[min_point][i]); // Get copy of items
    }
    if(check_time == 0){
        for(int i = 0; i < this->k; i++)
          old_leaders_size[i] = team_euclidian_new[i].size() + team_euclidian_new[i].size(); // Create first old to be sure bigger the 10%
    }
    int sum = 0;
    for(int j = 0; j < this->k; j++){
      double *new_leaders = AveragePoint(team_euclidian_new[j]);
      int a = old_leaders_size[j] - team_euclidian_new[j].size();
      sum = sum + team_euclidian_new[j].size();
      if(old_leaders_size[j]*0.1 < abs(a)){ // If there is more than 10% change
        check_flag = 0;
        memcpy(centrals.leaders[j], new_leaders, this->dimensions*sizeof(double));
        centrals.teams_euclidian[j].clear();
        for(int i = 0; i < this->how_many; i++){
          POINTDISTANCE pd;
          pd.distance = EuclidianDistance(centrals.leaders[j], this->all_p[i], this->dimensions);
          pd.point = i;
          centrals.teams_euclidian[j].push_back(pd);
        }
      } else {
        check_flag++;
      }
      old_leaders_size[j] = team_euclidian_new[j].size();
      delete[] new_leaders;
    }
    check_time++;
  } while(check_time < LLOYDSLIMIT && check_flag < this->k);

  for (int i= 0; i < this->k; i++)
    centrals.teams_euclidian[i] = team_euclidian_new[i]; // Get final results to return them

  return centrals;
}

TEAMSCENTRAL CENTRALCLUSTERING::Reverse(int method, int L_Probes, int R, int K, int N, int M){ // method is 1 for LSH or 2 for Hypercube

  vector<POINTDISTANCE> flags_to_change[this->k];
  double dist;
  double dist2 = DBL_MAX;
  for (int i = 0; i < this->k; i++){
    for (int j = i + 1; j < this->k; j++){
      dist = EuclidianDistance(centrals.leaders[i], centrals.leaders[j], this->dimensions);
      if (dist < dist2)
        dist2 = dist;
    }
  }
  int Range = (int)(dist2 / 2);
  int total = 0;
  int result_size = Range*16; // Do range search 4 times
  double ** new_p;

  if (this->how_many <= 0 ){
    return centrals;
  }

  new_p = new double*[this->how_many]; // Make copy of points
  for(int i = 0; i < this->how_many; i++)
     new_p[i] = new double[this->dimensions];

  for (int i = 0; i < this->how_many; i++)
    memcpy(new_p[i], this->all_p[i], this->dimensions*sizeof(double));

  int new_how_many = this->how_many;
  for (int i = 0; i < this->k; i++)
      centrals.teams_euclidian[i].clear(); // Ensure it is empty

  while((total < this->how_many) && (Range < result_size)){
    int size = 0;
    total = 0;
    result_size = 0;

    if ((this->k <= 0 ) || (this->how_many <= 0)){
      return centrals;
    }
    vector<POINTDISTANCE> results[this->k];
    if (method == 1){
      LSH lsh(new_p, new_how_many, N, this->dimensions, K, R, L_Probes);
      for (int i = 0; i < this->k; i++){ // Do range search
        lsh.RangeSearch(centrals.leaders[i], Range).swap(results[i]);
        sort(results[i].begin(), results[i].end(), PDCompareFunctionPoint());
      }

      for (int i = 0; i < this->k; i++){
        for (int p = 0; p < results[i].size(); p++){
          int min = i;
          for (int j = i + 1; j < this->k; j++){
            if (count_if(results[j].begin(), results[j].end(), PDCompareFunctionEqualPoint(results[i][p].point)) != 0){ // If point found in multiple centrals
              dist = EuclidianDistance(all_p[results[i][p].point], centrals.leaders[j], this->dimensions);
              dist2 = EuclidianDistance(all_p[results[i][p].point], centrals.leaders[i], this->dimensions);
              if (dist < dist2) // Enter point to the central it is closest to
                min = j;
            }
          }
          flags_to_change[min].push_back(results[i][p]);
          centrals.teams_euclidian[min].push_back(results[i][p]);
        }
      }
      for (int i = 0; i < this->k; i++){
        size += flags_to_change[i].size();
        total += flags_to_change[i].size();
        result_size += results[i].size();
        flags_to_change[i].clear(); // Update counters and clear array
      }
    } else if (method == 2){ // Similarly here
      HYPERCUBE hypercube(new_p, new_how_many, N, M, this->dimensions, K, R, L_Probes);
      for (int i = 0; i < this->k; i++){
        hypercube.RangeSearch(centrals.leaders[i], Range).swap(results[i]);
        sort(results[i].begin(), results[i].end(), PDCompareFunctionPoint());
      }

      for (int i = 0; i < this->k; i++){
        for (int p = 0; p < results[i].size(); p++){
          int min = i;
          for (int j = i + 1; j < this->k; j++){
            if (count_if(results[j].begin(), results[j].end(), PDCompareFunctionEqualPoint(results[i][p].point)) != 0){
              dist = EuclidianDistance(all_p[results[i][p].point], centrals.leaders[j], this->dimensions);
              dist2 = EuclidianDistance(all_p[results[i][p].point], centrals.leaders[i], this->dimensions);
              if (dist < dist2)
                min = j;
            }
          }
          centrals.teams_euclidian[min].push_back(results[i][p]);
          flags_to_change[min].push_back(results[i][p]);
        }
      }

      for (int i = 0; i < this->k; i++){
        size += flags_to_change[i].size();
        total += flags_to_change[i].size();
        result_size += results[i].size();
        flags_to_change[i].clear();
      }
    }

    for(int i = 0; i < new_how_many; i++)
        delete[] new_p[i];

    delete[] new_p; // Delete points created
    new_how_many = new_how_many - size;
    if (new_how_many <= 0){
      return centrals;
    }
    new_p = new double*[new_how_many];
    for (int i = 0; i < new_how_many; i++)
      new_p[i] = new double[this->dimensions];
    int r = 0;
    for(int i = 0; i < this->how_many; i++){
      int check = 0;
      for (int j = 0; j < this->k; j++){
        if (count_if(results[j].begin(), results[j].end(), PDCompareFunctionEqualPoint(i)) != 0) // If i was found results[j]
          check++;

        if((check == this->k) && (r < (sizeof(new_p) /sizeof(double **))))
          memcpy(new_p[r++], this->all_p[i], this->dimensions*sizeof(double)); // Store point
      }
    }
    Range = Range*2;
  }

  for(int i = 0; i < new_how_many; i++)
    delete[] new_p[i];

  delete[] new_p;

  if (total < this->how_many){ // If not all points have been stored in a cluster, store them
    for(int i = 0; i < this->how_many; i++){
      int flag = -1;
      for (int j = 0; j < this->k; j++){
        if (count_if(centrals.teams_euclidian[j].begin(), centrals.teams_euclidian[j].end(), PDCompareFunctionEqualPoint(i)) != 0){
          flag = j;
          break;
        }
      }
      if (flag == -1){
        dist = EuclidianDistance(centrals.leaders[0], all_p[i], this->dimensions);
        int min = 0;
        for (int x = 1; x < this->k; x++){
          dist2 = EuclidianDistance(all_p[i], centrals.leaders[x], this->dimensions);
          if (dist > dist2){
            dist = dist2;
            min = x;
          }
        }
        POINTDISTANCE pd;
        pd.point = i;
        pd.distance = dist;
        centrals.teams_euclidian[min].push_back(pd);
      }
    }
  }
  return centrals;
}

double CENTRALCLUSTERING::Silhouette(TEAMSCENTRAL cluster, FILE * outputfile){ // To find silhouettes and average
  double results;
  double s[this->how_many]; // Find all silhouettes
  int counter = 0;
  for (int j = 0; j < this->k; j++){

    if (cluster.teams_euclidian->size() > j){
    for (int i = 0; i < cluster.teams_euclidian[j].size(); i++){
      double a = AverageDistance(cluster.teams_euclidian[j], i);
      double dist = (double)INT_MAX;
      int minr = 0;
      for (int r = 0; r < this->k; r++){
        if (j == r)
          continue;
        if ((cluster.leaders.size() > r) && (cluster.teams_euclidian[j].size() > i)){
          if (cluster.teams_euclidian[j][i].point < (sizeof(all_p)/ sizeof(double*))){
          double curr = EuclidianDistance(all_p[cluster.teams_euclidian[j][i].point], cluster.leaders[r], this->dimensions);
          if (curr < dist){
            dist = curr;
            minr = r;
          }
        }}
      }
      double b = AverageDistance(cluster.teams_euclidian[minr], -1);
      if (max(a, b) != 0)
        s[counter++] = (b - a)/max(a, b); // Based on slides
      if (counter > 0)
        fprintf(outputfile, "%f, ", s[counter - 1]);
    }}
  }
  results = AverageSilhouette(s); // Return average of all silhouettes (stotal)
  return results;
}