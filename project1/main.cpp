#include "header.h"
#include "common.h"
#ifdef LSHH
#include "LSH.h"
#include "BruteForce.h"
#endif
#ifdef CUBEH
#include "Hypercube.h"
#include "BruteForce.h"
#endif
#ifdef CLUSTERH
#include "CentralClustering.h"
#include "LSH.h"
#include "Hypercube.h"
#endif
#define QUERYNUM 10 // how many images to take from query file

using namespace std;

int main(int argc, char *argv[]) {

  string inputfilename;
  string queryfilename;
  string outputfilename;
  string clusterfilename;
  string method = (string)"Classic"; // Set as default option
  int k = 4;
  int N = 1;
  int L = 5;
  int M = 10;
  int probes = 2;
  int complete = 0;
  int R = 10000; // default values

  for (int i = 1; i < argc - 1; i++){ // Get all values from command line

    if ((argv[i] == (string)"-d") || (argv[i] == (string)"-i")) {
      inputfilename = argv[i + 1];
    
    } else if (argv[i] == (string)"-c") {
      clusterfilename = argv[i + 1];

    } else if (argv[i] == (string)"-q") {
      queryfilename = argv[i + 1];

    } else if (argv[i] == (string)"-k") {
      k = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-L") {
      L = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-M") {
      M = atoi(argv[i + 1]);
    
    } else if (argv[i] == (string)"-probes") {
      probes = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-o") {
      outputfilename = argv[i + 1];

    } else if (argv[i] == (string)"-N") {
      N = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-R") {
      R = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-complete") {
      complete = 1;

    } else if (argv[i] == (string)"-m") {
      method = argv[i + 1];
    }
  }

  int flag = 0; // For repetition
  do {
    if (inputfilename.empty()){ // If not given using command line, user must give paths
      cout << "Enter input file path:" << endl;
      cin >> inputfilename;
    }
    
    if (queryfilename.empty()){
      cout << "Enter query file path:" << endl;
      cin >> queryfilename;
    }

    if (outputfilename.empty()){
      cout << "Enter output file path:" << endl;
      cin >> outputfilename;
    }

    FILE *inputfile = fopen(inputfilename.c_str(), "r"); // Open input file
    if (!inputfile) {
      cout << "Can't open input file." << endl;
      return -1;
    }

    int images, rows, columns, info[4];
    for (int j = 0; j < 4; j++){

      double byte[4];
      for (int i = 0; i < 4 ; i++)
        byte[i] = fgetc(inputfile);
      info[j] = byte[3] + byte[2]*256 + byte[1]*256*256 + byte[0]*256*256*256; // Turn hex into decimal
    }

    images = info[1]; // Number of images
    rows = info[2];
    columns = info[3];
    int dimensions = rows*columns;
    double *v[images];
    for (int i = 0; i < images; i++)
      v[i] = (double*)malloc(dimensions*sizeof(double));
    int c, count;

    for (int i = 0; i < images; i++){ // Fill each table with rows*columns pixels
      count = 0;
      do {
        c = fgetc(inputfile);
        v[i][count] = c; // Insert pixel into table
      } while ((++count < rows*columns) && (c != EOF)); // Ensure file hasn't ended
    }

    FILE *queryfile = fopen(queryfilename.c_str(), "r"); // Open query file
    if (!queryfile) {
      cout << "Can't open query file." << endl;
      return -1;
    }

    double *q[QUERYNUM]; // Create as an array for each query
    for (int i = 0; i < QUERYNUM; i++)
      q[i] = (double*)malloc(dimensions*sizeof(double));

    for (int i = 0; i < QUERYNUM; i++){ // Fill each array with rows*columns pixels
      count = 0;
      do {
        c = fgetc(queryfile);
        q[i][count] = c; // Insert pixel into array
      } while ((++count < rows*columns) && (c != EOF)); // Ensure file hasn't ended
    }

    if (outputfilename.empty())
      outputfilename = (string)"output.txt";
    FILE *outputfile = fopen(outputfilename.c_str(), "w"); // Open output file

    vector <POINTDISTANCE> nn[QUERYNUM];
    vector <POINTDISTANCE> rs[QUERYNUM];
    #if defined(LSHH) || defined(CUBEH)
      vector <POINTDISTANCE> brutenn[QUERYNUM]; // To store information
      vector <POINTDISTANCE> bruters[QUERYNUM];

      chrono::steady_clock::time_point start_brute = chrono::steady_clock::now();
      using clock = chrono::system_clock;
      using sec = chrono::duration<double>;

      BRUTEFORCE bruteforce(v, images, N, dimensions, R);
      const auto start_brutenn = clock::now(); // Start clock for NearestNeighbors
      for (int i = 0; i < QUERYNUM; i++)
        bruteforce.NNearestNeighbors(q[i]).swap(brutenn[i]);

      const sec total_time_brute = clock::now() - start_brutenn; // Stop the clock
      for (int i = 0; i < QUERYNUM; i++)
        bruteforce.RangeSearch(q[i]).swap(bruters[i]);
    #endif

    #ifdef LSHH
      LSH lsh(v, images, N, dimensions, k, R, L); // Create LSH object

      chrono::steady_clock::time_point start_nn = chrono::steady_clock::now();
      using clock = chrono::system_clock;
      using sec = chrono::duration<double>;
      const auto start_NN = clock::now();

      for (int i = 0; i < QUERYNUM; i++)
        lsh.NNearestNeighbors(q[i]).swap(nn[i]);
    
      const sec total_time_NN = clock::now() - start_NN;

      for (int i = 0; i < QUERYNUM; i++)
        lsh.RangeSearch(q[i], 0).swap(rs[i]);

      for (int i = 0; i < QUERYNUM; i++){ // Start printing output to output file

        fprintf(outputfile, "Query %d\n", i + 1);

        for(int j = 0; j < N; j++){
          fprintf(outputfile, "Nearest neighbor-%d: %d\ndistanceLSH: %f\n", j + 1, nn[i][j].point, nn[i][j].distance);
          fprintf(outputfile, "distanceTrue: %f\n", brutenn[i][j].distance);
        }
        fprintf(outputfile, "R-near neighbors:\n");
        for(int k = 0; k < rs[i].size(); k++) 
          fprintf(outputfile, "%d, True: %d\n", rs[i][k].point, bruters[i][k].point);
      }
      fprintf(outputfile, "tLSH: %f\n", total_time_NN.count()); // Total time for ALL queries
      fprintf(outputfile, "tTrue: %f\n", total_time_brute.count());
    #endif

    #ifdef CUBEH
      HYPERCUBE hypercube(v, images, N, M, dimensions, k, R, probes); // Similar to LSH
      chrono::steady_clock::time_point start_nn = chrono::steady_clock::now();
      using clock = chrono::system_clock;
      using sec = chrono::duration<double>;
      const auto start_NN = clock::now();
      for (int i = 0; i < QUERYNUM; i++)
        hypercube.NNearestNeighbors(q[i]).swap(nn[i]);
      const sec total_time_NN = clock::now() - start_NN;
      
      for (int i = 0; i < QUERYNUM; i++)
        hypercube.RangeSearch(q[i], 0).swap(rs[i]); 

      for (int i = 0; i < QUERYNUM; i++){

        fprintf(outputfile, "Query %d\n", i + 1);

        for(int j = 0; j < N; j++){
          fprintf(outputfile, "Nearest neighbor-%d: %d\ndistanceHypercube: %f\n", j + 1, nn[i][j].point, nn[i][j].distance);
          fprintf(outputfile, "distanceTrue: %f\n", brutenn[i][j].distance);
        }
        fprintf(outputfile, "R-near neighbors:\n");
        for(int k = 0; k < rs[i].size(); k++) 
          fprintf(outputfile, "%d, True: %d\n", rs[i][k].point, bruters[i][k].point);
      }
      fprintf(outputfile, "tHypercube: %f\n", total_time_NN.count());
      fprintf(outputfile, "tTrue: %f\n", total_time_brute.count());
    #endif

    #ifdef CLUSTERH
      ifstream clusterfile;
      clusterfile.open(clusterfilename.c_str()); // To open cluster.conf file
      int clusters, vectorht, vectorhf, maxM, hypercube_dimensions, num_probes;
      string value[6];
      for (int i = 0; i < 6; i++){ // Get info from the file
        string line;
        getline(clusterfile, line);
        string current;
        stringstream getwords(line);
        getwords >> current;
        getwords >> current; 
        if (current != (string)"//")
          value[i] = current;
      }

      clusters = atoi(value[0].c_str()); // Set defaults if file had missing options
      if (clusters == 10)
        clusters = 10;
      vectorht = atoi(value[1].c_str());
      if (vectorht == 0)
        vectorht = 3;
      vectorhf = atoi(value[2].c_str());
      if (vectorhf == 0)
        vectorhf = 4;
      maxM = atoi(value[3].c_str());
      if (maxM == 0)
        maxM = 10;
      hypercube_dimensions = atoi(value[4].c_str());
      if (hypercube_dimensions == 0)
        hypercube_dimensions = 3;
      num_probes = atoi(value[5].c_str());
      if (num_probes = 0)
        num_probes = 2; 

      CENTRALCLUSTERING cluster(v, images, dimensions, clusters);
      TEAMSCENTRAL centrals;
      chrono::steady_clock::time_point start_cl = chrono::steady_clock::now();
      using clock = chrono::system_clock;
      using sec = chrono::duration<double>;
      const auto start_CL = clock::now();
      if (method == (string)"Classic"){
        centrals = cluster.Classic();
      }
      else if (method == (string)"LSH"){
        centrals = cluster.Reverse(1, vectorht, R, vectorhf, N, maxM);
      }
      else if (method == (string)"Hypercube"){
        centrals = cluster.Reverse(2, num_probes, R, hypercube_dimensions, N, maxM);
      }
      const sec total_time_CL = clock::now() - start_CL;

      fprintf(outputfile, "Algorithm: %s\n", method.c_str()); // Print to output file
      for (int i = 0; i < clusters; i++){
        fprintf(outputfile, "CLUSTER-%d {size: %d, centroid:", i + 1, (int)((centrals.leaders).size()));
        for (int j = 0; j < dimensions; j++)
          fprintf(outputfile, " %f", centrals.leaders[i][j]);
        if (complete == 1){
          fprintf(outputfile, ", ");
          for (int j = 0; j < (centrals.teams_euclidian[i].size() - 1); j++)
            fprintf(outputfile, "%d, ", centrals.teams_euclidian[i][j].point);
          fprintf(outputfile, "%d", centrals.teams_euclidian[i][centrals.teams_euclidian[i].size() - 1].point);
        }
        fprintf(outputfile, "}\n");
      }
      fprintf(outputfile, "clustering_time: %f\n", total_time_CL.count());
      fprintf(outputfile,"Silhouette: [");
      double stotal = cluster.Silhouette(centrals, outputfile);
      fprintf(outputfile, "%f]\n", stotal);
      clusterfile.close();
      delete[] centrals.teams_euclidian;
    #endif

    for (int i = 0; i < images; i++) // Free and close files
      free(v[i]);
    for (int i = 0; i < QUERYNUM; i++)
      free(q[i]);
    fclose(queryfile);
    fclose(inputfile);
    fclose(outputfile);

    cout << "Rerun program? If yes, type 1. If no, type 0." << endl; // Option to rerun program
    cin >> flag;
    if (flag == 1){
      cout << "Change input file? If yes, give new name. If no, type 0." << endl;
      string answer = (string)"0";
      cin >> answer;
      if (answer != (string)"0")
        inputfilename = answer;
      cout << "Change query file? If yes, give new name. If no, type 0." << endl;
      cin >> answer;
      if (answer != (string)"0")
        queryfilename = answer;
      cout << "Change output file? If yes, give new name. If no, type 0." << endl;
      cin >> answer;
      if (answer != (string)"0")
        outputfilename = answer;
    }
  } while (flag == 1);
  return 1;
}