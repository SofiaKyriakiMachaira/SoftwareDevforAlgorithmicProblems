#include "header.h"
#include "common.h"
#include "GNNS.h"
#include "MRNG.h"
#include "BruteForce.h"

#define QUERYNUM 10 // how many images to take from query file

using namespace std;

int main(int argc, char *argv[]) {

  string inputfilename;
  string queryfilename;
  string outputfilename;
  int k = 50;
  int N = 1;
  int l = 20;
  int R = 1;
  int E = 30;
  int m = 1; // default values

  for (int i = 1; i < argc - 1; i++){ // Get all values from command line

    if (argv[i] == (string)"-d") {
      inputfilename = argv[i + 1];

    } else if (argv[i] == (string)"-q") {
      queryfilename = argv[i + 1];

    } else if (argv[i] == (string)"-k") {
      k = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-E") {
      E = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-R") {
      R = atoi(argv[i + 1]);
    
    } else if (argv[i] == (string)"-N") {
      N = atoi(argv[i + 1]);
    
    } else if (argv[i] == (string)"-l") {
      l = atoi(argv[i + 1]);
    
    } else if (argv[i] == (string)"-m") {
      m = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-o") {
      outputfilename = argv[i + 1];
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
    vector <POINTDISTANCE> brutenn[QUERYNUM]; // To store information

    chrono::steady_clock::time_point start_brute = chrono::steady_clock::now();
    using clock = chrono::system_clock;
    using sec = chrono::duration<double>;
    BRUTEFORCE bruteforce(v, images, N, dimensions, R);
    const auto start_brutenn = clock::now(); // Start clock for NearestNeighbors
    for (int i = 0; i < QUERYNUM; i++)
      bruteforce.NNearestNeighbors(q[i]).swap(brutenn[i]);

    const sec total_time_brute = clock::now() - start_brutenn; // Stop the clock

    if (m == 1){
      GNNS gnns(v, images, dimensions, k, N, R, E); // Create GNNS object
      chrono::steady_clock::time_point start_nn = chrono::steady_clock::now();
      using clock = chrono::system_clock;
      using sec = chrono::duration<double>;
      const auto start_NN = clock::now();
      for (int i = 0; i < QUERYNUM; i++)
        gnns.GraphNNSearch(q[i]).swap(nn[i]);
      const sec total_time_NN = clock::now() - start_NN;

      fprintf(outputfile, "GNNS Results\n");
      double max = 0.0;
      for (int i = 0; i < QUERYNUM; i++){ // Start printing output to output file

        fprintf(outputfile, "Query %d\n", i + 1);

        for(int j = 0; j < N; j++){
          double current = nn[i][j].distance / brutenn[i][j].distance;
          if (max < current)
            max = current;
          fprintf(outputfile, "Nearest neighbor-%d: %d\ndistanceApproximate: %f\n", j + 1, nn[i][j].point, nn[i][j].distance);
          fprintf(outputfile, "distanceTrue: %f\n", brutenn[i][j].distance);
        }
      }
      fprintf(outputfile, "tAverageApproximate: %f\n", total_time_NN.count()/QUERYNUM); // Total average time for queries
      fprintf(outputfile, "tAverageTrue: %f\n", total_time_brute.count()/QUERYNUM);
      fprintf(outputfile, "MAF: %f\n", max);

    } else if (m == 2) {
      MRNG searchongraph(v, images, dimensions, k, N, R, E, l); // Create GNNS object
      chrono::steady_clock::time_point start_nn = chrono::steady_clock::now();
      using clock = chrono::system_clock;
      using sec = chrono::duration<double>;
      const auto start_NN = clock::now();
      for (int i = 0; i < QUERYNUM; i++)
        searchongraph.SearchGraph(q[i]).swap(nn[i]);
      const sec total_time_NN = clock::now() - start_NN;

      fprintf(outputfile, "MRNG Results\n");
      double max = 0.0;
      for (int i = 0; i < QUERYNUM; i++){ // Start printing output to output file

        fprintf(outputfile, "Query %d\n", i + 1);

        for(int j = 0; j < N; j++){
          double current = nn[i][j].distance / brutenn[i][j].distance;
          if (max < current)
            max = current;

          fprintf(outputfile, "Nearest neighbor-%d: %d\ndistanceApproximate: %f\n", j + 1, nn[i][j].point, nn[i][j].distance);
          fprintf(outputfile, "distanceTrue: %f\n", brutenn[i][j].distance);
        }
      }
      fprintf(outputfile, "tAverageApproximate: %f\n", total_time_NN.count()/QUERYNUM); // Total average time for queries
      fprintf(outputfile, "tAverageTrue: %f\n", total_time_brute.count()/QUERYNUM);
      fprintf(outputfile, "MAF: %f\n", max);
    }

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