#include "header.h"
#include "common.h"
#include "GNNS.h"
#include "MRNG.h"
#include "BruteForce.h"
#include "CentralClustering.h"
#include "Hypercube.h"
#define QUERYNUM 10 // how many images to take from query file

using namespace std;

double ApproximationFactor(double *p, double *q, int dimensions){
  double dist = 0;
  for (int i = 1; i < dimensions; i++)
    dist += pow(p[i] - q[i], 2);
  return sqrt(dist);
}

int main(int argc, char *argv[]) {
  
  string inputfilename;
  string queryfilename;
  string outputfilename;
  string clusterfilename = (string)"cluster.conf";
  string neuralinputfilename;
  string neuralqueryfilename;
  string neuraloriginalinput; // Input file of latent_dim returned to original dimensions
  string neuraloriginalquery; // Query file of latent_dim returned to original dimensions
  string method = (string)"Classic"; // Set as default option
  int k = 20;
  int N = 1;
  int L = 5;
  int M = 10;
  int probes = 2;
  int complete = 0;
  int R = 10000;
  int K = 50; // for gnns/mrng
  int l = 20;
  int r = 1; // for gnns/mrng
  int E = 30;
  int m = 1; // default values for every algorithm 
  
  // READ EVERYTHING FROM COMMAND LINE
  for (int i = 1; i < argc - 1; i++){

    if (argv[i] == (string)"-dp"){
      neuraloriginalinput = argv[i + 1];

    }if (argv[i] == (string)"-dq"){
      neuraloriginalquery = argv[i + 1]; 

    } if (argv[i] == (string)"-od"){
      neuralinputfilename = argv[i + 1];
    
    } if (argv[i] == (string)"-oq"){
      neuralqueryfilename = argv[i + 1];

    } if (argv[i] == (string)"-d") {
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

    } else if (argv[i] == (string)"-c") {
      clusterfilename = argv[i + 1];

    } else if (argv[i] == (string)"-K") {
      K = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-L") {
      L = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-M") {
      M = atoi(argv[i + 1]);
    
    } else if (argv[i] == (string)"-probes") {
      probes = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-r") {
      r = atoi(argv[i + 1]);

    } else if (argv[i] == (string)"-complete") {
      complete = 1;

    } else if (argv[i] == (string)"-m") {
      method = argv[i + 1];
    }
  }

  if (inputfilename.empty()){ // If not given using command line, user must give paths
    cout << "Enter input file path:" << endl;
    cin >> inputfilename;
  }
    
  if (queryfilename.empty()){
    cout << "Enter query file path:" << endl;
    cin >> queryfilename;
  }

  FILE *inputfile = fopen(inputfilename.c_str(), "r"); // Open input file
  if (!inputfile) {
    cout << "Can't open input file." << endl;
    return -1;
  }

  FILE *queryfile = fopen(queryfilename.c_str(), "r"); // Open query file
  if (!queryfile) {
    cout << "Can't open query file." << endl;
    return -1;
  }

  if (outputfilename.empty())
    outputfilename = (string)"output.txt";
  FILE *outputfile = fopen(outputfilename.c_str(), "w"); // Open output file

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

  if (neuralinputfilename.empty()){ // If input file of Neural Network not given
    cout << "Enter input file path from Neural Net:" << endl;
    cin >> neuralinputfilename;
  }

  FILE *neuralinputfile = fopen(neuralinputfilename.c_str(), "r"); // Open input file of Neural Network
  if (!neuralinputfile) {
    cout << "Can't open input file from Neural Net." << endl;
    return -1;
  }
  
  // WITH ORIGINAL DIMENSIONS //
  if (neuraloriginalinput.empty()){ // If input file of Neural Network not given
    cout << "Enter input file path from Neural Net with original dimensions:" << endl;
    cin >> neuraloriginalinput;
  }
  // WITH ORIGINAL DIMENSIONS //
  FILE *neuraloriginalp = fopen(neuraloriginalinput.c_str(), "r"); // Open input file of Neural Network
  if (!neuraloriginalp) {
    cout << "Can't open input file from Neural Net with original dimensions." << endl;
    return -1;
  }

  int imagesN, rowsN, columnsN, infoN[4]; // Read NeuralNet files
  for (int j = 0; j < 4; j++){
    char x[1000];
    if (fscanf(neuralinputfile, " %1000s", x) == 1)
      infoN[j] = atof(x);
    char y[1000];
    int dump[4];
    if (fscanf(neuraloriginalp, " %1000s", y) == 1)
      dump[j] = atof(x);
  }

  imagesN = infoN[1]; // Number of images
  rowsN = infoN[2];
  columnsN = infoN[3];
  int latent_dim = rowsN*columnsN;
  double *vN[imagesN];
  double *vNOG[imagesN]; // WITH ORIGINAL DIMENSIONS

  for (int i = 0; i < imagesN; i++){
    vN[i] = (double*)malloc(latent_dim*sizeof(double));
    vNOG[i] = (double*)malloc(dimensions*sizeof(double));
  }

  for (int i = 0; i < imagesN; i++){ // Fill each table with rows*columns pixels
    count = 0;
    char x[1000];
    int flag = 0;
    do {
      if (fscanf(neuralinputfile, "%1000s", x) == 1){
        vN[i][count] = atof(x);
      } else {
        flag = 1;
      }
    } while ((++count < latent_dim) && (flag == 0));  // Ensure file hasn't ended

    count = 0;
    char y[1000];
    flag = 0;
    do {
      if (fscanf(neuraloriginalp, "%1000s", y) == 1){
        vNOG[i][count] = atof(y);
      } else {
        flag = 1;
      }
    } while ((++count < dimensions) && (flag == 0));  // Ensure file hasn't ended
  }

  if (neuralqueryfilename.empty()){ // If query file name of neural network not given
    cout << "Enter query file path from Neural Net:" << endl;
    cin >> neuralqueryfilename;
  }

  FILE *neuralqueryfile = fopen(neuralqueryfilename.c_str(), "r"); // Open query file of neural network
  if (!neuralqueryfile) {
    cout << "Can't open query file from Neural Net." << endl;
    return -1;
  }
  // WITH ORIGINAL DIMENSIONS
  if (neuraloriginalquery.empty()){ // If query file name of neural network not given
    cout << "Enter query file path from Neural Net with original dimensions:" << endl;
    cin >> neuraloriginalquery;
  }
  // WITH ORIGINAL DIMENSIONS
  FILE *neuraloriginalq = fopen(neuraloriginalquery.c_str(), "r"); // Open query file of neural network
  if (!neuraloriginalq) {
    cout << "Can't open query file from Neural Net with original dimensions." << endl;
    return -1;
  }

  double *qN[QUERYNUM]; // Create as an array for each query
  double *qNOG[QUERYNUM]; // WITH ORIGINAL DIMENSIONS
  for (int i = 0; i < QUERYNUM; i++){
    qN[i] = (double*)malloc(latent_dim*sizeof(double));
    qNOG[i] = (double*)malloc(dimensions*sizeof(double));
  }

  for (int i = 0; i < QUERYNUM; i++){ // Fill each array with rows*columns pixels
    count = 0;
    char x[1000];
    int flag = 0;
    do {
      if (fscanf(neuralqueryfile, "%1000s", x) == 1){
        qN[i][count] = atof(x);
      } else {
        flag = 1;
      }
    } while ((++count < latent_dim) && (flag == 0)); // Ensure file hasn't ended

    count = 0;
    flag = 0;
    do {
      if (fscanf(neuraloriginalq, "%1000s", x) == 1){
        qNOG[i][count] = atof(x);
      } else {
        flag = 1;
      }
    } while ((++count < dimensions) && (flag == 0)); // Ensure file hasn't ended
  }
  
  vector <POINTDISTANCE> brutenn[QUERYNUM]; // To store information
  vector <POINTDISTANCE> bruters[QUERYNUM];
  vector <POINTDISTANCE> brutennNN[QUERYNUM]; // NeuralNet results
  vector <POINTDISTANCE> brutersNN[QUERYNUM];

  chrono::steady_clock::time_point start_brute = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  BRUTEFORCE bruteforce(v, imagesN, N, dimensions, R); // Give NeuralNet number of images for fair comparison

  const auto start_brutenn = clock::now(); // Start clock for NearestNeighbors
  for (int i = 0; i < QUERYNUM; i++)
    bruteforce.NNearestNeighbors(q[i]).swap(brutenn[i]);
  const sec total_time_brute = clock::now() - start_brutenn; // Stop the clock

  const auto start_bruters = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    bruteforce.RangeSearch(q[i]).swap(bruters[i]);
  const sec total_time_bruters = clock::now() - start_bruters;

  chrono::steady_clock::time_point start_bruteNN = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  BRUTEFORCE bruteforceNN(vN, imagesN, N, latent_dim, R);
  
  const auto start_brutennNN = clock::now(); // Start clock for NearestNeighbors
  for (int i = 0; i < QUERYNUM; i++)
    bruteforceNN.NNearestNeighbors(qN[i]).swap(brutennNN[i]);
  const sec total_time_bruteNN = clock::now() - start_brutennNN; // Stop the clock

  const auto start_brutersNN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    bruteforceNN.RangeSearch(qN[i]).swap(brutersNN[i]);
  const sec total_time_brutersNN = clock::now() - start_brutersNN;

  vector <POINTDISTANCE> nn[QUERYNUM];
  vector <POINTDISTANCE> rs[QUERYNUM];
  vector <POINTDISTANCE> nnNN[QUERYNUM];
  vector <POINTDISTANCE> rsNN[QUERYNUM];
  // FOR QUESTION B
  LSH lsh(v, imagesN, N, dimensions, k, R, L); // Create LSH object

  chrono::steady_clock::time_point start_nn = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_NN = clock::now();

  for (int i = 0; i < QUERYNUM; i++)
    lsh.NNearestNeighbors(q[i], 0).swap(nn[i]);
    
  const sec total_time_NN = clock::now() - start_NN;

  const auto start_RS = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    lsh.RangeSearch(q[i], 0).swap(rs[i]);
  const sec total_time_RS = clock::now() - start_RS;
 
  LSH lshN(vN, imagesN, N, latent_dim, k, R, L); // Create LSH object

  chrono::steady_clock::time_point start_nnNN = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;

  const auto start_NNNN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    lshN.NNearestNeighbors(qN[i], 0).swap(nnNN[i]);
  const sec total_time_NNNN = clock::now() - start_NNNN;

  const auto start_RSNN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    lshN.RangeSearch(qN[i], 0).swap(rsNN[i]);
  const sec total_time_RSNN = clock::now() - start_RSNN;

  double MAF = 0;
  double MAFRS = 0;
  double MAFNN = 0;
  double MAFRSNN = 0;

  for (int i = 0; i < QUERYNUM; i++){ // Start printing output to output file

    for(int j = 0; j < N; j++){
      MAF += nn[i][j].distance / brutenn[i][j].distance;
      MAFRS += rs[i][j].distance / bruters[i][j].distance;
      MAFNN += ApproximationFactor(vNOG[nnNN[i][j].point], qN[i], dimensions) / brutenn[i][j].distance;
      MAFRSNN += ApproximationFactor(vNOG[rsNN[i][j].point], qN[i], dimensions) / bruters[i][j].distance;
    }
  }

  fprintf(outputfile, "tTrue: %f\n", total_time_brute.count() / QUERYNUM);
  fprintf(outputfile, "tTrue NeuralNet: %f\n", total_time_bruteNN.count() / QUERYNUM);

  fprintf(outputfile, "tTrue Range Search: %f\n", total_time_bruters.count() / QUERYNUM);
  fprintf(outputfile, "tTrue Range Search NeuralNet: %f\n", total_time_brutersNN.count() / QUERYNUM);

  fprintf(outputfile, "tLSH: %f\n", total_time_NN.count() / QUERYNUM); // Total time for ALL queries
  fprintf(outputfile, "MAF: %f\n", MAF / (QUERYNUM*N));
  fprintf(outputfile, "tLSH NeuralNet: %f\n", total_time_NNNN.count() / QUERYNUM);
  fprintf(outputfile, "MAF NeuralNet: %f\n", MAFNN / (QUERYNUM*N));

  fprintf(outputfile, "tLSH Range Search: %f\n", total_time_RS.count() / QUERYNUM);
  fprintf(outputfile, "MAF Range Search: %f\n", MAFRS / (QUERYNUM*N));
  fprintf(outputfile, "tLSH Range Search NeuralNet: %f\n", total_time_RSNN.count() / QUERYNUM);
  fprintf(outputfile, "MAF Range Search NeuralNet: %f\n", MAFRSNN / (QUERYNUM*N));
  
  // Reset variables after every algorithm 
  MAF = 0;
  MAFRS = 0;
  MAFNN = 0;
  MAFRSNN = 0;

  HYPERCUBE hypercube(v, imagesN, N, M, dimensions, k, R, probes); // Similar to LSH
  chrono::steady_clock::time_point start_hype_nn = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_hype_NN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    hypercube.NNearestNeighbors(q[i], 0).swap(nn[i]);
  const sec total_time_hype_NN = clock::now() - start_hype_NN;
      
  const auto start_hype_RS = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    hypercube.RangeSearch(q[i], 0).swap(rs[i]); 
  const sec total_time_hype_RS = clock::now() - start_hype_RS;

  HYPERCUBE hypeNN(vN, imagesN, N, M, latent_dim, k, R, probes); // Similar to LSH
  chrono::steady_clock::time_point start_hype_nnNN = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_hype_NNNN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    hypeNN.NNearestNeighbors(qN[i], 0).swap(nnNN[i]);
  const sec total_time_hype_NNNN = clock::now() - start_hype_NNNN;
      
  const auto start_hype_RSNN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    hypeNN.RangeSearch(qN[i], 0).swap(rsNN[i]); 
  const sec total_time_hype_RSNN = clock::now() - start_hype_RSNN;

  for (int i = 0; i < QUERYNUM; i++){ // Start printing output to output file

    for(int j = 0; j < N; j++){
      MAF += nn[i][j].distance / brutenn[i][j].distance;
      MAFRS += rs[i][j].distance / bruters[i][j].distance;
      MAFNN += ApproximationFactor(vNOG[nnNN[i][j].point], qN[i], dimensions) / brutenn[i][j].distance;
      MAFRSNN += ApproximationFactor(vNOG[rsNN[i][j].point], qN[i], dimensions) / bruters[i][j].distance;
    }
  }

  fprintf(outputfile, "tHypercube: %f\n", total_time_hype_NN.count() / QUERYNUM);
  fprintf(outputfile, "MAF: %f\n", MAF / (QUERYNUM*N));
  fprintf(outputfile, "tHypercube NeuralNet: %f\n", total_time_hype_NNNN.count() / QUERYNUM);
  fprintf(outputfile, "MAF NeuralNet: %f\n", MAFNN / (QUERYNUM*N));

  fprintf(outputfile, "tHypercube Range Search: %f\n", total_time_hype_RS.count() / QUERYNUM);
  fprintf(outputfile, "MAF Range Search: %f\n", MAFRS / (QUERYNUM*N));
  fprintf(outputfile, "tHypercube Range Search NeuralNet: %f\n", total_time_hype_NNNN.count() / QUERYNUM);
  fprintf(outputfile, "MAF Range Search NeuralNet: %f\n", MAFRSNN / (QUERYNUM*N));
  
  MAF = 0;
  MAFRS = 0;
  MAFNN = 0;
  MAFRSNN = 0;

  GNNS gnns(v, imagesN, dimensions, K, N, r, E); // Create GNNS object
  chrono::steady_clock::time_point start_gnns_nn = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_gnns_NN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    gnns.GraphNNSearch(q[i]).swap(nn[i]);
  const sec total_time_gnns_NN = clock::now() - start_gnns_NN;

  GNNS gnnsNN(vN, imagesN, latent_dim, K, N, r, E); // Create GNNS object
  chrono::steady_clock::time_point start_gnns_nnNN = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_gnns_NNNN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    gnnsNN.GraphNNSearch(qN[i]).swap(nnNN[i]);
  const sec total_time_gnns_NNNN = clock::now() - start_gnns_NNNN;
  for (int i = 0; i < QUERYNUM; i++){ // Start printing output to output file

    for(int j = 0; j < N; j++){
      MAF += nn[i][j].distance / brutenn[i][j].distance;
      MAFNN += ApproximationFactor(vNOG[nnNN[i][j].point], qN[i], dimensions) / brutenn[i][j].distance;
    }
  }

  fprintf(outputfile, "tGNNS: %f\n", total_time_gnns_NN.count() / QUERYNUM);
  fprintf(outputfile, "MAF: %f\n", MAF / (QUERYNUM*N));
  fprintf(outputfile, "tGNNS NeuralNet: %f\n", total_time_gnns_NNNN.count() / QUERYNUM);
  fprintf(outputfile, "MAF NeuralNet: %f\n", MAFNN / (QUERYNUM*N));

  MAF = 0;
  MAFNN = 0;
  
  MRNG mrng(v, imagesN, dimensions, K, N, r, E, l); // Create MRNG object
  chrono::steady_clock::time_point start_mrng_nn = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_mrng_NN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    mrng.SearchGraph(q[i]).swap(nn[i]);
  const sec total_time_mrng_NN = clock::now() - start_mrng_NN;

  MRNG mrngNN(vN, imagesN, latent_dim, K, N, r, E, l); // Create MRNG object
  chrono::steady_clock::time_point start_mrng_nnNN = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_mrng_NNNN = clock::now();
  for (int i = 0; i < QUERYNUM; i++)
    mrngNN.SearchGraph(qN[i]).swap(nnNN[i]);
  const sec total_time_mrng_NNNN = clock::now() - start_mrng_NNNN;

  for (int i = 0; i < QUERYNUM; i++){ // Start printing output to output file

    for(int j = 0; j < N; j++){
      MAF += nn[i][j].distance / brutenn[i][j].distance;
      MAFNN += ApproximationFactor(vNOG[nnNN[i][j].point], qN[i], dimensions) / brutenn[i][j].distance;
    }
  }
  
  fprintf(outputfile, "tMRNG: %f\n", total_time_mrng_NN.count() / QUERYNUM);
  fprintf(outputfile, "MAF: %f\n", MAF / (QUERYNUM*N));
  fprintf(outputfile, "tMNRG NeuralNet: %f\n", total_time_mrng_NNNN.count() / QUERYNUM);
  fprintf(outputfile, "MAF NeuralNet: %f\n", MAFNN / (QUERYNUM*N));
  
  // FOR QUESTION C
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
    if (current != (string)"//"){
      value[i] = current;
    } else {
      value[i] = (string)"0";
    }
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
 
  CENTRALCLUSTERING cluster(v, imagesN, dimensions, clusters);
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

  fprintf(outputfile, "Clustering_time: %f\n", total_time_CL.count() / QUERYNUM);
  fprintf(outputfile,"Silhouette: [");
  double stotal = cluster.Silhouette(centrals, outputfile);
  fprintf(outputfile, "%f]\n", stotal);
  clusterfile.close();
  delete[] centrals.teams_euclidian; 
 
  CENTRALCLUSTERING clusterNN(vN, imagesN, latent_dim, clusters);
  TEAMSCENTRAL centralsNN;
  chrono::steady_clock::time_point start_clNN = chrono::steady_clock::now();
  using clock = chrono::system_clock;
  using sec = chrono::duration<double>;
  const auto start_CLNN = clock::now();
  if (method == (string)"Classic"){
    centralsNN = clusterNN.Classic();
  }
  else if (method == (string)"LSH"){
    centralsNN = clusterNN.Reverse(1, vectorht, R, vectorhf, N, maxM);
  }
  else if (method == (string)"Hypercube"){
    centralsNN = clusterNN.Reverse(2, num_probes, R, hypercube_dimensions, N, maxM);
  }
  const sec total_time_CLNN = clock::now() - start_CLNN;

  fprintf(outputfile, "Clustering_time NeuralNet: %f\n", total_time_CLNN.count() / QUERYNUM);
  fprintf(outputfile,"Silhouette NeuralNet: [");
  double stotalNN = clusterNN.Silhouette(centralsNN, outputfile);
  fprintf(outputfile, "%f]\n", stotalNN);
  delete[] centralsNN.teams_euclidian; 

  for (int i = 0; i < images; i++) // Free and close files
    free(v[i]);

  for (int i = 0; i < imagesN; i++){
    free(vN[i]);
    free(vNOG[i]);
  }

  for (int i = 0; i < QUERYNUM; i++){
    free(q[i]);
    free(qN[i]);
    free(qNOG[i]);
  }

  fclose(queryfile);
  fclose(inputfile);
  fclose(neuralinputfile);
  fclose(neuralqueryfile);
  fclose(neuraloriginalp);
  fclose(neuraloriginalq);
  fclose(outputfile);
}