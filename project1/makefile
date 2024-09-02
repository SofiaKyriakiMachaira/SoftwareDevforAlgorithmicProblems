CC  =  g++ -std=c++11

CFLAGS  = -g -c
CFLAGSD = -g -c -D
CHEAD = header.h common.h BruteForce.h
all: lsh cube cluster

lsh: LSH.o BruteForce.o mainlsh.o 
	$(CC) -o lsh LSH.o BruteForce.o main.o 

cube: Hypercube.o BruteForce.o maincube.o
	$(CC) -o cube Hypercube.o BruteForce.o main.o

cluster: LSH.o Hypercube.o CentralClustering.o maincluster.o
	$(CC) -o cluster CentralClustering.o LSH.o Hypercube.o main.o

LSH.o: LSH.cpp
	$(CC) $(CFLAGS) LSH.cpp common.h LSH.h

Hypercube.o: Hypercube.cpp
	$(CC) $(CFLAGS) Hypercube.cpp common.h Hypercube.h

CentralClustering.o: CentralClustering.cpp
	$(CC) $(CFLAGS) CentralClustering.cpp common.h CentralClustering.h LSH.h Hypercube.h

BruteForce.o: BruteForce.cpp
	$(CC) $(CFLAGS) BruteForce.cpp common.h BruteForce.h

mainlsh.o: main.cpp
	$(CC) $(CFLAGSD) LSHH main.cpp $(CHEAD) LSH.h

maincube.o: main.cpp 
	$(CC) $(CFLAGSD) CUBEH main.cpp $(CHEAD) Hypercube.h

maincluster.o: main.cpp
	$(CC) $(CFLAGSD) CLUSTERH main.cpp header.h common.h CentralClustering.h LSH.h Hypercube.h

.PHONY: clean
clean:
	rm -rf *o lsh cube cluster common.h.gch header.h.gch Hypercube.h.gch LSH.h.gch CentralClustering.h.gch BruteForce.h.gch