CC  =  g++ -std=c++11

CFLAGS  = -g -O2 -c
CFLAGSD = -g -O2 -c -D
CHEAD = header.h common.h BruteForce.h
all: graph_search

LSH.o: LSH.cpp
	$(CC) $(CFLAGS) LSH.cpp common.h LSH.h

Hypercube.o: Hypercube.cpp
	$(CC) $(CFLAGS) Hypercube.cpp common.h Hypercube.h

BruteForce.o: BruteForce.cpp
	$(CC) $(CFLAGS) BruteForce.cpp common.h BruteForce.h

graph_searchcube: Hypercube.o BruteForce.o GNNScube.o MRNG.o main2.o
	$(CC) -o graph_search Hypercube.o BruteForce.o GNNScube.o MRNG.o main2.o

graph_search: LSH.o BruteForce.o GNNS.o MRNG.o main2.o
	$(CC) -o graph_search LSH.o BruteForce.o GNNS.o MRNG.o main2.o

GNNS.o: GNNS.cpp
	$(CC) $(CFLAGS) GNNS.cpp common.h GNNS.h LSH.h

GNNScube.o: GNNS.cpp
	$(CC) $(CFLAGS) GNNS.cpp common.h GNNS.h Hypercube.h

MRNG.o: MRNG.cpp
	$(CC) $(CFLAGS) MRNG.cpp common.h MRNG.h

main2.o: main2.cpp
	$(CC) $(CFLAGS) main2.cpp header.h common.h BruteForce.h GNNS.h MRNG.h

.PHONY: clean
clean:
	rm -rf *o graph_search graph_searchcube common.h.gch header.h.gch Hypercube.h.gch LSH.h.gch BruteForce.h.gch GNNS.h.gch MRNG.h.gch