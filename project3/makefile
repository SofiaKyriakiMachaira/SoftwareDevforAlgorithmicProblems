CC  =  g++ -std=c++11

CFLAGS  = -g -O2 -c
CFLAGSD = -g -O2 -c -D
CHEAD = header.h common.h BruteForce.h
all: reduced

LSH.o: LSH.cpp
	$(CC) $(CFLAGS) LSH.cpp common.h LSH.h

Hypercube.o: Hypercube.cpp
	$(CC) $(CFLAGS) Hypercube.cpp common.h Hypercube.h

CentralClustering.o: CentralClustering.cpp
	$(CC) $(CFLAGS) CentralClustering.cpp common.h CentralClustering.h LSH.h Hypercube.h

BruteForce.o: BruteForce.cpp
	$(CC) $(CFLAGS) BruteForce.cpp common.h BruteForce.h

reduced: Hypercube.o LSH.o BruteForce.o GNNS.o MRNG.o CentralClustering.o main3.o
	$(CC) -o reduced Hypercube.o LSH.o BruteForce.o GNNS.o MRNG.o CentralClustering.o main3.o

GNNS.o: GNNS.cpp
	$(CC) $(CFLAGS) GNNS.cpp common.h GNNS.h Hypercube.h LSH.h

MRNG.o: MRNG.cpp
	$(CC) $(CFLAGS) MRNG.cpp common.h MRNG.h

main3.o: main3.cpp
	$(CC) $(CFLAGS) main3.cpp header.h common.h BruteForce.h GNNS.h MRNG.h Hypercube.h CentralClustering.h

.PHONY: clean
clean:
	rm -rf *o reduced common.h.gch header.h.gch Hypercube.h.gch LSH.h.gch CentralClustering.h.gch BruteForce.h.gch GNNS.h.gch MRNG.h.gch