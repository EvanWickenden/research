

CXX = clang++
CC = clang++
CXXFLAGS = -g -Wall -std=c++11

run: lattice.o path-integral.o

run.o: path-integral.h

lattice.o: lattice.h

path-integral.o: lattice.h path-integral.h

.PHONY: clean

clean:
	rm -f *.o path-integral 

