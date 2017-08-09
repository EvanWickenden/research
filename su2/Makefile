
CXX = clang++
CC = clang++

CXXFLAGS = -g -Wall -std=c++1z -Ofast
LDFLAGS = -g

.PHONY: default clean

default: simulation2

simulation2: su2.o data-set.o

simulation2.o: lattice/lattice2.hpp su2.h path-integral.h

path-integral2.o: lattice/lattice2.hpp path-integral.h su2.h



simulation: lattice.o su2.o path-integral.o data-set.o

simulation.o: lattice.h path-integral.h data-set.h

lattice.o: su2.h lattice.h

su2.o: su2.h 

path-integral.o: lattice.h su2.h path-integral.h

data-set.o: data-set.h


clean: 
	rm -f *.o simulation
	rm -rf *.dSYM
