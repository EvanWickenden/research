
CXX = clang++
CC = clang++

CXXFLAGS = -g -Wall -std=c++1z -Ofast
LDFLAGS = -g

.PHONY: default clean


default: simulation

simulation: data-set.o su2.o

simulation.o: path-integral.hpp lattice.hpp su2.h data-set.h callback/callback.h observables.hpp

data-set.o:

su2.o: matrix/matrix.hpp



clean: 
	rm -f *.o simulation
	rm -rf *.dSYM
