
CXX = clang++

CXXFLAGS = -Wall -std=c++1z
LDFLAGS = 

.PHONY: default clean

testers = su2-tester lattice-tester path-integral-tester

default: $(testers)

su2-tester: ../su2.o

lattice-tester: ../su2.o

matrix-tester:

path-integral-tester: ../su2.o

#plaquette-tester: ../su2.o ../lattice.o

#path-integral-tester: ../path-integral.o ../su2.o ../lattice.o

#random-matrices-tester: ../path-integral.o ../su2.o ../lattice.o

clean:
	rm -f $(testers) *.o
	rm -rf *.dSYM
