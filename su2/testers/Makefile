
CXX = clang++

CXXFLAGS = -g -Wall -std=c++11
LDFLAGS = -g

.PHONY: default clean

testers = plaquette-tester su2-tester path-integral-tester random-matrices-tester

default: $(testers)

su2-tester: ../su2.o

plaquette-tester: ../su2.o ../lattice.o

path-integral-tester: ../path-integral.o ../su2.o ../lattice.o

random-matrices-tester: ../path-integral.o ../su2.o ../lattice.o

clean:
	rm -f $(testers) *.o
	rm -rf *.dSYM
