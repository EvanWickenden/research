

CXX = clang++
CC = clang++
CXXFLAGS = -g -Wall -std=c++11 -Ofast
LDFLAGS = -lpthread

generate-samples: lattice.o path-integral.o

generate-samples.o: path-integral.h

lattice.o: lattice.h

path-integral.o: lattice.h path-integral.h

.PHONY: clean all

clean:
	rm -f *.o generate-samples vgcore.*

all: clean generate-samples
