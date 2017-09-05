
CXX = clang++
CC = clang++
CXXFLAGS = -Wall -Werror -std=c++11 -Ofast
LDFLAGS = -lm

.PHONY: default clean all

default: simulation


simulation: path-integral.o data.o

simulation.o: path-integral.h warray.h data.h observables.h

path-integral.o: path-integral.h

data.o: data.h




clean:
	rm -f *.o simulation

all: clean default
