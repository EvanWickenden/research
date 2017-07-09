

CXX = clang++
CC = clang++
CXXFLAGS = -g -Wall -std=c++11 #-Ofast
LDFLAGS = -g -lpthread 


.PHONY: default
default: simulation 


simulation: path-integral.o monitor.o process.o ratio.o

simulation.o: path-integral.h monitor.h warray.h log.h process.h

path-integral.o: warray.h path-integral.h log.h ratio.h

process.o: process.h

analyze-samples: process.o monitor.o

analyze-samples.o: process.h log.h 

monitor.o: monitor.h

ratio.o: ratio.h

.PHONY: clean all

clean:
	rm -f *.o simulation analyze-samples vgcore.*

all: clean default
