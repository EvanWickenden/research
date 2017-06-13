

CXX = clang++
CC = clang++
CXXFLAGS = -Wpedantic -Wall -std=c++11 -Ofast
LDFLAGS = -lpthread 


.PHONY: default
default: simulation #analyze-samples


simulation: path-integral.o monitor.o process.o

simulation.o: path-integral.h monitor.h warray.h log.h process.h

path-integral.o: warray.h path-integral.h log.h

process.o: process.h

analyze-samples: process.o monitor.o

analyze-samples.o: process.h log.h 

monitor.o: monitor.h

.PHONY: clean all

clean:
	rm -f *.o simulation analyze-samples vgcore.*

all: clean default
