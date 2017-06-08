

CXX = clang++
CC = clang++
CXXFLAGS = -I/home/evan/lib/ -g -Wall -std=c++11 -Ofast
LDFLAGS = -L/home/evan/lib/ -lpthread -lmonitor


.PHONY: default
default: generate-samples analyze-samples

generate-samples: lattice.o path-integral.o ../monitor/monitor.o

generate-samples.o: path-integral.h settings.h log.h

lattice.o: lattice.h 

path-integral.o: lattice.h path-integral.h settings.h log.h

process.o: process.h

analyze-samples: process.o ../monitor/monitor.o

analyze-samples.o: process.h settings.h log.h 

.PHONY: clean all

clean:
	rm -f *.o generate-samples analyze-samples vgcore.*

all: clean generate-samples
