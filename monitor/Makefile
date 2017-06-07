
CC = clang++
CXX = clang++
CXXFLAGS = -Wall -std=c++11 -Wpedantic -Ofast
LDFLAGS = 

.PHONY: install install_h install_a

install: install_h install_a

install_h: monitor.h
	ln -f monitor.h ~/lib/monitor.h 

install_a: libmonitor.a
	ln -f libmonitor.a ~/lib/libmonitor.a 

libmonitor.a: monitor.o
	ar crf libmonitor.a monitor.o
	ranlib libmonitor.a

monitor.o: monitor.h



