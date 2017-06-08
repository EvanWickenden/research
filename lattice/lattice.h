
/* a simple wrapper around an array of doubles  */
#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <iostream>

struct Lattice
{
	double *lattice; 
	int N;

	Lattice(int N) : lattice(new double[N]), N(N) {}
	
	~Lattice() { delete[] lattice; }

	/* automatically impose periodic boundary condition */
	double& operator[](int i) { return lattice[i % N]; }

	const double& operator[](int i) const { return lattice[i % N]; }
	
};

std::ostream& operator<<(std::ostream& s, const Lattice& lattice);

#endif
