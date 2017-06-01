
/* a simple wrapper around an array of doubles  */
#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <iostream>

struct Lattice
{
	double *lattice; /* public while in development */
	int N;

	Lattice(int N) : lattice(new double[N]), N(N) {}
	
	~Lattice() { delete[] lattice; }

	double& operator[](int i) { return lattice[i]; }

	const double& operator[](int i) const { return lattice[i]; }
	
};

std::ostream& operator<<(std::ostream& s, const Lattice& lattice);

#endif
