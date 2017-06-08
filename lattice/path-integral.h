#ifndef __PATH_INTEGRAL_H__
#define __PATH_INTEGRAL_H__

#include <iostream>
#include <math.h>
#include <random>
#include "lattice.h"
 
typedef double (*Lagrangian)(const double& q, const double& qdot); 
typedef void (*Observable)(const Lattice& lattice, void *arg);

class PathIntegral
{
	std::default_random_engine generator;
	std::normal_distribution<double> delta;
	std::uniform_int_distribution<int> index;
	std::uniform_real_distribution<double> accept;

	int N; /* reduntant; also appears in lattice struct */
	double tau;
	Lattice lattice; /* operator[] automatically imposes % N */

	public:

	long progress;

	PathIntegral(int N, double tau) :
		generator((int) &N),
		delta(0, 1),
		index(1, N - 1),
		accept(0, 1),
		N(N),
		tau(tau),
		lattice(N),
		progress(0)
	{ }

	void populate_lattice(double start, double end);

	void run(int nruns, Lagrangian lagrangian, Observable observable, void *arg);
};


#endif
