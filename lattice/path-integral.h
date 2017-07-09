#ifndef __PATH_INTEGRAL_H__
#define __PATH_INTEGRAL_H__

#include <iostream>
#include <math.h>
#include <random>
#include "lattice.h"
 
typedef double (*Lagrangian)(double q, double qdot); 
typedef void (*Observable)(const Lattice& lattice, void *arg);

class PathIntegral
{
	std::default_random_engine generator;
	std::normal_distribution<double> delta;
	std::uniform_int_distribution<int> index;
	std::uniform_real_distribution<double> accept;

	int N; /* reduntant; also appears in lattice struct */
	int progress;
	double tau;
	Lattice lattice;
	Lagrangian lagrangian;

	public:

	PathIntegral(int N, double tau, Lagrangian lagrangian) :
		generator((int) &N),
		delta(0, 1),
		index(1, N - 2),
		accept(0, 1),
		N(N),
		progress(0),
		tau(tau),
		lattice(N),
		lagrangian(lagrangian)
	{ }

	void populate_lattice(double start, double end);

	void run(int nruns, Observable observable, void *arg);

	int get_progress();
};


#endif
