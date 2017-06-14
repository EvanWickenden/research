#ifndef __PATH_INTEGRAL_H__
#define __PATH_INTEGRAL_H__

#include <iostream>
#include <math.h>
#include <random>

#include "warray.h"
#include "monitor.h"
#include "ratio.h"

typedef WArray<double> Lattice;

typedef double (*Lagrangian)(double q, double qdot); 
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

	Ratio acceptance_ratio;
	Monitor monitor;

	PathIntegral(int N, double tau); 

	void populate_lattice(double start, double end);

	void run(int nruns, Lagrangian lagrangian, Observable observable, void *arg);
};


#endif
