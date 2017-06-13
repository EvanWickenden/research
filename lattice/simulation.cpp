

/* parameters: */

/* lattice dimensions */
/* time step */
/* number of iterations */

#include <math.h>
#include <iostream>
#include <cstdio>

#include "path-integral.h"
#include "monitor.h"
#include "warray.h"
#include "process.h"


/* Euclidean lagrangian */

double _omega;
double harmonic_oscillator(double q, double qdot)
{
	return (pow(qdot, 2) + _omega * pow(qdot, 2));
}

/* operator */
int nelts;
int k;
void inline x1_x2(const Lattice& lattice, double *products)
{
	static int i = lattice.N / 3;
	static int j = 2 * i;
	nelts++;
	products[k++] = lattice[i] * lattice[j];
}


FILE *out;
void simulate(int lattice_dimensions, double time_step, long nr_iterations, double omega)
{
	_omega = omega;
	nelts = 0;
	k = 0;

	/* observable only called once per lattice_dimensions iterations */
	WArray<double> products(nr_iterations);

	/* minimize memory usage */
	{
		PathIntegral p(lattice_dimensions, time_step);
		p.populate_lattice(0, 0);
		p.monitor.start();
		p.run(nr_iterations, (Lagrangian) &harmonic_oscillator, (Observable) &x1_x2, products.data);
		p.monitor.end();
	}

	products.N = nelts;

	double mean, tau_exp, tau_int;

	Process process(products);
	process.monitor.start();
	mean = process.mean();
	process.unnormalized_autocorrelation_function();
	tau_exp = process.exponential_autocorrelation_time();
	tau_int = process.integrated_autocorrelation_time(10, 4);
	process.monitor.end();

	fprintf(out, "simulation parameters: lattice dimensions %d, time step %lf, MC steps %ld, omega %lf\n",
		lattice_dimensions,
		time_step,
		nr_iterations,
		omega);

	fprintf(out, "mean = %lf, tau_exp = %lf, tau_int = %lf\n", mean, tau_exp, tau_int);
}


#define NELTS(_array)   ( sizeof (_array) / sizeof (*_array) )

int main()
{
	double omegas[] = { 0, 1, 2, 4, 6, 8, 10 };
	int dimensions[] = { 3, 6, 9, 12, 15, 20 };
	long nruns = 500000;
	double time_step = 0.1;

	out = stdout;

	int i, j;
	for (i = 0; i < NELTS(omegas); i++)
	{
		for (j = 0; j < NELTS(dimensions); j++)
		{
			simulate(dimensions[j], time_step, nruns, omegas[i]);
		}
	}
}


