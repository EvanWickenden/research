

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

/* lattice element parameters */

int _start_index = 0;
int difference = 10;

void inline x1_x2(const Lattice& lattice, double *products)
{
	nelts++;
	products[k++] = lattice[_start_index] * lattice[_start_index + difference];
}


#define less_than_point_2(_tau_exp, _nr_iterations) \
	({ 	float ratio = _tau_exp / (float) _nr_iterations; \
		(0 <= ratio && ratio < 0.2) ? (int) (ratio * _nr_iterations) : (int) (0.2 * _nr_iterations) ; })


FILE *out;
void simulate(int lattice_dimensions, double time_step, long nr_iterations, double omega, int start_index)
{
	_omega = omega;
	_start_index = start_index;
	nelts = 0;
	k = 0;

	/* observable only called once per lattice_dimensions iterations */
	WArray<double> products(nr_iterations);

	/* minimize memory usage */
	{
		PathIntegral p(lattice_dimensions, time_step);
		p.acceptance_ratio.denominator = nr_iterations;
		p.populate_lattice(0, 0);
		p.monitor.start();
		p.run(nr_iterations, (Lagrangian) &harmonic_oscillator, (Observable) &x1_x2, products.data);
		p.monitor.end();

		LOG("acceptance ratio: %lf", p.acceptance_ratio());
//		fprintf(out, "acceptance ratio: %lf", p.acceptance_ratio());
	}

	products.N = nelts;

	double mean, tau_exp, tau_int;

	Process process(products);
	process.monitor.start();
	process.mean(); 
	process.unnormalized_autocorrelation_function();
	process.monitor.end();

///* DELETE ME */
//
//	LOG(out, "simulation parameters: lattice dimensions %d, time step %lf, MC steps %ld, omega %lf\n",
//		lattice_dimensions,
//		time_step,
//		nr_iterations,
//		omega);
//
//	LOG(out, "mean = %lf; start_index %d", mean, _start_index);
//
//	return;
//
//
///* __END__ */

	tau_exp = 4 * process.exponential_autocorrelation_time();

	/* trim off 4 tau_exp, or 20% of sample, whichever is fewest */
	int new_start = less_than_point_2(tau_exp, process.n);

	process.data += new_start;
	process.n -= new_start;

	/* worth recalculating autocorrelation function? */
	tau_int = process.integrated_autocorrelation_time(10, 4);
	/* recalculate mean */
	mean = process.mean();


	fprintf(out, "simulation parameters: lattice dimensions %d, time step %lf, MC steps %ld, omega %lf\n",
		lattice_dimensions,
		time_step,
		nr_iterations,
		omega);

	fprintf(out, "mean = %lf, tau_exp = %lf, nr truncated = %d, tau_int = %lf\n", mean, tau_exp, new_start, tau_int);
}


#define NELTS(_array)   ( sizeof (_array) / sizeof (*_array) )

int main()
{
	long nruns = 100000;
	double time_step = 20;
	out = stderr
	int i, j;

	/* element offset tester */

	int start_index[] = { 0, 5, 10, 15, 40, 60, 90 };

	for (i = 0; i < NELTS(start_index); i++)
	{
		simulate(100, time_step, nruns, 5, start_index[i]);
	}

	return 0;

	/* old */


	double omegas[] = { 0, 2, 5, 10, 15 };
	int dimensions[] = { 3, 6, 9, 15 };



	for (i = 0; i < NELTS(omegas); i++)
	{
		for (j = 0; j < NELTS(dimensions); j++)
		{
//			simulate(dimensions[j], time_step, nruns, omegas[i]);
		}
	}
}


