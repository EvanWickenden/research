

/* parameters: */

/* lattice dimensions */
/* time step */
/* number of iterations */

#include <math.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "path-integral.h"
#include "monitor.h"
#include "warray.h"
#include "process.h"


/* Euclidean lagrangian */

double omega;
double harmonic_oscillator(double q, double qdot)
{
	return (pow(qdot, 2) + omega * pow(qdot, 2));
}


WArray<double> *p2, *p3, *p4, *p5;

/* operator */
int nelts;
int k;
int start;
int offset; 
int cut_off;
void inline x1_x2(const Lattice& lattice, double *products)
{
	if (cut_off++ < 2000) return;
	nelts++;
	products[k] = lattice[start] * lattice[start + 1];
//	p2->data[k] = lattice[start] * lattice[start + 2];
//	p3->data[k] = lattice[start] * lattice[start + 3];
//	p4->data[k] = lattice[start] * lattice[start + 4];
//	p5->data[k] = lattice[start] * lattice[start];

	static int t = 0;
	if (t++ % lattice.N == 0)
	{
//		LOG("product %lf, difference %lf", products[k], lattice[start + 1] - lattice[start]);
//		LOG("%lf, %lf, %lf, %lf, %lf", products[k], p2->data[k],  p3->data[k],  p3->data[k],  p4->data[k]); 
	}

	k++;
}


#define less_than_point_2(_tau_exp, _nr_iterations) \
	({ 	float ratio = _tau_exp / (float) _nr_iterations; \
		(0 <= ratio && ratio < 0.2) ? (int) (ratio * _nr_iterations) : (int) (0.2 * _nr_iterations) ; })


FILE *out;
void simulate(int lattice_dimensions, double time_step, long nr_iterations, double _omega, int _start, int _offset,
double random_distribution_width, double initial_pos)
{
	omega = _omega;
	start = _start;
	offset = _offset;
	cut_off = 0;
	nelts = 0;
	k = 0;

//	int elts = nr_iterations / lattice_dimensions;
	int elts = nr_iterations;

	/* observable only called once per lattice_dimensions iterations */
	WArray<double> products(elts);
	WArray<double> products2(elts);
	WArray<double> products3(elts);
	WArray<double> products4(elts);
	WArray<double> products5(elts);

	p2 = &products2;	
	p3 = &products3;	
	p4 = &products4;	
	p5 = &products5;	


	/* minimize memory usage */
	{
		PathIntegral p(lattice_dimensions, time_step, random_distribution_width, rand());
		p.acceptance_ratio.denominator = nr_iterations;
		p.populate_lattice(initial_pos, initial_pos);
		p.monitor.start();
		p.run(nr_iterations, (Lagrangian) &harmonic_oscillator, (Observable) &x1_x2, products.data);
		p.monitor.end();

		LOG("acceptance ratio: %lf", p.acceptance_ratio());
		LOG("increases %d, decreases %d, ratio %lf; acceptance_ration * icreases", p.up_down.numerator, p.up_down.denominator, p.up_down());
//		fprintf(out, "acceptance ratio: %lf", p.acceptance_ratio());
	}

	products.N = nelts;
	products2.N = nelts;
	products3.N = nelts;
	products4.N = nelts;
	products5.N = nelts;

	double mean, tau_exp, tau_int;

	Process pr2(*p2);
	Process pr3(*p3);
	Process pr4(*p4);
	Process pr5(*p5);

	double m2, m3, m4, m5;
	
	m2 = pr2.mean();
	m3 = pr3.mean();
	m4 = pr4.mean();
	m5 = pr5.mean();

	Process process(products);
	mean = process.mean(); 

	LOG("mean given offset: 0 - %lf, 1 - %lf, 2 - %lf, 3 - %lf, 4 - %lf", m5, mean, m2, m3, m4);

	return;


	process.monitor.start();
	process.unnormalized_autocorrelation_function();
	tau_exp = process.exponential_autocorrelation_time();
	process.monitor.end();


	/* trim off 4 tau_exp, or 20% of sample, whichever is fewest */
//	int new_start = less_than_point_2(tau_exp, process.n);


//	process.data += new_start;
//	process.n -= new_start;

	/* worth recalculating autocorrelation function? */
	tau_int = process.integrated_autocorrelation_time(10, 4);
	/* recalculate mean */
//	mean = process.mean();

	float variance = (2 * tau_int * process.C[0]) / nelts;

	fprintf(out, "simulation parameters: lattice dimensions %d, time step %lf, samples generated %d, omega %lf, start_index %d, offset %d\n",
		lattice_dimensions,
		time_step,
		nelts,
		omega,
		start,
		offset);


	fprintf(out, "mean = %lf, tau_exp = %lf, tau_int = %lf, variance = %f\n", mean, tau_exp, tau_int, variance);
}


#define NELTS(_array)   ( sizeof (_array) / sizeof (*_array) )

//void simulate(int lattice_dimensions, double time_step, long nr_iterations, double _omega, int _start, int _offset,
//double random_distribution_width) 

int main()
{

	double omegas[] = { 100 };
//	double omegas[] = { 5 };
	int dimensions[] = { 120 };
//	int starts[] = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
	int starts[] = { 10 };
//	int offsets[] = { 1, 2, 3, 5, 7, 9 };
//	int offsets[] = { 1, 1, 1, 2, 2, 2 };
	int offsets[] = { 1, 5, 10 };
	long nruns = 1200000;
//	long nruns = 120;
	double time_step = 1 / (float) 10;

	
	srand(((int) &time_step) + time(NULL));

	out = stdout;

//	simulate(2, 0.1, 1000000, 1, 0, 1, 0.05);

	double ps = 0.0; /* ps; probability distribution width */
	for (; ps < 1000; ps += 10)
	{
		simulate(2, 0.01, 100000, 1, 0, 1, 0.1, ps);
	}


//	simulate(100, time_step, nruns, omegas[0], 40, 1);
//	simulate(100, time_step, nruns, omegas[0], 40, 1);
//	simulate(100, time_step, nruns, omegas[0], 40, 1);
//	simulate(100, time_step, nruns, omegas[0], 40, 1);

	/*
	int i, j, k, l;
	for (i = 0; i < NELTS(omegas); i++)
	{
		for (j = 0; j < NELTS(dimensions); j++)
		{
			for (k = 0; k < NELTS(starts); k++)
			{
				for (l = 1; l < 25; l+=2)
				{
					simulate(dimensions[j], time_step, nruns, omegas[i], starts[k], l);
				}
			}
			
		}
		LOG("__________________");
	}
	*/
}


