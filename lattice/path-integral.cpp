
#include "path-integral.h"

#include <time.h>

PathIntegral::PathIntegral(int N, double tau, double random_dist_width, int seed) :
	generator(seed), /* random seed */
	delta(0, random_dist_width),
	index(0, N - 1),
	accept(0, 1),
	N(N),
	tau(tau),
	lattice(N),
	acceptance_ratio(0, N),
	up_down(0, 0)
{}


void PathIntegral::populate_lattice(double start, double end)
{
	int i;
	lattice[0] = start;
	lattice[N - 1] = end;
	for (i = 1; i < N-1; i++)
	{
//		lattice[i] = delta(generator);
		//cold start
		lattice[i] = 0;

		//hot start
		lattice[i] = 100;
	}
}

/* note that result should be multiplied by tau (ommitted for computational efficiency) */
#define _local_lagrangian(_x1, _x2, _lagrangian, _tau) 	\
	(_lagrangian(_x1, (_x2 - _x1) / _tau))

#define local_lagrangian(_lattice, _lagrangian, _index, _tau)   \
	(_local_lagrangian(_lattice[_index], _lattice[_index + 1], _lagrangian, _tau))
	

void PathIntegral::run(int nruns, Lagrangian lagrangian, Observable observable, void *arg)
{
	int i;
	int k = 0;
//	double *action = new double[N];
	WArray<double> action(N);

	double sum_delta_s = 0.0;

	monitor.prime("generate samples", nruns);

	for (i = 0; i < N; i++)
	{
		action[i] = local_lagrangian(lattice, lagrangian, i, tau);
	}

	Ratio delta_direction(0, 0);

	for (i = 0; i < nruns; i++)
	{
		monitor++;
		/* propose; conditionally accept; measure observable */
		int j = index(generator); 

		double temp = delta(generator);

//		LOG("proposed change %lf", temp);

		double new_pos = lattice[j] + temp;

		double l1 = _local_lagrangian(lattice[j - 1], new_pos, lagrangian, tau);
		double l2 = _local_lagrangian(new_pos, lattice[j +1], lagrangian, tau);

//		LOG("i = %d; -1 = %lf, 0 = %lf, 1 = %lf; delta = %lf", j, lattice[j - 1], lattice[j], lattice[j + 1], temp);

		double delta_action = l1 + l2 - action[j - 1] - action[j];

		if (temp < 0) delta_direction++;
		else ++delta_direction;

		/* up if action is increased; down if decreased */
		if (delta_action < 0) up_down++; 
		else ++up_down;
	
//		LOG("change in action %lf", delta_action);

		double e = exp( - delta_action * tau);

		sum_delta_s += e;

//		LOG("average e^-S %lf", sum_delta_s / (float) i);

		/* metropolis rule */
		if (delta_action < 0
			|| accept(generator) < e)
		{
			lattice[j] = new_pos;
			action[j - 1] = l1;
			action[j] = l2;
			++acceptance_ratio; 
		}

		/* sample, more or less, once per entire path change */
//		k = (k+1) % N;
//		if (k == 0)
		{
			observable(lattice, arg);
		}
	}

	sum_delta_s /= nruns;
	LOG("ratio direction of change %lf", delta_direction());
//	LOG("average change in S: %lf", sum_delta_s);
}




