
#include "path-integral.h"


PathIntegral::PathIntegral(int N, double tau) :
	generator((int) &N), /* random seed */
	delta(0, 1),
	index(1, N),
	accept(0, 1),
	N(N),
	tau(tau),
	lattice(N),
	acceptance_ratio(0, N)
{}


void PathIntegral::populate_lattice(double start, double end)
{
	int i;
//	lattice[0] = start;
//	lattice[N - 1] = end;
	for (i = 0; i < N; i++)
	{
		lattice[i] = delta(generator);
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
	WArray<double> action(N);

	monitor.prime("generate samples", nruns);

	for (i = 0; i < N; i++)
	{
		action[i] = local_lagrangian(lattice, lagrangian, i, tau);
	}

	for (i = 0; i < nruns; i++)
	{
		monitor++;
		/* propose; conditionally accept; measure observable */
		int j = index(generator); 
		double new_pos = lattice[j] + delta(generator);

		double l1 = _local_lagrangian(lattice[j - 1], new_pos, lagrangian, tau);
		double l2 = _local_lagrangian(new_pos, lattice[j +1], lagrangian, tau);

		double delta_action = l1 + l2 - action[(j-1) % N] - action[j];

		/* average delta S should be 1 */

		/* metropolis rule */
		if (delta_action < 0
			|| accept(generator) < exp( - delta_action * tau ))
		{
			lattice[j] = new_pos;
			action[(j - 1) % N] = l1;
			action[j] = l2;
			++acceptance_ratio; 
		}

		/* sample, more or less, once per entire path change */
		static int k = 0;
		if (k++ % N)
		{
			observable(lattice, arg);
		}
	}
}

