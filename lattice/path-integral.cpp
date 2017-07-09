
/* first draft: procedural; specific */
#include "path-integral.h"


int PathIntegral::get_progress()
{
	return progress;
}


void PathIntegral::populate_lattice(double start, double end)
{
	int i;
	lattice[0] = start;
	lattice[N - 1] = end;
	for (i = 1; i < N - 1; i++)
	{
		lattice[i] = delta(generator);
	}
}

/* note that result should be multiplied by tau (ommitted for computational efficiency) */
#define _local_lagrangian(_x1, _x2, _lagrangian, _tau) 	\
	(_lagrangian(_x1, (_x2 - _x1) / _tau))

#define local_lagrangian(_lattice, _lagrangian, _index, _tau)   \
	(_local_lagrangian(_lattice[_index], _lattice[_index + 1], _lagrangian, _tau))
	

void PathIntegral::run(int nruns, Observable observable, void *arg)
{
	int i;
	double *action = new double[N];

	for (i = 0; i < N - 1; i++)
	{
		action[i] = local_lagrangian(lattice, lagrangian, i, tau);
	}

	for (i = 0; i < nruns; i++)
	{
		/* propose; conditionally accept; measure observable */
		int j = index(generator); 
		double new_pos = lattice[j] + delta(generator);

		double l1 = _local_lagrangian(lattice[j - 1], new_pos, lagrangian, tau);
		double l2 = _local_lagrangian(new_pos, lattice[j +1], lagrangian, tau);

		double delta_action = l1 + l2 - action[j-1] - action[j];

		progress++;

		/* metropolis rule */
		if (delta_action < 0
			|| accept(generator) < exp( - delta_action * tau ))
		{
			lattice[j] = new_pos;
			action[j - 1] = l1;
			action[j] = l2;
		}

		static int k = 0;
		if (k++ % 1000 == 0)
		{
			observable(lattice, arg);
		}
	}

	delete[] action;
}

