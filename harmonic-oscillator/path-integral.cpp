
#include <math.h>
#include <random>
#include <iostream>

#include "path-integral.h"


/* accepts a pre-populated lattice */
void path_integral(Lattice& lattice, long steps, double delta_x_width, DeltaAction& DS, Observable& observable)
{
	int N = lattice.N;

	std::random_device ran;
	std::mt19937 gen(ran());

	std::uniform_int_distribution<int> _index(0, N - 1);
	std::uniform_real_distribution<double> _accept(0, 1);
	std::normal_distribution<double> _delta_x(0, delta_x_width);

	int index;
	double Delta_x;
	double Delta_S = 0;

	int acceptance_ratio = 0;

	long i;
	for (i = 0; i < steps; i++)
	{
		index = _index(gen);
		Delta_x = _delta_x(gen);
		Delta_S = DS(lattice, Delta_x, index);
		if (_accept(gen) < exp(-Delta_S))
		{
			++acceptance_ratio;
			lattice[index] += Delta_x;
		}

		if (i % N == 0)
		{
			observable(lattice);
		}
	}

	std::cout << "acceptance ratio: " << acceptance_ratio / (float) steps << std::endl;
}
