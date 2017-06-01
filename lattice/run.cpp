
#include "path-integral.h"
#include "lattice.h"
#include <math.h>
#include <iostream>

double lagrangian(double q, double qdot)
{
	return (pow(qdot, 2) + pow(qdot, 2));
}

void observable(const Lattice& lattice, void *arg)
{

}

int main()
{
	const int nruns = 10000000000;
	const double tau = .01;

	PathIntegral p(nruns, tau, &lagrangian);

	p.populate_lattice(0, 0);


}
