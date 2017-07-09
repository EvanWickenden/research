
#include "path-integral.h"
#include "lattice.h"
#include <math.h>
#include <iostream>
#include <cstdio>
#include <pthread.h>
#include <unistd.h>

#include "monitor.h"

#include "settings.h" /* nruns; tau; N */



double harmonic_oscillator(const double& q, const double& qdot)
{
	return (pow(qdot, 2) + pow(qdot, 2));
}

void inline x1_x2(const Lattice& lattice, FILE *arg)
{
	static int i = lattice.N / 3;
	static int j = 2 * i;

	//fprintf(arg, "%lf, %lf\n", lattice[i], lattice[j]);
	double product = lattice[i] * lattice[j];

	fwrite(&product, sizeof (double), 1, arg);

}

pthread_t progress;

int main()
{
	PathIntegral p(N, tau, &harmonic_oscillator);

	p.populate_lattice(0, 0);

	FILE *fp = fopen("HO-x1_x2.csv", "w");
	if (fp == NULL)
	{
		fprintf(stderr, "fopen() failed\n");
		exit(1);
	}


	Monitor monitor;

	monitor.start("generate observable samples", &p.progress, &nruns);

	p.run(nruns, (Observable) &x1_x2, fp);

	monitor.end();

	std::cout << "\ndone :)" << std::endl;

	fclose(fp);
}


