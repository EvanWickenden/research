
#include "path-integral.h"
#include "lattice.h"
#include <math.h>
#include <iostream>
#include <cstdio>
#include <pthread.h>
#include <unistd.h>


const int N = 100;
const double tau = .01;
const int nruns = 10000000;

void *monitor_progress(PathIntegral* p)
{
	int progress;
	double percent;
	while (1)
	{
		progress = p->get_progress();
		percent =  progress / (float) nruns;
		fprintf(stderr, "\r  %d / %d :  %.4f", progress, nruns, percent);

		sleep(1);
	}
}


double harmonic_oscillator(const double& q, const double& qdot)
{
	return (pow(qdot, 2) + pow(qdot, 2));
}

void inline x1_x2(const Lattice& lattice, FILE *arg)
{
	static int i = lattice.N / 3;
	static int j = 2 * i;

	fprintf(arg, "%lf, %lf\n", lattice[i], lattice[j]);
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

	pthread_create(&progress, NULL, (void *(*)(void *)) monitor_progress, &p);

	p.run(nruns, (Observable) &x1_x2, fp);

	pthread_cancel(progress);

	char zeros[100] = { ' ' };

	std::cout << "\ndone :)" << zeros << std::endl;

	fclose(fp);
}


