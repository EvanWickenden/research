
#include <pthread.h>
#include "process.h"
#include <cstdio>


/**
 * 	\mu = 1/n \sum_i=1^n f_i
 */
double Process::mean()
{
	int i = 0;
	double sum = 0;
	progress = 0;
	for (; i < n; i++)
	{
		progress++;
		sum += data[i];
	}
	return mu = sum / n;
}


/**
 * C(t) = \frac 1 {n-|t|}  \sum_{i=1}^{n-|t|}  (f_i - \mu)(f_{i+|t|} - \mu)
 */

void Process::unnormalized_autocorrelation_function()
{
	int i = 0;
	progress = 0;
	fprintf(stderr, "\r%ld", progress);
	for (;i < n; i++)
	{
		/* this fucking optimization! */
		if (++progress == 0)
			fprintf(stderr, "\r%ld", progress);

		double sum = 0;
		int j = 0;
		for (; j < n-i; j++)
		{
			sum += (data[j] - mu) * (data[j+i] - mu);
		}	
		C[i] = sum / (n-i); 

	}
}

void Process::record_C(char *filename)
{
	FILE *fp = fopen(filename, "w");
	fwrite(C, sizeof (*C), n, fp);
}


/**
 *	\tau_exp = \limsup_{t -> \infty}  \frac {- t} {ln(C_{ff}(t)}
 */

/* overestimate by replacing limsup with max */
double Process::exponential_autocorrelation_time()
{
	double tau = 0;
	double temp;
	int i = 0;
	progress = 0;
	for (; i < n; i++)
	{
		progress++;
		temp = - i / std::log(C[i]);
		if (temp > tau) tau = temp;
	}
	return tau;
}


/**
 *	\tau_int  =  1/2  \sum_{t = -M}^{M} \rho(t)  =  1/2 + \sum_{t=1}^M \rho(t)
 *  \rho(t) = C_ff(t) / C_ff(0);
 */

/**
 *  iterative process: increase M until M >= 5 tau_exp(M)
 *	might be sensitive to properly chosen initial value for M
 */

double Process::integrated_autocorrelation_time(float c)
{
	int i;
	double tau = 0.5;
	progress = 0;
	for (i = 1; i < M; i++)
	{
		progress++;
		tau += C[i] / C[0];
	}

	while (M <= c*tau && M < n)
	{
		progress++;
		tau += C[M] / C[0];
		M++;
	}

	return tau;
}

