
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
//	monitor.prime("mean:", n);

	for (; i < n; i++)
	{
//		monitor++;
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
	monitor.prime("unnormalized autocorrelation function:", n);
	
	for (; i < n; i++)
	{
		monitor++;

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

double Process::exponential_autocorrelation_time()
{
	double tau = 1;
	double temp, temp2;
	int i = 0;

	double log_C_0 = std::log(C[0]);


	/* limit supremum */

//	LOG("\nexponential autocorrelation time");

	// -i / std::log(C[i])

increasing:
	
	for (; i < n; i++)
	{
		if (C[i] == 0) continue;
		temp = - i / (std::log(C[i]) - log_C_0);
		if (temp > tau) tau = temp;
		else 
		{
//			LOG("local max %lf, C[i] = %lf; C[0] = %lf", tau, C[i], C[0]);
			temp2 = temp;
			goto decreasing;
		}
	}

decreasing:

	for (; i < n; i++)
	{
		if (C[i] == 0) continue;
		if (temp2 < (temp = - i / (std::log(C[i]) - log_C_0)))
		{
			tau = temp;
			goto increasing;
		}
		temp2 = temp;
	}

//	LOG("done");
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
double Process::integrated_autocorrelation_time(long _M, float c) 
{
	int i;
	double tau = 0.5;
	M = _M;
//	monitor.prime("inegrated autocorrelation time:", M);

	for (i = 1; i < M; i++)
	{
//		monitor++;
		tau += C[i] / C[0];
	}

	while (M <= c*tau && M < n)
	{
		tau += C[M] / C[0];
		M++;
//		monitor++;
//		monitor.N++;
	}

	return tau;
}

