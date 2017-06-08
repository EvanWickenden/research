#ifndef __PROCESS_H__
#define __PROCESS_H__

#include <math.h>


struct Process
{
	double *data; 		/* data to analyze */
	int n;				/* dimensions of data array */
	double mu; 			/* mean */
	double *C; 			/* unnormalized autocorrelation function */
	long M;				/* autocorrelation window */
	long progress;

	/* calculate the unnormalized autocorrelation function over the domain provided */
	Process(double *data, int n, int M = 100) : 
		data(data), 
		n(n), 
		mu(0), 
		C(new double[n]),
		M(M), /* sensible initial value depends on data set */
		progress(0)
	{
		mean();
	}

	~Process() { delete[] C; }

	long& get_progress() { return progress; }

	double mean();
	void unnormalized_autocorrelation_function();
	void record_C(char *name);
	double exponential_autocorrelation_time();
	double integrated_autocorrelation_time(float c = 5.0);
};



#endif
