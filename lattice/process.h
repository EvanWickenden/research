#ifndef __PROCESS_H__
#define __PROCESS_H__

#include <math.h>

#include "monitor.h"
#include "warray.h"

struct Process
{
	double *data; 		/* data to analyze */
	int &n;				/* dimensions of data array */
	double mu; 			/* mean */
	double *C; 			/* unnormalized autocorrelation function */
	long M;				/* autocorrelation window */

	Monitor monitor;

	Process(WArray<double> &arr) : 
		data(arr.data), 
		n(arr.N), 
		mu(0), 
		C(new double[arr.N]),
		M(1) 
	{}

	~Process() { delete[] C; }

	double mean();
	void unnormalized_autocorrelation_function();
	void record_C(char *name);
	double exponential_autocorrelation_time();
	double integrated_autocorrelation_time(long _M, float c);
};



#endif
