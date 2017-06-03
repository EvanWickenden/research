#ifndef __PROCESS_H__
#define __PROCESS_H__

#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <vector>

class Process
{
	std::vector<float>& data;
	int n; 				/* total number of samples */
	int tau_exp; 		/* exponential autocorrelation time; treat as starting index */
	double tau_int; 	/* integrated autocorrelation time */
	double mu; 			/* mean */
	double var;			/* variance */
	double* C[2]; 		/* unnormalized autocorrelation function */
	double M;			/* windowing cutoff */

	public:

		Process(std::vector<float>& data, int n, double tau_exp = 0) : data(data), n(n), tau_exp(tau_exp) {}
};



#endif
