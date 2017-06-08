
#include <stdio.h>

#include "process.h"
#include "settings.h"
#include "log.h"

#include "monitor.h"


int main()
{
	double *data = new double[nruns / N]; /* a large array */

	FILE* f = fopen("HO-x1_x2.csv", "r");
	if (f == NULL) DIE("fopen() failed");

	int ret;
	//if ((ret = fread(data, sizeof (double), nruns / N, f)) < nruns)
	//	DIE("fread() failed; %d", ret);
	int n = fread(data, sizeof (double), nruns / N, f);

	Process process(data, n, 5000);

	long _n = n;

	Monitor monitor;
	monitor.start("mean:", &process.progress, &_n);

	double mu = process.mean();

	LOG("what the fuck?? %p", &process.progress);

	monitor.restart("unnormalzide autocorrelation function:", &process.progress, &_n);

	process.unnormalized_autocorrelation_function();

	monitor.end();

//	process.record_C("unnormalized_autocorrelation_function");


//	double tau_exp = process.exponential_autocorrelation_time();

	double tau_int = process.integrated_autocorrelation_time(1);



//	monitor.start(process.progress, process.M);

	tau_int = process.integrated_autocorrelation_time(100);


	LOG("tau_int = %lf", tau_int);

//	LOG("mean = %f; tau_exp = %f; tau_int = %f", mu, tau_exp, tau_int);
}
