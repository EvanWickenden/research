#ifndef __MONITOR_H__
#define __MONITOR_H__

#include <cstdio>
#include <pthread.h>
#include <semaphore.h>

class Monitor
{
	/* inelegant but seemingly necessary */
	struct ThreadArgs
	{
		const long &i;
		const long &N;
		ThreadArgs(const long &i, const long &N);
	} *thread_args;

	FILE *fp;
	pthread_t thread;

	static void* execute(Monitor *m);

	public:

	Monitor(FILE *fp);

	void start(const long &i, const long &N); 
	void restart(const long &i, const long &N); 
	void end();
};

#endif
