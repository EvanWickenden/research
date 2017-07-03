#ifndef __MONITOR_H__
#define __MONITOR_H__

#include <cstdio>
#include <pthread.h>
#include <semaphore.h>

struct Monitor
{
	const char *msg;
	volatile long N;
	volatile long i;

	pthread_t thread;

	static void* execute(Monitor *m);

	Monitor();

	void prime(const char *msg, long N_initial);
	Monitor& operator++(int);

	void start(); 
	void end();
};

#endif
