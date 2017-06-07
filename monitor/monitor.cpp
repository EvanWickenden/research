
#include <unistd.h>
#include <pthread.h>

#include "monitor.h"

Monitor::ThreadArgs::ThreadArgs(const long &i, const long &N) :
	i(i),
	N(N)
{}


Monitor::Monitor(FILE *fp) :
	fp(fp)
{}

void* Monitor::execute(Monitor *m)
{
	while (1)
	{
		fprintf(m->fp, "\r %ld / %ld;  %.4f %% complete", 100 * m->thread_args->i, m->thread_args->N, (m->thread_args->i / (float) m->thread_args->N));
		fflush(m->fp);
		sleep(1);
	}
}

void Monitor::start(const long &i, const long &N) 
{
	thread_args = new ThreadArgs(i, N);
	pthread_create(&thread, NULL, (void *(*)(void *)) execute, this);
}

void Monitor::restart(const long &i, const long &N) 
{
	end();
	start(i, N);
}

void Monitor::end()
{
	fprintf(fp, "\n");
	pthread_cancel(thread);
	pthread_join(thread, NULL);
	delete thread_args;
}


