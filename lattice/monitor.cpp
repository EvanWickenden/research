
#include <unistd.h>
#include <pthread.h>
#include <cstdlib>

#include "monitor.h"



void *Monitor::execute(Monitor *m)
{
	while (1)
	{
		fprintf(stderr, "\r%s  %ld / %ld; %.4f %%", m->msg, m->i, m->N, 100 * (m->i / (float) (m->N)));
		sleep(1);
	}
}

void Monitor::prime(const char *msg, long N_initial)
{
	this->msg = msg;
	this->i = 0;
	this->N = N_initial;
	fprintf(stderr, "\n");
}


Monitor& Monitor::operator++(int null)
{
	i++;
	return *this;
}

void Monitor::start()
{
	pthread_create(&thread, NULL, (void *(*)(void *)) Monitor::execute, this);
}


void Monitor::end()
{
	fprintf(stderr, "\n");
	pthread_cancel(thread);
	pthread_join(thread, NULL);
}
