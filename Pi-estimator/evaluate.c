#include <pthread.h>
#include <semaphore.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>

struct coord
{
    float x;
    float y;
};

#define NCOORDS 400
#define NTHREADS 10

#define LOG(msg, ...)   fprintf(stderr, "log: " msg "\n", ##__VA_ARGS__)


#ifdef __DEBUG__
#define DEBUG(msg, ...)     fprintf(stderr, "debug: " msg "\n", ##__VA_ARGS__)
#else
#define DEBUG(msg, ...)     {};
#endif


struct coord coords[NCOORDS];

pthread_t workers[NTHREADS];

sem_t empty;
sem_t full;
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

int start;
int end;

long long n_runs = 0;
long long n_accepts = 0;
pthread_mutex_t stats = PTHREAD_MUTEX_INITIALIZER;
    

#define generate_random() (rand() / (float) RAND_MAX);


void *execute(void *num)
{
    struct coord coord;
    int do_accept;

    int n = (int) num;

    LOG("thread %d: entering execute", n);

produce:

    DEBUG("produce");

    coord.x = generate_random();
    coord.y = generate_random();

    if (sem_trywait(&empty) != 0) 
    {
        if (errno == EAGAIN) 
        {
            LOG("switch to consume; %d", n);
            goto consume;
        }
        else goto end;
    }

    pthread_mutex_lock(&lock); /* spinlock? */
    coords[end] = coord;
    end = (end + 1) % NCOORDS;
    pthread_mutex_unlock(&lock);
    sem_post(&full);

    goto produce;


consume:

    DEBUG("consume");

    if (sem_trywait(&full) != 0) 
    {
        if (errno == EAGAIN) 
        {
            LOG("switch to produce; %d", n);
            goto produce;
        }
        else goto end;
    }

    pthread_mutex_lock(&lock);
    start = (start + 1) % NCOORDS;
    coord = coords[start];
    pthread_mutex_unlock(&lock);
    sem_post(&empty);

    /* accept or reject; add to statistics */
    // PI
    //do_accept = (coord.x * coord.x + coord.y * coord.y < 1.0);

    // x^10 - 1
    float y = powf(coord.x + 1, 10) - 1;
    do_accept = (coord.y * pow(2,10) < y);

    pthread_mutex_lock(&stats);
    n_runs++;
    n_accepts += do_accept;
    pthread_mutex_unlock(&stats);

    goto consume;

end: 
    LOG("exeting execute");
    return NULL;
}




int main()
{   
    int i;
    srand((unsigned int) &i);

    sem_init(&empty, 
            0, /* use among threads */
            NCOORDS);
    sem_init(&full,
            0,
            0);

    end = start = 0;

    sigset_t sigset;
    sigemptyset(&sigset);
    sigaddset(&sigset, SIGINT);
    pthread_sigmask(SIG_BLOCK, &sigset, NULL);

    for (i = 0; i < NTHREADS; i++)
    {
        pthread_create(workers + i, NULL, execute, NULL);
    }

    sigset_t pending;
    while (sigpending(&pending) == 0)
    {
        if (sigismember(&pending, SIGINT))
        {
            printf("Estimate for pi: N_TOTAL = %lld; N_ACCEPTED = %lld; integral ~ %f \n", 
                    n_runs, 
                    n_accepts,
                    pow(2,10) * n_accepts / (float) n_runs); 
        }
    }
}


