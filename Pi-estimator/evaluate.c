#include <pthread.h>
#include <semaphore.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

struct coord
{
    float x;
    float y;
};

#define NCOORDS 400
#define NTHREADS 16

#define LOG(msg, ...)   fprintf(stderr, "log: " msg "\n", ##__VA_ARGS__)

struct coord coords[NCOORDS];

pthread_t workers[NTHREADS];

sem_t empty;
sem_t full;
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

struct coord *start;
struct coord *end;

long long n_runs;
long long n_accepts;
pthread_mutex_t stats = PTHREAD_MUTEX_INITIALIZER;


#define generate_random() (rand() / (float) RAND_MAX);


void *execute(void *null)
{
    struct coord *index;
    struct coord coord;
    int do_accept;

    LOG("entering execute");

produce:
    coord.x = generate_random();
    coord.y = generate_random();

    if (sem_trywait(&empty) != 0) 
    {
        if (errno == EAGAIN) goto consume;
        goto end;
    }
    pthread_mutex_lock(&lock); /* spinlock? */
    *++end = coord;
    pthread_mutex_unlock(&lock);
    sem_post(&full);

    goto produce;


consume:
    if (sem_trywait(&full) != 0) 
    {
        if (errno == EAGAIN) goto produce;
        goto end;
    }

    pthread_mutex_lock(&lock);
    coord = *++start;
    pthread_mutex_unlock(&lock);
    sem_post(&empty);

    /* accept or reject; add to statistics */
    do_accept = (coord.x * coord.x + coord.y * coord.y < 1.0);

    pthread_mutex_lock(&stats);
    n_runs++;
    n_accepts += do_accept;
    pthread_mutex_unlock(&stats);

    goto consume;

    /* terminate execution synchronously by destroying semaphores */
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

    end = start = coords;

    for (i = 0; i < NTHREADS; i++)
        pthread_create(workers + i, NULL, execute, NULL);

    sleep(10);
    
    sem_destroy(&empty);
    sem_destroy(&full);

    for (i = 0; i < NTHREADS; i++)
        pthread_join(workers[i], NULL); 

    printf("Estimate for pi: N_TOTAL = %lld; N_ACCEPTED = %lld;\n", 
            n_runs, 
            n_accepts); 
}


