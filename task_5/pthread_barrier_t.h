#ifdef __APPLE__

#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>

typedef int pthread_barrierattr_t;
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    unsigned int count;
    unsigned int tripCount;
    unsigned int generation;
} pthread_barrier_t;


int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count)
{
    if (attr != 0)
        return EINVAL;

    if(count == 0)
        return EINVAL;

    int error = pthread_mutex_init(&barrier->mutex, 0);
    if (error != 0)
        return error;
    error = pthread_cond_init(&barrier->cond, 0);
    if (error != 0)
    {
        pthread_mutex_destroy(&barrier->mutex);
        return error;
    }
    barrier->tripCount = count;
    barrier->count = 0;
    barrier->generation = 0;

    return 0;
}

int pthread_barrier_destroy(pthread_barrier_t *barrier)
{
    int error = pthread_cond_destroy(&barrier->cond);
    if (error != 0)
        return error;
    return pthread_mutex_destroy(&barrier->mutex);
}

int pthread_barrier_wait(pthread_barrier_t *barrier)
{
    int error = pthread_mutex_lock(&barrier->mutex);
    if (error != 0)
        return error;
    unsigned int generation = barrier->generation;
    ++(barrier->count);
    if(barrier->count >= barrier->tripCount)
    {
        barrier->count = 0;
        barrier->generation++;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    }
    else
    {
        // Ложное пробуждение не должно выпускать поток из текущего поколения барьера.
        while (generation == barrier->generation)
        {
            error = pthread_cond_wait(&barrier->cond, &barrier->mutex);
            if (error != 0)
            {
                pthread_mutex_unlock(&barrier->mutex);
                return error;
            }
        }
        pthread_mutex_unlock(&barrier->mutex);
        return 0;
    }
}

#endif // PTHREAD_BARRIER_H_
#endif // __APPLE__
