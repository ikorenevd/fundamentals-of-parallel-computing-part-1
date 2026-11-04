#ifndef HEADERS_AND_DEFINES_H
#define HEADERS_AND_DEFINES_H

// сами потоки
#include <pthread.h>
// (инвалидное яблоко не имеет реализацию барьеров)
#ifdef __APPLE__
    #include "./pthread_barrier_t.h"
#endif

// машинный эпсилон
#define MACHINE_EPS __DBL_EPSILON__

// для распределения тредов по ядрам (инвалидное яблоко не позволяет этого сделать)
#ifdef __linux__
    #include <sched.h>
    #include <sys/sysinfo.h>
#endif

// стандартные хедеры
#include <cstdio>
#include <cfloat>
#include <algorithm>
#include <time.h>
#include <new>
#include <cstring>
#include <cmath>

#include <sys/time.h>


#endif