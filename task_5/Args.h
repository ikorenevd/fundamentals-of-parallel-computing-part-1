#ifndef ARGS_H
#define ARGS_H

#include <pthread.h>

#ifdef __APPLE__
    #include "./pthread_barrier_t.h"
#endif

class Args
{
public:
    // общие данные системы
    double* coeff_matrix    = 0;
    double* right_hand_side = 0;
    double* solution        = 0;
    int*    block_col_perm  = 0; // перестановки матричных столбцов
    double* pivot_block_inverse = 0;

    int matrix_size       = 0;
    int block_size        = 0;
    int full_blocks_count = 0;
    int remainder_size    = 0;
    int blocks_count      = 0;
    double* original_matrix_norm = 0;

    // локальные блоки - для каждого потока свои - своя память
    double* pivot_block_row = 0;
    double* block1 = 0;
    double* block2 = 0;
    double* block3 = 0;
    int*    perm   = 0;
    double* buffer = 0;

    // своя память необходимая для поиска лучшего обратного блока
    double best_block_norm = 0.;
    int    best_block_i = 0;
    int    best_block_j = 0;

    // информация о потоках и ошибках
    int threads_count = 0;
    int thread_number = 0;
    bool mem_alloc = 0;

    // инфорация о времени
    double t1 = 0.;
    double t2 = 0.;

    // информация о запуске
    int r           = 0;
    int s           = 0;
    char* file_name = 0;

    // баррьер
    pthread_barrier_t* barrier;

    Args() = default;
    Args(const Args&) = delete;
    Args(Args&&) = delete;
    Args& operator=(const Args&) = delete;
    Args& operator=(Args&&) = delete;

    ~Args() = default;
};


#endif
