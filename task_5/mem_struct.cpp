#include "./mem_struct.h"

void GENERAL_MEMORY::set(int m_size, int b_size, int th_count, int _r, int _s, char* _file_name)
{
    matrix_size = m_size;
    block_size = b_size;
    remainder_size = matrix_size % block_size;
    full_block_count = matrix_size / block_size;
    block_count = full_block_count + (remainder_size > 0);
    threads_count = th_count;
    s = _s;
    r = _r;
    file_name = _file_name;
}

bool GENERAL_MEMORY::create_and_init_barrier()
{
    barrier = new pthread_barrier_t;

    if (!barrier)
        return 0;

    if (pthread_barrier_init(barrier, 0, threads_count) != 0)
        return 0;

    return 1;
}

bool GENERAL_MEMORY::allocate()
{
    coeff_matrix = new double[matrix_size * matrix_size];
    if (coeff_matrix == 0)
        return 0;

    right_hand_side = new double[matrix_size];
    if (right_hand_side == 0)
        return 0;
    
    solution = new double[matrix_size];
    if (solution == 0)
        return 0;
    
    block_col_perm = new int[full_block_count];
    if (block_col_perm == 0)
        return 0;
    
    inverse_pivot_block = new double[block_size * block_size];
    if (inverse_pivot_block == 0)
        return 0;
    
    return 1;
}

void GENERAL_MEMORY::free()
{
    if (coeff_matrix)
        delete[] coeff_matrix;
    if (right_hand_side)
        delete[] right_hand_side;
    if (solution)
        delete[] solution;
    if (block_col_perm)
        delete[] block_col_perm;
    if (inverse_pivot_block)
        delete[] inverse_pivot_block;
    if (barrier)
        delete barrier;

    coeff_matrix = 0;
    right_hand_side = 0;
    solution = 0;
    block_col_perm = 0;
    inverse_pivot_block = 0;
}

GENERAL_MEMORY::~GENERAL_MEMORY()
{
    if (barrier)
        pthread_barrier_destroy(barrier);
    free();
}

bool THREAD_MEMORY::allocate(int matrix_size, int block_size)
{
    pivot_row = new double[block_size * matrix_size];
    if (pivot_row == 0)
        return 0;

    inverse_pivot_block_candidate = new double[block_size * block_size];
    if (inverse_pivot_block_candidate == 0)
        return 0;

    current_block = new double[block_size * block_size];
    if (current_block == 0)
        return 0;

    best_inverse_pivot_block = new double[block_size * block_size];
    if (best_inverse_pivot_block == 0)
        return 0;
    
    perm = new int[block_size];
    if (perm == 0)
        return 0;
    
    buffer = new double[block_size];
    if (buffer == 0)
        return 0;

    return 1;
}

void THREAD_MEMORY::free()
{
    if (pivot_row)
        delete[] pivot_row;
    if (inverse_pivot_block_candidate)
        delete[] inverse_pivot_block_candidate;
    if (current_block)
        delete[] current_block;
    if (best_inverse_pivot_block)
        delete[] best_inverse_pivot_block;
    if (perm)
        delete[] perm;
    if (buffer)
        delete[] buffer;
}

THREAD_MEMORY::~THREAD_MEMORY()
{
    free();
}