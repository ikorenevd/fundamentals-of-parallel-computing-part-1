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
    try
    {
        barrier = new pthread_barrier_t;
    }
    catch (const std::bad_alloc&)
    {
        return 0;
    }

    if (pthread_barrier_init(barrier, 0, threads_count) != 0)
    {
        delete barrier;
        barrier = 0;
        return 0;
    }

    return 1;
}

bool GENERAL_MEMORY::allocate()
{
    try
    {
        coeff_matrix = new double[(size_t)matrix_size * (size_t)matrix_size];
        right_hand_side = new double[matrix_size];
        solution = new double[matrix_size];
        block_col_perm = new int[full_block_count];
        inverse_pivot_block = new double[(size_t)block_size * (size_t)block_size];
    }
    catch (const std::bad_alloc&)
    {
        return 0;
    }
    
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
    try
    {
        pivot_row = new double[(size_t)block_size * (size_t)matrix_size];
        inverse_pivot_block_candidate = new double[(size_t)block_size * (size_t)block_size];
        current_block = new double[(size_t)block_size * (size_t)block_size];
        best_inverse_pivot_block = new double[(size_t)block_size * (size_t)block_size];
        perm = new int[block_size];
        buffer = new double[block_size];
    }
    catch (const std::bad_alloc&)
    {
        return 0;
    }

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
