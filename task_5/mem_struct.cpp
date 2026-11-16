#include "./mem_struct.h"

void SharedMemory::set(int matrix_order, int block_order, int number_of_threads, int output_limit, int formula, const char* input_file_name)
{
    matrix_size = matrix_order;
    block_size = block_order;
    remainder_size = matrix_size % block_size;
    full_block_count = matrix_size / block_size;
    block_count = full_block_count + (remainder_size > 0);
    thread_count = number_of_threads;
    matrix_formula = formula;
    print_limit = output_limit;
    file_name = input_file_name;
}

bool SharedMemory::create_and_init_barrier()
{
    try
    {
        barrier = new pthread_barrier_t;
    }
    catch (const std::bad_alloc&)
    {
        return 0;
    }

    if (pthread_barrier_init(barrier, 0, thread_count) != 0)
    {
        delete barrier;
        barrier = 0;
        return 0;
    }

    return 1;
}

bool SharedMemory::allocate()
{
    try
    {
        coeff_matrix = new double[(size_t)matrix_size * (size_t)matrix_size];
        right_hand_side = new double[matrix_size];
        solution = new double[matrix_size];
        block_column_order = new int[full_block_count];
        inverse_pivot_block = new double[(size_t)block_size * (size_t)block_size];
    }
    catch (const std::bad_alloc&)
    {
        return 0;
    }

    return 1;
}

void SharedMemory::free()
{
    if (coeff_matrix)
        delete[] coeff_matrix;
    if (right_hand_side)
        delete[] right_hand_side;
    if (solution)
        delete[] solution;
    if (block_column_order)
        delete[] block_column_order;
    if (inverse_pivot_block)
        delete[] inverse_pivot_block;
    if (barrier)
    {
        pthread_barrier_destroy(barrier);
        delete barrier;
    }

    coeff_matrix = 0;
    right_hand_side = 0;
    solution = 0;
    block_column_order = 0;
    inverse_pivot_block = 0;
    barrier = 0;
}

SharedMemory::~SharedMemory()
{
    pthread_mutex_destroy(&launch_mutex);
    free();
}

bool ThreadMemory::allocate(int matrix_size, int block_size)
{
    try
    {
        pivot_row_tail = new double[(size_t)block_size * (size_t)matrix_size];
        inverse_candidate = new double[(size_t)block_size * (size_t)block_size];
        working_block = new double[(size_t)block_size * (size_t)block_size];
        best_inverse_block = new double[(size_t)block_size * (size_t)block_size];
        column_order = new int[block_size];
        vector_product = new double[block_size];
    }
    catch (const std::bad_alloc&)
    {
        return 0;
    }

    return 1;
}

void ThreadMemory::free()
{
    if (pivot_row_tail)
        delete[] pivot_row_tail;
    if (inverse_candidate)
        delete[] inverse_candidate;
    if (working_block)
        delete[] working_block;
    if (best_inverse_block)
        delete[] best_inverse_block;
    if (column_order)
        delete[] column_order;
    if (vector_product)
        delete[] vector_product;

    pivot_row_tail = 0;
    inverse_candidate = 0;
    working_block = 0;
    best_inverse_block = 0;
    column_order = 0;
    vector_product = 0;
}

ThreadMemory::~ThreadMemory()
{
    free();
}
