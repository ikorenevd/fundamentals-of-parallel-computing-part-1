#include "./thread_function.h"
#include "./headers_and_defines.h"
#include "./linear_system_init_functions.h"
#include "./linear_system_solver.h"
#include "./matrix_io.h"
#include "./residuals.h"
#include "./thread_args.h"

static double get_full_time(void)
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1e6;
}

// double get_cpu_time(void)
// {
//     struct timespec buf;
//     if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &buf) != 0)
//         return -1.;
//     return buf.tv_sec + buf.tv_nsec / 1e9;
// }

static void set_thread_affinity([[maybe_unused]] int thread_index)
{
    // закреплеяем тред за ядром

    #ifdef __linux__
        cpu_set_t cpu;
        CPU_ZERO(&cpu);
        int cpu_count = get_nprocs();
        int cpu_index = cpu_count - 1 - (thread_index % cpu_count);
        CPU_SET(cpu_index, &cpu);
        pthread_t thread = pthread_self();
        pthread_setaffinity_np(thread, sizeof(cpu), &cpu);
    #endif
}

void* thread_function(void* ptr)
{
    ThreadArgs* arg = static_cast<ThreadArgs*>(ptr);

    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int print_limit = arg->shared->print_limit;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    pthread_barrier_t* barrier = arg->shared->barrier;
    const int thread_index = arg->local.thread_index;

    pthread_mutex_lock(&arg->shared->launch_mutex);
    bool launch_succeeded = arg->shared->launch_succeeded;
    pthread_mutex_unlock(&arg->shared->launch_mutex);

    if (!launch_succeeded)
        return 0;

    set_thread_affinity(thread_index);
    initialize_thread_memory(arg);

    // инициализация данных
    if (!initialize_system(arg, 1))
        return 0;

    pthread_barrier_wait(barrier);

    if (thread_index == 0)
    {
        print_matrix(matrix_size, matrix_size, block_size, coeff_matrix, print_limit, "Coefficient matrix\n");
        print_matrix(1, matrix_size, block_size, right_hand_side, print_limit, "Right hand side:\n");
    }

    if (thread_index == 0)
        arg->shared->solve_time = get_full_time();

    pthread_barrier_wait(barrier);

    if (!forward_elimination(arg))
    {
        if (thread_index == 0)
        {
            arg->shared->solve_time = get_full_time() - arg->shared->solve_time;
            std::fprintf(stderr, "Error: couldn't find pivot block\n");
        }
        return 0;
    }

    if (!solve_remainder_block(arg))
    {
        if (thread_index == 0)
        {
            arg->shared->solve_time = get_full_time() - arg->shared->solve_time;
            std::fprintf(stdout, "Error: couldn't solve small system\n");
        }
        return 0;
    }

    backward_substitution(arg);

    if (thread_index == 0)
        arg->shared->solve_time = get_full_time() - arg->shared->solve_time;

    pthread_barrier_wait(barrier);

    // реинициализация данных
    if (!initialize_system(arg, 0))
        return 0;

    if (thread_index == 0)
        arg->shared->residual_time = get_full_time();

    pthread_barrier_wait(barrier);

    compute_residuals(arg);

    pthread_barrier_wait(barrier);

    if (thread_index == 0)
        arg->shared->residual_time = get_full_time() - arg->shared->residual_time;

    return 0;
}

// if (thread_index == 0)
// {
//     std::fprintf(stdout, "Alpha: %d\n", full_block_count + 1);
//     print_matrix(matrix_size, matrix_size, block_size, coeff_matrix, print_limit, "Coefficient matrix\n");
//     print_matrix(1, matrix_size, block_size, right_hand_side, print_limit, "Right hand side:\n");
// }
