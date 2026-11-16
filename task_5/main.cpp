#include "./headers_and_defines.h"
#include "./floating_point_traps.h"
#include "./matrix_io.h"
#include "./thread_args.h"
#include "./thread_function.h"

static bool parse_arguments(int argc, char** argv, SharedMemory& shared)
{
    int matrix_size = 0;
    int block_size = 0;
    int print_limit = 0; // не более чего выводим
    int matrix_formula = 0; // тип функции
    int thread_count = 0;
    char* file_name = nullptr;

    // считываем входные данные
    if (!((argc == 6 || argc == 7)
        && (std::sscanf(argv[1], "%d", &matrix_size) == 1)
        && (std::sscanf(argv[2], "%d", &block_size) == 1)
        && (std::sscanf(argv[3], "%d", &thread_count) == 1)
        && (std::sscanf(argv[4], "%d", &print_limit) == 1)
        && (std::sscanf(argv[5], "%d", &matrix_formula) == 1)))
    {
        std::printf("Usage: %s n m p r s [filename]\n", argv[0]);
        return 0;
    }

    if (argc == 7)
        file_name = argv[6];

    // проверяем на корректность
    if (matrix_size < 1
        || block_size < 1
        || print_limit < 0
        || thread_count < 1
        || (matrix_formula > 4)
        || (matrix_formula == 0 && file_name == nullptr)
        || (argc == 6 && matrix_formula < 0)
        || (argc == 7 && matrix_formula != 0))
    {
        std::fprintf(stderr, "Wrong arguments\n");
        return 0;
    }

    shared.set(matrix_size, block_size, thread_count, print_limit, matrix_formula, file_name);
    return 1;
}

static bool allocate_thread_memory(SharedMemory& shared, ThreadArgs* args)
{
    const int thread_count = shared.thread_count;

    // инициализируем аргументы и выделяем локальные данные потоков
    for (int thread_index = 0; thread_index < thread_count; thread_index++)
    {
        args[thread_index].shared = &shared;
        args[thread_index].local.thread_index = thread_index;
        if (!args[thread_index].local.allocate(shared.matrix_size, shared.block_size))
            return 0;
    }

    return 1;
}

static bool run_threads(SharedMemory& shared, ThreadArgs* args, pthread_t* threads)
{
    const int thread_count = shared.thread_count;

    // работаем == запускаем потоки
    pthread_mutex_lock(&shared.launch_mutex);
    int started = 1;
    for (; started < thread_count; started++)
        if (pthread_create(threads + started, 0, thread_function, args + started) != 0)
            break;

    shared.launch_succeeded = started == thread_count;
    pthread_mutex_unlock(&shared.launch_mutex);

    if (shared.launch_succeeded)
        thread_function(args + 0);

    for (int thread_index = 1; thread_index < started; thread_index++)
        pthread_join(*(threads + thread_index), 0);

    return shared.launch_succeeded;
}

static void print_results(const char* program_name, const SharedMemory& shared, const ThreadArgs* args)
{
    const int matrix_size = shared.matrix_size;
    const int block_size = shared.block_size;
    const int thread_count = shared.thread_count;
    const int print_limit = shared.print_limit;
    const int matrix_formula = shared.matrix_formula;
    int task_number = 11;
    double relative_residual = -1.;
    double relative_error = -1.;

    // решение закончено
    bool solved = 1;
    for (int thread_index = 0; thread_index < thread_count; thread_index++)
        solved = solved && args[thread_index].local.succeeded;

    if (solved)
    {
        print_matrix(1, matrix_size, block_size, shared.solution, print_limit, "Solution x:\n");

        double residual = 0.;
        double rhs_norm = 0.;
        double solution_error = 0.;

        for (int thread_index = 0; thread_index < thread_count; thread_index++)
        {
            residual += args[thread_index].local.residual;
            rhs_norm += args[thread_index].local.rhs_norm;
            solution_error += args[thread_index].local.solution_error;
        }

        if (rhs_norm > 0.)
            relative_residual = residual / rhs_norm;
        else
            relative_residual = (residual <= 0.) ? 0. : HUGE_VAL;

        relative_error = solution_error / (matrix_size / 2 + matrix_size % 2);
    }

    // конечный вывод
    std::printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", program_name, task_number, relative_residual, relative_error, shared.solve_time, shared.residual_time, matrix_formula, matrix_size, block_size, thread_count);
}

int main(int argc, char** argv)
{
    if (!configure_floating_point_traps())
    {
        std::fprintf(stderr, "error: cannot enable floating-point traps\n");
        return 1;
    }

    SharedMemory shared;
    ThreadArgs*  args = 0;
    pthread_t*   threads = 0;

    if (!parse_arguments(argc, argv, shared))
        return 1;

    // аргументы потоков и потоки
    try
    {
        args = new ThreadArgs[shared.thread_count];
        threads = new pthread_t[shared.thread_count];
    }
    catch (const std::bad_alloc&)
    {
        std::fprintf(stderr, "Error: couldn't allocate memory for args or threads\n");

        if (args)
            delete[] args;
        if (threads)
            delete[] threads;

        return 1;
    }

    // общие данные потоков
    if (!shared.allocate())
    {
        std::fprintf(stderr, "Error: couldn't allocate memory for general data\n");

        delete[] args;
        delete[] threads;

        return 1;
    }

    // барьер
    if (!shared.create_and_init_barrier())
    {
        std::fprintf(stderr, "Error: couldn't create or allocate barrier\n");

        delete[] args;
        delete[] threads;

        return 1;
    }

    if (!allocate_thread_memory(shared, args))
    {
        std::fprintf(stderr, "Error: couldn't allocate threads own memory\n");

        delete[] args;
        delete[] threads;

        return 1;
    }

    if (!run_threads(shared, args, threads))
    {
        std::fprintf(stderr, "Error: couldn't create threads\n");

        delete[] args;
        delete[] threads;

        return 1;
    }

    print_results(argv[0], shared, args);

    if (args)
        delete[] args;
    if (threads)
        delete[] threads;

    return 0;
}
