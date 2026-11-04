#include "./headers_and_defines.h"

#include "./mem_struct.h"
#include "./matrix_io.h"
#include "./Args.h"
#include "./thread_function.h"

#define ENABLE_FPE 0

#if ENABLE_FPE
#include <cfenv>

namespace
{
    bool enable_fpe()
    {
        const int exceptions = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
        if (std::feclearexcept(exceptions) != 0)
            return false;
        #if defined(__linux__)
            return ::feenableexcept(exceptions) != -1;
        #elif defined(__APPLE__)
            std::fenv_t env;
            if (std::fegetenv(&env) != 0)
                return false;
            #if defined(__aarch64__) || defined(__arm64__)
                // ARM FPCR exception-enable bits are shifted by 8 from FE_* flags.
                const auto traps = static_cast<decltype(env.__fpcr)>(exceptions) << 8;
                env.__fpcr |= traps;
                if (std::fesetenv(&env) != 0 || std::fegetenv(&env) != 0)
                    return false;
                return (env.__fpcr & traps) == traps;
            #elif defined(__x86_64__) || defined(__i386__)
                // Unmask exceptions in both the x87 and SSE control registers.
                env.__control = static_cast<decltype(env.__control)>(env.__control & ~exceptions);
                env.__mxcsr &= ~(static_cast<unsigned int>(exceptions) << 7);
                if (std::fesetenv(&env) != 0 || std::fegetenv(&env) != 0)
                    return false;
                return (env.__control & exceptions) == 0
                    && (env.__mxcsr & (static_cast<unsigned int>(exceptions) << 7)) == 0;
            #else
                #error Unsupported macOS architecture for floating-point traps
            #endif
        #else
            #error Floating-point traps require Linux or macOS
        #endif
    }
}
#endif


int main(int argc, char** argv)
{
    #if ENABLE_FPE
        if (!enable_fpe())
        {
            std::fprintf(stderr, "error: cannot enable floating-point traps\n");
            return 1;
        }
    #endif

    int task = 11;
    int matrix_size = 0;
    int block_size = 0;
    int r = 0; // не более чего выводим
    int s = 0; // тип функции
    int threads_count = 0;
    double r1 = -1.;
    double r2 = -1.;
    char* file_name = nullptr;

    GENERAL_MEMORY general_data;
    Args*          args = 0;
    pthread_t*     threads = 0;

    // считываем входные данные
    if (!((argc == 6 || argc == 7)
        && (std::sscanf(argv[1], "%d", &matrix_size) == 1)
        && (std::sscanf(argv[2], "%d", &block_size) == 1)
        && (std::sscanf(argv[3], "%d", &threads_count) == 1)
        && (std::sscanf(argv[4], "%d", &r) == 1)
        && (std::sscanf(argv[5], "%d", &s) == 1)))
    {
        std::printf("Usage: %s n m p r s [filename]\n", argv[0]);
        return 1;
    }

    if (argc == 7)
        file_name = argv[6];

    // проверяем на корректность
    if (matrix_size < 1
        || block_size < 1
        || r < 0 
        || threads_count < 1
        || (s > 4)
        || (s == 0 && file_name == nullptr)
        || (argc == 6 && s < 0)
        || (argc == 7 && s != 0))
    {
        std::fprintf(stderr, "Wrong arguments\n");
        return 1;
    }

    // аргументы потоков и потоки 
    try
    {
        args = new Args[threads_count];
        threads = new pthread_t[threads_count];
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
    general_data.set(matrix_size, block_size, threads_count, r, s, file_name);
    if (!general_data.allocate())
    {
        std::fprintf(stderr, "Error: couldn't allocate memory for general data\n");
        
        delete[] args;
        delete[] threads;
        
        return 1;
    }

    // барьер
    if (!general_data.create_and_init_barrier())
    {
        std::fprintf(stderr, "Error: couldn't create or allocate barrier\n");

        delete[] args;
        delete[] threads;

        return 1;
    }

    // инициализируем аргументы и выделяем локальные данные потоков
    for (int thread_num = 0; thread_num < threads_count; thread_num++)
    {
        args[thread_num].gen = &general_data;
        args[thread_num].own.thread_number = thread_num;
        if (!args[thread_num].own.allocate(matrix_size, block_size))
        {
            std::fprintf(stderr, "Error: couldn't allocate threads own memory\n");

            delete[] args;
            delete[] threads;

            return 1;
        }
    }

    // работаем == запускаем потоки
    for (int th_num = 1; th_num < threads_count; th_num++)
    {
        if (pthread_create(threads + th_num, 0, thread_function, args + th_num) != 0)
        {
            std::fprintf(stderr, "Error: couldn't create threads\n");

            delete[] args;
            delete[] threads;

            return 1;
        }
    }
    thread_function(args + 0);

    for (int th_num = 1; th_num < threads_count; th_num++)
        pthread_join(*(threads + th_num), 0);

    // решение закончено
    bool solved = 1;
    for (int th_num = 0; th_num < threads_count; th_num++)
        solved = solved && args[th_num].own.okay;

    if (solved)
    {
        print_matrix(1, matrix_size, block_size, general_data.solution, r, "Solution x:\n");

        double residual = 0.;
        double norm_b = 0.;
        double error = 0.;

        for (int th_num = 0; th_num < threads_count; th_num++)
        {
            residual += args[th_num].own.residual;
            norm_b += args[th_num].own.norm_b;
            error += args[th_num].own.error;
        }

        if (norm_b > 0.)
            r1 = residual / norm_b;
        else
            r1 = (residual <= 0.) ? 0. : HUGE_VAL;

        r2 = error / (matrix_size / 2 + matrix_size % 2);
    }

    // конечный вывод
    std::printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], task, r1, r2, general_data.t1, general_data.t2, s, matrix_size, block_size, threads_count);

    if (args)
        delete[] args;
    if (threads)
        delete[] threads;

    return 0;
}
