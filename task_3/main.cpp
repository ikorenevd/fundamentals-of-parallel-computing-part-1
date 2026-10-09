#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE 1
#endif

#define ENABLE_FPE 0

#include <cstdio>
#include <ctime>
#ifdef __linux__
#include <sched.h>
#endif

#include "./LinearSystem.h"

#if ENABLE_FPE
#include <cfenv>
// нейрослопокод
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
    int n = 0, m = 0, r = 0, s = 0;
    double t1 = 0., t2 = 0., r1 = -1., r2 = -1.;
    char* file_name = nullptr;
    bool solved = 0;

    if (!((argc == 5 || argc == 6)
        && (std::sscanf(argv[1], "%d", &n) == 1)
        && (std::sscanf(argv[2], "%d", &m) == 1)
        && (std::sscanf(argv[3], "%d", &r) == 1)
        && (std::sscanf(argv[4], "%d", &s) == 1)))
    {
        std::printf("Usage: %s n m r s\n", argv[0]);
        return 1;
    }

    if (argc == 6)
    {
        file_name = argv[5];
    }

    if ((n < 1 || m < 1 || r < 0)
        || (s > 4)
        || (s == 0 && file_name == nullptr)
        || (argc == 5 && s < 0)
        || (argc == 6 && s != 0)
    )
    {
        std::fprintf(stderr, "Wrong input\n");
        return 1;
    }

    int block_size = m > n ? n : m;

#ifdef __linux__
    // привязка к последнему доступному ядру
    cpu_set_t available_cpus;
    CPU_ZERO(&available_cpus);
    if (sched_getaffinity(0, sizeof(available_cpus), &available_cpus) == -1)
    {
        std::perror("sched_getaffinity");
        return 1;
    }

    int last_cpu = CPU_SETSIZE - 1;
    while (last_cpu >= 0 && !CPU_ISSET(last_cpu, &available_cpus))
    {
        last_cpu--;
    }
    if (last_cpu < 0)
    {
        std::fprintf(stderr, "error: no available CPUs\n");
        return 1;
    }

    cpu_set_t cpu_mask;
    CPU_ZERO(&cpu_mask);
    CPU_SET(last_cpu, &cpu_mask);
    if (sched_setaffinity(0, sizeof(cpu_mask), &cpu_mask) == -1)
    {
        std::perror("sched_setaffinity");
        return 1;
    }
#endif

    // создание системы
    LinearSystem system(n, block_size);

    // выделение памяти
    if (!system.memory_alloc())
    {
        std::fprintf(stderr, "error: allocating memory\n");
        return 1;
    }

    // инициализация
    if (!system.init(s, file_name))
    {
        std::fprintf(stderr, "error: when inizialization\n");
        return 1;
    }

    // вывод матрицы и правой части(второй вывод не нужен)
    std::printf("Matrix A:\n");
    system.print_matrix(r);

    // решаем
    std::clock_t start = std::clock();
    solved = system.solve() == 0;
    t1 = (double)(std::clock() - start) / CLOCKS_PER_SEC;

    // считаем невязки
    if (solved)
    {
        std::printf("Solution x:\n");
        system.print_solution(r);
        // восстанавливаем A и b для невязки
        if (!system.init(s, file_name))
        {
            std::fprintf(stderr, "error: restoring matrix and rhs\n");
            return 1;
        }

        start = std::clock();
        system.compute_residuals(r1, r2);
        t2 = (double)(std::clock() - start) / CLOCKS_PER_SEC;
    }

    // конечный вывод
    std::printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], task, r1, r2, t1, t2, s, n, m);

    // память освобождается в деструкторе
    return 0;
}
