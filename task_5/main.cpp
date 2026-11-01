#include <cstdio>
#include <ctime>

#include <algorithm>

#include <pthread.h>
#ifdef __APPLE__
    #include "./pthread_barrier_t.h"
#endif

#ifdef __linux__
    #include <pthread_barrier.h>
#endif
#include "./LinearSystem.h"
#include "./Args.h"

void* thread_function(void* a)
{
    Args* args = (Args*)a;

    double* coeff_matrix    = args->coeff_matrix;
    double* right_hand_side = args->right_hand_side;
    double* solution        = args->solution;
    int*    block_col_perm  = args->block_col_perm; // перестановки матричных столбцов
    double* pivot_block_inverse = args->pivot_block_inverse;

    int matrix_size      = args->matrix_size;
    int block_size       = args->block_size;
    int full_block_count = args->full_block_count;
    int remainder_size   = args->remainder_size;
    int block_count      = args->block_count;
    double* original_matrix_norm = args->original_matrix_norm;

    int threads_count = args->threads_count;
    int thread_number = args->thread_number;

    pthread_barrier_t* barrier = args->barrier;

    std::fprintf(stdout, "hello from %d thread num out of %d amount\n", thread_number, threads_count);

    // привязка к процессору
    #ifdef __linux__
        cpu_set_t cpu;
        CPU_ZERO(&cpu);
        int n_cpus = get_nprocs();
        int cpu_id = n_cpus - 1 - (thread_number % n_cpus);
        CPU_SET(cpu_id, &cpu);
        pthread_t tid = pthread_self();
        pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
    #endif

    // инициализируем свою память(ту часть памяти coeff_matrix, right_hand_side, solution, с которой будет работать этот поток)
    // идем по блочным строкам этого потока
    for (int i = thread_number; i < full_block_count; i += threads_count)
    {
        std::memset(coeff_matrix + i * block_size * matrix_size, 0, block_size * matrix_size * sizeof(double));
        std::memset(right_hand_side + i * block_size * 1, 0, block_size * 1 * sizeof(double));
        std::memset(solution + i * block_size * 1, 0, block_size * 1 * sizeof(double));
    }
    
    
    std::fprintf(stdout, "after hello from %d thread num out of %d amount\n", thread_number, threads_count);
    
    // выделяем локальную память с проверкой
    
    // барьер
    pthread_barrier_wait(barrier);
    
    // основной поток проверяет все ли выделено, если нет, то нужно выйти
    return 0;
}

int main(int argc, char** argv)
{
    int task = 11;
    int matrix_size = 0,
        block_size = 0,
        r = 0,
        s = 0,
        threads_count = 0,
        full_block_count = 0,
        remainder_size = 0,
        block_count = 0;
    double t1 = 0., t2 = 0., r1 = -1., r2 = -1.;
    char* file_name = nullptr;
    // bool solved = 0;
    // общие данные системы - глобальная память
    double* coeff_matrix    = 0;
    double* right_hand_side = 0;
    double* solution        = 0;
    int*    perm            = 0; // перестановки матричных столбцов
    double* pivot_block_inverse = 0;

    // потоки и барьер
    pthread_t* pid             = 0;
    pthread_barrier_t* barrier = 0;
    Args* args                 = 0;

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
    if ((matrix_size < 1 || block_size < 1 || r < 0 || threads_count < 1)
        || (s > 4)
        || (s == 0 && file_name == nullptr)
        || (argc == 6 && s < 0)
        || (argc == 7 && s != 0))
    {
        std::fprintf(stderr, "Wrong input\n");
        return 1;
    }

    // создание системы
    block_size       = std::min(matrix_size, block_size);
    full_block_count = matrix_size / block_size;
    remainder_size   = matrix_size % block_size;
    block_count      = full_block_count + (remainder_size != 0);

    if (threads_count > block_count)
        threads_count = block_count;

    // выделение памяти (нужно только глобальную)
    {
        coeff_matrix    = new double[matrix_size * matrix_size];
        right_hand_side = new double[matrix_size];
        solution        = new double[matrix_size];
        perm            = new int[block_count];
        pivot_block_inverse = new double[block_size * block_size];

        if (!coeff_matrix || !right_hand_side || !solution || !perm || !pivot_block_inverse)
        {
            if (coeff_matrix)
                delete[] coeff_matrix;
            if (right_hand_side)
                delete[] right_hand_side;
            if (solution)
                delete[] solution;
            if(perm)
                delete[] perm;
            if (pivot_block_inverse)
                delete[] pivot_block_inverse;

            std::fprintf(stderr, "error: memory allocation\n");
            return 1;
        }

        pid     = new pthread_t[threads_count - 1];
        args    = new Args[threads_count];
        barrier = new pthread_barrier_t;

        if (!pid || !args || !barrier)
        {
            if (pid)
                delete[] pid;
            if (args)
                delete[] args;
            if (barrier)
                delete barrier;

            std::fprintf(stderr, "error: memory allocation\n");
            return 1;
        }

        int k;
        for (k = 0; k < threads_count; k++)
        {
            args[k].coeff_matrix    = coeff_matrix;
            args[k].right_hand_side = right_hand_side;
            args[k].solution        = solution;
            args[k].perm            = perm;
            
            args[k].matrix_size = matrix_size;
            args[k].block_size = block_size;
            args[k].full_block_count = full_block_count;
            args[k].remainder_size = remainder_size;
            args[k].block_count = block_count;

            args[k].pivot_block_row = new double[matrix_size * block_size];
            args[k].block1          = new double[block_size * block_size];
            args[k].block2          = new double[block_size * block_size];
            args[k].block3          = new double[block_size * block_size];
            args[k].perm            = new int[block_size];
            args[k].buffer          = new double[block_size];

            args[k].threads_count = threads_count;
            args[k].thread_number = k;

            args[k].barrier = barrier;

            if (!args[k].pivot_block_row ||
                !args[k].block1 ||
                !args[k].block2 ||
                !args[k].block3 ||
                !args[k].perm   ||
                !args[k].buffer)
            {
                for (int v = 0; v <= k; v++)
                {
                    if (args[v].pivot_block_row)
                        delete[] args[v].pivot_block_row;
                    if (args[v].block1)
                        delete[] args[v].block1;
                    if (args[v].block2)
                        delete[] args[v].block2;
                    if (args[v].block3)
                        delete[] args[v].block3;
                    if (args[v].perm)
                        delete[] args[v].perm;
                    if (args[v].buffer)
                        delete[] args[v].buffer;
                }
                
                std::fprintf(stderr, "error: mem allocation\n");
                return 1;
            }
        }
    }

    // инициализируем барьер
    if(pthread_barrier_init(barrier, 0, threads_count) != 0)
    {
        // todo обработка ошибки и отчиска памяти
        return 1;
    }

    // запускаем треды
    for (int i = 0; i < threads_count - 1; i++)
    {
        int flag = pthread_create(pid + i, 0, thread_function, args + i + 1);

        if (flag != 0)
        {
            std::fprintf(stderr, "error: could not create threads\n");

            // todo: осводобить глобальную и локальную память

            return 1;
        }
    }
    thread_function(args + 0);

    // ждем их окончания
    for (int i = 0; i < threads_count - 1; i++)
    {
        int flag = pthread_join(*(pid + i), 0);
        
        if (flag != 0)
        {
            std::fprintf(stderr, "error: could not join thread\n");
            // todo: осводобить глобальную и локальную память
        
            return 1;
        }
    }


    // конечный вывод
    std::printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], task, r1, r2, t1, t2, s, matrix_size, block_size, threads_count);
    // todo: осовбождение оставшейся памяти
    
    return 0;
}
