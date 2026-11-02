#include <cstdio>
#include <cstring>
#include <ctime>
#include <new>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <climits>
#include <utility>

#include <pthread.h>
#ifdef __APPLE__
    #include "./pthread_barrier_t.h"
#endif

#ifdef __linux__
    #include <sched.h>
    #include <sys/sysinfo.h>
#endif
#include "./LinearSystem.h"
#include "./Args.h"
#include "./matrix_io.h"

inline int matrix_offset(int n, int m, int row, int col)
{
    int i = row / m;
    int j = col / m;

    int wr = std::min(m, n - i * m);
    int wc = std::min(m, n - j * m);

    return (i * n + j * wr) * m + (row % m) * wc + col % m;
}

inline double f(int n, int s, int i, int j)
{
    int t;
    switch (s)
    {
    case 1:
        t = i > j ? i : j;
        return ((double)(n - t + 1));
        break;
    case 2:
        t = i > j ? i : j;
        return (double)t;
        break;
    case 3:
        t = i - j;
        return (double)(t >= 0 ? t : -t);
        break;
    case 4:
        return (double)1 / (double)(i + j - 1);
        break;
    }
    return 0;
}

bool init_coeff_matrix_from_file(double* coeff_matrix, int matrix_size, int block_size, char* file_name)
{
    if (file_name == nullptr)
        return 0;

    std::FILE* file = std::fopen(file_name, "r");
    if (file == nullptr)
        return 0;

    for (int i = 0; i < matrix_size; i++)
        for (int j = 0; j < matrix_size; j++)
        {
            double* value = coeff_matrix + matrix_offset(matrix_size, block_size, i, j);

            if (std::fscanf(file, "%lf", value) != 1 || std::isnan(*value) || std::isinf(*value))
            {
                std::fclose(file);
                return 0;
            }
        }

    char c;
    bool flag = std::fscanf(file, " %c", &c) == EOF && !std::ferror(file);

    std::fclose(file);
    return flag;
}

bool init_coeff_matrix_from_formula(double* coeff_matrix, int matrix_size, int block_size, int blocks_count, int s, int threads_count, int thread_number)
{
    // будет идти по блокам и заполнять в каждом блоке
    for (int blocked_row = thread_number; blocked_row < blocks_count; blocked_row += threads_count)
    {
        double* block = coeff_matrix + blocked_row * block_size * matrix_size;
        int height = std::min(block_size, matrix_size - blocked_row * block_size);

        for (int blocked_col = 0; blocked_col < blocks_count; blocked_col++)
        {
            int width = std::min(block_size, matrix_size - blocked_col * block_size);

            for (int q = 0; q < height; q++)
                for (int p = 0; p < width; p++)
                    *block++ = f(matrix_size, s, blocked_row * block_size + q + 1, blocked_col * block_size + p + 1);
        }
    }
    return 1;
}

bool init_right_hand_size(double* coeff_matrix, double* right_hand_side, int matrix_size, int block_size, int blocks_count, int threads_count, int thread_number)
{
    // будет идти по блокам и заполнять в каждом блоке
    for (int blocked_row = thread_number; blocked_row < blocks_count; blocked_row += threads_count)
    {
        double* block = coeff_matrix + blocked_row * block_size * matrix_size;
        double* rhs_chunk = right_hand_side + blocked_row * block_size;
        int height = std::min(block_size, matrix_size - blocked_row * block_size);
        for (int u = 0; u < height; u++)
            *(rhs_chunk + u) = 0.;
        for (int blocked_col = 0; blocked_col < blocks_count; blocked_col++)
        {
            int width = std::min(block_size, matrix_size - blocked_col * block_size);

            for (int q = 0; q < height; q++)
            {
                double sum = 0;
                for (int p = 0; p < width; p++)
                {
                    int j = blocked_col * block_size + p + 1;
                    if (j % 2 == 1)
                    {
                        double value = *block;
              
                        if (!std::isfinite(value))
                            return 0;
 
                        if ((sum > 0. && value > 0.) || (sum < 0. && value < 0.))
                        {
                            double abs_sum   = std::fabs(sum);
                            double abs_value = std::fabs(value);
                            if (std::min(abs_sum, abs_value) > DBL_MAX - std::max(abs_sum, abs_value))
                                return 0;
                        }

                        sum += value;
                    }
                    block++;
                }
                *(rhs_chunk + q) += sum;
            }
        }
    }
    return 1;
}

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
    int full_blocks_count = args->full_blocks_count;
    int remainder_size   = args->remainder_size;
    int blocks_count      = args->blocks_count;
    double* original_matrix_norm = args->original_matrix_norm;

    int threads_count = args->threads_count;
    int thread_number = args->thread_number;

    pthread_barrier_t* barrier = args->barrier;

    // std::fprintf(stdout, "hello from %d thread num out of %d amount\n", thread_number, threads_count);

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
    // идем по блочным строкам этого потока, включая неполную последнюю строку
    // pthread_barrier_wait(barrier);

    for (int blocked_row = thread_number; blocked_row < blocks_count; blocked_row += threads_count)
    {
        // std::fprintf(stderr, "thread_number: %d, blocked_row: %d, blocks_count: %d\n", thread_number, blocked_row, blocks_count);
        int height = std::min(block_size, matrix_size - blocked_row * block_size);
        std::memset(coeff_matrix + blocked_row * block_size * matrix_size, thread_number + 1, height * matrix_size * sizeof(double));
        std::memset(right_hand_side + blocked_row * block_size, 0, height * sizeof(double));
        std::memset(solution + blocked_row * block_size, 0, height * sizeof(double));
    }

    if (thread_number == 0)
        for (int i = 0; i < full_blocks_count; i++)
            block_col_perm[i] = i;
        
    pthread_barrier_wait(barrier);
    // выделяем локальную память с проверкой
    args->pivot_block_row = new double[matrix_size * block_size];
    args->block1          = new double[block_size * block_size];
    args->block2          = new double[block_size * block_size];
    args->block3          = new double[block_size * block_size];
    args->perm            = new int[block_size];
    args->buffer          = new double[block_size];
    args->mem_alloc = args->pivot_block_row && args->block1 && args->block2
        && args->block3 && args->perm && args->buffer;
    
    // барьер
    pthread_barrier_wait(barrier);
    
    // основной поток проверяет все ли выделено, если нет, то нужно выйти
    Args* all_args = args - thread_number;
    if (thread_number == 0)
    {
        for (int i = 1; i < threads_count; i++)
            if (!all_args[i].mem_alloc)
                all_args[0].mem_alloc = 0;

        if (!all_args[0].mem_alloc)
            std::fprintf(stderr, "error: local memory allocation\n");
    }

    // все потоки должны получить результат проверки до продолжения работы
    pthread_barrier_wait(barrier);
    if (!all_args[0].mem_alloc)
        return 0;

    pthread_barrier_wait(barrier);
    // инициализация матрицы
    // todo: добавить флаг корректной инициализации и обработать ошибки (в особенности чтения из файла)
    if (args->s == 0)
    {
        if (thread_number == 0)
            init_coeff_matrix_from_file(coeff_matrix, matrix_size, block_size, args->file_name);
    }
    
    if (args->s != 0)
        init_coeff_matrix_from_formula(coeff_matrix, matrix_size, block_size, blocks_count, args->s, threads_count, thread_number);

    pthread_barrier_wait(barrier);

    // инициализация правой части
    // todo: проверка на корректность
    init_right_hand_size(coeff_matrix, right_hand_side, matrix_size, block_size, blocks_count, threads_count, thread_number);

    if (thread_number == 0)
    {
        print_matrix(matrix_size, matrix_size, block_size, coeff_matrix, args->r, "Coefficient Matrix:\n");
        print_matrix(1, matrix_size, block_size, right_hand_side, args->r, "Right Hand Side:\n");
    }
    
    pthread_barrier_wait(barrier);


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
        full_blocks_count = 0,
        remainder_size = 0,
        blocks_count = 0;
    double t1 = 0., t2 = 0., r1 = -1., r2 = -1.;
    char* file_name = nullptr;
    // bool solved = 0;
    // общие данные системы - глобальная память
    double* coeff_matrix    = 0;
    double* right_hand_side = 0;
    double* solution        = 0;
    int*    block_col_perm      = 0; // перестановки матричных столбцов
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
    full_blocks_count = matrix_size / block_size;
    remainder_size   = matrix_size % block_size;
    blocks_count      = full_blocks_count + (remainder_size != 0);

    // if (threads_count > blocks_count)
    //     threads_count = blocks_count;

    // выделение памяти (нужно только глобальную)
    {
        coeff_matrix    = new double[matrix_size * matrix_size];
        right_hand_side = new double[matrix_size];
        solution        = new double[matrix_size];
        block_col_perm  = new int[full_blocks_count];
        pivot_block_inverse = new double[block_size * block_size];

        if (!coeff_matrix || !right_hand_side || !solution || !block_col_perm || !pivot_block_inverse)
        {
            if (coeff_matrix)
                delete[] coeff_matrix;
            if (right_hand_side)
                delete[] right_hand_side;
            if (solution)
                delete[] solution;
            if(block_col_perm)
                delete[] block_col_perm;
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
            args[k].block_col_perm  = block_col_perm;
            args[k].pivot_block_inverse = pivot_block_inverse;
            
            args[k].matrix_size = matrix_size;
            args[k].block_size = block_size;
            args[k].full_blocks_count = full_blocks_count;
            args[k].remainder_size = remainder_size;
            args[k].blocks_count = blocks_count;

            args[k].threads_count = threads_count;
            args[k].thread_number = k;

            args[k].r         = r;
            args[k].s         = s;
            args[k].file_name = file_name;

            args[k].barrier = barrier;
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

    int status = args[0].mem_alloc ? 0 : 1;

    // конечный вывод
    if (status == 0)
        std::printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], task, r1, r2, t1, t2, s, matrix_size, block_size, threads_count);

    // освобождаем память после завершения всех потоков, в том числе при ошибке выделения
    if (args)
    {
        for (int i = 0; i < threads_count; i++)
        {
            if (args[i].pivot_block_row)
                delete[] args[i].pivot_block_row;
            if (args[i].block1)
                delete[] args[i].block1;
            if (args[i].block2)
                delete[] args[i].block2;
            if (args[i].block3)
                delete[] args[i].block3;
            if (args[i].perm)
                delete[] args[i].perm;
            if (args[i].buffer)
                delete[] args[i].buffer;
        }
    }
    if (barrier)
    {
        pthread_barrier_destroy(barrier);
        delete barrier;
    }
    if (pid)
        delete[] pid;
    if (args)
        delete[] args;
    if (coeff_matrix)
        delete[] coeff_matrix;
    if (right_hand_side)
        delete[] right_hand_side;
    if (solution)
        delete[] solution;
    if (block_col_perm)
        delete[] block_col_perm;
    if (pivot_block_inverse)
        delete[] pivot_block_inverse;

    return 0;
}
