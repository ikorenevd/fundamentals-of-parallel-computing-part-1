#ifndef MEM_STRUCT_H
#define MEM_STRUCT_H

#include "./headers_and_defines.h"

class SharedMemory
{
    public:
        double solve_time = 0.; // время решения
        double residual_time = 0.; // время подсчета невязок

        // параметры системы
        int matrix_size = 0;
        int block_size = 0;
        int remainder_size = 0;
        int full_block_count = 0;
        int block_count = 0;

        // параметры запуска
        int thread_count = 0;
        int matrix_formula = 0;
        int print_limit = 0;
        const char* file_name = 0;

        // сама система
        double* coeff_matrix = 0;
        double* right_hand_side = 0;
        double* solution = 0;

        // перестановки блочных столбцов для восстановления
        int* block_column_order = 0;

        // данные связанные с поиском главного блока и какой поток его содержит
        int    pivot_block_row  = -1;
        int    pivot_block_column  = -1;

        // обратная матрица для главного блока
        double* inverse_pivot_block  = 0;

        // барьер
        pthread_barrier_t* barrier = 0;
        pthread_mutex_t launch_mutex = PTHREAD_MUTEX_INITIALIZER;
        bool launch_succeeded = 0;

        // информация, которая считается до запуска и не изменяется во время работы
        double original_matrix_norm = 0.; // todo: будет считаться во время инициализации правой части
    public:
        // запуск
        void set(int matrix_order, int block_order, int number_of_threads, int output_limit, int formula, const char* input_file_name);
        bool create_and_init_barrier();

        // работа с памятью
        bool allocate();
        void free();

        // по православному
        SharedMemory() = default;
        SharedMemory(const SharedMemory&) = delete;
        SharedMemory(SharedMemory&&) = delete;
        SharedMemory& operator=(const SharedMemory&) = delete;
        SharedMemory& operator=(SharedMemory&&) = delete;

        ~SharedMemory();
};

class ThreadMemory
{
    public:
        int thread_index = 0;
        double* pivot_row_tail = 0;
        double candidate_norm = -1.;
        int    candidate_row = -1;
        int    candidate_column = -1;
        double* inverse_candidate = 0;
        double* working_block = 0;
        double* best_inverse_block = 0;
        int* column_order = 0; // перестановки для поиска обратного блока
        double* vector_product = 0; // нужен для умножения матрицы на вектор (в обратном ходе)

        double residual = 0.;
        double rhs_norm = 0.;
        double solution_error = 0.;

        bool succeeded = 0;

        bool allocate(int matrix_size, int block_size);
        void free();

        ThreadMemory() = default;
        ThreadMemory(const ThreadMemory&) = delete;
        ThreadMemory(ThreadMemory&&) = delete;
        ThreadMemory& operator=(const ThreadMemory&) = delete;
        ThreadMemory& operator=(ThreadMemory&&) = delete;

        ~ThreadMemory();
};


#endif
