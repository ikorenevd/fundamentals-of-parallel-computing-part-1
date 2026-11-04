#ifndef MEM_STRUCT_H
#define MEM_STRUCT_H

#include "./headers_and_defines.h"

// todo: использовать умные указатели на GENERAL_MEMORY
class GENERAL_MEMORY
{
    public:
        bool solved = 0;
        double t1 = 0.; // время решения
        double t2 = 0.; // время подсчета невязок

        // параметры системы
        int matrix_size = 0;
        int block_size = 0;
        int remainder_size = 0;
        int full_block_count = 0;
        int block_count = 0;

        // параметры запуска
        int threads_count = 0;
        int s = 0;
        int r = 0;
        char* file_name = 0;

        // сама система
        double* coeff_matrix = 0;
        double* right_hand_side = 0;
        double* solution = 0;

        // перестановки блочных столбцов для восстановления
        int* block_col_perm = 0;

        // данные связанные с поиском главного блока и какой поток его содержит
        double pivot_block_norm = -1;
        int    pivot_block_row  = -1;
        int    pivot_block_col  = -1;
        int    pivot_thread_num = -1;

        // обратная матрица для главного блока
        double* inverse_pivot_block  = 0;

        // барьер
        pthread_barrier_t* barrier = 0;

        // информация, которая считается до запуска и не изменяется во время работы
        double original_matrix_norm = 0.; // todo: будет считаться во время инициализации правой части
    public:
        // запуск
        void set(int m_size, int b_size, int th_count, int _r, int _s, char* _file_name);
        bool create_and_init_barrier();

        // работа с памятью
        bool allocate();
        void free();

        // по православному
        GENERAL_MEMORY() = default;
        GENERAL_MEMORY(const GENERAL_MEMORY&) = delete;
        GENERAL_MEMORY(GENERAL_MEMORY&&) = delete;
        GENERAL_MEMORY& operator=(const GENERAL_MEMORY&) = delete;
        GENERAL_MEMORY& operator=(GENERAL_MEMORY&&) = delete;

        ~GENERAL_MEMORY();
};

class THREAD_MEMORY
{
    public:
        int thread_number = 0;
        double* pivot_row = 0;
        double pivot_block_candidate_norm = -1.;
        int    pivot_block_candidate_row = -1;
        int    pivot_block_candidate_col = -1;
        double* inverse_pivot_block_candidate = 0;
        double* current_block = 0;
        double* best_inverse_pivot_block = 0;
        int* perm = 0; // перестановки для поиска обратного блока
        double* buffer = 0; // нужен для умножения матрицы на вектор (в обратном ходе)

        double residual = 0.;
        double norm_b = 0.;
        double error = 0.;

        bool okay = 0;

        bool allocate(int matrix_size, int block_size);
        void free();

        THREAD_MEMORY() = default;
        THREAD_MEMORY(const THREAD_MEMORY&) = delete;
        THREAD_MEMORY(THREAD_MEMORY&&) = delete;
        THREAD_MEMORY& operator=(const THREAD_MEMORY&) = delete;
        THREAD_MEMORY& operator=(THREAD_MEMORY&&) = delete;
        
        ~THREAD_MEMORY();
};


#endif