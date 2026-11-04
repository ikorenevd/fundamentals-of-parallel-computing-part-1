#include "./thread_function.h"

#include "./headers_and_defines.h"
#include "./helper.h"
#include "./matrix_io.h"
#include "./linear_system_init_functions.h"
#include "./linear_system_solver.h"
#include "./Args.h"

void* thread_function(void* ptr)
{
    Args* arg = (Args*)ptr;

    // global
    const int g_matrix_size = arg->gen->matrix_size;
    const int g_block_size = arg->gen->block_size;
    const int g_remainder_size = arg->gen->remainder_size;
    const int g_full_block_count = arg->gen->full_block_count;
    const int g_blocks_count = arg->gen->block_count;

    const int g_threads_count = arg->gen->threads_count;
    const int g_s = arg->gen->s;
    const int g_r = arg->gen->r;
    const char* g_file_name = arg->gen->file_name;

    double* g_coeff_matrix = arg->gen->coeff_matrix;
    double* g_right_hand_side = arg->gen->right_hand_side;
    double* g_solution = arg->gen->solution;
    
    int* g_block_col_perm = arg->gen->block_col_perm;

    pthread_barrier_t* g_barrier = arg->gen->barrier;
    
    // own
    const int o_thread_number = arg->own.thread_number;
    double* o_pivot_row = arg->own.pivot_row;
    double* o_inverse_pivot_block_candidate = arg->own.inverse_pivot_block_candidate;
    double* o_current_block = arg->own.current_block;
    double* o_best_inverse_pivot_block = arg->own.best_inverse_pivot_block;
    int* o_perm = arg->own.perm;
    double* o_buffer = arg->own.buffer;

    // закреплеяем тред за ядром
    {
        #ifdef __linux__
            cpu_set_t cpu;
            CPU_ZERO(&cpu);
            int n_cpus = get_nprocs();
            int cpu_id = n_cpus - 1 - (o_thread_number % n_cpus);
            CPU_SET(cpu_id, &cpu);
            pthread_t tid = pthread_self();
            pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
        #endif
    }
    
    // инициализируем свою память(ту часть памяти coeff_matrix, right_hand_side, solution, с которой будет работать этот поток)
    // и собственную память потока
    // идем по блочным строкам этого потока, включая неполную последнюю строку
    // pthread_barrier_wait(barrier);
    {
        std::memset(o_pivot_row, 0, g_block_size * g_matrix_size * sizeof(double));
        std::memset(o_inverse_pivot_block_candidate, 0, g_block_size * g_block_size * sizeof(double));
        std::memset(o_current_block, 0, g_block_size * g_block_size * sizeof(double));
        std::memset(o_best_inverse_pivot_block, 0, g_block_size * g_block_size * sizeof(double));
        std::memset(o_perm, 0, g_block_size * sizeof(int));
        std::memset(o_buffer, 0, g_block_size * sizeof(double));

        for (int blocked_row = o_thread_number; blocked_row < g_blocks_count; blocked_row += g_threads_count)
        {
            int height = std::min(g_block_size, g_matrix_size - blocked_row * g_block_size);
            std::memset(g_coeff_matrix + blocked_row * g_block_size * g_matrix_size, 0, height * g_matrix_size * sizeof(double));
            std::memset(g_right_hand_side + blocked_row * g_block_size, 0, height * sizeof(double));
            std::memset(g_solution + blocked_row * g_block_size, 0, height * sizeof(double));
        }

        if (o_thread_number == 0)
            for (int i = 0; i < g_full_block_count; i++)
                g_block_col_perm[i] = i;
    }

    pthread_barrier_wait(g_barrier);

    // инициализация данных
    {
        // инициализация матрицы коэффициентов и проверка все ли хорошо
        // arg->own.okay = 0;
        if (g_s == 0)
        {
            if (o_thread_number == 0)
                arg->own.okay = init_coeff_matrix_from_file(g_coeff_matrix, g_matrix_size, g_block_size, g_file_name);
            else
                arg->own.okay = 1;
        }
        else
            arg->own.okay = init_coeff_matrix_from_formula(g_coeff_matrix, g_matrix_size, g_block_size, g_blocks_count, g_s, g_threads_count, o_thread_number);

        pthread_barrier_wait(g_barrier);

        if (o_thread_number == 0)
        {
            bool everything_okay = arg->own.okay;
            for (int thread_num = 1; thread_num < g_threads_count; thread_num++)
                everything_okay = everything_okay && (arg + thread_num)->own.okay;
            
            for (int thread_num = 0; thread_num < g_threads_count; thread_num++)
                (arg + thread_num)->own.okay = everything_okay;
        }

        pthread_barrier_wait(g_barrier);

        if (!arg->own.okay && o_thread_number == 0)
            std::fprintf(stderr, "Error: couldn't init coefficient matrix\n");

        if (!arg->own.okay)
            return 0;

        // инициализация правой части
        pthread_barrier_wait(g_barrier);

        arg->own.okay = init_right_hand_size(g_coeff_matrix, g_right_hand_side, g_matrix_size, g_block_size, g_blocks_count, g_threads_count, o_thread_number);
        
        pthread_barrier_wait(g_barrier);

        if (o_thread_number == 0)
        {
            bool everything_okay = arg->own.okay;
            for (int thread_num = 1; thread_num < g_threads_count; thread_num++)
                everything_okay = everything_okay && (arg + thread_num)->own.okay;
            
            for (int thread_num = 0; thread_num < g_threads_count; thread_num++)
                (arg + thread_num)->own.okay = everything_okay;

            arg->gen->original_matrix_norm = get_matrix_norm(g_coeff_matrix, g_matrix_size, g_block_size, g_blocks_count);
        }

        pthread_barrier_wait(g_barrier);

        if (!arg->own.okay && o_thread_number == 0)
            std::fprintf(stderr, "Error: couldn't init right hand side\n");

        if (!arg->own.okay)
            return 0;
    }

    pthread_barrier_wait(g_barrier);

    if (o_thread_number == 0)
    {
        print_matrix(g_matrix_size, g_matrix_size, g_block_size, g_coeff_matrix, g_r, "Coefficient matrix\n");
        print_matrix(1, g_matrix_size, g_block_size, g_right_hand_side, g_r, "Right hand side:\n");
    }

    double g_original_matrix_norm = arg->gen->original_matrix_norm;
    double eps = (MACHINE_EPS * MACHINE_EPS) * g_original_matrix_norm;

    if (o_thread_number == 0)
        arg->gen->t1 = get_full_time();

    pthread_barrier_wait(g_barrier);

    // прямой ход
    // для квадратных блоков
    for (int alpha = 0; alpha < g_full_block_count; alpha++)
    {
        // каждый поток считаем ищет кандидата на главных блок
        arg->own.okay = find_block_pivot(arg, alpha);

        // синхранизируемся после него
        pthread_barrier_wait(g_barrier);
        
        // нулевой поток выбирает главных блок
        if (o_thread_number == 0)
        {
            double best_norm = -1.;
            int best_thread = -1;
            // ищем первый поток с найденным кандидатом
            for (int thread_number = 0; thread_number < g_threads_count; thread_number++)
                if ((arg + thread_number)->own.okay)
                {
                    best_thread = thread_number;
                    break;
                }

            // если нашли
            if (best_thread != -1)
            {
                best_norm = (arg + best_thread)->own.pivot_block_candidate_norm;

                for (int thread_number = best_thread + 1; thread_number < g_threads_count; thread_number++)
                    if ((arg + thread_number)->own.okay && 
                        best_norm > (arg + thread_number)->own.pivot_block_candidate_norm)
                    {
                        best_thread = thread_number;
                        best_norm = (arg + thread_number)->own.pivot_block_candidate_norm;
                    }

                for (int thread_number = 0; thread_number < g_threads_count; thread_number++)
                    (arg + thread_number)->own.okay = 1;

                arg->gen->pivot_block_norm = best_norm;
                arg->gen->pivot_block_row = (arg + best_thread)->own.pivot_block_candidate_row;
                arg->gen->pivot_block_col = (arg + best_thread)->own.pivot_block_candidate_col;
                arg->gen->pivot_thread_num = best_thread;

                std::memcpy(arg->gen->inverse_pivot_block, (arg + best_thread)->own.best_inverse_pivot_block, g_block_size * g_block_size * sizeof(double));

                // сразу меняем блочные строки 
                if (alpha != arg->gen->pivot_block_row)
                {
                    std::swap_ranges(
                        g_coeff_matrix + alpha * g_block_size * g_matrix_size,
                        g_coeff_matrix + (alpha + 1) * g_block_size * g_matrix_size,
                        g_coeff_matrix + arg->gen->pivot_block_row * g_block_size * g_matrix_size);
                    std::swap_ranges(
                        g_right_hand_side + alpha * g_block_size,
                        g_right_hand_side + (alpha + 1) * g_block_size, 
                        g_right_hand_side + arg->gen->pivot_block_row * g_block_size
                    );
                }
            }
            // не нашли, зануляем флаг, чтобы выйти после синхронизации
            else
            {
                for (int thread_number = 0; thread_number < g_threads_count; thread_number++)
                    (arg + thread_number)->own.okay = 0;
            }
        }
    
        pthread_barrier_wait(g_barrier);

        // не нашли
        if (!arg->own.okay)
        {
            if (o_thread_number == 0)
            {
                arg->gen->t1 = get_full_time() - arg->gen->t1;
                std::fprintf(stderr, "Error: couldn't find pivot block\n");
            }
            return 0;
        }

        // меняем строки (зависит от потока)
        swap_block_rows(arg, alpha, arg->gen->pivot_block_col);

        // так как строки и лок, столбцы уже поменяны, то можно уже начать обновлять строку
        if (o_thread_number == (alpha % g_threads_count))
        {
            // копируем обратный главный блок в собственную память
            std::memcpy(o_current_block, arg->gen->inverse_pivot_block, g_block_size * g_block_size * sizeof(double));

            // на место главного блока ставим единичную матрицу
            double* pivot_block = get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, alpha, alpha);
            std::memset(pivot_block, 0, g_block_size * g_block_size * sizeof(double));
            for (int i = 0; i < g_block_size; i++)
                pivot_block[i * g_block_size + i] = 1.;

            // обновляем все блоки справа от него
            for (int j = alpha + 1; j < g_blocks_count; j++)
            {
                matrix_multiplication(
                    o_current_block,
                    get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, alpha, j),
                    o_best_inverse_pivot_block,
                    g_block_size,
                    g_block_size,
                    std::min(g_block_size, g_matrix_size - j * g_block_size)
                );
                set_block(g_coeff_matrix, g_matrix_size, g_block_size, alpha, j, o_best_inverse_pivot_block);
            }

            // чанк правой части
            {
                matrix_multiplication(o_current_block, g_right_hand_side + alpha * g_block_size, o_best_inverse_pivot_block, g_block_size, g_block_size, 1);
                std::memcpy(g_right_hand_side + alpha * g_block_size, o_best_inverse_pivot_block, g_block_size * sizeof(double));
            }
            
            // Единичный диагональный блок задаёт нижнюю границу её масштаба.
            {
                double* pivot_row = get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, alpha, alpha);
                int pivot_row_size = g_block_size * (g_matrix_size - alpha * g_block_size);
                double pivot_row_max = 0.;
                
                for (int t = g_block_size * g_block_size; t < pivot_row_size; t++) // ищем максимальный по модулю элемент в оставшейся строке
                    pivot_row_max = std::max(pivot_row_max, std::fabs(pivot_row[t]));
                
                double pivot_row_eps = (MACHINE_EPS * MACHINE_EPS) * pivot_row_max;
                
                for (int t = g_block_size * g_block_size; t < pivot_row_size; t++)
                    if (std::fabs(pivot_row[t]) < pivot_row_eps)
                        pivot_row[t] = 0.;
            }
        }

        pthread_barrier_wait(g_barrier);

        int first = alpha + 1;
        int start_row = first + (o_thread_number - first % g_threads_count + g_threads_count) % g_threads_count;

        // копируем только блоки справа от главного и только если есть свои нижние строки
        if (start_row < g_blocks_count)
        {
            std::memcpy(
                o_pivot_row,
                get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, alpha, first),
                (size_t)g_block_size * (g_matrix_size - first * g_block_size) * sizeof(double)
            );
            std::memcpy(o_current_block, g_right_hand_side + alpha * g_block_size, g_block_size * sizeof(double));
        }

        pthread_barrier_wait(g_barrier);

        // обновляем строки 
        for (int i = start_row; i < g_blocks_count; i += g_threads_count)
        {
            int height = std::min(g_block_size, g_matrix_size - i * g_block_size);

            double* i_alpha = get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, i, alpha);

            for (int t = 0; t < height * g_block_size; t++)
                if (std::fabs(i_alpha[t]) < eps)
                    i_alpha[t] = 0.;
            
            for (int j = alpha + 1; j < g_blocks_count; j++)
            {
                matrix_multiplication_subtract(
                    i_alpha,
                    o_pivot_row + (j - first) * g_block_size * g_block_size,
                    get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, i, j), 
                    height, 
                    g_block_size, 
                    std::min(g_block_size, g_matrix_size - j * g_block_size)
                );
            }

            matrix_multiplication(
                i_alpha,
                o_current_block,
                // g_right_hand_side + alpha * g_block_size,
                o_buffer,
                height,
                g_block_size,
                1
            );

            for (int t = 0; t < height; t++)
                g_right_hand_side[i * g_block_size + t] -= o_buffer[t];

            // std::memset(i_alpha, 0, sizeof(double) * height * g_block_size); //  можно не занулять, так как эти данные трогать не будем
        }
    }

    pthread_barrier_wait(g_barrier);

    arg->own.okay = 1;
    // если есть неквадратные блоки
    if (g_remainder_size != 0 && (g_full_block_count % g_threads_count == o_thread_number))
    {
        std::memcpy(o_best_inverse_pivot_block, g_right_hand_side + g_full_block_count * g_block_size, sizeof(double) * g_remainder_size);
        double* sm_coeff_matrix = get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, g_full_block_count, g_full_block_count);
        arg->own.okay = naive_gauss_solve_block(sm_coeff_matrix, g_remainder_size, o_best_inverse_pivot_block, o_perm, g_remainder_size * g_original_matrix_norm * MACHINE_EPS);
        std::memcpy(g_right_hand_side + g_full_block_count * g_block_size, o_best_inverse_pivot_block, sizeof(double) * g_remainder_size);
    }

    pthread_barrier_wait(g_barrier);

    for (int thread_num = 0; thread_num < g_threads_count; thread_num++)
        if (!(arg - o_thread_number + thread_num)->own.okay)
        {
            if (o_thread_number == 0)
            {
                arg->gen->t1 = get_full_time() - arg->gen->t1;
                std::fprintf(stdout, "Error: couldn't solve small system\n");
            }
            return 0;
        }

    // обратный ход
    pthread_barrier_wait(g_barrier);


    // записываем последний блок решения
    int last_row = g_blocks_count - 1;
    int last_height = std::min(g_block_size, g_matrix_size - last_row * g_block_size);

    if (o_thread_number == last_row % g_threads_count)
        std::memcpy(
            g_solution + last_row * g_block_size,
            g_right_hand_side + last_row * g_block_size,
            last_height * sizeof(double)
        );

    pthread_barrier_wait(g_barrier);

    for (int alpha = last_row; alpha > 0; alpha--)
    {
        int width = std::min(g_block_size, g_matrix_size - alpha * g_block_size);

        // копируем готовый блок решения в собственную память
        std::memcpy(
            o_current_block,
            g_solution + alpha * g_block_size,
            width * sizeof(double)
        );

        // каждый поток обновляет свои строки выше alpha
        for (int i = o_thread_number; i < alpha; i += g_threads_count)
        {
            matrix_multiplication(
                get_block_address(g_coeff_matrix, g_matrix_size, g_block_size, i, alpha),
                o_current_block,
                o_buffer,
                g_block_size,
                width,
                1
            );

            for (int t = 0; t < g_block_size; t++)
                g_right_hand_side[i * g_block_size + t] -= o_buffer[t];
        }

        // владелец следующей строки записывает готовый блок решения
        if (o_thread_number == (alpha - 1) % g_threads_count)
            std::memcpy(
                g_solution + (alpha - 1) * g_block_size,
                g_right_hand_side + (alpha - 1) * g_block_size,
                g_block_size * sizeof(double)
            );

        pthread_barrier_wait(g_barrier);
    }

    // восстанавливаем исходный порядок неизвестных
    if (o_thread_number == 0)
        for (int i = 0; i < g_full_block_count; i++)
            std::memcpy(g_solution + g_block_col_perm[i] * g_block_size, g_right_hand_side + i * g_block_size, g_block_size * sizeof(double));

    pthread_barrier_wait(g_barrier);

    if (o_thread_number == 0)
    {
        arg->gen->t1 = get_full_time() - arg->gen->t1;
        arg->gen->solved = 1;
    }

    pthread_barrier_wait(g_barrier);
    
    // реинициализация данных
    {
        // инициализация матрицы коэффициентов и проверка все ли хорошо
        if (g_s == 0)
        {
            if (o_thread_number == 0)
                arg->own.okay = init_coeff_matrix_from_file(g_coeff_matrix, g_matrix_size, g_block_size, g_file_name);
            else
                arg->own.okay = 1;
        }
        else
            arg->own.okay = init_coeff_matrix_from_formula(g_coeff_matrix, g_matrix_size, g_block_size, g_blocks_count, g_s, g_threads_count, o_thread_number);

        pthread_barrier_wait(g_barrier);

        if (o_thread_number == 0)
        {
            bool everything_okay = arg->own.okay;
            for (int thread_num = 1; thread_num < g_threads_count; thread_num++)
                everything_okay = everything_okay && (arg + thread_num)->own.okay;
            
            for (int thread_num = 0; thread_num < g_threads_count; thread_num++)
                (arg + thread_num)->own.okay = everything_okay;
        }

        pthread_barrier_wait(g_barrier);

        if (!arg->own.okay && o_thread_number == 0)
            std::fprintf(stderr, "Error: couldn't init coefficient matrix\n");

        if (!arg->own.okay)
            return 0;

        // инициализация правой части
        pthread_barrier_wait(g_barrier);

        arg->own.okay = init_right_hand_size(g_coeff_matrix, g_right_hand_side, g_matrix_size, g_block_size, g_blocks_count, g_threads_count, o_thread_number);
        
        pthread_barrier_wait(g_barrier);

        if (o_thread_number == 0)
        {
            bool everything_okay = arg->own.okay;
            for (int thread_num = 1; thread_num < g_threads_count; thread_num++)
                everything_okay = everything_okay && (arg + thread_num)->own.okay;
            
            for (int thread_num = 0; thread_num < g_threads_count; thread_num++)
                (arg + thread_num)->own.okay = everything_okay;
        }

        pthread_barrier_wait(g_barrier);

        if (!arg->own.okay && o_thread_number == 0)
            std::fprintf(stderr, "Error: couldn't init right hand side\n");

        if (!arg->own.okay)
            return 0;
    }

    if (o_thread_number == 0)
        arg->gen->t2 = get_full_time();

    pthread_barrier_wait(g_barrier);

    thread_compute_residuals(arg);

    pthread_barrier_wait(g_barrier);

    if (o_thread_number == 0)
        arg->gen->t2 = get_full_time() - arg->gen->t2;

    return 0;
}

void* thread_compute_residuals(void* ptr)
{
    Args* arg = (Args*)ptr;

    const int g_matrix_size = arg->gen->matrix_size;
    const int g_block_size = arg->gen->block_size;
    const int g_blocks_count = arg->gen->block_count;
    const int g_threads_count = arg->gen->threads_count;

    const int o_thread_number = arg->own.thread_number;

    const double* g_coeff_matrix = arg->gen->coeff_matrix;
    const double* g_right_hand_side = arg->gen->right_hand_side;
    const double* g_solution = arg->gen->solution;

    double residual = 0.;
    double norm_b = 0.;
    double error = 0.;

    for (int block_row = o_thread_number; block_row < g_blocks_count; block_row += g_threads_count)
    {
        int height = std::min(g_block_size, g_matrix_size - block_row * g_block_size);

        for (int local_row = 0; local_row < height; local_row++)
        {
            int row = block_row * g_block_size + local_row;
            double ax = 0.;

            for (int block_col = 0; block_col < g_blocks_count; block_col++)
            {
                int col = block_col * g_block_size;
                int width = std::min(g_block_size, g_matrix_size - col);

                int offset = block_row * g_block_size * g_matrix_size
                    + block_col * g_block_size * height
                    + local_row * width;

                for (int j = 0; j < width; j++)
                    ax += g_coeff_matrix[offset + j] * g_solution[col + j];
            }

            residual += std::fabs(ax - g_right_hand_side[row]);
            norm_b += std::fabs(g_right_hand_side[row]);
            error += std::fabs(g_solution[row] - (row % 2 == 0 ? 1. : 0.));
        }
    }

    arg->own.residual = residual;
    arg->own.norm_b = norm_b;
    arg->own.error = error;

    return 0;
}

// if (o_thread_number == 0)
// {
//     std::fprintf(stdout, "Alpha: %d\n", g_full_block_count + 1);
//     print_matrix(g_matrix_size, g_matrix_size, g_block_size, g_coeff_matrix, g_r, "Coefficient matrix\n");
//     print_matrix(1, g_matrix_size, g_block_size, g_right_hand_side, g_r, "Right hand side:\n");
// }
