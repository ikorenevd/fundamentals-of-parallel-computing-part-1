#include "./linear_system_solver.h"
#include "./block_gauss.h"
#include "./block_matrix.h"
#include "./headers_and_defines.h"
#include "./matrix_products.h"
#include "./pivot_search.h"
#include "./thread_args.h"

static void normalize_pivot_row(ThreadArgs* arg, int alpha)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int block_count = arg->shared->block_count;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    double* working_block = arg->local.working_block;
    double* best_inverse_block = arg->local.best_inverse_block;

    // копируем обратный главный блок в собственную память
    std::memcpy(working_block, arg->shared->inverse_pivot_block, block_size * block_size * sizeof(double));

    // на место главного блока ставим единичную матрицу
    double* pivot_block = get_block_address(coeff_matrix, matrix_size, block_size, alpha, alpha);
    std::memset(pivot_block, 0, block_size * block_size * sizeof(double));
    for (int i = 0; i < block_size; i++)
        pivot_block[i * block_size + i] = 1.;

    // обновляем все блоки справа от него
    for (int j = alpha + 1; j < block_count; j++)
    {
        matrix_multiplication(
            working_block,
            get_block_address(coeff_matrix, matrix_size, block_size, alpha, j),
            best_inverse_block,
            block_size,
            block_size,
            std::min(block_size, matrix_size - j * block_size)
        );
        set_block(coeff_matrix, matrix_size, block_size, alpha, j, best_inverse_block);
    }

    // чанк правой части
    {
        matrix_multiplication(working_block, right_hand_side + alpha * block_size, best_inverse_block, block_size, block_size, 1);
        std::memcpy(right_hand_side + alpha * block_size, best_inverse_block, block_size * sizeof(double));
    }

    // Единичный диагональный блок задаёт нижнюю границу её масштаба.
    {
        double* pivot_row_tail = get_block_address(coeff_matrix, matrix_size, block_size, alpha, alpha);
        int pivot_row_size = block_size * (matrix_size - alpha * block_size);
        double pivot_row_max = 0.;

        for (int t = block_size * block_size; t < pivot_row_size; t++) // ищем максимальный по модулю элемент в оставшейся строке
            pivot_row_max = std::max(pivot_row_max, std::fabs(pivot_row_tail[t]));

        double pivot_row_eps = (MACHINE_EPS * MACHINE_EPS) * pivot_row_max;

        for (int t = block_size * block_size; t < pivot_row_size; t++)
            if (std::fabs(pivot_row_tail[t]) < pivot_row_eps)
                pivot_row_tail[t] = 0.;
    }
}

static void copy_pivot_row(ThreadArgs* arg, int alpha, int first)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    double* pivot_row_tail = arg->local.pivot_row_tail;
    double* working_block = arg->local.working_block;

    std::memcpy(
        pivot_row_tail,
        get_block_address(coeff_matrix, matrix_size, block_size, alpha, first),
        (size_t)block_size * (matrix_size - first * block_size) * sizeof(double)
    );
    std::memcpy(working_block, right_hand_side + alpha * block_size, block_size * sizeof(double));
}

static void eliminate_block_row(ThreadArgs* arg, int alpha, int block_row, double eps)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int block_count = arg->shared->block_count;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    double* pivot_row_tail = arg->local.pivot_row_tail;
    double* working_block = arg->local.working_block;
    double* vector_product = arg->local.vector_product;
    int first = alpha + 1;

    int height = std::min(block_size, matrix_size - block_row * block_size);

    double* multiplier = get_block_address(coeff_matrix, matrix_size, block_size, block_row, alpha);

    for (int entry = 0; entry < height * block_size; entry++)
        if (std::fabs(multiplier[entry]) < eps)
            multiplier[entry] = 0.;

    for (int block_column = alpha + 1; block_column < block_count; block_column++)
    {
        matrix_multiplication_subtract(
            multiplier,
            pivot_row_tail + (block_column - first) * block_size * block_size,
            get_block_address(coeff_matrix, matrix_size, block_size, block_row, block_column),
            height,
            block_size,
            std::min(block_size, matrix_size - block_column * block_size)
        );
    }

    matrix_multiplication(
        multiplier,
        working_block,
        // right_hand_side + alpha * block_size,
        vector_product,
        height,
        block_size,
        1
    );

    for (int entry = 0; entry < height; entry++)
        right_hand_side[block_row * block_size + entry] -= vector_product[entry];

    // std::memset(multiplier, 0, sizeof(double) * height * block_size); //  можно не занулять, так как эти данные трогать не будем
}

bool forward_elimination(ThreadArgs* arg)
{
    const int full_block_count = arg->shared->full_block_count;
    const int block_count = arg->shared->block_count;
    const int thread_count = arg->shared->thread_count;
    pthread_barrier_t* barrier = arg->shared->barrier;
    const int thread_index = arg->local.thread_index;
    const double original_matrix_norm = arg->shared->original_matrix_norm;
    double eps = (MACHINE_EPS * MACHINE_EPS) * original_matrix_norm;

    // прямой ход
    // для квадратных блоков
    for (int alpha = 0; alpha < full_block_count; alpha++)
    {
        // каждый поток считаем ищет кандидата на главных блок
        arg->local.succeeded = find_block_pivot(arg, alpha);

        // синхранизируемся после него
        pthread_barrier_wait(barrier);

        // нулевой поток выбирает главных блок
        if (thread_index == 0)
            select_pivot(arg, alpha);

        pthread_barrier_wait(barrier);

        // не нашли
        if (!arg->local.succeeded)
            return 0;

        // меняем строки (зависит от потока)
        swap_pivot_columns(arg, alpha, arg->shared->pivot_block_column);

        // так как строки и лок, столбцы уже поменяны, то можно уже начать обновлять строку
        if (thread_index == (alpha % thread_count))
            normalize_pivot_row(arg, alpha);

        pthread_barrier_wait(barrier);

        int first = alpha + 1;
        int start_row = first + (thread_index - first % thread_count + thread_count) % thread_count;

        // копируем только блоки справа от главного и только если есть свои нижние строки
        if (start_row < block_count)
            copy_pivot_row(arg, alpha, first);

        pthread_barrier_wait(barrier);

        // обновляем строки
        for (int block_row = start_row; block_row < block_count; block_row += thread_count)
            eliminate_block_row(arg, alpha, block_row, eps);
    }

    return 1;
}

bool solve_remainder_block(ThreadArgs* arg)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int remainder_size = arg->shared->remainder_size;
    const int full_block_count = arg->shared->full_block_count;
    const int thread_count = arg->shared->thread_count;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    pthread_barrier_t* barrier = arg->shared->barrier;
    const int thread_index = arg->local.thread_index;
    double* best_inverse_block = arg->local.best_inverse_block;
    int* column_order = arg->local.column_order;
    const double original_matrix_norm = arg->shared->original_matrix_norm;

    pthread_barrier_wait(barrier);

    arg->local.succeeded = 1;
    // если есть неквадратные блоки
    if (remainder_size != 0 && (full_block_count % thread_count == thread_index))
    {
        std::memcpy(best_inverse_block, right_hand_side + full_block_count * block_size, sizeof(double) * remainder_size);
        double* diagonal_block = get_block_address(coeff_matrix, matrix_size, block_size, full_block_count, full_block_count);
        arg->local.succeeded = naive_gauss_solve_block(diagonal_block, remainder_size, best_inverse_block, column_order, remainder_size * original_matrix_norm * MACHINE_EPS);
        std::memcpy(right_hand_side + full_block_count * block_size, best_inverse_block, sizeof(double) * remainder_size);
    }

    pthread_barrier_wait(barrier);

    for (int other_thread_index = 0; other_thread_index < thread_count; other_thread_index++)
        if (!(arg - thread_index + other_thread_index)->local.succeeded)
            return 0;

    return 1;
}

void backward_substitution(ThreadArgs* arg)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int full_block_count = arg->shared->full_block_count;
    const int block_count = arg->shared->block_count;
    const int thread_count = arg->shared->thread_count;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    double* solution = arg->shared->solution;
    int* block_column_order = arg->shared->block_column_order;
    pthread_barrier_t* barrier = arg->shared->barrier;
    const int thread_index = arg->local.thread_index;
    double* working_block = arg->local.working_block;
    double* vector_product = arg->local.vector_product;

    // обратный ход
    pthread_barrier_wait(barrier);


    // записываем последний блок решения
    int last_row = block_count - 1;
    int last_height = std::min(block_size, matrix_size - last_row * block_size);

    if (thread_index == last_row % thread_count)
        std::memcpy(
            solution + last_row * block_size,
            right_hand_side + last_row * block_size,
            last_height * sizeof(double)
        );

    pthread_barrier_wait(barrier);

    for (int alpha = last_row; alpha > 0; alpha--)
    {
        int width = std::min(block_size, matrix_size - alpha * block_size);

        // копируем готовый блок решения в собственную память
        std::memcpy(
            working_block,
            solution + alpha * block_size,
            width * sizeof(double)
        );

        // каждый поток обновляет свои строки выше alpha
        for (int i = thread_index; i < alpha; i += thread_count)
        {
            matrix_multiplication(
                get_block_address(coeff_matrix, matrix_size, block_size, i, alpha),
                working_block,
                vector_product,
                block_size,
                width,
                1
            );

            for (int t = 0; t < block_size; t++)
                right_hand_side[i * block_size + t] -= vector_product[t];
        }

        // владелец следующей строки записывает готовый блок решения
        if (thread_index == (alpha - 1) % thread_count)
            std::memcpy(
                solution + (alpha - 1) * block_size,
                right_hand_side + (alpha - 1) * block_size,
                block_size * sizeof(double)
            );

        pthread_barrier_wait(barrier);
    }

    // восстанавливаем исходный порядок неизвестных
    if (thread_index == 0)
        for (int i = 0; i < full_block_count; i++)
            std::memcpy(solution + block_column_order[i] * block_size, right_hand_side + i * block_size, block_size * sizeof(double));

    pthread_barrier_wait(barrier);
}
