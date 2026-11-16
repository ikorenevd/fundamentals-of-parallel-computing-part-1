#include "./pivot_search.h"
#include "./block_gauss.h"
#include "./block_matrix.h"
#include "./headers_and_defines.h"
#include "./thread_args.h"

bool find_block_pivot(ThreadArgs* arg, int alpha)
{
    double* coeff_matrix = arg->shared->coeff_matrix;

    int matrix_size = arg->shared->matrix_size;
    int block_size = arg->shared->block_size;
    int full_block_count = arg->shared->full_block_count;

    int thread_count = arg->shared->thread_count;
    int thread_index = arg->local.thread_index;

    double* working_block = arg->local.working_block;
    double* inverse_candidate = arg->local.inverse_candidate;
    int* column_order = arg->local.column_order;
    double* best_inverse = arg->local.best_inverse_block;

    int best_row = -1;
    int best_column = -1;
    double best_norm = HUGE_VAL;
    double pivot_eps = (matrix_size * MACHINE_EPS) * arg->shared->original_matrix_norm;

    int start_row = alpha + (thread_index - alpha % thread_count + thread_count) % thread_count;
    for (int block_row = start_row; block_row < full_block_count; block_row += thread_count)
    {
        for (int block_column = alpha; block_column < full_block_count; block_column++)
        {
            get_block(coeff_matrix, matrix_size, block_size, block_row, block_column, working_block);
            if (!invert_block(working_block, inverse_candidate, block_size, column_order, pivot_eps))
                continue;

            double inverse_norm = calculate_block_norm(inverse_candidate, block_size, 0);
            if (!std::isfinite(inverse_norm))
                continue;

            // Первый подходящий блок задаёт начальную норму для сравнения.
            if (best_row == -1 || inverse_norm < best_norm)
            {
                best_row = block_row;
                best_column = block_column;
                best_norm = inverse_norm;
                std::swap(best_inverse, inverse_candidate);
            }
        }
    }

    if (best_row == -1)
        return 0;

    double eps = MACHINE_EPS * MACHINE_EPS * (best_norm);
    for (int entry = 0; entry < block_size * block_size; entry++)
    {
        double value = best_inverse[entry];
        arg->local.best_inverse_block[entry] = std::fabs(value) < eps ? 0. : value;
    }

    arg->local.candidate_row = best_row;
    arg->local.candidate_column = best_column;
    arg->local.candidate_norm = best_norm;

    return 1;
}

void select_pivot(ThreadArgs* arg, int alpha)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int thread_count = arg->shared->thread_count;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;

    double best_norm = -1.;
    int best_thread = -1;
    // ищем первый поток с найденным кандидатом
    for (int thread_index = 0; thread_index < thread_count; thread_index++)
        if ((arg + thread_index)->local.succeeded)
        {
            best_thread = thread_index;
            break;
        }

    // если нашли
    if (best_thread != -1)
    {
        best_norm = (arg + best_thread)->local.candidate_norm;

        for (int thread_index = best_thread + 1; thread_index < thread_count; thread_index++)
            if ((arg + thread_index)->local.succeeded &&
                best_norm > (arg + thread_index)->local.candidate_norm)
            {
                best_thread = thread_index;
                best_norm = (arg + thread_index)->local.candidate_norm;
            }

        for (int thread_index = 0; thread_index < thread_count; thread_index++)
            (arg + thread_index)->local.succeeded = 1;

        arg->shared->pivot_block_row = (arg + best_thread)->local.candidate_row;
        arg->shared->pivot_block_column = (arg + best_thread)->local.candidate_column;

        std::memcpy(arg->shared->inverse_pivot_block, (arg + best_thread)->local.best_inverse_block, block_size * block_size * sizeof(double));

        // сразу меняем блочные строки
        if (alpha != arg->shared->pivot_block_row)
        {
            std::swap_ranges(
                coeff_matrix + alpha * block_size * matrix_size,
                coeff_matrix + (alpha + 1) * block_size * matrix_size,
                coeff_matrix + arg->shared->pivot_block_row * block_size * matrix_size);
            std::swap_ranges(
                right_hand_side + alpha * block_size,
                right_hand_side + (alpha + 1) * block_size,
                right_hand_side + arg->shared->pivot_block_row * block_size
            );
        }
    }
    // не нашли, зануляем флаг, чтобы выйти после синхронизации
    else
    {
        for (int thread_index = 0; thread_index < thread_count; thread_index++)
            (arg + thread_index)->local.succeeded = 0;
    }
}

void swap_pivot_columns(ThreadArgs* arg, int alpha, int pivot_column)
{
    if (alpha == pivot_column)
        return;

    double* coeff_matrix = arg->shared->coeff_matrix;
    int thread_count = arg->shared->thread_count;
    int thread_index = arg->local.thread_index;
    int matrix_size = arg->shared->matrix_size;
    int block_size = arg->shared->block_size;
    int block_count = arg->shared->block_count;

    // меняем блоки в матрице коэффициентов
    for (int block_row = thread_index; block_row < block_count; block_row += thread_count)
    {
        int height = std::min(block_size, matrix_size - block_row * block_size);
        double* tail = coeff_matrix + block_row * matrix_size * block_size;
        std::swap_ranges(
            tail + alpha * height * block_size,
            tail + (alpha + 1) * height * block_size,
            tail + pivot_column * height * block_size
        );
    }

    // меняем перестановки
    if (thread_index == 0)
        std::swap(arg->shared->block_column_order[alpha], arg->shared->block_column_order[pivot_column]);
}
