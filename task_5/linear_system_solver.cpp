#include "./linear_system_solver.h"

#include "./helper.h"
#include "./headers_and_defines.h"
#include "./Args.h"

bool find_block_pivot(Args* arg, int alpha)
{
    double* coeff_matrix = arg->gen->coeff_matrix;

    int matrix_size = arg->gen->matrix_size;
    int block_size = arg->gen->block_size;
    int full_blocks_count = arg->gen->full_block_count;

    int threads_count = arg->gen->threads_count;
    int thread_number = arg->own.thread_number;

    double* current_block = arg->own.current_block;
    double* inverse_candidate = arg->own.inverse_pivot_block_candidate;
    int* perm = arg->own.perm;
    double* inverse_best = arg->own.best_inverse_pivot_block;

    int best_i = -1;
    int best_j = -1;
    double best_norm = HUGE_VAL;
    double pivot_eps = (matrix_size * MACHINE_EPS) * arg->gen->original_matrix_norm;

    int start_row = alpha + (thread_number - alpha % threads_count + threads_count) % threads_count;
    for (int i = start_row; i < full_blocks_count; i += threads_count)
    {
        for (int j = alpha; j < full_blocks_count; j++)
        {
            get_block(coeff_matrix, matrix_size, block_size, i, j, current_block);
            if (!invert_block(current_block, inverse_candidate, block_size, perm, pivot_eps))
                continue;

            double inverse_norm = calculate_block_norm(inverse_candidate, block_size, 0);
            if (!std::isfinite(inverse_norm))
                continue;

            // Первый подходящий блок задаёт начальную норму для сравнения.
            if (best_i == -1 || inverse_norm < best_norm)
            {
                best_i = i;
                best_j = j;
                best_norm = inverse_norm;
                std::swap(inverse_best, inverse_candidate);
            }
        }
    }

    if (best_i == -1)
        return 0;

    double eps = MACHINE_EPS * MACHINE_EPS * (best_norm);
    for (int t = 0; t < block_size * block_size; t++)
    {
        double value = inverse_best[t];
        arg->own.best_inverse_pivot_block[t] = std::fabs(value) < eps ? 0. : value;
    }

    arg->own.pivot_block_candidate_row = best_i;
    arg->own.pivot_block_candidate_col = best_j;
    arg->own.pivot_block_candidate_norm = best_norm;

    return 1;
}

void swap_block_rows(Args* arg, int i, int j)
{
    if (i == j)
        return;
    
    double* coeff_matrix = arg->gen->coeff_matrix;
    int threads_count = arg->gen->threads_count;
    int thread_number = arg->own.thread_number;
    int matrix_size = arg->gen->matrix_size;
    int block_size = arg->gen->block_size;
    int block_count = arg->gen->block_count;

    // меняем блоки в матрице коэффициентов
    for (int row = thread_number; row < block_count; row += threads_count)
    {
        int height = std::min(block_size, matrix_size - row * block_size);
        double* tail = coeff_matrix + row * matrix_size * block_size;
        std::swap_ranges(
            tail + i * height * block_size,
            tail + (i + 1) * height * block_size,
            tail + j * height * block_size
        );
    }

    // меняем перестановки
    if (thread_number == 0)
        std::swap(arg->gen->block_col_perm[i], arg->gen->block_col_perm[j]);
}