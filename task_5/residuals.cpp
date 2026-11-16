#include "./residuals.h"
#include "./headers_and_defines.h"
#include "./thread_args.h"

void compute_residuals(ThreadArgs* arg)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int block_count = arg->shared->block_count;
    const int thread_count = arg->shared->thread_count;

    const int thread_index = arg->local.thread_index;

    const double* coeff_matrix = arg->shared->coeff_matrix;
    const double* right_hand_side = arg->shared->right_hand_side;
    const double* solution = arg->shared->solution;

    double residual = 0.;
    double rhs_norm = 0.;
    double solution_error = 0.;

    for (int block_row = thread_index; block_row < block_count; block_row += thread_count)
    {
        int height = std::min(block_size, matrix_size - block_row * block_size);

        for (int local_row = 0; local_row < height; local_row++)
        {
            int row = block_row * block_size + local_row;
            double product = 0.;

            for (int block_column = 0; block_column < block_count; block_column++)
            {
                int column = block_column * block_size;
                int width = std::min(block_size, matrix_size - column);

                int offset = block_row * block_size * matrix_size
                    + block_column * block_size * height
                    + local_row * width;

                for (int entry = 0; entry < width; entry++)
                    product += coeff_matrix[offset + entry] * solution[column + entry];
            }

            residual += std::fabs(product - right_hand_side[row]);
            rhs_norm += std::fabs(right_hand_side[row]);
            solution_error += std::fabs(solution[row] - (row % 2 == 0 ? 1. : 0.));
        }
    }

    arg->local.residual = residual;
    arg->local.rhs_norm = rhs_norm;
    arg->local.solution_error = solution_error;
}
