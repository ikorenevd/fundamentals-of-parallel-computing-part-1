#include "./residuals.h"

#include <cmath>
#include <algorithm>

// todo: заместо HUGE_VAL переписать с поправкой на вычиселия машинного эпсиолна
void LinearSystem::compute_residuals(double& r1, double& r2) const
{
    r1 = -1; r2 = -1;

    double residual = 0.;
    double norm_b   = 0.;
    double error    = 0.;

    // обход строк в блочном хранении
    for (int row = 0; row < matrix_size; row++)
    {
        int block_row = row / block_size;
        int height    = std::min(block_size, matrix_size - block_row * block_size);
        int local_row = row % block_size;
        double ax     = 0.;

        for (int col = 0, block_col = 0; col < matrix_size; block_col++)
        {
            int width  = std::min(block_size, matrix_size - col);
            int offset = block_row * block_size * matrix_size + block_col * block_size * height + local_row * width;
            for (int j = 0; j < width; j++)
                ax += coeff_matrix[offset + j] * solution[col + j];

            col += width;
        }

        residual += std::fabs(ax - right_hand_side[row]);
        norm_b   += std::fabs(right_hand_side[row]);
        error    += std::fabs(solution[row] - (row % 2 == 0 ? 1. : 0.));
    }

    r1 = norm_b > 0. ? residual / norm_b : (residual <= 0. ? 0. : HUGE_VAL);
    r2 = error / (matrix_size / 2 + matrix_size % 2);
}
