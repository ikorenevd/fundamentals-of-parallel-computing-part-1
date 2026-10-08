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
    for (int row = 0; row < n; row++)
    {
        int block_row = row / m;
        int height    = std::min(m, n - block_row * m);
        int local_row = row % m;
        double ax     = 0.;

        for (int col = 0, block_col = 0; col < n; block_col++)
        {
            int width  = std::min(m, n - col);
            int offset = block_row * m * n + block_col * m * height + local_row * width;
            for (int j = 0; j < width; j++)
                ax += A[offset + j] * solution[col + j];

            col += width;
        }

        residual += std::fabs(ax - b[row]);
        norm_b   += std::fabs(b[row]);
        error    += std::fabs(solution[row] - (row % 2 == 0 ? 1. : 0.));
    }

    r1 = norm_b > 0. ? residual / norm_b : (residual <= 0. ? 0. : HUGE_VAL);
    r2 = error / (n / 2 + n % 2);
}
