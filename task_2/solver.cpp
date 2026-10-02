#include "./solver.h"
#include "./LinearSystem.h"

#include <math.h>

double calculate_matrix_norm(int n, int m, const double* A, double* workspace)
{
    if (n < 1 || m < 1 || A == nullptr || workspace == nullptr)
        return HUGE_VAL; // todo: по-человечески

    for (int col = 0; col < n; col++)
        workspace[col] = 0.;

    const double* block = A;
    if (m == 1)
    {
        // при m = 1 матрица хранится по строкам
        for (int row = 0; row < n; row++)
        {
            for (int col = 0; col < n; col++)
                workspace[col] += fabs(block[col]);

            block += n;
        }
    }
    else
    {
        for (int row = 0; row < n; )
        {
            int height = m < n - row ? m : n - row;

            for (int col = 0; col < n; )
            {
                int width = m < n - col ? m : n - col;

                for (int p = 0; p < height; p++)
                {
                    for (int q = 0; q < width; q++)
                        workspace[col + q] += fabs(block[q]);

                    block += width;
                }

                col += width;
            }

            row += height;
        }
    }

    double norm = 0.;

    for (int col = 0; col < n; col++)
    {
        if (!isfinite(workspace[col]))
            return workspace[col];

        if (workspace[col] > norm)
            norm = workspace[col];
    }

    return norm;
}

int solve(int n, int m, double* A, double* b, double* x, double* workspace)
{
    if (n < 1 || m < 1 || A == nullptr || b == nullptr
        || x == nullptr || workspace == nullptr)
        return 1;

    double matrix_norm = calculate_matrix_norm(n, m, A, workspace);

    if (!isfinite(matrix_norm))
        return 1;

    // todo: заменить заглушку на что-то осмысленное потом
    for (int i = 0; i < n; i++)
        x[i] = i % 2 == 0 ? 1. : 0.;

    return 0;
}

int LinearSystem::solve()
{
    return ::solve(n, m, A, b, x, norm_workspace);
}
