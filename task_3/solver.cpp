#include "./solver.h"
#include "./block_ops.h"
#include "./LinearSystem.h"

#include <cmath>
#include <algorithm>
#include <cstring>

bool LinearSystem::finding_block_pivot(int alpha, int& pivot_i, int& pivot_j, double* ws_block1, double* ws_block2, double* ws_vector_norms, double* block_inverse)
{
    pivot_i = -1;
    pivot_j = -1;

    double pivot_norm = 0.;

    for (int i = alpha; i < k; i++)
    {
        for (int j = alpha; j < k; j++)
        {
            // ws_block1 хранит кандидата, ws_block2 — его обратную матрицу.
            get_block(i, j, ws_block1);
            if (!invert_block(ws_block1, ws_block2, m, w_perm))
                continue;

            double inverse_norm = calculate_block_norm(ws_block2, m, ws_vector_norms);
            if (!std::isfinite(inverse_norm))
                continue;

            // Первый подходящий блок задаёт начальную норму для сравнения.
            if (pivot_i == -1 || inverse_norm < pivot_norm)
            {
                pivot_i = i;
                pivot_j = j;
                pivot_norm = inverse_norm;
                std::memcpy(block_inverse, ws_block2, sizeof(double) * m * m);
            }
        }
    }

    if (pivot_i == -1)
        return 0;

    return 1;
}

bool LinearSystem::blocked_matrix_to_triangular()
{
    // используем block1 для обратного блока, block2 и block3 - для работы, block4 - тоже(в начале для вычисления норма, потом для свободных членов)
    for (int alpha = 0; alpha < k; alpha++)
    {
        int pivot_i, pivot_j;
        if (!finding_block_pivot(alpha, pivot_i, pivot_j, block2, block3, block4, block1))
            return 0;
        
        swap_blocked_rows(alpha, pivot_i, alpha);
        swap_blocked_columns(alpha, pivot_j);

        {
            // на место главного блока ставим единичный
            std::memset(block2, 0, sizeof(double) * m * m);
            for (int i = 0; i < m; i++)
                block2[i * m + i] = 1.;
            set_block(alpha, alpha, block2);

            // обновляем все блоки справа от него
            for (int j = alpha + 1; j < k + (l != 0); j++)
            {
                int width = std::min(m, n - j * m);
                get_block(alpha, j, block2);
                matrix_multiplication(block1, block2, block3, m, m, width);
                set_block(alpha, j, block3);
            }

            // столбец свободных коэффициентов
            matrix_multiplication(block1, b + m * alpha, block2, m, m, 1);
            std::memcpy(b + m * alpha, block2, sizeof(double) * m);
        }

        {
            // обновляем нижние блоки и свободные члены(все блоки и квадратные, и прямоугольные)
            for (int i = alpha + 1; i < k + (l != 0); i++)
            {
                int height = std::min(m, n - i * m);
                get_block(i, alpha, block1); // block1 = A_{i, \alpha}

                for (int j = alpha + 1; j < k + (l != 0); j++)
                {
                    int width = std::min(m, n - j * m);
                    get_block(alpha, j, block2); // block2 = A_{\alpha, j}
                    get_block(i, j, block4);
                    matrix_multiplication_subtract(block1, block2, block4, height, m, width);
                    set_block(i, j, block4);
                }

                matrix_multiplication(block1, b + alpha * m, block2, height, m, 1);

                for (int t = 0; t < height; t++)
                    b[i * m + t] -= block2[t];
                
                std::memset(block1, 0, sizeof(double) * height * m);
                set_block(i, alpha, block1);
            }
        }
    }

    // Последний блок решаем один раз после исключения всех полных блоков.
    if (l != 0)
    {
        get_block(k, k, block1);
        std::memcpy(block2, b + k * m, sizeof(double) * l);
        if (!naive_gauss_solve_block(block1, l, block2, w_perm))
            return 0;
        std::memcpy(b + k * m, block2, sizeof(double) * l);
    }

    return 1;
}

bool LinearSystem::triangular_blocked_to_solution()
{
    if (l != 0)
        std::memcpy(solution + m * k, b + m * k, sizeof(double) * l);

    for (int i = k - 1; i >= 0; i--)
    {
        for (int row = 0; row < m; row++)
            solution[perm[i] * m + row] = b[i * m + row];

        for (int j = i + 1; j < k + (l != 0); j++)
        {
            int width     = std::min(m, n - j * m);
            double* block = A + i * m * n + j * m * m;
            int offset    = (j < k ? perm[j] : k) * m;

            for (int row = 0; row < m; row++)
            {
                double sum = 0.;
                for (int col = 0; col < width; col++)
                    sum += block[row * width + col] * solution[offset + col];
                solution[perm[i] * m + row] -= sum;
            }
        }
    }

    return 1;
}

int LinearSystem::solve()
{
    bool flag = 0;

    flag = blocked_matrix_to_triangular();
    if (!flag)
        return -1;

    flag = triangular_blocked_to_solution();
    if(!flag)
        return -1;

    return 0;
}
