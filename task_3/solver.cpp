#include "./solver.h"
#include "./block_ops.h"
#include "./LinearSystem.h"

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <cstring>

bool LinearSystem::finding_block_pivot(int alpha, int& pivot_i, int& pivot_j, double* ws_block1, double* ws_block2, double* ws_vector_norms, double* block_inverse)
{
    pivot_i = -1;
    pivot_j = -1;

    double pivot_norm = 0.;
    double pivot_eps = (n * MACHINE_EPS) * original_matrix_norm;
    double* inverse_candidate = ws_block2;
    double* inverse_best = block_inverse;

    for (int i = alpha; i < k; i++)
    {
        for (int j = alpha; j < k; j++)
        {
            get_block(i, j, ws_block1);
            if (!invert_block(ws_block1, inverse_candidate, m, w_perm, pivot_eps))
                continue;

            double inverse_norm = calculate_block_norm(inverse_candidate, m, ws_vector_norms);
            if (!std::isfinite(inverse_norm))
                continue;

            // Первый подходящий блок задаёт начальную норму для сравнения.
            if (pivot_i == -1 || inverse_norm < pivot_norm)
            {
                pivot_i = i;
                pivot_j = j;
                pivot_norm = inverse_norm;
                std::swap(inverse_best, inverse_candidate);
            }
        }
    }

    if (pivot_i == -1)
        return 0;

    double eps = MACHINE_EPS * MACHINE_EPS * pivot_norm;
    for (int t = 0; t < m * m; t++)
    {
        double value = inverse_best[t];
        block_inverse[t] = std::fabs(value) < eps ? 0. : value;
    }
    return 1;
}

bool LinearSystem::blocked_matrix_to_triangular()
{
    double column_eps = (MACHINE_EPS * MACHINE_EPS) * original_matrix_norm;
    // используем block1 для обратного блока, block2 и block3 - для работы, block4 - тоже(в начале для вычисления норма, потом для свободных членов)
    // block1 == pivot_inverse 
    for (int alpha = 0; alpha < k; alpha++)
    {
        int pivot_i, pivot_j;
        if (!finding_block_pivot(alpha, pivot_i, pivot_j, block2, block3, block4, block1))
            return 0;
        
        swap_blocked_rows(alpha, pivot_i, alpha);
        swap_blocked_columns(alpha, pivot_j);

        {
            // на место главного блока ставим единичный
            double* pivot_block = get_block_address(alpha, alpha);
            std::memset(pivot_block, 0, sizeof(double) * m * m);
            for (int i = 0; i < m; i++)
                pivot_block[i * m + i] = 1.;

            // обновляем все блоки справа от него
            for (int j = alpha + 1; j < block_count; j++)
            {
                matrix_multiplication(block1, get_block_address(alpha, j), block3, m, m, std::min(m, n - j * m));
                set_block(alpha, j, block3);
            }

            // столбец свободных коэффициентов
            matrix_multiplication(block1, b + m * alpha, block2, m, m, 1);
            std::memcpy(b + m * alpha, block2, sizeof(double) * m);
        }

        {
            // Единичный диагональный блок задаёт нижнюю границу её масштаба.
            double* pivot_row = get_block_address(alpha, alpha);
            int pivot_row_size = m * (n - alpha * m);
            double pivot_row_max = 0.;
            for (int t = m * m; t < pivot_row_size; t++) // ищем максимальный по модулю элемент в оставшейся строке
                pivot_row_max = std::max(pivot_row_max, std::fabs(pivot_row[t]));
            double pivot_row_eps = (MACHINE_EPS * MACHINE_EPS) * pivot_row_max;
            for (int t = m * m; t < pivot_row_size; t++)
                if (std::fabs(pivot_row[t]) < pivot_row_eps)
                    pivot_row[t] = 0.;
        }

        {
            // обновляем нижние блоки и свободные члены(все блоки и квадратные, и прямоугольные)
            for (int i = alpha + 1; i < block_count; i++)
            {
                int height = std::min(m, n - i * m);
                double* i_alpha = get_block_address(i, alpha);
                for (int t = 0; t < height * m; t++)
                    if (std::fabs(i_alpha[t]) < column_eps)
                        i_alpha[t] = 0.;

                for (int j = alpha + 1; j < block_count; j++)
                {
                    matrix_multiplication_subtract(i_alpha, get_block_address(alpha, j), get_block_address(i, j), height, m, std::min(m, n - j * m));
                }

                matrix_multiplication(i_alpha, b + alpha * m, block2, height, m, 1);

                for (int t = 0; t < height; t++)
                    b[i * m + t] -= block2[t];
                
                // std::memset(i_alpha, 0, sizeof(double) * height * m);
            }
        }
    }

    // Последний блок решаем один раз после исключения всех полных блоков.
    if (l != 0)
    {
        get_block(k, k, block1);
        std::memcpy(block2, b + k * m, sizeof(double) * l);
        if (!naive_gauss_solve_block(block1, l, block2, w_perm, l * original_matrix_norm * MACHINE_EPS))
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

        for (int j = i + 1; j < block_count; j++)
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
