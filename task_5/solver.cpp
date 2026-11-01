#include "./solver.h"
#include "./block_ops.h"
#include "./LinearSystem.h"

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <cstring>

bool LinearSystem::finding_block_pivot(int alpha, int& pivot_i, int& pivot_j, double* ws_block1, double* ws_block2, double* block_inverse)
{
    pivot_i = -1;
    pivot_j = -1;

    double pivot_norm = 0.;
    double pivot_eps = (matrix_size * MACHINE_EPS) * original_matrix_norm;
    double* inverse_candidate = ws_block2;
    double* inverse_best = block_inverse;

    for (int i = alpha; i < full_block_count; i++)
    {
        for (int j = alpha; j < full_block_count; j++)
        {
            get_block(i, j, ws_block1);
            if (!invert_block(ws_block1, inverse_candidate, block_size, w_perm, pivot_eps))
                continue;

            double inverse_norm = calculate_block_norm(inverse_candidate, block_size);
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
    for (int t = 0; t < block_size * block_size; t++)
    {
        double value = inverse_best[t];
        block_inverse[t] = std::fabs(value) < eps ? 0. : value;
    }
    return 1;
}

bool LinearSystem::blocked_matrix_to_triangular()
{
    double column_eps = (MACHINE_EPS * MACHINE_EPS) * original_matrix_norm;
    // используем block1 для обратного блока, block2 и block3 — для работы
    // block1 == pivot_inverse 
    for (int alpha = 0; alpha < full_block_count; alpha++)
    {
        int pivot_i, pivot_j;
        if (!finding_block_pivot(alpha, pivot_i, pivot_j, block2, block3, block1))
            return 0;
        
        swap_blocked_rows(alpha, pivot_i, alpha);
        swap_blocked_columns(alpha, pivot_j);

        {
            // на место главного блока ставим единичный
            double* pivot_block = get_block_address(alpha, alpha);
            std::memset(pivot_block, 0, sizeof(double) * block_size * block_size);
            for (int i = 0; i < block_size; i++)
                pivot_block[i * block_size + i] = 1.;

            // обновляем все блоки справа от него
            for (int j = alpha + 1; j < block_count; j++)
            {
                matrix_multiplication(block1, get_block_address(alpha, j), block3, block_size, block_size, std::min(block_size, matrix_size - j * block_size));
                set_block(alpha, j, block3);
            }

            // столбец свободных коэффициентов
            matrix_multiplication(block1, b + block_size * alpha, block2, block_size, block_size, 1);
            std::memcpy(b + block_size * alpha, block2, sizeof(double) * block_size);
        }

        {
            // Единичный диагональный блок задаёт нижнюю границу её масштаба.
            double* pivot_row = get_block_address(alpha, alpha);
            int pivot_row_size = block_size * (matrix_size - alpha * block_size);
            double pivot_row_max = 0.;
            for (int t = block_size * block_size; t < pivot_row_size; t++) // ищем максимальный по модулю элемент в оставшейся строке
                pivot_row_max = std::max(pivot_row_max, std::fabs(pivot_row[t]));
            double pivot_row_eps = (MACHINE_EPS * MACHINE_EPS) * pivot_row_max;
            for (int t = block_size * block_size; t < pivot_row_size; t++)
                if (std::fabs(pivot_row[t]) < pivot_row_eps)
                    pivot_row[t] = 0.;
        }

        {
            // обновляем нижние блоки и свободные члены(все блоки и квадратные, и прямоугольные)
            for (int i = alpha + 1; i < block_count; i++)
            {
                int height = std::min(block_size, matrix_size - i * block_size);
                double* i_alpha = get_block_address(i, alpha);
                for (int t = 0; t < height * block_size; t++)
                    if (std::fabs(i_alpha[t]) < column_eps)
                        i_alpha[t] = 0.;

                for (int j = alpha + 1; j < block_count; j++)
                {
                    matrix_multiplication_subtract(i_alpha, get_block_address(alpha, j), get_block_address(i, j), height, block_size, std::min(block_size, matrix_size - j * block_size));
                }

                matrix_multiplication(i_alpha, b + alpha * block_size, block2, height, block_size, 1);

                for (int t = 0; t < height; t++)
                    b[i * block_size + t] -= block2[t];
                
                // std::memset(i_alpha, 0, sizeof(double) * height * block_size);
            }
        }
    }

    // Последний блок решаем один раз после исключения всех полных блоков.
    if (remainder_size != 0)
    {
        get_block(full_block_count, full_block_count, block1);
        std::memcpy(block2, b + full_block_count * block_size, sizeof(double) * remainder_size);
        if (!naive_gauss_solve_block(block1, remainder_size, block2, w_perm, remainder_size * original_matrix_norm * MACHINE_EPS))
            return 0;
        std::memcpy(b + full_block_count * block_size, block2, sizeof(double) * remainder_size);
    }

    return 1;
}

bool LinearSystem::triangular_blocked_to_solution()
{
    if (remainder_size != 0)
        std::memcpy(solution + block_size * full_block_count, b + block_size * full_block_count, sizeof(double) * remainder_size);

    for (int i = full_block_count - 1; i >= 0; i--)
    {
        for (int row = 0; row < block_size; row++)
            solution[perm[i] * block_size + row] = b[i * block_size + row];

        for (int j = i + 1; j < block_count; j++)
        {
            int width     = std::min(block_size, matrix_size - j * block_size);
            double* block = A + i * block_size * matrix_size + j * block_size * block_size;
            int offset    = (j < full_block_count ? perm[j] : full_block_count) * block_size;

            for (int row = 0; row < block_size; row++)
            {
                double sum = 0.;
                for (int col = 0; col < width; col++)
                    sum += block[row * width + col] * solution[offset + col];
                solution[perm[i] * block_size + row] -= sum;
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
