#include "./block_ops.h"

#include <cmath>
#include <algorithm>
#include <cstring>

double calculate_block_norm(const double* block, int m, double* vector_norms)
{
    if (block == nullptr || m < 1 || vector_norms == nullptr)
        return HUGE_VAL;

    double max = 0.;
    std::memset(vector_norms, 0, sizeof(double) * m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            vector_norms[i] += std::fabs(block[i * m + j]);

    for (int i = 0; i < m; i++)
        if (vector_norms[i] > max)
            max = vector_norms[i];

    return max;
}

// A = f × l, B = l × k, result = f × k.
void matrix_multiplication(const double* A, const double* B, double* result, int f, int l, int k)
{
    int i = 0;
    for (; i + 2 < f; i += 3)
    {
        const double* a0 = A + i * l;
        const double* a1 = a0 + l;
        const double* a2 = a1 + l;

        double* r0 = result + i * k;
        double* r1 = r0 + k;
        double* r2 = r1 + k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c00 = 0., c01 = 0., c02 = 0.;
            double c10 = 0., c11 = 0., c12 = 0.;
            double c20 = 0., c21 = 0., c22 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;

                double x0 = a0[u];
                double x1 = a1[u];
                double x2 = a2[u];

                double y0 = b[0];
                double y1 = b[1];
                double y2 = b[2];

                c00 += x0 * y0; c10 += x1 * y0; c20 += x2 * y0;
                c01 += x0 * y1; c11 += x1 * y1; c21 += x2 * y1;
                c02 += x0 * y2; c12 += x1 * y2; c22 += x2 * y2;
            }

            r0[j]     = c00;
            r0[j + 1] = c01;
            r0[j + 2] = c02;

            r1[j]     = c10;
            r1[j + 1] = c11;
            r1[j + 2] = c12;

            r2[j]     = c20;
            r2[j + 1] = c21;
            r2[j + 2] = c22;
        }

        // Остаточные столбцы 3x1.
        for (; j < k; j++)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                double b = B[u * k + j];
                c0 += a0[u] * b;
                c1 += a1[u] * b;
                c2 += a2[u] * b;
            }

            r0[j] = c0;
            r1[j] = c1;
            r2[j] = c2;
        }
    }

    // Остаточные строки
    for (; i < f; i++)
    {
        const double* a = A + i * l;
        double* r       = result + i * k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;
                double x = a[u];

                c0 += x * b[0];
                c1 += x * b[1];
                c2 += x * b[2];
            }

            r[j]     = c0;
            r[j + 1] = c1;
            r[j + 2] = c2;
        }

        // Остаточные элементы 1x1
        for (; j < k; j++)
        {
            double sum = 0.;
            for (int u = 0; u < l; u++)
                sum += a[u] * B[u * k + j];
            r[j] = sum;
        }
    }
}

// A = f × l, B = l × k, result -= A × B; result = f × k.
void matrix_multiplication_subtract(const double* A, const double* B, double* result, int f, int l, int k)
{
    int i = 0;
    for (; i + 2 < f; i += 3)
    {
        const double* a0 = A + i * l;
        const double* a1 = a0 + l;
        const double* a2 = a1 + l;

        double* r0 = result + i * k;
        double* r1 = r0 + k;
        double* r2 = r1 + k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c00 = 0., c01 = 0., c02 = 0.;
            double c10 = 0., c11 = 0., c12 = 0.;
            double c20 = 0., c21 = 0., c22 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;

                double x0 = a0[u];
                double x1 = a1[u];
                double x2 = a2[u];

                double y0 = b[0];
                double y1 = b[1];
                double y2 = b[2];

                c00 += x0 * y0; c10 += x1 * y0; c20 += x2 * y0;
                c01 += x0 * y1; c11 += x1 * y1; c21 += x2 * y1;
                c02 += x0 * y2; c12 += x1 * y2; c22 += x2 * y2;
            }

            r0[j]     -= c00;
            r0[j + 1] -= c01;
            r0[j + 2] -= c02;

            r1[j]     -= c10;
            r1[j + 1] -= c11;
            r1[j + 2] -= c12;

            r2[j]     -= c20;
            r2[j + 1] -= c21;
            r2[j + 2] -= c22;
        }

        // Остаточные столбцы 3x1.
        for (; j < k; j++)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                double b = B[u * k + j];
                c0 += a0[u] * b;
                c1 += a1[u] * b;
                c2 += a2[u] * b;
            }

            r0[j] -= c0;
            r1[j] -= c1;
            r2[j] -= c2;
        }
    }

    // Остаточные строки
    for (; i < f; i++)
    {
        const double* a = A + i * l;
        double* r       = result + i * k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;
                double x = a[u];

                c0 += x * b[0];
                c1 += x * b[1];
                c2 += x * b[2];
            }

            r[j]     -= c0;
            r[j + 1] -= c1;
            r[j + 2] -= c2;
        }

        // Остаточные элементы 1x1
        for (; j < k; j++)
        {
            double sum = 0.;
            for (int u = 0; u < l; u++)
                sum += a[u] * B[u * k + j];
            r[j] -= sum;
        }
    }
}

// 1 если нашли успешно и матрица обратима, 0 если обратное
// inverse и block должны быть выделены m * m
// block_perm m
// при поиске обратной block изменяется, его нужно выделять с помощью get_block
bool invert_block(double* block, double* inverse, int m, int* block_perm)
{
    // делаем inverse единичной
    std::memset(inverse, 0, sizeof(double) * m * m);
    for (int i = 0; i < m; i++)
    {
        inverse[i * m + i] = 1.;
        // перестановки столбцов
        block_perm[i] = i;
    }

    for (int alpha = 0; alpha < m; alpha++)
    {
        int pivot_i      = alpha;
        int pivot_j      = alpha;
        double abs_pivot = 0.;

        // главный элемент
        for (int i = alpha; i < m; i++)
            for (int j = alpha; j < m; j++)
                if(std::fabs(block[i * m + j]) > abs_pivot)
                {
                    abs_pivot = std::fabs(block[i * m + j]);
                    pivot_i   = i;
                    pivot_j   = j;
                }

        if (abs_pivot < MACHINE_EPS)
            return 0;

        // меняем строки местами
        if (pivot_i != alpha)
            for (int j = 0; j < m; j++)
            {
                std::swap(block[alpha * m + j], block[pivot_i * m + j]);
                std::swap(inverse[alpha * m + j], inverse[pivot_i * m + j]);
            }

        // меняем столбцы
        if (pivot_j != alpha)
        {
            for (int i = 0; i < m; i++)
                std::swap(block[i * m + alpha], block[i * m + pivot_j]);
            // индексы столбцов
            std::swap(block_perm[alpha], block_perm[pivot_j]);
        }

        // обнуляем
        for (int i = alpha + 1; i < m; i++)
        {
            double coeff = block[i * m + alpha] / block[alpha * m + alpha];

            for (int j = alpha + 1; j < m; j++)
                block[i * m + j]   -= coeff * block[alpha * m + j];
            
            for (int j = 0; j < m; j++)
                inverse[i * m + j] -= coeff * inverse[alpha * m + j];
            
            block[i * m + alpha] = 0.;
        }

    }

    for (int i = m - 1; i >= 0; i--)
    {
        for (int k = i + 1; k < m; k++)
        {
            double coeff = block[i * m + k];

            for (int j = 0; j < m; j++)
                inverse[i * m + j] -= coeff * inverse[k * m + j];
        }

        double r = 1. / block[i * m + i];
        for (int j = 0; j < m; j++)
            inverse[i * m + j] *= r;
    }

    for (int i = 0; i < m; )
    {
        int k = block_perm[i];
        if (k == i)
        {
            i++;
            continue;
        }

        for (int j = 0; j < m; j++)
            std::swap(inverse[i * m + j], inverse[k * m + j]);
        std::swap(block_perm[i], block_perm[k]);
    }

    return 1;
}

bool naive_gauss_solve_block(double* block, int m, double* b, int* block_perm)
{
    int alpha = 0;

    for (alpha = 0; alpha < m; alpha++)
    {
        double abs_pivot = 0.;
        int pivot_i = 0, pivot_j = 0;

        // ищем наибольший по модулю элемент
        for (int i = alpha; i < m; i++)
            for (int j = alpha; j < m; j++)
                if (std::fabs(block[i * m + j]) > abs_pivot)
                {
                    abs_pivot = std::fabs(block[i * m + j]);
                    pivot_i   = i;
                    pivot_j   = j;
                }

        // todo: нужно ли делать относительный эпсион(домножать на норму матрицы)
        // если главный элемент не подходит, то выходим
        if (abs_pivot < MACHINE_EPS)
            return 0;

        // физически меняем строки alpha и pivot_i
        for (int j = alpha; j < m; j++)
            std::swap(block[pivot_i * m + j], block[alpha * m + j]);
        std::swap(b[pivot_i], b[alpha]);

        // физически меняем столбцы alpha и pivot_j
        for (int i = 0; i < m; i++)
            std::swap(block[i * m + alpha], block[i * m + pivot_j]);
        std::swap(block_perm[alpha], block_perm[pivot_j]);

        // зануляем
        for (int i = alpha + 1; i < m; i++)
        {
            double coeff = block[i * m + alpha] / block[alpha * m + alpha];

            for (int j = alpha + 1; j < m; j++)
                block[i * m + j] -= coeff * block[alpha * m + j];

            b[i]                -= coeff * b[alpha];
            block[i * m + alpha] = 0.;
        }
    }

    for (int i = m - 1; i >= 0; i--)
    {
        double diagonal = block[i * m + i];
        if (std::fabs(diagonal) < MACHINE_EPS)
            return 0;

        double rhs = b[i];
        for (int j = i + 1; j < m; j++)
            rhs -= block[i * m + j] * b[j];

        double value = rhs / diagonal;
        if (!std::isfinite(value))
            return 0;

        b[i] = value;
    }

    // Восстанавливаем исходный порядок неизвестных.
    for (int i = 0; i < m;)
    {
        int k = block_perm[i];
        if (k == i)
        {
            i++;
            continue;
        }

        std::swap(b[i], b[k]);
        std::swap(block_perm[i], block_perm[k]);
    }

    return 1;
}
