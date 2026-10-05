#include "./solver.h"
#include "./LinearSystem.h"

#include <math.h>
#include <algorithm>
#include <string.h>

#define MACHINE_EPS 1e-16

double calculate_block_norm(const double* block, int m)
{
    if (block == nullptr || m < 1)
        return HUGE_VAL;

    double norm = 0.;
    for (int j = 0; j < m; j++)
    {
        double sum = 0.;
        for (int i = 0; i < m; i++)
            sum += fabs(block[i * m + j]);

        if (!isfinite(sum))
            return HUGE_VAL;

        norm = std::max(norm, sum);
    }
    return norm;
}

// A = f x l, B = l x k, A * B = f x k, все матрицы хранятся по строчно
void matrix_multiplication(const double* A, const double *B, double* result, int f, int l, int k)
{
    memset(result, 0, sizeof(double) * f * k);
    for (int i = 0; i < f; i++)
        for (int j = 0; j < k; j++)
            for (int u = 0; u < l; u++)
                result[i * k + j] += A[i * l + u] * B[u * k + j];
}

// 1 если нашли успешно и матрица обратима, 0 если обратное
// inverse и block должны быть выделены m * m
// block_perm m
// при поиске обратной block изменяется, его нужно выделять с помощью get_block
// при невозможности посчитать обратную получаем abs_inverse_det = 0
bool invert_block(double* block, double* inverse, int m, int* block_perm, double& abs_inverse_det)
{
    abs_inverse_det = 0.;

    double determinant = 1.;

    // делаем inverse единичной
    for (int i = 0; i < m * m; i++)
        inverse[i] = (i % (m + 1) == 0) ? 1. : 0.;
    
    // перестановки столбцов
    for (int i = 0; i < m; i++)
        block_perm[i] = i;
    
    for (int alpha = 0; alpha < m; alpha++)
    {
        int pivot_i = alpha;
        int pivot_j = alpha;
        double abs_pivot = 0.;

        // главный элемент
        for (int i = alpha; i < m; i++)
            for (int j = alpha; j < m; j++)
                if(fabs(block[i * m + j]) > abs_pivot)
                {
                    abs_pivot = fabs(block[i * m + j]);
                    pivot_i = i;
                    pivot_j = j;
                }

        if (abs_pivot < MACHINE_EPS)
            return 0;

        // важен только знак определителя
        // if (pivot_i != alpha)
        //     determinant = -determinant;
        // if (pivot_j != alpha)
        //     determinant = -determinant;

        // меняем строки местами
        for (int j = 0; j < m; j++)
        {
            double temp = block[alpha * m + j];
            block[alpha * m + j] = block[pivot_i * m + j];
            block[pivot_i * m + j] = temp;

            temp = inverse[alpha * m + j];
            inverse[alpha * m + j] = inverse[pivot_i * m + j];
            inverse[pivot_i * m + j] = temp;
        }

        // меняем столбцы
        for (int i = 0; i < m; i++)
        {
            double temp = block[i * m + alpha];
            block[i * m + alpha] = block[i * m + pivot_j];
            block[i * m + pivot_j] = temp;
        }

        // индексы столбцов
        int temp_int = block_perm[alpha];
        block_perm[alpha] =  block_perm[pivot_j];
        block_perm[pivot_j] = temp_int;

        // обнуляем
        for (int i = alpha + 1; i < m; i++)
        {
            double coeff = block[i * m + alpha] / block[alpha * m + alpha];

            for (int j = alpha + 1; j < m; j++)
                block[i * m + j] -= coeff * block[alpha * m + j];
            
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
        determinant *= r;
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
        {
            double temp = inverse[i * m + j];
            inverse[i * m + j] = inverse[k * m + j];
            inverse[k * m + j] = temp;
        }

        block_perm[i] = block_perm[k];
        block_perm[k] = k;
    }

    abs_inverse_det = fabs(determinant);
    return 1;
}

bool LinearSystem::finding_block_pivot(int alpha, int& pivot_i, int& pivot_j, double* block_inverse)
{
    pivot_i = -1;
    pivot_j = -1;

    double pivot_norm = 0.;
    double abs_inverse_det = 0.;

    for (int i = alpha; i < k; i++)
    {
        for (int j = alpha; j < k; j++)
        {
            // block2 хранит кандидата, block3 — его обратную матрицу.
            get_block(i, j, block2);
            if (!invert_block(block2, block3, m, w_perm, abs_inverse_det))
                continue;

            const double inverse_norm = calculate_block_norm(block3, m);
            if (!isfinite(inverse_norm))
                continue;

            // Первый подходящий блок задаёт начальную норму для сравнения.
            if (pivot_i == -1 || inverse_norm < pivot_norm)
            {
                pivot_i = i;
                pivot_j = j;
                pivot_norm = inverse_norm;
                std::copy(block3, block3 + m * m, block_inverse);
            }
        }
    }

    if (pivot_i == -1)
        return 0;

    return 1;
}

// Построчное хранение: m == n или m == 1.
bool LinearSystem::naive_full_matrix_to_triangular()
{
    int alpha = 0;

    for (alpha = 0; alpha < n; alpha++)
    {
        double abs_pivot = 0., temp = 0.;
        int pivot_i = 0, pivot_j = 0;
        // ищем наибольший по модулю элемент
        for (int i = alpha; i < n; i++)
            for (int j = alpha; j < n; j++)
                if (fabs(A[i * n + j]) > abs_pivot)
                {
                    abs_pivot = fabs(A[i * n + j]);
                    pivot_i = i;
                    pivot_j = j;
                }

        // todo: нужно ли делать относительный эпсион(домножать на норму матрицы)
        // если главный элемент не подходит, то выходим
        if (!isfinite(abs_pivot) || abs_pivot < MACHINE_EPS)
            return 0;

        // физически меняем строки alpha и pivot_i
        for (int j = alpha; j < n; j++)
        {
            temp = A[pivot_i * n + j];
            A[pivot_i * n + j] = A[alpha * n + j];
            A[alpha * n + j] = temp;
        }
        temp = b[pivot_i];
        b[pivot_i] = b[alpha];
        b[alpha] = temp;

        // физически меняем столбцы alpha и pivot_j
        for (int i = 0; i < n; i++)
        {
            temp = A[i * n + alpha];
            A[i * n + alpha] = A[i * n + pivot_j];
            A[i * n + pivot_j] = temp;
        }
        int t = perm[alpha];
        perm[alpha] = perm[pivot_j];
        perm[pivot_j] = t;


        // зануляем
        for (int i = alpha + 1; i < n; i++)
        {
            double coeff = A[i * n + alpha] / A[alpha * n + alpha];
            if (!isfinite(coeff))
                return 0;

            for (int j = alpha + 1; j < n; j++)
            {
                A[i * n + j] -= coeff * A[alpha * n + j];
                if (!isfinite(A[i * n + j]))
                    return 0;
            }
            b[i] -= coeff * b[alpha];
            if (!isfinite(b[i]))
                return 0;

            A[i * n + alpha] = 0.;
        }
    }

    return 1;
}

// todo: какая-то из проверок лишняя
bool LinearSystem::naive_traingular_solution()
{
    for (int i = n - 1; i >= 0; i--)
    {
        double diagonal = A[i * n + i];
        if (!isfinite(diagonal) || fabs(diagonal) < MACHINE_EPS)
            return 0;

        double rhs = b[i];
        for (int j = i + 1; j < n; j++)
            rhs -= A[i * n + j] * solution[perm[j]];

        if (!isfinite(rhs))
            return 0;

        double value = rhs / diagonal;
        if (!isfinite(value))
            return 0;

        solution[perm[i]] = value;
    }

    return 1;
}

bool LinearSystem::blocked_matrix_to_triangular()
{
    // используем block1 для обратного блока
    for (int alpha = 0; alpha < k; alpha++)
    {
        int pivot_i, pivot_j;
        if (!finding_block_pivot(alpha, pivot_i, pivot_j, block1))
            return 0;
        
        swap_blocked_rows(alpha, pivot_i);
        swap_blocked_columns(alpha, pivot_j);

        // todo: сделай единой формулой
        {
            // block2 = identiy matrix
            for (int i = 0; i < m * m; i++)
                block2[i] = (i % (m + 1) == 0) ? 1. : 0.;
            set_block(alpha, alpha, block2);
            
            // изменяем квадратные блоки
            for (int j = alpha + 1; j < k; j++)
            {
                get_block(alpha, j, block2);
                matrix_multiplication(block1, block2, block3, m, m, m);
                set_block(alpha, j, block3);
            }
            
            // неквадратный
            if (l != 0)
            {
                get_block(alpha, k, block2);
                matrix_multiplication(block1, block2, block3, m, m, l);
                set_block(alpha, k, block3);
            }
            
            // блок свободных членов
            for (int i = 0; i < m; i++)
            {
                double sum = 0.;
                for (int col = 0; col < m; col++)
                    sum += block1[i * m + col] * b[alpha * m + col];

                if (!isfinite(sum))
                    return 0;

                block3[i] = sum;
            }

            std::copy(block3, block3 + m, b + alpha * m);
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
                    
                    matrix_multiplication(block1, block2, block3, height, m, width);
                    
                    get_block(i, j, block4);
                    
                    for (int t = 0; t < width * height; t++)
                    {
                        block4[t] -= block3[t];
                        if (!isfinite(block4[t]))
                            return 0;
                    }
                    
                    set_block(i, j, block4);
                }

                matrix_multiplication(block1, b + alpha * m, block2, height, m, 1);

                for (int t = 0; t < height; t++)
                {
                    b[i * m + t] -= block2[t];
                    if (!isfinite(b[i * m + t]))
                        return 0;
                }
                
                
                memset(block1, 0, sizeof(double) * height * m);
                set_block(i, alpha, block1);
            }
        }

    }

    // Последний блок решаем один раз после исключения всех полных блоков.
    if (l != 0)
    {
        double abs_inverse_det = 0.;
        get_block(k, k, block1);

        if (!invert_block(block1, block2, l, w_perm, abs_inverse_det))
            return 0;

        matrix_multiplication(block2, b + k * m, block1, l, l, 1);
        std::copy(block1, block1 + l, b + k * m);
    }

    return 1;
}

bool LinearSystem::triangular_blocked_to_solution()
{
    if (l != 0)
        std::copy(b + m * k, b + m * k + l, solution + m * k);

    for (int i = k - 1; i >= 0; i--)
    {
        for (int row = 0; row < m; row++)
            solution[perm[i] * m + row] = b[i * m + row];

        for (int j = i + 1; j < k + (l != 0); j++)
        {
            int width = std::min(m, n - j * m);
            double* block = A + i * m * n + j * m * m;
            int offset = (j < k ? perm[j] : k) * m;

            for (int row = 0; row < m; row++)
            {
                double sum = 0.;
                for (int col = 0; col < width; col++)
                    sum += block[row * width + col] * solution[offset + col];

                if (!isfinite(sum))
                    return 0;

                solution[perm[i] * m + row] -= sum;
                if (!isfinite(solution[perm[i] * m + row]))
                    return 0;
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
