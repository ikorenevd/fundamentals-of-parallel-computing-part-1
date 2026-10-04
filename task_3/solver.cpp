#include "./solver.h"
#include "./LinearSystem.h"

#include <math.h>

#define MACHINE_EPS 1e-16

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

// Построчное хранение: m == n или m == 1.
bool LinearSystem::naive_full_matrix_to_triangular()
{
    if (n < 1 || m < 1 || (m != n && m != 1)
        || A == nullptr || b == nullptr || perm == nullptr)
        return 0;

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
        if (abs_pivot < MACHINE_EPS)
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
            for (int j = alpha + 1; j < n; j++)
                A[i * n + j] -= coeff * A[alpha * n + j];
            b[i] -= coeff * b[alpha];
            A[i * n + alpha] = 0.;
        }
    }

    return 1;
}

// todo: поработать с делением нуля
bool LinearSystem::naive_traingular_solution()
{
    if (n < 1 || m < 1 || (m != n && m != 1)
        || A == nullptr || b == nullptr || solution == nullptr || perm == nullptr)
        return 0;

    for (int i = n - 1; i >= 0; i--)
    {
        double diagonal = A[i * n + i];

        double rhs = b[i];
        for (int j = i + 1; j < n; j++)
            rhs -= A[i * n + j] * solution[perm[j]];

        double value = rhs / diagonal;

        solution[perm[i]] = value;
    }

    return 1;
}

int LinearSystem::solve()
{
    bool flag = 0;

    if (n < 1 || m < 1 || A == nullptr || b == nullptr || solution == nullptr || norm_workspace == nullptr || perm == nullptr)
        return -1;

    double matrix_norm = calculate_matrix_norm(n, m, A, norm_workspace);

    if (!isfinite(matrix_norm))
        return -1;

    flag = naive_full_matrix_to_triangular();
    if (!flag)
        return -1;

    flag = naive_traingular_solution();
    if(!flag)
        return -1;

    return 0;
}
