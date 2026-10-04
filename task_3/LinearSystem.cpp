#include "./LinearSystem.h"
#include "./matrix_io.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

namespace
{
    inline int min(int a, int b)
    {
        return (a > b ? b : a);
    }

    // inline int max(int a, int b)
    // {
    //     return (a > b ? a : b);
    // }

    inline double f(int n, int s, int i, int j)
    {
        int t;
        switch (s)
        {
        case 1:
            t = i > j ? i : j;
            return ((double)(n - t + 1));
            break;
        case 2:
            t = i > j ? i : j;
            return (double)t;
            break;
        case 3:
            t = i - j;
            return (double)(t >= 0 ? t : -t);
            break;
        case 4:
            return (double)1 / (double)(i + j - 1);
            break;
        }
        return 0;
    }

    inline int matrix_offset(int n, int m, int row, int col)
    {
        int i = row / m;
        int j = col / m;

        int wr = min(m, n - i * m);
        int wc = min(m, n - j * m);

        return (i * n + j * wr) * m + (row % m) * wc + col % m;
    }

}

LinearSystem::LinearSystem(int _n, int _m) : n(_n), m(_m)
{
}

LinearSystem::~LinearSystem()
{
    free_memory();
}

bool LinearSystem::memory_alloc()
{
    A = (double*)malloc(n * n * sizeof(double));
    if (A == nullptr)
        return 0;
    memset(A, 0, sizeof(double) * n * n);

    b = (double*)malloc(n * sizeof(double));
    if (b == nullptr)
        return 0;
    memset(b, 0, sizeof(double) * n);

    solution = (double*)malloc(sizeof(double) * n);
    if (solution == nullptr)
        return 0;
    memset(solution, 0, sizeof(double) * n);

    norm_workspace = (double*)malloc(sizeof(double) * n);
    if (norm_workspace == nullptr)
        return 0;
    memset(norm_workspace, 0, sizeof(double) * n);

    perm = (int*)malloc(sizeof(int) * n);
    if (perm == nullptr)
        return 0;
    memset(perm, 0, sizeof(int) * n);

    return 1;
}

void LinearSystem::free_memory()
{
    free(A);
    free(b);
    free(solution);
    free(norm_workspace);
    free(perm);
    A = nullptr;
    b = nullptr;
    solution = nullptr;
    norm_workspace = nullptr;
    perm = nullptr;
}

bool LinearSystem::init(int s, char* file_name)
{
    bool flag = 0;
    
    if (s == 0)
        flag = init_matrix_from_file(file_name);
    else
        flag = init_matrix_from_formula(s);

    if (!flag)
        return 0;

    flag = init_rhs();

    if (!flag)
        return 0;
    
    flag = init_perm();

    return flag;
}

bool LinearSystem::init_matrix_from_formula(int s)
{
    double* block = A;
    for (int row = 0; row < n; )
    {
        int height = min(m, n - row);

        for (int col = 0; col < n; )
        {
            int width = min(m, n - col);

            for (int p = 0; p < height; p++)
                for (int q = 0; q < width; q++)
                    *block++ = f(n, s, row + p + 1, col + q + 1);

            col += width;
        }
        row += height;
    }
    return 1;
}

// todo: remove matrix_offset
bool LinearSystem::init_matrix_from_file(char* file_name)
{
    if (file_name == nullptr)
        return 0;

    FILE* file = fopen(file_name, "r");
    if (file == nullptr)
        return 0;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            double* value = A + matrix_offset(n, m, i, j);

            if (fscanf(file, "%lf", value) != 1 || isnan(*value))
            {
                fclose(file);
                return 0;
            }
        }

    char c;
    bool flag = fscanf(file, " %c", &c) == EOF && !ferror(file);

    fclose(file);
    return flag;
}

bool LinearSystem::init_rhs()
{
    if (b == nullptr)
        return 0;

    const double* block = A;

    for (int row = 0; row < n; )
    {
        int height = min(m, n - row);

        for (int p = 0; p < height; p++)
            b[row + p] = 0.;

        for (int col = 0; col < n; )
        {
            int width = min(m, n - col);

            for (int p = 0; p < height; p++)
            {
                double sum = b[row + p];

                // столбцы с чётным глобальным индексом
    
                for (int q = col % 2; q < width; q += 2)
                    sum += block[q];

                b[row + p] = sum;
                block += width;
            }

            col += width;
        }

        for (int p = 0; p < height; p++)
            if (!isfinite(b[row + p]))
                return 0;

        row += height;
    }
    return 1;
}

bool LinearSystem::init_perm()
{
    for (int j = 0; j < n; j++)
        perm[j] = j;
    return 1;
}

void LinearSystem::print_matrix(int r) const
{
    ::print_matrix(n, n, m, A, r);
}

void LinearSystem::print_rhs(int r) const
{
    ::print_matrix(1, n, m, b, r);
}

void LinearSystem::print_solution(int r) const
{
    ::print_matrix(1, n, m, solution, r);
}

void LinearSystem::get_block(int i, int j, double* dest) const
{
    int height = min(m, n - i * m);
    int width = min(m, n - j * m);
    int offset = i * m * n + j * m * height;

    memcpy(dest, A + offset, height * width * sizeof(double));
}

void LinearSystem::set_block(int i, int j, const double* src)
{
    int height = min(m, n - i * m);
    int width = min(m, n - j * m);
    int offset = i * m * n + j * m * height;

    memcpy(A + offset, src, height * width * sizeof(double));
}

// todo: заместо HUGE_VAL переписать с поправкой на вычиселия машинного эпсиолна
void LinearSystem::compute_residuals(double& r1, double& r2) const
{
    r1 = -1; r2 = -1;
    if (n < 1 || m < 1 || A == nullptr || b == nullptr || solution == nullptr)
        return;

    double residual = 0.;
    double norm_b = 0.;
    double error = 0.;

    // обход строк в блочном хранении
    for (int row = 0; row < n; row++)
    {
        int block_row = row / m;
        int height = min(m, n - block_row * m);
        int local_row = row % m;
        double ax = 0.;

        for (int col = 0, block_col = 0; col < n; block_col++)
        {
            int width = min(m, n - col);
            int offset = block_row * m * n + block_col * m * height + local_row * width;

            for (int j = 0; j < width; j++)
                ax += A[offset + j] * solution[col + j];

            col += width;
        }

        residual += fabs(ax - b[row]);
        norm_b += fabs(b[row]);
        error += fabs(solution[row] - (row % 2 == 0 ? 1. : 0.));
    }

    r1 = norm_b > 0. ? residual / norm_b : (residual <= 0. ? 0. : HUGE_VAL);

    r2 = error / (n / 2 + n % 2);
}
