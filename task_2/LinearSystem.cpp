#include "./LinearSystem.h"

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

    x = (double*)malloc(sizeof(double) * n);
    if (x == nullptr)
        return 0;
    memset(x, 0, sizeof(double) * n);

    return 1;
}

void LinearSystem::free_memory()
{
    free(A);
    free(b);
    free(x);
    A = nullptr;
    b = nullptr;
    x = nullptr;
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

    return flag;
}

bool LinearSystem::init_matrix_from_formula(int s)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            *(A + matrix_offset(n, m, i, j)) = f(n, s, i + 1, j + 1);
    return 1;
}

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

    char extra;
    bool complete = fscanf(file, " %c", &extra) == EOF && !ferror(file);
    fclose(file);
    return complete;
}

bool LinearSystem::init_rhs()
{
    if (b == nullptr)
        return 0;

    for (int i = 0; i < n; i++)
    {
        double sum = 0.;

        for (int k = 0; k < n; k += 2)
            sum += *(A + matrix_offset(n, m, i, k));

        *(b + i) = sum;

        if (isinf(*(b + i)))
            return 0;
    }


    return 1;
}

// todo: удалить после написания нормального решения
bool LinearSystem::init_solution()
{
    if (x == nullptr)
        return 0;

    for (int i = 0; i < n; i ++)
        x[i] = (i % 2 == 0) ? 1. : 0.;

    return 0;
}

void LinearSystem::print_matrix(int r) const
{
    int s = min(r, n);

    for (int i = 0; i < s; i++)
    {
        for (int j = 0; j < s; j++)
            printf(" %10.3e", *(A + matrix_offset(n, m, i, j)));
        printf("\n");
    }
}

void LinearSystem::print_rhs(int r) const
{
    int s = min(r, n);

    for (int i = 0; i < s; i++)
        printf(" %10.3e", *(b + i));
    
    printf("\n");
}

void LinearSystem::print_solution(int r) const
{
    int s = min(r, n);

    for (int i = 0; i < s; i++)
        printf(" %10.3e", *(x + i));
    
    printf("\n");
}


// todo когда нормально потребуют
bool LinearSystem::solve()
{
    if (isinf(matrix_norm) || isnan(matrix_norm))
        return 0;

    init_solution();

    return 1;
}

void LinearSystem::get_block(int i, int j, double* dest) const
{
    const int height = min(m, n - i * m);
    const int width = min(m, n - j * m);
    const size_t offset = (size_t)i * m * n + (size_t)j * m * height;
    memcpy(dest, A + offset, (size_t)height * width * sizeof(double));
}

void LinearSystem::compute_residuals(double& r1, double& r2) const
{
    r1 = -1; r2 = -1;
    if (n < 1 || m < 1 || A == nullptr || b == nullptr || x == nullptr)
        return;

    const size_t block_size = (size_t)min(m, n);
    double* work = (double*)malloc((block_size * block_size + 2 * block_size)
                                  * sizeof(double));
    if (work == nullptr)
        return;
    double* work_a = work;
    double* work_x = work_a + block_size * block_size;
    double* ax = work_x + block_size;

    double residual = 0.;
    double norm_b = 0.;
    double error = 0.;

    for (int row = 0, block_row = 0; row < n; ++block_row)
    {
        const int height = min(m, n - row);
        memset(ax, 0, height * sizeof(double));

        for (int col = 0, block_col = 0; col < n; ++block_col)
        {
            const int width = min(m, n - col);
            get_block(block_row, block_col, work_a);
            memcpy(work_x, x + col, width * sizeof(double));
            const double* block = work_a;
            for (int i = 0; i < height; ++i)
            {
                double sum = ax[i];
                for (int j = 0; j < width; ++j)
                    sum += block[j] * work_x[j];
                ax[i] = sum;
                block += width;
            }
            col += width;
        }

        for (int i = 0; i < height; ++i)
        {
            const int index = row + i;
            residual += fabs(ax[i] - b[index]);
            norm_b += fabs(b[index]);
            error += fabs(x[index] - (index % 2 == 0 ? 1. : 0.));
        }
        row += height;
    }

    free(work);
    r1 = norm_b > 0. ? residual / norm_b
        : (residual <= 0. ? 0. : HUGE_VAL);
    r2 = error / (n / 2 + n % 2);
}

void LinearSystem::calculate_matrix_norm()
{
    double sum = 0.;
    if (n < 1)
        return;


    for (int i = 0; i < n; i++)
        sum += fabs(*(A + matrix_offset(n, m, i, 0)));
    

    for (int j = 1; j < n; j++)
    {
        sum = 0.f;

        for (int i = 0; i < n; i++)
            sum += fabs(*(A + matrix_offset(n, m, i, j)));
        
        if (sum > matrix_norm)
            matrix_norm = sum;
    }
}
