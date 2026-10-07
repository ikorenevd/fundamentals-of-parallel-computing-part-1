#include "./LinearSystem.h"
#include "./matrix_io.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <climits>
#include <utility>

namespace
{
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

        int wr = std::min(m, n - i * m);
        int wc = std::min(m, n - j * m);

        return (i * n + j * wr) * m + (row % m) * wc + col % m;
    }

}

LinearSystem::LinearSystem(int _n, int _m)
{
    if (_n < 1 || _m < 1 || _m > _n)
        return;

    n = _n;
    m = _m;
    k = n / m;
    l = n % m;
}

LinearSystem::~LinearSystem()
{
    free_memory();
}

bool LinearSystem::memory_alloc()
{
    A = (double*)std::malloc(n * n * sizeof(double));
    if (A == nullptr)
        return 0;
    std::memset(A, 0, sizeof(double) * n * n);

    b = (double*)std::malloc(n * sizeof(double));
    if (b == nullptr)
        return 0;
    std::memset(b, 0, sizeof(double) * n);

    solution = (double*)std::malloc(sizeof(double) * n);
    if (solution == nullptr)
        return 0;
    std::memset(solution, 0, sizeof(double) * n);

    perm = (int*)std::malloc(sizeof(int) * k);
    if (perm == nullptr)
        return 0;
    std::memset(perm, 0, sizeof(int) * k);

    w_perm = (int*)std::malloc(sizeof(int) * m);
    if (w_perm == nullptr)
        return 0;

    block1 = (double*)std::malloc(sizeof(double) * m * m);
    if (block1 == nullptr)
        return 0;
    
    block2 = (double*)std::malloc(sizeof(double) * m * m);
    if (block2 == nullptr)
        return 0;

    block3 = (double*)std::malloc(sizeof(double) * m * m);
    if (block3 == nullptr)
        return 0;

    block4 = (double*)std::malloc(sizeof(double) * m * m);
    if (block4 == nullptr)
        return 0;

    return 1;
}

void LinearSystem::free_memory()
{
    std::free(A);
    std::free(b);
    std::free(solution);
    std::free(perm);
    std::free(w_perm);
    std::free(block1);
    std::free(block2);
    std::free(block3);
    std::free(block4);

    A = nullptr;
    b = nullptr;
    solution = nullptr;
    perm = nullptr;
    w_perm = nullptr;
    block1 = nullptr;
    block2 = nullptr;
    block3 = nullptr;
    block4 = nullptr;
}

bool LinearSystem::init(int s, char* file_name)
{
    if (n < 1 || m < 1 || m > n || n > INT_MAX / n || s < 0 || s > 4
        || (s == 0 && file_name == nullptr)
        || A == nullptr || b == nullptr || solution == nullptr || perm == nullptr
        || w_perm == nullptr || block1 == nullptr || block2 == nullptr
        || block3 == nullptr || block4 == nullptr)
        return 0;

    bool flag = 0;
    
    if (s == 0)
        flag = init_matrix_from_file(file_name);
    else
        flag = init_matrix_from_formula(s);

    if (!flag)
        return 0;

    if (!std::isfinite(get_matrix_norm()))
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
        int height = std::min(m, n - row);

        for (int col = 0; col < n; )
        {
            int width = std::min(m, n - col);

            for (int p = 0; p < height; p++)
                for (int q = 0; q < width; q++)
                    *block++ = f(n, s, row + p + 1, col + q + 1);

            col += width;
        }
        row += height;
    }
    return 1;
}

bool LinearSystem::init_matrix_from_file(char* file_name)
{
    if (file_name == nullptr)
        return 0;

    std::FILE* file = std::fopen(file_name, "r");
    if (file == nullptr)
        return 0;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            double* value = A + matrix_offset(n, m, i, j);

            if (std::fscanf(file, "%lf", value) != 1 || std::isnan(*value))
            {
                std::fclose(file);
                return 0;
            }
        }

    char c;
    bool flag = std::fscanf(file, " %c", &c) == EOF && !std::ferror(file);

    std::fclose(file);
    return flag;
}

bool LinearSystem::init_rhs()
{
    if (b == nullptr)
        return 0;

    const double* block = A;

    for (int row = 0; row < n; )
    {
        int height = std::min(m, n - row);

        for (int p = 0; p < height; p++)
            b[row + p] = 0.;

        for (int col = 0; col < n; )
        {
            int width = std::min(m, n - col);

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

        row += height;
    }
    return 1;
}

bool LinearSystem::init_perm()
{
    if (m < 1 || n < m || perm == nullptr)
        return 0;
    for (int j = 0; j < k; j++)
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
    int height = std::min(m, n - i * m);
    int width  = std::min(m, n - j * m);
    int offset = i * m * n + j * m * height;

    std::memcpy(dest, A + offset, height * width * sizeof(double));
}

void LinearSystem::set_block(int i, int j, const double* src)
{
    int height = std::min(m, n - i * m);
    int width  = std::min(m, n - j * m);
    int offset = i * m * n + j * m * height;

    std::memcpy(A + offset, src, height * width * sizeof(double));
}

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
        double ax    = 0.;

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


// todo: переписать, суммируя в блоках
double LinearSystem::get_matrix_norm() const
{
    double norm     = 0.;
    int block_count = k + (l != 0);

    for (int i = 0; i < block_count; i++)
    {
        int height = std::min(m, n - i * m);

        for (int row = 0; row < height; row++)
        {
            double sum = 0.;

            for (int j = 0; j < block_count; j++)
            {
                int width           = std::min(m, n - j * m);
                const double* block = A + i * m * n + j * m * height;

                for (int col = 0; col < width; col++)
                    sum += std::fabs(block[row * width + col]);
            }

            if (sum > norm)
                norm = sum;
        }
    }

    return norm;
}

// Перестановка двух полных блочных строк.
void LinearSystem::swap_blocked_rows(int i, int j, int first_column)
{
    if (i == j)
        return;

    std::swap_ranges(A + i * n * m + first_column * m * m,
                     A + (i + 1) * n * m,
                     A + j * n * m + first_column * m * m);
    std::swap_ranges(b + i * m, b + (i + 1) * m, b + j * m);
}

void LinearSystem::swap_blocked_columns(int i, int j)
{
    if (i == j)
        return;
    
    // меняем квадратные блоки
    for (int row = 0; row < k; row++)
        std::swap_ranges(A + row * n * m + i * m * m,
                         A + row * n * m + (i + 1) * m * m,
                         A + row * n * m + j * m * m);

    int l = n % m;
    if (l != 0)
    {
        double* tail = A + (k) * m * n;
        std::swap_ranges(tail + i * m * l,
                         tail + (i + 1) * m * l,
                         tail + j * m * l);
    }

    std::swap(perm[i], perm[j]);
}
