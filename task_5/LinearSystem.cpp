#include "./LinearSystem.h"
#include "./matrix_io.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <climits>
#include <utility>

namespace
{
    inline int matrix_offset(int n, int m, int row, int col)
    {
        int i = row / m;
        int j = col / m;

        int wr = std::min(m, n - i * m);
        int wc = std::min(m, n - j * m);

        return (i * n + j * wr) * m + (row % m) * wc + col % m;
    }

}

LinearSystem::LinearSystem(int _matrix_size, int _block_size)
{
    matrix_size      = _matrix_size;
    block_size       = std::min(_matrix_size, _block_size);
    full_block_count = matrix_size / block_size;
    remainder_size   = matrix_size % block_size;
    block_count      = full_block_count + (remainder_size != 0);
}

LinearSystem::~LinearSystem()
{
    free_memory();
}

bool LinearSystem::memory_alloc()
{
    coeff_matrix = (double*)std::malloc(matrix_size * matrix_size * sizeof(double));
    if (coeff_matrix == nullptr)
        return 0;
    std::memset(coeff_matrix, 0, sizeof(double) * matrix_size * matrix_size);

    right_hand_side = (double*)std::malloc(matrix_size * sizeof(double));
    if (right_hand_side == nullptr)
        return 0;
    std::memset(right_hand_side, 0, sizeof(double) * matrix_size);

    solution = (double*)std::malloc(sizeof(double) * matrix_size);
    if (solution == nullptr)
        return 0;
    std::memset(solution, 0, sizeof(double) * matrix_size);

    perm = (int*)std::malloc(sizeof(int) * full_block_count);
    if (perm == nullptr)
        return 0;
    std::memset(perm, 0, sizeof(int) * full_block_count);

    w_perm = (int*)std::malloc(sizeof(int) * block_size);
    if (w_perm == nullptr)
        return 0;

    block1 = (double*)std::malloc(sizeof(double) * block_size * block_size);
    if (block1 == nullptr)
        return 0;
    
    block2 = (double*)std::malloc(sizeof(double) * block_size * block_size);
    if (block2 == nullptr)
        return 0;

    block3 = (double*)std::malloc(sizeof(double) * block_size * block_size);
    if (block3 == nullptr)
        return 0;

    return 1;
}

void LinearSystem::free_memory()
{
    std::free(coeff_matrix);
    std::free(right_hand_side);
    std::free(solution);
    std::free(perm);
    std::free(w_perm);
    std::free(block1);
    std::free(block2);
    std::free(block3);

    coeff_matrix = nullptr;
    right_hand_side = nullptr;
    solution = nullptr;
    perm = nullptr;
    w_perm = nullptr;
    block1 = nullptr;
    block2 = nullptr;
    block3 = nullptr;
}

bool LinearSystem::init(int s, char* file_name)
{
    if (matrix_size < 1 || block_size < 1 || block_size > matrix_size || matrix_size > INT_MAX / matrix_size || s < 0 || s > 4
        || (s == 0 && file_name == nullptr)
        || coeff_matrix == nullptr || right_hand_side == nullptr || solution == nullptr || perm == nullptr
        || w_perm == nullptr || block1 == nullptr || block2 == nullptr || block3 == nullptr)
        return 0;

    bool flag = 0;
    
    // if (s == 0)
    //     flag = init_matrix_from_file(file_name);
    // else
    //     flag = init_matrix_from_formula(s);

    if (!flag)
        return 0;

    original_matrix_norm = get_matrix_norm();
    if (!std::isfinite(original_matrix_norm))
        return 0;

    // flag = init_rhs();

    if (!flag)
        return 0;
    
    flag = init_perm();

    return flag;
}

bool LinearSystem::init_perm()
{
    if (block_size < 1 || matrix_size < block_size || perm == nullptr)
        return 0;
    for (int j = 0; j < full_block_count; j++)
        perm[j] = j;
    return 1;
}

void LinearSystem::print_matrix(int r) const
{
    ::print_matrix(matrix_size, matrix_size, block_size, coeff_matrix, r);
}

void LinearSystem::print_rhs(int r) const
{
    ::print_matrix(1, matrix_size, block_size, right_hand_side, r);
}

void LinearSystem::print_solution(int r) const
{
    ::print_matrix(1, matrix_size, block_size, solution, r);
}

void LinearSystem::get_block(int i, int j, double* dest) const
{
    int height = std::min(block_size, matrix_size - i * block_size);
    int width  = std::min(block_size, matrix_size - j * block_size);
    int offset = i * block_size * matrix_size + j * block_size * height;

    std::memcpy(dest, coeff_matrix + offset, height * width * sizeof(double));
}

void LinearSystem::set_block(int i, int j, const double* src)
{
    int height = std::min(block_size, matrix_size - i * block_size);
    int width  = std::min(block_size, matrix_size - j * block_size);
    int offset = i * block_size * matrix_size + j * block_size * height;

    std::memcpy(coeff_matrix + offset, src, height * width * sizeof(double));
}

double LinearSystem::get_matrix_norm() const
{
    double norm = 0.;

    for (int i = 0; i < block_count; i++)
    {
        int height = std::min(block_size, matrix_size - i * block_size);

        for (int row = 0; row < height; row++)
        {
            double sum = 0.;

            for (int j = 0; j < block_count; j++)
            {
                int width           = std::min(block_size, matrix_size - j * block_size);
                const double* block = coeff_matrix + i * block_size * matrix_size + j * block_size * height;
                for (int col = 0; col < width; col++)
                {
                    double value = std::fabs(block[row * width + col]);
                    if (std::min(sum, value) > DBL_MAX - std::max(sum, value))
                        return HUGE_VAL; // return infinity
                    sum += value;
                }
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

    std::swap_ranges(coeff_matrix + i * matrix_size * block_size + first_column * block_size * block_size,
                     coeff_matrix + (i + 1) * matrix_size * block_size,
                     coeff_matrix + j * matrix_size * block_size + first_column * block_size * block_size);
    std::swap_ranges(right_hand_side + i * block_size,
                     right_hand_side + (i + 1) * block_size, 
                     right_hand_side + j * block_size);
}

void LinearSystem::swap_blocked_columns(int i, int j)
{
    if (i == j)
        return;
    
    // меняем квадратные блоки
    for (int row = 0; row < full_block_count; row++)
        std::swap_ranges(coeff_matrix + row * matrix_size * block_size + i * block_size * block_size,
                         coeff_matrix + row * matrix_size * block_size + (i + 1) * block_size * block_size,
                         coeff_matrix + row * matrix_size * block_size + j * block_size * block_size);

    if (remainder_size != 0)
    {
        double* tail = coeff_matrix + (full_block_count) * block_size * matrix_size;
        std::swap_ranges(tail + i * block_size * remainder_size,
                         tail + (i + 1) * block_size * remainder_size,
                         tail + j * block_size * remainder_size);
    }

    std::swap(perm[i], perm[j]);
}

double* LinearSystem::get_block_address(int i, int j) const
{
    int height = std::min(block_size, matrix_size - i * block_size);
    return (coeff_matrix + i * block_size * matrix_size + j * block_size * height);
}
