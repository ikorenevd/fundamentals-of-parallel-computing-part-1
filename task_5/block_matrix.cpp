#include "./block_matrix.h"
#include "./headers_and_defines.h"

int matrix_offset(int matrix_size, int block_size, int row, int column)
{
    int block_row = row / block_size;
    int block_column = column / block_size;

    int height = std::min(block_size, matrix_size - block_row * block_size);
    int width = std::min(block_size, matrix_size - block_column * block_size);

    return (block_row * matrix_size + block_column * height) * block_size + (row % block_size) * width + column % block_size;
}

void get_block(const double* coeff_matrix, int matrix_size, int block_size, int block_row, int block_column, double* destination)
{
    int height = std::min(block_size, matrix_size - block_row * block_size);
    int width  = std::min(block_size, matrix_size - block_column * block_size);
    int offset = block_row * block_size * matrix_size + block_column * block_size * height;

    std::memcpy(destination, coeff_matrix + offset, height * width * sizeof(double));
}

double get_matrix_norm(const double* coeff_matrix, int matrix_size, int block_size, int block_count)
{
    double norm = 0.;

    for (int block_column = 0; block_column < block_count; block_column++)
    {
        int width = std::min(block_size, matrix_size - block_column * block_size);

        for (int column = 0; column < width; column++)
        {
            double sum = 0.;

            for (int block_row = 0; block_row < block_count; block_row++)
            {
                int height          = std::min(block_size, matrix_size - block_row * block_size);
                const double* block = coeff_matrix + block_row * block_size * matrix_size + block_column * block_size * height;
                for (int row = 0; row < height; row++)
                {
                    double value = std::fabs(block[row * width + column]);
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

double* get_block_address(double* coeff_matrix, int matrix_size, int block_size, int block_row, int block_column)
{
    int height = std::min(block_size, matrix_size - block_row * block_size);
    return (coeff_matrix + block_row * block_size * matrix_size + block_column * block_size * height);
}

void set_block(double* coeff_matrix, int matrix_size, int block_size, int block_row, int block_column, const double* source)
{
    int height = std::min(block_size, matrix_size - block_row * block_size);
    int width = std::min(block_size, matrix_size - block_column * block_size);
    int offset = block_row * block_size * matrix_size + block_column * block_size * height;

    memcpy(coeff_matrix + offset, source, height * width * sizeof(double));
}
