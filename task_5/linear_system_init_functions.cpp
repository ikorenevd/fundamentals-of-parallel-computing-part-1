#include "./linear_system_init_functions.h"

#include "./headers_and_defines.h"
#include "./helper.h"

bool init_coeff_matrix_from_file(double* coeff_matrix, int matrix_size, int block_size, const char* file_name)
{
    if (file_name == nullptr)
        return 0;

    std::FILE* file = std::fopen(file_name, "r");
    if (file == nullptr)
        return 0;

    for (int i = 0; i < matrix_size; i++)
        for (int j = 0; j < matrix_size; j++)
        {
            double* value = coeff_matrix + matrix_offset(matrix_size, block_size, i, j);

            if (std::fscanf(file, "%lf", value) != 1 || std::isnan(*value) || std::isinf(*value))
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

bool init_coeff_matrix_from_formula(double* coeff_matrix, int matrix_size, int block_size, int blocks_count, int s, int threads_count, int thread_number)
{
    // будет идти по блокам и заполнять в каждом блоке
    for (int blocked_row = thread_number; blocked_row < blocks_count; blocked_row += threads_count)
    {
        double* block = coeff_matrix + blocked_row * block_size * matrix_size;
        int height = std::min(block_size, matrix_size - blocked_row * block_size);

        for (int blocked_col = 0; blocked_col < blocks_count; blocked_col++)
        {
            int width = std::min(block_size, matrix_size - blocked_col * block_size);

            for (int q = 0; q < height; q++)
                for (int p = 0; p < width; p++)
                    *block++ = f(matrix_size, s, blocked_row * block_size + q + 1, blocked_col * block_size + p + 1);
        }
    }
    return 1;
}

bool init_right_hand_size(double* coeff_matrix, double* right_hand_side, int matrix_size, int block_size, int blocks_count, int threads_count, int thread_number)
{
    // будет идти по блокам и заполнять в каждом блоке
    for (int blocked_row = thread_number; blocked_row < blocks_count; blocked_row += threads_count)
    {
        double* block = coeff_matrix + blocked_row * block_size * matrix_size;
        double* rhs_chunk = right_hand_side + blocked_row * block_size;
        int height = std::min(block_size, matrix_size - blocked_row * block_size);
        for (int u = 0; u < height; u++)
            *(rhs_chunk + u) = 0.;
        for (int blocked_col = 0; blocked_col < blocks_count; blocked_col++)
        {
            int width = std::min(block_size, matrix_size - blocked_col * block_size);

            for (int q = 0; q < height; q++)
            {
                double sum = 0;
                for (int p = 0; p < width; p++)
                {
                    int j = blocked_col * block_size + p + 1;
                    if (j % 2 == 1)
                    {
                        double value = *block;
              
                        if (!std::isfinite(value))
                        {
                            std::fprintf(stderr, "adsfaf\n");
                            return 0;
                        }
 
                        if ((sum > 0. && value > 0.) || (sum < 0. && value < 0.))
                        {
                            double abs_sum   = std::fabs(sum);
                            double abs_value = std::fabs(value);
                            if (std::min(abs_sum, abs_value) > DBL_MAX - std::max(abs_sum, abs_value))
                                return 0;
                        }

                        sum += value;
                    }
                    block++;
                }
    
                if (!std::isfinite(*(rhs_chunk + q) + sum))
                    return 0;
    
                *(rhs_chunk + q) += sum;
            }
        }
    }
    return 1;
}