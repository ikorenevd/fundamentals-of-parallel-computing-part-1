#include "./matrix_init.h"

#include "./headers_and_defines.h"

namespace
{
    bool init_matrix_from_file(double* matrix, int matrix_size, const char* file_name)
    {
        if (matrix == 0 || file_name == 0)
            return 0;

        std::FILE* file = std::fopen(file_name, "r");
        if (file == 0)
        {
            std::fprintf(stderr, "Error: couldn't open file\n");
            return 0;
        }

        for (int i = 0; i < matrix_size; i++)
            for (int j = 0; j < matrix_size; j++)
            {
                double* value = matrix + i * matrix_size + j;
                if (std::fscanf(file, "%lf", value) != 1 || std::isnan(*value) || std::isinf(*value))
                {
                    std::fclose(file);
                    return 0;
                }
            }

        char extra;
        bool valid = std::fscanf(file, " %c", &extra) == EOF && !std::ferror(file);
        std::fclose(file);
        
        if (!valid)
            std::fprintf(stderr, "Error: while reading elements\n");
        return valid;
    }

    // "настоящие" значения в матрица, т.е. от 1 до n
    double formula_entry(int matrix_size, int matrix_formula, int row, int column)
    {
        int value;
        
        switch (matrix_formula)
        {
        case 1:
            value = matrix_size - std::max(row, column) + 1;
            return ((double)(value));
            break;
        case 2:
            if (row == column)
                return 2.;
            if (std::abs(row - column) == 1)
                return -1.;
            return 0.;
            break;
        case 3:
            if (row == column && row < matrix_size)
                return 1.;
            if (column == matrix_size)
                return (double)row;
            if (row == matrix_size)
                return (double)column;
            return 0.;
            break;
        case 4:
            return (double)1 / (double)(row + column - 1);
            break;
        }
        return 0;
    }

    bool init_matrix_from_formula(double* matrix, int matrix_size, int matrix_formula)
    {
        if (matrix == 0)
            return 0;

        for (int i = 0; i < matrix_size; i++)
            for (int j = 0; j < matrix_size; j++)
                matrix[i * matrix_size + j] = formula_entry(matrix_size, matrix_formula, i + 1, j + 1);

        return 1;
    }
}

bool init_matrix(double* matrix, int matrix_size, int matrix_formula, const char* file_name)
{
    if (matrix == 0)
        return 0;
    
    if (matrix_formula == 0 && file_name == 0)
        return 0;
    
    if (matrix_formula == 0)
        return init_matrix_from_file(matrix, matrix_size, file_name);
    else
        return init_matrix_from_formula(matrix, matrix_size, matrix_formula);

    return 0;
}
