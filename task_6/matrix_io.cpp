#include "./matrix_io.h"

#include "./headers_and_defines.h"

void print_matrix(const double* matrix, int height, int width, int print_size, const char* pre_message)
{
    if (height < 1 || width < 1 || print_size < 1 || matrix == 0)
        return;

    int printed_rows = std::min(print_size, height);
    int printed_cols = std::min(print_size, width);
    
    if ((pre_message != 0) && (printed_rows != 0 || printed_cols != 0))
        std::fprintf(stdout, "%s", pre_message);

    for (int row = 0; row < printed_rows; row++)
    {
        for (int col = 0; col < printed_cols; col++)
            std::fprintf(stdout, " %10.3e", matrix[row * width + col]);
        std::fprintf(stdout, "\n");   
    }
}