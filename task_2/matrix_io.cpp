#include "./matrix_io.h"

#include <stdio.h>

void print_matrix(int rows, int cols, int m, const double* data, int r)
{
    if (rows < 1 || cols < 1 || m < 1 || r < 1 || data == nullptr)
        return;

    int printed_rows = r < rows ? r : rows;
    int printed_cols = r < cols ? r : cols;
    
    for (int row = 0; row < printed_rows; row++)
    {
        int row_start = row / m * m;
        int height = m < rows - row_start ? m : rows - row_start;

        for (int col = 0; col < printed_cols; )
        {
            int width = m < cols - col ? m : cols - col;
            int count = width < printed_cols - col ? width : printed_cols - col;
            const double* values = data + row_start * cols
                + col * height + (row - row_start) * width;

            for (int q = 0; q < count; q++)
                printf(" %10.3e", values[q]);

            col += count;
        }
        printf("\n");
    }
}
