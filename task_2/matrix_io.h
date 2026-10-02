#ifndef MATRIX_IO_H
#define MATRIX_IO_H

// вывод до r строк и столбцов блочной матрицы rows × cols
// блоки до m × m, внутри — по строкам; вектор задаётся как 1 × n
void print_matrix(int rows, int cols, int m, const double* data, int r);

#endif
