#ifndef MATRIX_IO_H
#define MATRIX_IO_H

// вывод до print_size строк и столбцов блочной матрицы height × width
// вектор задаётся как 1 × n
// pre_message сообщение, которое печается перед матрицей (без \n) 
void print_matrix(const double* matrix, int height, int width, int print_size, const char* pre_message = 0);

#endif