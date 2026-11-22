#ifndef MATRIX_INIT_H
#define MATRIX_INIT_H

// status заполняется только если что-то пошло не так
// -1 не смогли открыть файл
// -2 не смогли прочитать элемент или дочитать файл(+после элемента есть элементы)
// -3 матрица не симметрична
bool init_matrix(double* matrix, int matrix_size, int matrix_formula, const char* file_name, int& status);

#endif
