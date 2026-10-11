#ifndef BLOCKS_OPS_H
#define BLOCKS_OPS_H

#define MACHINE_EPS __DBL_EPSILON__

// first исключает уже обработанные строки и столбцы из оценки нормы.
double calculate_block_norm(const double* block, int m, double* vector_norms, int first = 0);
void matrix_multiplication(const double* A, const double* B, double* result, int f, int l, int k);
void matrix_multiplication_subtract(const double* A, const double* B, double* result, int f, int l, int k);
// pivot_eps сохраняет исходный масштаб системы на всех шагах исключения.
bool invert_block(double* block, double* inverse, int m, int* block_perm, double pivot_eps);
bool naive_gauss_solve_block(double* block, int m, double* b, int* block_perm, double pivot_eps);

#endif
