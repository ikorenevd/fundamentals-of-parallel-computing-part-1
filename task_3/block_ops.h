#ifndef BLOCKS_OPS_H
#define BLOCKS_OPS_H

#define MACHINE_EPS 1e-16

double calculate_block_norm(const double* block, int m, double* vector_norms);
// eps задаёт порог зануления результата, для свободных членов оставляем 0.
void matrix_multiplication(const double* A, const double* B, double* result, int f, int l, int k, double eps = 0.);
void matrix_multiplication_subtract(const double* A, const double* B, double* result, int f, int l, int k, double eps = 0.);
bool invert_block(double* block, double* inverse, int m, int* block_perm);
bool naive_gauss_solve_block(double* block, int m, double* b, int* block_perm);

#endif
