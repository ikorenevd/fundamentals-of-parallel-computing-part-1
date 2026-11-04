#ifndef HELPER_H
#define HELPER_H

int matrix_offset(int n, int m, int row, int col);
double f(int n, int s, int i, int j);
void get_block(const double* coeff_matrix, int matrix_size, int block_size, int i, int j, double* dest);
// 1 если нашли успешно и матрица обратима, 0 если обратное
// inverse и block должны быть выделены m * m
// block_perm m
// при поиске обратной block изменяется, его нужно выделять с помощью get_block
double calculate_block_norm(const double* block, int m, int first);
bool invert_block(double* block, double* inverse, int m, int* block_perm, double pivot_eps);
double find_block_max(double* block, int m, int alpha, int& pivot_i, int& pivot_j);
double get_matrix_norm(const double* coeff_matrix, int matrix_size, int block_size, int block_count);

double* get_block_address(double* coeff_matrix, int matrix_size, int block_size, int i, int j);

void matrix_multiplication(const double* A, const double* B, double* result, int f, int l, int k);

void set_block(double* coeff_matrix, int matrix_size, int block_size, int i, int j, const double* src);
void matrix_multiplication_subtract(const double* A, const double* B, double* result, int f, int l, int k);

bool naive_gauss_solve_block(double* block, int m, double* b, int* block_perm, double pivot_eps);

#endif