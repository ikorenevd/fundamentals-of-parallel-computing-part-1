#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H

int matrix_offset(int matrix_size, int block_size, int row, int column);
void get_block(const double* coeff_matrix, int matrix_size, int block_size, int block_row, int block_column, double* destination);
double get_matrix_norm(const double* coeff_matrix, int matrix_size, int block_size, int block_count);
double* get_block_address(double* coeff_matrix, int matrix_size, int block_size, int block_row, int block_column);
void set_block(double* coeff_matrix, int matrix_size, int block_size, int block_row, int block_column, const double* source);

#endif
