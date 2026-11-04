#ifndef LINEAR_SYSTEM_INIT_FUNCTIONS_H
#define LINEAR_SYSTEM_INIT_FUNCTIONS_H

bool init_coeff_matrix_from_file(double* coeff_matrix, int matrix_size, int block_size, const char* file_name);

bool init_coeff_matrix_from_formula(double* coeff_matrix, int matrix_size, int block_size, int blocks_count, int s, int threads_count, int thread_number);

bool init_right_hand_size(double* coeff_matrix, double* right_hand_side, int matrix_size, int block_size, int blocks_count, int threads_count, int thread_number);

#endif