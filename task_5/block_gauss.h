#ifndef BLOCK_GAUSS_H
#define BLOCK_GAUSS_H

double calculate_block_norm(const double* block, int size, int first);
// 1 если нашли успешно и матрица обратима, 0 если обратное
// inverse и block должны быть выделены size * size
// column_order size
// при поиске обратной block изменяется, его нужно выделять с помощью get_block
bool invert_block(double* block, double* inverse, int size, int* column_order, double pivot_eps);
bool naive_gauss_solve_block(double* block, int size, double* right_hand_side, int* column_order, double pivot_eps);

#endif
