#ifndef SOLVER_H
#define SOLVER_H

// Максимальная сумма модулей элементов столбца квадратного блока (норма 1).
// Блок хранится построчно; некорректные данные дают HUGE_VAL.
double calculate_block_norm(const double* block, int m, double* vector_norms);

#endif
