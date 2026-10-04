#ifndef SOLVER_H
#define SOLVER_H

// A хранится по блокам, внутри блока — по строкам
// workspace: n элементов, память выделяет вызывающий код отдельно от A, b, x
// норма = максимальная сумма модулей по столбцам
// todo(?): обойсись без workspace
double calculate_matrix_norm(int n, int m, const double* A, double* workspace);

// заглушка: x = (1, 0, 1, ...)
int solve(int n, int m, double* A, double* b, double* x, double* workspace);

#endif
