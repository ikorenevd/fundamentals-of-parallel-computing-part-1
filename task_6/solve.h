#ifndef SOLVE_H
#define SOLVE_H

// возвращает 0, если все хорошо
// отрицательные, если плохо
int solve(double* matrix, double* eigenvalues, int matrix_size, double& t1, double& t2);

// todo: добавить проверку на деление
void compute_residuals(double* matrix, double* eigenvalues, int matrix_size, double& res1, double& res2);

double compute_matrix_norm(double* matrix, int matrix_size);

#endif
