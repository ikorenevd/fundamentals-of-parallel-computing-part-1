#ifndef SOLVE_H
#define SOLVE_H

// возвращает 0, если все хорошо
// отрицательные, если плохо
int solve(double* matrix, double* eigenvalues, int matrix_size, double eps, int& iterations_count, double& t1, double& t2, double* diagonal, double* sub_diag);
void compute_residuals(const double* matrix, const double* eigenvalues, int matrix_size, double& res1, double& res2);

// подсчет infinity norm
double compute_matrix_norm(const double* matrix, int matrix_size);

#endif
