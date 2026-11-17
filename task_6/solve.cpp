#include "./solve.h"

#include "./headers_and_defines.h"

namespace
{
    bool matrix_to_almost_triangular(double* matrix, int matrix_size, double& t1)
    {
        if (matrix == 0)
            return 0;
    
        clock_t start = clock();

        // приводим к почти треугольному виду

        t1 = (double)(clock() - start) / CLOCKS_PER_SEC;
    }
}

int solve(double* matrix, double* eigenvalues, int matrix_size, double& t1, double& t2)
{
    return 1;
}

void compute_residuals(const double* matrix, const double* eigenvalues, int matrix_size, double& res1, double& res2)
{
    res1 = res2 = 0;
    double sum_of_sq_eigenvalues = 0;
    double matrix_norm = compute_matrix_norm(matrix, matrix_size);

    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
            // res2 += matrix[i * matrix_size + j] * matrix[j * matrix_size + i];
            res2 += matrix[i * matrix_size + j] * matrix[i * matrix_size + j];

        res1 += matrix[i * matrix_size + i] - eigenvalues[i];
        sum_of_sq_eigenvalues += eigenvalues[i] * eigenvalues[i];
    }


    res1 = std::fabs(res1);
    res1 /= matrix_norm;

    res2 = std::fabs(std::sqrt(res2) - std::sqrt(sum_of_sq_eigenvalues)) / matrix_norm;
}


// todo: узнать какую норму считаем
// пока реализована евклидова норма в /mathhb{R}^{n^2}
double compute_matrix_norm(const double* matrix, int matrix_size)
{
    double sum = 0.;

    for (int i = 0; i < matrix_size * matrix_size; i++)
        sum += matrix[i];

    return std::sqrt(sum);
}