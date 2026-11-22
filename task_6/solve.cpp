#include "./solve.h"

#include "./headers_and_defines.h"

// локальные функции
namespace
{
    bool matrix_to_almost_triangular(double* matrix, int matrix_size, double eps, double& t1)
    {
        t1 = 0.;
        if (matrix == 0 || matrix_size <= 0 || eps < 0.)
            return 0;

        clock_t start = clock();

        // приводим к почти треугольному симметричному виду, т.е. к трехдиагональному
        for (int step = 0; step + 2 < matrix_size; step++)
        {
            // print_matrix(matrix, matrix_size, matrix_size, matrix_size, "Matrix\n");

            int p = step + 1;
            double* row_p = matrix + p * matrix_size;
            double x = row_p[step]; // a_{p, step}

            for (int q = p + 1; q < matrix_size; q++)
            {
                double* row_q = matrix + q * matrix_size;
                double y = row_q[step]; // a_{q, step}
                if (std::fabs(y) <= 0.)
                    continue;

                double scale = std::max(std::fabs(x), std::fabs(y)),
                       u = x / scale,
                       v = y / scale,
                       length = std::sqrt(u * u + v * v),
                       c = u / length,
                       s = v / length;

                x = scale * length;
                row_q[step] = 0.;

                // обновляем пары элементов в нижнем треугольнике
                for (int j = p + 1; j < q; j++)
                {
                    double* row_j = matrix + j * matrix_size;
                    double jp = row_j[p],
                           qj = row_q[j];
                    row_j[p] = c * jp + s * qj;
                    row_q[j] = -s * jp + c * qj;
                }

                for (int j = q + 1; j < matrix_size; j++)
                {
                    double* row_j = matrix + j * matrix_size;
                    double jp = row_j[p],
                           jq = row_j[q];
                    row_j[p] = c * jp + s * jq;
                    row_j[q] = -s * jp + c * jq;
                }

                // выполняем двустороннее вращение блока на строках p и q и на столбцах p и q
                double pp = row_p[p],
                       qp = row_q[p],
                       qq = row_q[q],
                       cc = c * c,
                       ss = s * s,
                       cs = c * s,
                       cross = (2 * cs) * qp;

                // просто явное перемножение G A G^t
                row_p[p] = cc * pp + ss * qq + cross;
                row_q[p] = (cc - ss) * qp + (cs * qq - cs * pp);
                row_q[q] = ss * pp + cc * qq - cross;
            }

            row_p[step] = x;
        }

        // верхний треугольник можно восстановить по симметрии, но для LU нужны только диагональ и элементы ниже
        
        t1 = (double)(clock() - start) / CLOCKS_PER_SEC;
        
        return 1;
    }

    // n_-(lambda) = S(T - lambda * I), где S число перемен знака ведущих миноров
    // include_lambda задает правый предел n_-(lambda + 0) для правой границы отрезка
    int n_minus(const double* diagonal, const double* sub_diag, int matrix_size, double lambda, bool include_lambda)
    {
        double pivot_min = 4. * DBL_MIN; // матрица предварительно масштабирована
        int n_lambda = 0;
        double correction = 0.; // l_{i,i-1} * u_{i-1,i}

        for (int i = 0; i < matrix_size; i++)
        {
            double l_ii = diagonal[i] - lambda - correction;
            double l_next_i = i + 1 < matrix_size ? sub_diag[i] : 0.;

            if (std::fabs(l_ii) <= 0.)
            {
                if (std::fabs(l_next_i) > 0.)
                {
                    n_lambda++;
                    i++;
                }
                else if (include_lambda)
                    n_lambda++;

                correction = 0.;
            }
            else
            {
                // delta_i = delta_{i-1} * l_ii
                if (l_ii < 0.)
                    n_lambda++;

                l_ii = std::copysign(std::max(std::fabs(l_ii), pivot_min), l_ii); // абсолютное значение первого со знаком второго
                double u_i_next = l_next_i / l_ii;
                correction = l_next_i * u_i_next;
            }
        }

        return n_lambda;
    }

    // находим собственные значения на отрезке [a_0, b_0] в порядке возрастания
    int bisect_eigenvalues(const double* diagonal, const double* sub_diag, double* eigenvalues, int matrix_size,
                          double a_0, double b_0, double eps, int& iterations_count)
    {
        int n_a_0 = n_minus(diagonal, sub_diag, matrix_size, a_0, 0);
        int n_b_0 = n_minus(diagonal, sub_diag, matrix_size, b_0, 1);
        int n_a = n_a_0;
        double a = a_0;

        while (n_a < n_b_0)
        {
            int k = n_a + 1; // номер lambda_k в упорядоченном спектре
            double b = b_0;
            int n_b = n_b_0;

            while (b - a > eps)
            {
                iterations_count++;

                double c = a + (b - a) / 2.;

                if (c <= a || c >= b)
                    break;

                int n_c = n_minus(diagonal, sub_diag, matrix_size, c, 0);
                if (n_c < k)
                    a = c;
                else
                {
                    b = c;
                    n_b = std::min(n_c, n_b);
                }
            }

            double lambda_k = a + (b - a) / 2.;
            for (int j = n_a; j < n_b; j++)
                eigenvalues[j - n_a_0] = lambda_k;

            n_a = n_b;
            a = b;
        }

        return n_a - n_a_0;
    }

    // границы Гершгорина для matrix - origin * E
    // void gershgorin_bounds(const double* matrix, int matrix_size, double origin, double& left, double& right)
    // {
    //     left = right = 0.;
    //     double previous = 0.;
    //     for (int i = 0; i < matrix_size; i++)
    //     {
    //         double next = i + 1 < matrix_size ? std::fabs(matrix[(i + 1) * matrix_size + i]) : 0.;
    //         double radius = previous + next;
    //         double diagonal = matrix[i * matrix_size + i] - origin;

    //         left = std::min(left, diagonal - radius);
    //         right = std::max(right, diagonal + radius);
    //         previous = next;
    //     }
    // }
}

int solve(double* matrix, double* eigenvalues, int matrix_size, double eps, int& iterations_count, double& t1, double& t2, double* diagonal, double* sub_diag)
{
    t1 = t2 = 0.;
    iterations_count = 0;
    if (matrix == 0 || eigenvalues == 0 || matrix_size <= 0 || !std::isfinite(eps) || eps < 0. || diagonal == 0 || sub_diag == 0)
        return -1;

    double norm_A = compute_matrix_norm(matrix, matrix_size);
    if (!std::isfinite(norm_A))
        return -1;

    // норма матрицы буквально ноль, можно считать, что все собств. значения 0
    if (norm_A <= 0.)
    {
        for (int i = 0; i < matrix_size; i++)
            eigenvalues[i] = 0.;
        return 0;
    }

    // масштабируем матрицу, чтобы не было работы fpe
    // int exponent = 0;
    // std::frexp(norm_A, &exponent);
    // double scale = std::ldexp(1., exponent - 1);
    // for (int i = 0; i < matrix_size * matrix_size; i++)
    //     matrix[i] /= scale;

    // приводим к трехдиагональному виду
    if (!matrix_to_almost_triangular(matrix, matrix_size, eps, t1))
        return -1;

    clock_t start = clock();

    for (int i = 0; i < matrix_size; i++)
    {
        diagonal[i] = matrix[i * matrix_size + i];
        sub_diag[i] = i + 1 < matrix_size ? matrix[(i + 1) * matrix_size + i] : 0.;
    }

    // идея: сдвигаем матрицу так, чтобы весь спект находился в отрезке между -1 и 1
    // вычисляем границы относительно первого диагонального элемента
    // double origin = matrix[0];
    // double left = 0., right = 0.;
    // gershgorin_bounds(matrix, matrix_size, origin, left, right);
    // double center = (left + right) / 2.;
    // double half_width = (right - left) / 2.;
    // if (half_width <= 0.)
    // {
    //     for (int i = 0; i < matrix_size; i++)
    //         eigenvalues[i] = origin * scale;

    //     t2 = (double)(clock() - start) / CLOCKS_PER_SEC;
    //     return 0;
    // }

    // храним сдвиг и остаток округления
    // double mu = origin + center;
    // double rounded_center = mu - origin;
    // double mu_error = (origin - (mu - rounded_center)) + (center - rounded_center);

    // // масштабируем саму матрицу
    // for (int i = 0; i < matrix_size; i++)
    // {
    //     double* row = matrix + i * matrix_size;
    //     row[i] = ((row[i] - mu) - mu_error) / half_width;
    //     if (i > 0)
    //         row[i - 1] /= half_width;
    // }

    // double scaled_eps = std::min(eps, norm_A) / scale;
    // scaled_eps = std::min(scaled_eps, (b_0 - a_0) * half_width) / half_width;
    double b_0 = (1 + (double)matrix_size * DBL_EPSILON) * norm_A; // границы поиска
    double a_0 = -b_0;
    double scaled_eps = eps * norm_A;
    int found = bisect_eigenvalues(diagonal, sub_diag, eigenvalues, matrix_size, a_0, b_0, scaled_eps, iterations_count);

    t2 = (double)(clock() - start) / CLOCKS_PER_SEC;

    if (found != matrix_size)
        return -1;

    // можно еще обрубать собственные значения, кооторые больше по модулю max нормыыыы
    // приводим собтвенные значения к оригинальному виду
    // for (int i = 0; i < matrix_size; i++)
    //     eigenvalues[i] = (mu + (half_width * eigenvalues[i] + mu_error)) * scale;

    return 0;
}

void compute_residuals(const double* matrix, const double* eigenvalues, int matrix_size, double& res1, double& res2)
{
    res1 = res2 = 0;
    double matrix_norm = compute_matrix_norm(matrix, matrix_size);
    if (matrix_norm <= 0. || matrix == 0)
        return;

    double matrix_squares = 0.;
    double eigenvalue_squares = 0.;
    for (int i = 0; i < matrix_size; i++)
    {
        const double* row = matrix + i * matrix_size;
        for (int j = 0; j < matrix_size; j++)
        {
            double value = row[j] / matrix_norm;
            matrix_squares += value * value;
        }

        double eigenvalue = eigenvalues[i] / matrix_norm;
        res1 += row[i] / matrix_norm - eigenvalue;
        eigenvalue_squares += eigenvalue * eigenvalue;
    }

    res1 = std::fabs(res1);
    res2 = std::fabs(std::sqrt(matrix_squares) - std::sqrt(eigenvalue_squares));
}

double compute_matrix_norm(const double* matrix, int matrix_size)
{
    double norm = 0.;
    for (int i = 0; i < matrix_size; i++)
    {
        const double* row = matrix + i * matrix_size;
        double row_sum = 0.;
        for (int j = 0; j < matrix_size; j++)
            row_sum += std::fabs(row[j]);

        norm = std::max(norm, row_sum);
    }

    return norm;
}
