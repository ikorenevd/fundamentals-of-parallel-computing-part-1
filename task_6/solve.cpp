#include "./solve.h"

#include "./headers_and_defines.h"

// локальные функции
namespace
{
    inline void rotate_pair(double& x, double& y, double c, double s)
    {
        double old_x = x,
               old_y = y;

        x =  c * old_x + s * old_y;
        y = -s * old_x + c * old_y;
    }

    void finish_rows(double* matrix, int matrix_size, int p, const double* cosines, const double* sines)
    {
        double* row_p = matrix + p * matrix_size;
        int j = p + 1;
        for (; j + 3 < matrix_size; j += 4)
        {
            double* rows[4];
            double x[4];
            for (int r = 0; r < 4; r++)
            {
                rows[r] = matrix + (j + r) * matrix_size;
                x[r]   = row_p[j + r];

                double* row_r = rows[r];
                for (int q = j + r + 1; q < j + 4; q++)
                    rotate_pair(x[r], row_r[q], cosines[q], sines[q]);
            }

            for (int q = j + 4; q < matrix_size; q++)
            {
                double c = cosines[q], s = sines[q];

                for (int r = 0; r < 4; r++)
                    rotate_pair(x[r], rows[r][q], c, s);
            }

            // for (int r = 0; r < 4; r++)
            //     row_p[j + r] = x[r];
            row_p[j + 0] = x[0];
            row_p[j + 1] = x[1];
            row_p[j + 2] = x[2];
            row_p[j + 3] = x[3];
        }

        for (; j < matrix_size; j++)
        {
            double* row_j = matrix + j * matrix_size;
            double x = row_p[j];

            for (int q = j + 1; q < matrix_size; q++)
                rotate_pair(x, row_j[q], cosines[q], sines[q]);
            row_p[j] = x;
        }
    }

    // todo: можно применять вначале считать коэффициенты поворота, а только потом поворачивать сразу целые блоки
    // matrix — построчный массив matrix_size * matrix_size, используется и изменяется только верхний треугольник
    // eps — абс. эпсилон, значения по модулю меньше принимаем за 0 и пропускаем
    // cosines и sines — рабочие массивы длины matrix_size для коэффициентов вращений
    // t1 — процессорное время преобразования в секундах.
    // Возвращает 1 при успехе, 0 при нулевых указателях, matrix_size <= 0 или eps < 0.
    bool matrix_to_almost_triangular(double* matrix, int matrix_size, double eps, double& t1, double* cosines, double* sines)
    {
        t1 = 0.;
        if (matrix == 0 || cosines == 0 || sines == 0 || matrix_size <= 0 || eps < 0.)
            return 0;

        clock_t start = clock();

        // приводим к почти треугольному симметричному виду, т.е. к трехдиагональному
        for (int step = 0; step + 2 < matrix_size; step++)
        {
            int p = step + 1;
            double* row_step = matrix + step * matrix_size;
            double* row_p    = matrix +    p * matrix_size;
            double x = row_step[p]; // a_{step, p}

            bool dense = true;
            for (int q = p + 1; q < matrix_size; q++)
                if (std::fabs(row_step[q]) <= eps)
                {
                    dense = false;
                    break;
                }

            for (int q = p + 1; q < matrix_size; q++)
            {
                double* row_q = matrix + q * matrix_size;
                double y = row_step[q]; // a_{step, q}
                if (std::fabs(y) <= eps)
                    continue;

                double scale = std::max(std::fabs(x), std::fabs(y)),
                       u = x / scale,
                       v = y / scale,
                       length = std::sqrt(u * u + v * v),
                       c = u / length,
                       s = v / length;

                x = scale * length;
                row_step[q] = 0.;

                if (dense)
                {
                    cosines[q] = c;
                    sines[q]   = s;
                }
                else
                    for (int j = p + 1; j < q; j++)
                        rotate_pair(row_p[j], matrix[j * matrix_size + q], c, s);

                for (int j = q + 1; j < matrix_size; j++)
                {
                    double pj = row_p[j],
                           qj = row_q[j];

                    row_p[j] =  c * pj + s * qj;
                    row_q[j] = -s * pj + c * qj;
                }

                // выполняем двустороннее вращение блока на строках p и q и на столбцах p и q
                double pp = row_p[p],
                       pq = row_p[q],
                       qq = row_q[q],
                       cc = c * c,
                       ss = s * s,
                       cs = c * s,
                       cross = (2 * cs) * pq;

                // просто явное перемножение G A G^t
                row_p[p] = cc * pp + ss * qq + cross;
                row_p[q] = (cc - ss) * pq + (cs * qq - cs * pp);
                row_q[q] = ss * pp + cc * qq - cross;
            }

            if (dense)
                finish_rows(matrix, matrix_size, p, cosines, sines);

            row_step[p] = x;
        }

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


    // приводим к трехдиагональному виду
    if(!matrix_to_almost_triangular(matrix, matrix_size, eps * norm_A, t1, diagonal, sub_diag))
        return -1;

    clock_t start = clock();

    for (int i = 0; i < matrix_size; i++)
    {
        diagonal[i] = matrix[i * matrix_size + i];
        sub_diag[i] = i + 1 < matrix_size ? matrix[(i) * matrix_size + (i + 1)] : 0.;
    }

    double b_0 = (1. + (double)matrix_size * DBL_EPSILON) * norm_A; // границы поиска
    double a_0 = -b_0;
    double scaled_eps = eps * norm_A;
    int found = bisect_eigenvalues(diagonal, sub_diag, eigenvalues, matrix_size, a_0, b_0, scaled_eps, iterations_count);

    t2 = (double)(clock() - start) / CLOCKS_PER_SEC;

    if (found != matrix_size)
        return -1;

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
