#include "./matrix_products.h"

void matrix_multiplication(const double* A, const double* B, double* result, int f, int l, int k)
{
    int i = 0;
    for (; i + 2 < f; i += 3)
    {
        const double* a0 = A + i * l;
        const double* a1 = a0 + l;
        const double* a2 = a1 + l;

        double* r0 = result + i * k;
        double* r1 = r0 + k;
        double* r2 = r1 + k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c00 = 0., c01 = 0., c02 = 0.;
            double c10 = 0., c11 = 0., c12 = 0.;
            double c20 = 0., c21 = 0., c22 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;

                double x0 = a0[u];
                double x1 = a1[u];
                double x2 = a2[u];

                double y0 = b[0];
                double y1 = b[1];
                double y2 = b[2];

                c00 += x0 * y0; c10 += x1 * y0; c20 += x2 * y0;
                c01 += x0 * y1; c11 += x1 * y1; c21 += x2 * y1;
                c02 += x0 * y2; c12 += x1 * y2; c22 += x2 * y2;
            }

            r0[j]     = c00;
            r0[j + 1] = c01;
            r0[j + 2] = c02;

            r1[j]     = c10;
            r1[j + 1] = c11;
            r1[j + 2] = c12;

            r2[j]     = c20;
            r2[j + 1] = c21;
            r2[j + 2] = c22;
        }

        // Остаточные столбцы 3x1.
        for (; j < k; j++)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                double b = B[u * k + j];
                c0 += a0[u] * b;
                c1 += a1[u] * b;
                c2 += a2[u] * b;
            }

            r0[j] = c0;
            r1[j] = c1;
            r2[j] = c2;
        }
    }

    // Остаточные строки
    for (; i < f; i++)
    {
        const double* a = A + i * l;
        double* r       = result + i * k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;
                double x = a[u];

                c0 += x * b[0];
                c1 += x * b[1];
                c2 += x * b[2];
            }

            r[j]     = c0;
            r[j + 1] = c1;
            r[j + 2] = c2;
        }

        // Остаточные элементы 1x1
        for (; j < k; j++)
        {
            double sum = 0.;
            for (int u = 0; u < l; u++)
                sum += a[u] * B[u * k + j];
            r[j] = sum;
        }
    }
}

// A = f × l, B = l × k, result -= A × B; result = f × k.
void matrix_multiplication_subtract(const double* A, const double* B, double* result, int f, int l, int k)
{
    int i = 0;
    for (; i + 2 < f; i += 3)
    {
        const double* a0 = A + i * l;
        const double* a1 = a0 + l;
        const double* a2 = a1 + l;

        double* r0 = result + i * k;
        double* r1 = r0 + k;
        double* r2 = r1 + k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c00 = 0., c01 = 0., c02 = 0.;
            double c10 = 0., c11 = 0., c12 = 0.;
            double c20 = 0., c21 = 0., c22 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;

                double x0 = a0[u];
                double x1 = a1[u];
                double x2 = a2[u];

                double y0 = b[0];
                double y1 = b[1];
                double y2 = b[2];

                c00 += x0 * y0; c10 += x1 * y0; c20 += x2 * y0;
                c01 += x0 * y1; c11 += x1 * y1; c21 += x2 * y1;
                c02 += x0 * y2; c12 += x1 * y2; c22 += x2 * y2;
            }

            c00       = r0[j] - c00;
            r0[j]     = c00;
            c01       = r0[j + 1] - c01;
            r0[j + 1] = c01;
            c02       = r0[j + 2] - c02;
            r0[j + 2] = c02;

            c10       = r1[j] - c10;
            r1[j]     = c10;
            c11       = r1[j + 1] - c11;
            r1[j + 1] = c11;
            c12       = r1[j + 2] - c12;
            r1[j + 2] = c12;

            c20       = r2[j] - c20;
            r2[j]     = c20;
            c21       = r2[j + 1] - c21;
            r2[j + 1] = c21;
            c22       = r2[j + 2] - c22;
            r2[j + 2] = c22;
        }

        // Остаточные столбцы 3x1.
        for (; j < k; j++)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                double b = B[u * k + j];
                c0 += a0[u] * b;
                c1 += a1[u] * b;
                c2 += a2[u] * b;
            }

            c0    = r0[j] - c0;
            r0[j] = c0;
            c1    = r1[j] - c1;
            r1[j] = c1;
            c2    = r2[j] - c2;
            r2[j] = c2;
        }
    }

    // Остаточные строки
    for (; i < f; i++)
    {
        const double* a = A + i * l;
        double* r       = result + i * k;

        int j = 0;
        for (; j + 2 < k; j += 3)
        {
            double c0 = 0., c1 = 0., c2 = 0.;

            for (int u = 0; u < l; u++)
            {
                const double* b = B + u * k + j;
                double x = a[u];

                c0 += x * b[0];
                c1 += x * b[1];
                c2 += x * b[2];
            }

            c0       = r[j] - c0;
            r[j]     = c0;
            c1       = r[j + 1] - c1;
            r[j + 1] = c1;
            c2       = r[j + 2] - c2;
            r[j + 2] = c2;
        }

        // Остаточные элементы 1x1
        for (; j < k; j++)
        {
            double sum = 0.;
            for (int u = 0; u < l; u++)
                sum += a[u] * B[u * k + j];
            sum  = r[j] - sum;
            r[j] = sum;
        }
    }
}
