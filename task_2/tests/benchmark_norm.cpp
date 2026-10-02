#include "../solver.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>
#ifdef __linux__
#include <sched.h>
#endif

double baseline(int n, int m, const double* a, double*)
{
    double norm = 0.;
    for (int col = 0; col < n; col++)
    {
        double sum = 0.;
        for (int row = 0; row < n; row++)
        {
            int bi = row / m, bj = col / m;
            int h = std::min(m, n - bi * m);
            int w = std::min(m, n - bj * m);
            int offset = bi * m * n + bj * m * h
                + (row % m) * w + col % m;
            sum += std::fabs(a[offset]);
        }
        norm = std::max(norm, sum);
    }
    return norm;
}

double block_columns(int n, int m, const double* a, double* sums)
{
    double norm = 0.;
    for (int col = 0; col < n; )
    {
        int w = std::min(m, n - col);
        std::fill(sums, sums + w, 0.);
        for (int row = 0; row < n; )
        {
            int h = std::min(m, n - row);
            const double* block = a + row * n + col * h;
            for (int p = 0; p < h; p++)
                for (int q = 0; q < w; q++)
                    sums[q] += std::fabs(block[p * w + q]);
            row += h;
        }
        for (int q = 0; q < w; q++)
            norm = std::max(norm, sums[q]);
        col += w;
    }
    return norm;
}

int main()
{
#ifdef __linux__
    cpu_set_t cpus;
    CPU_ZERO(&cpus);
    if (sched_getaffinity(0, sizeof(cpus), &cpus) != 0)
        return 1;
    int cpu = CPU_SETSIZE - 1;
    while (cpu >= 0 && !CPU_ISSET(cpu, &cpus))
        cpu--;
    if (cpu < 0)
        return 1;
    CPU_ZERO(&cpus);
    CPU_SET(cpu, &cpus);
    if (sched_setaffinity(0, sizeof(cpus), &cpus) != 0)
        return 1;
#endif
    typedef double (*Norm)(int, int, const double*, double*);
    Norm variants[] = {baseline, block_columns, calculate_matrix_norm};
    std::puts("n,m,baseline_ms,block_columns_ms,storage_order_ms");
    for (int n : {512, 2048, 2051, 4096})
        for (int m : {1, 16, 90, 256, n})
        {
            std::vector<double> a(n * n), sums(n);
            for (int k = 0; k < n * n; k++)
                a[k] = (double)((int)(k % 257) - 128);
            double reference = baseline(n, m, a.data(), sums.data());
            for (Norm norm : variants)
                if (std::fabs(norm(n, m, a.data(), sums.data()) - reference) > 0.)
                    return 2;
            std::vector<double> times[3];
            for (int repeat = 0; repeat < 7; repeat++)
                for (int k = 0; k < 3; k++)
                {
                    int variant = (k + repeat) % 3;
                    auto start = std::chrono::steady_clock::now();
                    double value = variants[variant](n, m, a.data(), sums.data());
                    auto stop = std::chrono::steady_clock::now();
                    if (std::fabs(value - reference) > 0.)
                        return 2;
                    times[variant].push_back(std::chrono::duration<double, std::milli>(stop - start).count());
                }
            for (auto& time : times)
                std::sort(time.begin(), time.end());
            std::printf("%d,%d,%.4f,%.4f,%.4f\n", n, m, times[0][3], times[1][3], times[2][3]);
        }
}
