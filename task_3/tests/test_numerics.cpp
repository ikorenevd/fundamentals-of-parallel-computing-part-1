#include "../solver.h"
#include "../matrix_io.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

// упаковка по блокам без функции адресации из программы
static std::vector<double> pack(const std::vector<double>& a, int rows, int cols, int m)
{
    std::vector<double> result;
    for (int row = 0; row < rows; row += m)
        for (int col = 0; col < cols; col += m)
            for (int i = row; i < std::min(rows, row + m); i++)
                for (int j = col; j < std::min(cols, col + m); j++)
                    result.push_back(a[i * cols + j]);
    return result;
}

int main()
{
    for (int n = 1; n <= 19; n++)
        for (int m = 1; m <= n + 2; m++)
        {
            std::vector<double> a(n * n), b(n, 0.), x(n, -1.), scratch(n + 2, -17.);
            // несимметричные матрицы с целыми и дробными элементами
            for (int variant = 0; variant < 3; variant++)
            {
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        a[i * n + j] = variant == 0 ? 0.
                            : (17. * i - 9. * j + 3.) / (variant == 1 ? 1. : i + j + 1.);
                double expected = 0.;
                for (int j = 0; j < n; j++)
                {
                    double sum = 0.;
                    for (int i = 0; i < n; i++)
                        sum += std::fabs(a[i * n + j]);
                    expected = std::max(expected, sum);
                }
                std::vector<double> blocks = pack(a, n, n, m);
                const std::vector<double> original = blocks;
                double actual = calculate_matrix_norm(n, m, blocks.data());
                assert(std::isfinite(actual));
                assert(std::fabs(actual - expected) <= 1e-14 * std::max(1., expected));
                assert(std::fabs(scratch.front() + 17.) < 1e-15);
                assert(std::fabs(scratch.back() + 17.) < 1e-15);
                assert(blocks == original);
                assert(solve(n, m, blocks.data(), b.data(), x.data(), scratch.data() + 1) == 0);
                for (int i = 0; i < n; i++)
                    assert(std::fabs(x[i] - (i % 2 == 0 ? 1. : 0.)) < 1e-15);
            }
        }

    for (double invalid : {std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::max()})
        for (int m : {1, 2, 3})
        {
            std::vector<double> a(4, invalid), b(2, 0.), x(2, -7.), scratch(2);
            assert(!std::isfinite(calculate_matrix_norm(2, m, a.data())));
            assert(solve(2, m, a.data(), b.data(), x.data(), scratch.data()) != 0);
            for (double value : x)
                assert(std::fabs(value + 7.) < 1e-15);
        }

    // вывод и неполные блоки проверяются в test_cli.py
    for (int rows : {1, 3, 5})
        for (int cols : {1, 3, 5})
            for (int m : {1, 2, 4, 7})
            {
                std::vector<double> a(rows * cols);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        a[i * cols + j] = 10. * i + j + 1.;
                const std::vector<double> blocks = pack(a, rows, cols, m);
                for (int r : {0, 1, 3, 7})
                    print_matrix(rows, cols, m, blocks.data(), r);
            }
}
