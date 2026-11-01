#include <cstdio>
#include <ctime>

#include "./LinearSystem.h"

int main(int argc, char** argv)
{
    int task = 11;
    int n = 0, m = 0, r = 0, s = 0;
    double t1 = 0., t2 = 0., r1 = -1., r2 = -1.;
    char* file_name = nullptr;
    bool solved = 0;

    if (!((argc == 5 || argc == 6)
        && (std::sscanf(argv[1], "%d", &n) == 1)
        && (std::sscanf(argv[2], "%d", &m) == 1)
        && (std::sscanf(argv[3], "%d", &r) == 1)
        && (std::sscanf(argv[4], "%d", &s) == 1)))
    {
        std::printf("Usage: %s n m r s [filename]\n", argv[0]);
        return 1;
    }

    if (argc == 6)
        file_name = argv[5];

    // проверяем на корректность
    if ((n < 1 || m < 1 || r < 0)
        || (s > 4)
        || (s == 0 && file_name == nullptr)
        || (argc == 5 && s < 0)
        || (argc == 6 && s != 0))
    {
        std::fprintf(stderr, "Wrong input\n");
        return 1;
    }

    // создание системы
    LinearSystem system(n, m);

    // выделение памяти
    if (!system.memory_alloc())
    {
        std::fprintf(stderr, "error: allocating memory\n");
        return 1;
    }

    // инициализация
    if (!system.init(s, file_name))
    {
        std::fprintf(stderr, "error: when inizialization\n");
        return 1;
    }

    // вывод матрицы и правой части(второй вывод не нужен)
    std::printf("Matrix A:\n");
    system.print_matrix(r);

    // решаем
    std::clock_t start = std::clock();
    solved = system.solve() == 0;
    t1 = (double)(std::clock() - start) / CLOCKS_PER_SEC;

    // считаем невязки
    if (solved)
    {
        std::printf("Solution x:\n");
        system.print_solution(r);
        // восстанавливаем A и b для невязки
        if (!system.init(s, file_name))
        {
            std::fprintf(stderr, "error: restoring matrix and rhs\n");
            return 1;
        }

        start = std::clock();
        system.compute_residuals(r1, r2);
        t2 = (double)(std::clock() - start) / CLOCKS_PER_SEC;
    }

    // конечный вывод
    std::printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], task, r1, r2, t1, t2, s, n, m);

    // память освобождается в деструкторе
    return 0;
}
