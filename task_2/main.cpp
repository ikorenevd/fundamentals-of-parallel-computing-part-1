#include <stdio.h>
#include <time.h>

#include "./LinearSystem.h"

int main(int argc, char** argv)
{
    int task = 11;
    int n = 0, m = 0, r = 0, s = 0;
    double t1 = 0., t2 = 0., r1 = 0., r2 = 0.;
    char* file_name = nullptr;
    bool solved = 0;

    if (!((argc == 5 || argc == 6)
        && (sscanf(argv[1], "%d", &n) == 1)
        && (sscanf(argv[2], "%d", &m) == 1)
        && (sscanf(argv[3], "%d", &r) == 1)
        && (sscanf(argv[4], "%d", &s) == 1)))
    {
        printf("Usage: %s n m r s\n", argv[0]);
        return 1;
    }

    if (argc == 6)
    {
        file_name = argv[5];
    }

    if ((n < 1 || m < 1 || r < 0)
        || (s > 4)
        || (s == 0 && file_name == nullptr)
        || (argc == 5 && s < 0)
        || (argc == 6 && s != 0)
    )
    {
        fprintf(stderr, "Wrong input\n");
        return 1;
    }

    if (m > n)
        m = n;

    // todo: привязка к последнему ядро через <schudl.h>

    // создание системы
    LinearSystem system(n, m);

    // выделение памяти
    if (!system.memory_alloc())
    {
        fprintf(stderr, "error: allocating memory\n");
        return 1;
    }

    // инициализация
    if (!system.init(s, file_name))
    {
        fprintf(stderr, "error: when inizialization\n");
        return 1;
    }

    // выводим матрицу коэффициентов и столбец неизвестных
    system.print_matrix(r);
    system.print_rhs(r);

    // todo: проверка решено или возникли вырожденные блоки
    // решаем
    t1 = clock();
    solved = system.solve();
    t1 = (t1 - clock()) / CLOCKS_PER_SEC;

    system.print_solution(r);
    
    // считаем несвязки
    if (solved)
    {
        // Решение изменяет A и b; восстанавливаем их, сохраняя найденный x.
        if (!system.init(s, file_name))
        {
            fprintf(stderr, "error: restoring matrix and right-hand side\n");
            return 1;
        }

        t2 = clock();
        system.compute_residuals(r1, r2);
        t2 = (t2 - clock()) / CLOCKS_PER_SEC;
    }

    // конечный вывод
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], task, r1, r2, t1, t2, s, n, m);

    // память выделяется только в LinearSystem::init, освобождается в деструкторе 
    return 0;
}
