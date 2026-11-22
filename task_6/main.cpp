#include "./headers_and_defines.h"

#include "./matrix_init.h"
#include "./matrix_io.h"
#include "./solve.h"

#include "./floating_point_traps.h"

#ifdef __linux__
#include <sched.h>
#endif

int main(int argc, char** argv)
{
    // Метод бисекции нахождения всех собственных значений симметричной матрицы с приведением
    // ее к трехдиагональному виду методом вращений и вычислением числа перемен знака в последо-
    // вательности главных миноров с помощью LU-разложения.
    // const int task = 13;
    int matrix_size = 0;
    int print_size = 0;
    double eps = 0.;
    int matrix_formula = 0;
    int iterations_count = 0;
    const char* file_name = 0;

    double t1 = 0.;
    double t2 = 0.;
    double res1 = 0.;
    double res2 = 0.;

    double* matrix      = 0;
    double* eigenvalues = 0;

    // после привода к симметричному трехдиагонмальному виду для работы нужны будут только диагональ и элементы ниже
    double* diagonal = 0;
    double* sub_diag = 0;

    int status;
    
#ifdef ENABLE_FPE
    if (!configure_floating_point_traps())
    {
        std::fprintf(stderr, "Error: cannot enable fpe\n");
        return -1;
    }
#endif

    // считываем и проверяем данные на корректность
    if (!(argc == 5 || argc == 6))
    {
        std::fprintf(stdout, "Usage: %s n m eps k [filename]\n", argv[0]);
        return 0;
    }

    if (!((std::sscanf(argv[1], "%d",  &matrix_size)    == 1) &&
          (std::sscanf(argv[2], "%d",  &print_size)     == 1) &&
          (std::sscanf(argv[3], "%lf", &eps)            == 1) &&
          (std::sscanf(argv[4], "%d",  &matrix_formula) == 1)))
    {
        std::fprintf(stdout, "Error: cannot read arguments\n");
        return -1;   
    }

    if ((matrix_formula == 0) && (argc != 6))
    {
        std::fprintf(stdout, "Usage: %s n m eps k [filename]\n", argv[0]);
        return -1;
    }

    if ((matrix_formula != 0) && (argc != 5))
    {
        std::fprintf(stdout, "Usage: %s n m eps k [filename]\n", argv[0]);
        return -1;
    }

    if (matrix_formula == 0)
        file_name = argv[5];

    if (!(matrix_size > 0 && print_size >= 0 && std::isfinite(eps) && eps >= 0.
          && (matrix_formula >= 0 && matrix_formula <= 4)))
    {
        std::fprintf(stdout, "Error: wrong arguments\n");
        return -1;
    }

    // аллоцируем память
    try
    {
        matrix      = new double[matrix_size * matrix_size];
        eigenvalues = new double[matrix_size];
        diagonal    = new double[matrix_size];
        sub_diag    = new double[matrix_size];
    }
    catch(const std::bad_alloc&)
    {
        std::fprintf(stderr, "Error: couldn't allocate memory\n");
        if (!matrix)
            delete[] matrix;
        if (!eigenvalues)
            delete[] eigenvalues;
        if(!diagonal)
            delete[] diagonal;
        if(!sub_diag)
            delete[] sub_diag;
        return -1;
    }
    
    // инициализируем матрицу
    if (!init_matrix(matrix, matrix_size, matrix_formula, file_name, status))
    {
        switch (status)
        {
        case -1:
            std::fprintf(stderr, "Error: (matrix inizialization) - couldn't open the file: %s\n", file_name);
            break;
        case -2:
            std::fprintf(stderr, "Error: (matrix inizialization) - couldn't read elements in %s\n", file_name);
            break;
        case -3:
            std::fprintf(stderr, "Error: (matrix inizialization) - matrix in file %s is not symmetrical\n", file_name);
            break;
        default:
            std::fprintf(stderr, "Error: couldn't init matrix\n");
            break;
        }
        delete[] matrix;
        delete[] eigenvalues;
        delete[] diagonal;
        delete[] sub_diag;
        return -1;
    }

    print_matrix(matrix, matrix_size, matrix_size, print_size, "Matrix:\n");

#ifdef __linux__
    // привязка к последнему доступному ядру
    cpu_set_t available_cpus;
    CPU_ZERO(&available_cpus);
    if (sched_getaffinity(0, sizeof(available_cpus), &available_cpus) == -1)
    {
        perror("sched_getaffinity");
        return 1;
    }

    int last_cpu = CPU_SETSIZE - 1;
    while (last_cpu >= 0 && !CPU_ISSET(last_cpu, &available_cpus))
    {
        last_cpu--;
    }
    if (last_cpu < 0)
    {
        fprintf(stderr, "error: no available CPUs\n");
        return 1;
    }

    cpu_set_t cpu_mask;
    CPU_ZERO(&cpu_mask);
    CPU_SET(last_cpu, &cpu_mask);
    if (sched_setaffinity(0, sizeof(cpu_mask), &cpu_mask) == -1)
    {
        perror("sched_setaffinity");
        return 1;
    }
#endif

    // решаем задачу
    if (solve(matrix, eigenvalues, matrix_size, eps, iterations_count, t1, t2) != 0)
    {
        std::fprintf(stderr, "Error: couldn't find all eigenvalues\n");
        delete[] matrix;
        delete[] eigenvalues;
        delete[] diagonal;
        delete[] sub_diag;
        return -1;
    }

    print_matrix(eigenvalues, 1, matrix_size, print_size, "Eigenvalues:\n");

    // восстанавливаем изначальную матрицу для подсчета невязок
    if (!init_matrix(matrix, matrix_size, matrix_formula, file_name, status))
    {
        switch (status)
        {
        case -1:
            std::fprintf(stderr, "Error: (matrix reinizialization) - couldn't open the file: %s\n", file_name);
            break;
        case -2:
            std::fprintf(stderr, "Error: (matrix reinizialization) - couldn't read elements in %s\n", file_name);
            break;
        case -3:
            std::fprintf(stderr, "Error: (matrix reinizialization) - matrix in file %s is not symmetrical\n", file_name);
            break;
        default:
            std::fprintf(stderr, "Error: couldn't reinit matrix when computing residuals\n");
            break;
        }
        delete[] matrix;
        delete[] eigenvalues;
        delete[] diagonal;
        delete[] sub_diag;
        return 1;
    }

    compute_residuals(matrix, eigenvalues, matrix_size, res1, res2);

    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, iterations_count, iterations_count / matrix_size, t1, t2);

    delete[] matrix;
    delete[] eigenvalues;
    delete[] diagonal;
    delete[] sub_diag;

    return 0;
}
