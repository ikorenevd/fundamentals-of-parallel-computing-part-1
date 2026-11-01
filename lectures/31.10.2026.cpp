// args.h

#include <pthread.h>

// если используем барьер, то нужно положить сюда ссылку на него
class Args
{
    public:
        double *a = nullptr;
        double *b = nullptr;
        double *x = nullptr;

        int n = 0;
        int m = 0;
        int p = 0;
        int k = 0;
        double t = 0; // cpu time
        // pthread_barrier_t*  bar;
        char* name = 0;

        Args() = default;
        Args(const Args&) = delete;
        Args(Args&&) = default;
        Args& operator=(const Args&) = delete;
        Args& operator=(Args&&) = delete;

        ~Args() = default;
};

// main.cpp

void *thread_function(void* ptr)
{
    Args* ap = (Args*)ptr;

    // для удобвства
    double* a = ap->a;
    // и т.д.

    cpu_set_t cpu;

    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (k % n_cpus);
    CPU_SET(cpu_id, &cpu);

    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    // отображаем виртуальную память на физическу
    // нужно привязать (обратиться к своем памяти)
    int i;
    for (i = k * m; i < n; i += p * m)
    {
        int h = (i + m < n) ? m : i + m - n;

        memset(a + i * n, 0, h * n * sizeof(double)); // первое обращение к матрице, поэтому тут будет привязка страниц памяти к физической
        memset(b + i, 0, h * sizeof(double));
        // ...
    }

    // локальная память
    double* block1 = new double[m * m];
    if (name != nullptr)
    {
        double res = 0; // результат
        if (k == 0)
        {
            res = read_matrix(a, n, name);
        }
        reduce_sum(p, &res, 1); // из тетреди
        if (res < 0)
        {
            // особождаем свою память
        }
    }
    else
        init_matrix(a, n, s, k, p);
    
    init_b(a, n, b, k, p); // инит парвую часть тоже параллельно

    double t_full = get_full_time();
    double t = get_cpu_time();
    solve(a, b, x, n, m, k, p);
    t_full = get_full_time() - t_full;
    t = get_cpu_time() - t;

    ap->t = t;
    ap->t_full = t_full;

    if (k == 0)
    {
        // невязка
    }
    // освобождаем локальную память

    // ставим барьер
    reduce_sum(...);
    return 0;
}

void init_matrix(double *a, int n, int m, int s, int k, int p)
{
    int i, j, i2;

    for (i = k * m; i < n; i += p * n)
    {
        int h = (i + m < n) ? m : i + m - n;

        for (i2 = i; i2 < i + h; i2++)
            for (j = 0; j < n; j++)
                a[i2 + n + j] = f(i2, j, n, s);
    }

    // синхронизуем
    redue_sum(p);
}

// умножение тоже должно быть распараллеленым
// void matrix_multiplication(double )


int main(int argc, char** argv)
{
    double *a = 0, *b = 0, *x = 0;
    int n = 0, m = 0, k = 0, p = 0;
    Args *ap = 0;
    

    // чтение аргуменров из командной строки

    pthread_t* tid = 0;

    ap = new Args[p];
    tid = new pthread_t[p];
    // обработка ошибок
    if (ap == 0 || tid == 0)
    {
        // удалить у кого можем
    }

    a = new double[n * n];
    b = new double[n];
    x = new double[n];

    // обработка ошибок

    // барьер, если хотим
    // pthread_barrier_t bar;
    // pthread_barrier_init(&bar, 0, p);


    // инициализируем
    for (k = 0; k < p; k++)
    {
        ap[k].a = a;
        ap[k].b = b;
        //... и т д
        //ap[k].bar = &b;
    }

    // инициализируем матрицу
    // если из файла - из одного потока - остальные ждут
    // по формуле - многопоточно

    // начать считать время
    // сама работа
    for (k = 1; k < p; k++)
        if (pthread_create(tid + k, 0, thread_function, ap + k))
        {
        }
        else
        {
        }
    thread_function(ap + 0);
    // закончить

    pthread_create(tid, ...);

    for (k = 1; k < p; k++)
        pthread_join(tid[k], 0);
    
    // вывод

    return 0;
}