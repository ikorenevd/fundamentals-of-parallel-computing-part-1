#include "./linear_system_init_functions.h"
#include "./headers_and_defines.h"
#include "./block_matrix.h"
#include "./thread_args.h"

static bool init_coeff_matrix_from_file(double* coeff_matrix, int matrix_size, int block_size, const char* file_name)
{
    if (file_name == nullptr)
        return 0;

    std::FILE* file = std::fopen(file_name, "r");
    if (file == nullptr)
        return 0;

    for (int i = 0; i < matrix_size; i++)
        for (int column_index = 0; column_index < matrix_size; column_index++)
        {
            double* value = coeff_matrix + matrix_offset(matrix_size, block_size, i, column_index);

            if (std::fscanf(file, "%lf", value) != 1 || std::isnan(*value) || std::isinf(*value))
            {
                std::fclose(file);
                return 0;
            }
        }

    char extra;
    bool valid = std::fscanf(file, " %c", &extra) == EOF && !std::ferror(file);

    std::fclose(file);
    return valid;
}

static double formula_entry(int matrix_size, int matrix_formula, int row, int column)
{
    int value;
    switch (matrix_formula)
    {
    case 1:
        value = row > column ? row : column;
        return ((double)(matrix_size - value + 1));
        break;
    case 2:
        value = row > column ? row : column;
        return (double)value;
        break;
    case 3:
        value = row - column;
        return (double)(value >= 0 ? value : -value);
        break;
    case 4:
        return (double)1 / (double)(row + column - 1);
        break;
    }
    return 0;
}

static bool init_coeff_matrix_from_formula(double* coeff_matrix, int matrix_size, int block_size, int block_count, int matrix_formula, int thread_count, int thread_index)
{
    // будет идти по блокам и заполнять в каждом блоке
    for (int block_row = thread_index; block_row < block_count; block_row += thread_count)
    {
        double* block = coeff_matrix + block_row * block_size * matrix_size;
        int height = std::min(block_size, matrix_size - block_row * block_size);

        for (int block_column = 0; block_column < block_count; block_column++)
        {
            int width = std::min(block_size, matrix_size - block_column * block_size);

            for (int row = 0; row < height; row++)
                for (int column = 0; column < width; column++)
                    *block++ = formula_entry(matrix_size, matrix_formula, block_row * block_size + row + 1, block_column * block_size + column + 1);
        }
    }
    return 1;
}

static bool init_right_hand_side(double* coeff_matrix, double* right_hand_side, int matrix_size, int block_size, int block_count, int thread_count, int thread_index)
{
    // будет идти по блокам и заполнять в каждом блоке
    for (int block_row = thread_index; block_row < block_count; block_row += thread_count)
    {
        double* block = coeff_matrix + block_row * block_size * matrix_size;
        double* rhs_block = right_hand_side + block_row * block_size;
        int height = std::min(block_size, matrix_size - block_row * block_size);
        for (int row = 0; row < height; row++)
            *(rhs_block + row) = 0.;
        for (int block_column = 0; block_column < block_count; block_column++)
        {
            int width = std::min(block_size, matrix_size - block_column * block_size);

            for (int row = 0; row < height; row++)
            {
                double sum = 0;
                for (int column = 0; column < width; column++)
                {
                    int column_index = block_column * block_size + column + 1;
                    if (column_index % 2 == 1)
                    {
                        double value = *block;

                        if (!std::isfinite(value))
                        {
                            std::fprintf(stderr, "adsfaf\n");
                            return 0;
                        }

                        if ((sum > 0. && value > 0.) || (sum < 0. && value < 0.))
                        {
                            double abs_sum   = std::fabs(sum);
                            double abs_value = std::fabs(value);
                            if (std::min(abs_sum, abs_value) > DBL_MAX - std::max(abs_sum, abs_value))
                                return 0;
                        }

                        sum += value;
                    }
                    block++;
                }

                if (!std::isfinite(*(rhs_block + row) + sum))
                    return 0;

                *(rhs_block + row) += sum;
            }
        }
    }
    return 1;
}

void initialize_thread_memory(ThreadArgs* arg)
{
    // global
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int block_count = arg->shared->block_count;
    const int full_block_count = arg->shared->full_block_count;
    const int thread_count = arg->shared->thread_count;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    double* solution = arg->shared->solution;
    int* block_column_order = arg->shared->block_column_order;
    pthread_barrier_t* barrier = arg->shared->barrier;

    // own
    const int thread_index = arg->local.thread_index;
    double* pivot_row_tail = arg->local.pivot_row_tail;
    double* inverse_candidate = arg->local.inverse_candidate;
    double* working_block = arg->local.working_block;
    double* best_inverse_block = arg->local.best_inverse_block;
    int* column_order = arg->local.column_order;
    double* vector_product = arg->local.vector_product;

    // инициализируем свою память(ту часть памяти coeff_matrix, right_hand_side, solution, с которой будет работать этот поток)
    // и собственную память потока
    // идем по блочным строкам этого потока, включая неполную последнюю строку
    // pthread_barrier_wait(barrier);
    std::memset(pivot_row_tail, 0, (size_t)block_size * matrix_size * sizeof(double));
    std::memset(inverse_candidate, 0, (size_t)block_size * block_size * sizeof(double));
    std::memset(working_block, 0, (size_t)block_size * block_size * sizeof(double));
    std::memset(best_inverse_block, 0, (size_t)block_size * block_size * sizeof(double));
    std::memset(column_order, 0, block_size * sizeof(int));
    std::memset(vector_product, 0, block_size * sizeof(double));

    for (int block_row = thread_index; block_row < block_count; block_row += thread_count)
    {
        int height = std::min(block_size, matrix_size - block_row * block_size);
        std::memset(coeff_matrix + block_row * block_size * matrix_size, 0, height * matrix_size * sizeof(double));
        std::memset(right_hand_side + block_row * block_size, 0, height * sizeof(double));
        std::memset(solution + block_row * block_size, 0, height * sizeof(double));
    }

    if (thread_index == 0)
        for (int i = 0; i < full_block_count; i++)
            block_column_order[i] = i;

    pthread_barrier_wait(barrier);
}

static void collect_initialization_status(ThreadArgs* arg, int thread_count)
{
    bool everything_okay = arg->local.succeeded;

    for (int thread_index = 1; thread_index < thread_count; thread_index++)
        everything_okay = everything_okay && (arg + thread_index)->local.succeeded;

    for (int thread_index = 0; thread_index < thread_count; thread_index++)
        (arg + thread_index)->local.succeeded = everything_okay;
}

bool initialize_system(ThreadArgs* arg, bool calculate_norm)
{
    const int matrix_size = arg->shared->matrix_size;
    const int block_size = arg->shared->block_size;
    const int block_count = arg->shared->block_count;
    const int thread_count = arg->shared->thread_count;
    const int matrix_formula = arg->shared->matrix_formula;
    const char* file_name = arg->shared->file_name;
    double* coeff_matrix = arg->shared->coeff_matrix;
    double* right_hand_side = arg->shared->right_hand_side;
    pthread_barrier_t* barrier = arg->shared->barrier;
    const int thread_index = arg->local.thread_index;

    // инициализация матрицы коэффициентов и проверка все ли хорошо
    // arg->local.succeeded = 0;
    if (matrix_formula == 0)
    {
        if (thread_index == 0)
            arg->local.succeeded = init_coeff_matrix_from_file(coeff_matrix, matrix_size, block_size, file_name);
        else
            arg->local.succeeded = 1;
    }
    else
        arg->local.succeeded = init_coeff_matrix_from_formula(coeff_matrix, matrix_size, block_size, block_count, matrix_formula, thread_count, thread_index);

    pthread_barrier_wait(barrier);
    if (thread_index == 0)
        collect_initialization_status(arg, thread_count);

    pthread_barrier_wait(barrier);

    if (!arg->local.succeeded && thread_index == 0)
        std::fprintf(stderr, "Error: couldn't init coefficient matrix\n");

    if (!arg->local.succeeded)
        return 0;

    // инициализация правой части
    pthread_barrier_wait(barrier);

    arg->local.succeeded = init_right_hand_side(coeff_matrix, right_hand_side, matrix_size, block_size, block_count, thread_count, thread_index);

    pthread_barrier_wait(barrier);
    if (thread_index == 0)
    {
        collect_initialization_status(arg, thread_count);
        if (calculate_norm)
            arg->shared->original_matrix_norm = get_matrix_norm(coeff_matrix, matrix_size, block_size, block_count);
    }

    pthread_barrier_wait(barrier);

    if (!arg->local.succeeded && thread_index == 0)
        std::fprintf(stderr, "Error: couldn't init right hand side\n");

    if (!arg->local.succeeded)
        return 0;

    return 1;
}
