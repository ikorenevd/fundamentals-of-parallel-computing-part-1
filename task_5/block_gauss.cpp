#include "./block_gauss.h"
#include "./headers_and_defines.h"

static double find_block_max(double* block, int size, int alpha, int& pivot_row, int& pivot_column)
{
    double pivot_abs = 0.;

    for (int i = alpha; i < size; i++)
    {
        double* row_address = block + i * size;
        double row_max = 0.;
        for (int j = alpha; j < size; j++)
        {
            double abs = std::fabs(*(row_address + j));
            if (abs > row_max)
                row_max = abs;
        }

        if (row_max > pivot_abs)
        {
            pivot_abs = row_max;
            pivot_row = i;
            for (int j = alpha; j < size; j++)
            {
                if (std::fabs(*(row_address + j)) >= row_max)
                {
                    pivot_column = j;
                    break;
                }
            }
        }
    }

    return pivot_abs;
}

double calculate_block_norm(const double* block, int size, int first)
{
    if (block == nullptr || size < 1 || first < 0 || first >= size)
        return HUGE_VAL;
    double max = 0.;
    for (int i = first; i < size; i++)
    {
        double sum = 0.;
        for (int j = first; j < size; j++)
        {
            double value = std::fabs(block[i * size + j]);
            if (!std::isfinite(value) || std::min(sum, value) > DBL_MAX - std::max(sum, value))
                return HUGE_VAL;
            sum += value;
        }
        max = std::max(max, sum);
    }

    return max;
}

// 1 если нашли успешно и матрица обратима, 0 если обратное
// inverse и block должны быть выделены size * size
// column_order size
// при поиске обратной block изменяется, его нужно выделять с помощью get_block
bool invert_block(double* block, double* inverse, int size, int* column_order, double pivot_eps)
{
    if (!std::isfinite(pivot_eps) || pivot_eps < 0.)
        return 0;

    // делаем inverse единичной
    std::memset(inverse, 0, sizeof(double) * size * size);
    for (int i = 0; i < size; i++)
    {
        inverse[i * size + i] = 1.;
        // перестановки столбцов
        column_order[i] = i;
    }

    for (int alpha = 0; alpha < size; alpha++)
    {
        int pivot_row      = alpha;
        int pivot_column      = alpha;
        double pivot_abs = find_block_max(block, size, alpha, pivot_row, pivot_column);

        if (!std::isfinite(pivot_abs) || pivot_abs <= pivot_eps)
            return 0;
        double norm = calculate_block_norm(block, size, alpha);
        double eps = MACHINE_EPS * norm;
        if (!std::isfinite(eps))
            return 0;

        // меняем строки местами
        if (pivot_row != alpha)
            for (int j = 0; j < size; j++)
            {
                std::swap(block[alpha * size + j], block[pivot_row * size + j]);
                std::swap(inverse[alpha * size + j], inverse[pivot_row * size + j]);
            }

        // меняем столбцы
        if (pivot_column != alpha)
        {
            for (int i = 0; i < size; i++)
                std::swap(block[i * size + alpha], block[i * size + pivot_column]);
            // индексы столбцов
            std::swap(column_order[alpha], column_order[pivot_column]);
        }

        // обнуляем
        double row_eps = MACHINE_EPS * pivot_abs;
        for (int i = alpha + 1; i < size; i++)
        {
            if (std::fabs(block[i * size + alpha]) < row_eps)
            {
                block[i * size + alpha] = 0.;
                continue;
            }
            double multiplier = block[i * size + alpha] / block[alpha * size + alpha];

            for (int j = alpha + 1; j < size; j++)
            {
                double value = block[i * size + j] - multiplier * block[alpha * size + j];
                block[i * size + j] = std::fabs(value) < eps ? 0. : value;
            }

            for (int j = 0; j < size; j++)
                inverse[i * size + j] -= multiplier * inverse[alpha * size + j];

            block[i * size + alpha] = 0.;
        }

    }

    for (int i = size - 1; i >= 0; i--)
    {
        for (int other_row = i + 1; other_row < size; other_row++)
        {
            double multiplier = block[i * size + other_row];
            if (std::fabs(multiplier) < MACHINE_EPS * std::fabs(block[i * size + i]))
                continue;

            for (int j = 0; j < size; j++)
                inverse[i * size + j] -= multiplier * inverse[other_row * size + j];
        }

        double inverse_diagonal = 1. / block[i * size + i];
        for (int j = 0; j < size; j++)
            inverse[i * size + j] *= inverse_diagonal;
    }

    for (int i = 0; i < size; )
    {
        int other_row = column_order[i];
        if (other_row == i)
        {
            i++;
            continue;
        }

        for (int j = 0; j < size; j++)
            std::swap(inverse[i * size + j], inverse[other_row * size + j]);
        std::swap(column_order[i], column_order[other_row]);
    }

    return 1;
}

bool naive_gauss_solve_block(double* block, int size, double* right_hand_side, int* column_order, double pivot_eps)
{
    if (!std::isfinite(pivot_eps) || pivot_eps < 0.)
        return 0;
    double eps = 0.;
    for (int i = 0; i < size; i++)
        column_order[i] = i;
    int alpha = 0;

    for (alpha = 0; alpha < size; alpha++)
    {
        int pivot_row = 0, pivot_column = 0;
        double pivot_abs = find_block_max(block, size, alpha, pivot_row, pivot_column);

        // если главный элемент не подходит, то выходим
        if (!std::isfinite(pivot_abs) || pivot_abs < pivot_eps)
            return 0;
        double norm = calculate_block_norm(block, size, alpha);
        eps = (MACHINE_EPS) * norm;
        if (!std::isfinite(eps))
            return 0;

        // физически меняем строки alpha и pivot_row
        if (pivot_row != alpha)
        {
            for (int j = alpha; j < size; j++)
            std::swap(block[pivot_row * size + j], block[alpha * size + j]);
            std::swap(right_hand_side[pivot_row], right_hand_side[alpha]);
        }

        // физически меняем столбцы alpha и pivot_column
        if (pivot_column != alpha)
        {
            for (int i = 0; i < size; i++)
            std::swap(block[i * size + alpha], block[i * size + pivot_column]);
            std::swap(column_order[alpha], column_order[pivot_column]);
        }

        // зануляем
        double row_eps = MACHINE_EPS * pivot_abs;
        for (int i = alpha + 1; i < size; i++)
        {
            if (std::fabs(block[i * size + alpha]) < row_eps)
            {
                block[i * size + alpha] = 0.;
                continue;
            }
            double multiplier = block[i * size + alpha] / block[alpha * size + alpha];

            for (int j = alpha + 1; j < size; j++)
            {
                double value = block[i * size + j] - multiplier * block[alpha * size + j];
                block[i * size + j] = std::fabs(value) < eps ? 0. : value;
            }

            right_hand_side[i]                -= multiplier * right_hand_side[alpha];
            block[i * size + alpha] = 0.;
        }
    }

    // обратный ход
    for (int i = size - 1; i >= 0; i--)
    {
        double diagonal = block[i * size + i];
        if (std::fabs(diagonal) < pivot_eps)
            return 0;

        double rhs = right_hand_side[i];
        for (int j = i + 1; j < size; j++)
            rhs -= block[i * size + j] * right_hand_side[j];

        double value = rhs / diagonal;
        if (!std::isfinite(value))
            return 0;

        right_hand_side[i] = value;
        // block[i * size + i] = 1.;
    }

    // Восстанавливаем исходный порядок неизвестных.
    for (int i = 0; i < size;)
    {
        int other_row = column_order[i];
        if (other_row == i)
        {
            i++;
            continue;
        }

        std::swap(right_hand_side[i], right_hand_side[other_row]);
        std::swap(column_order[i], column_order[other_row]);
    }

    return 1;
}
