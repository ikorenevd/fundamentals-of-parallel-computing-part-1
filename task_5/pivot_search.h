#ifndef PIVOT_SEARCH_H
#define PIVOT_SEARCH_H

class ThreadArgs;

bool find_block_pivot(ThreadArgs* arg, int alpha);
void select_pivot(ThreadArgs* arg, int alpha);
void swap_pivot_columns(ThreadArgs* arg, int alpha, int pivot_column);

#endif
