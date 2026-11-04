#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

class Args;

bool find_block_pivot(Args* arg, int alpha);
void swap_block_rows(Args* arg, int i, int j);

#endif