#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

class ThreadArgs;

bool forward_elimination(ThreadArgs* arg);
bool solve_remainder_block(ThreadArgs* arg);
void backward_substitution(ThreadArgs* arg);

#endif
