#ifndef LINEAR_SYSTEM_INIT_FUNCTIONS_H
#define LINEAR_SYSTEM_INIT_FUNCTIONS_H

class ThreadArgs;

void initialize_thread_memory(ThreadArgs* arg);
bool initialize_system(ThreadArgs* arg, bool calculate_norm);

#endif
