#ifndef THREAD_ARGS_H
#define THREAD_ARGS_H

#include "./mem_struct.h"

class ThreadArgs
{
    public:
        SharedMemory* shared = 0;
        ThreadMemory local;

        ThreadArgs() = default;
        ThreadArgs(const ThreadArgs&) = delete;
        ThreadArgs(ThreadArgs&&) = delete;
        ThreadArgs& operator=(const ThreadArgs&) = delete;
        ThreadArgs& operator=(ThreadArgs&&) = delete;
        ~ThreadArgs() = default;
};

#endif
