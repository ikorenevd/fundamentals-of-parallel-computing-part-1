#ifndef ARGS_H
#define ARGS_H

#include "./mem_struct.h"

class Args
{
    public:
        GENERAL_MEMORY* gen;
        THREAD_MEMORY own;

        Args() = default;
        Args(const Args&) = delete;
        Args(Args&&) = delete;
        Args& operator=(const Args&) = delete;
        Args& operator=(Args&&) = delete;
        ~Args() = default;
};

#endif