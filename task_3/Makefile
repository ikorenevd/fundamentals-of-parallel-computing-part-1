UNAME_S := $(shell uname -s)

# Respect a compiler supplied through the environment or command line.
ifeq ($(origin CXX),default)
	ifeq ($(UNAME_S),Darwin)
		CXX = clang++
	else
		CXX = g++
	endif
endif

# Use flags supported by both GCC and Clang, including on Apple Silicon.
CXXFLAGS ?= -std=c++11 -O3 -fstack-protector-all -g -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat -Wwrite-strings -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual
MEMCHECK_ARGS ?= 40 2 4 3

# LeakSanitizer is unavailable with Apple's AddressSanitizer runtime.
ifeq ($(UNAME_S),Darwin)
	ASAN_OPTIONS ?= detect_leaks=0:halt_on_error=1
else
	ASAN_OPTIONS ?= detect_leaks=1:halt_on_error=1
endif

.PHONY: clean memcheck

a.out: main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp LinearSystem.h solver.h matrix_io.h Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) -r build a.out

memcheck: main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp LinearSystem.h solver.h matrix_io.h Makefile
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -O1 -fsanitize=address -fno-omit-frame-pointer main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp $(LDFLAGS) $(LDLIBS) -o build/memcheck
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck $(MEMCHECK_ARGS)
