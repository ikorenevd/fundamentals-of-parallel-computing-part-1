UNAME_S := $(shell uname -s)

ifeq ($(origin CXX),default)
	ifeq ($(UNAME_S),Darwin)
		CXX = clang++
	else
		CXX = g++
	endif
endif

COMMON_FLAGS = -O3 -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual

ifeq ($(UNAME_S),Darwin)
	# Apple Clang uses -Wformat; SSE math and suggest-attribute are GCC-specific.
	CXXFLAGS ?= -std=c++17 $(filter-out -Wformat=1,$(COMMON_FLAGS)) -Wformat
	FPFLAGS = -ffp-model=strict -ffp-contract=off
else
	CXXFLAGS ?= -std=c++17 $(COMMON_FLAGS) -mfpmath=sse -Wno-suggest-attribute=format
	FPFLAGS = -ftrapping-math -frounding-math -ffp-contract=off
endif

ifeq ($(UNAME_S),Darwin)
	ASAN_OPTIONS ?= detect_leaks=0:halt_on_error=1
else
	ASAN_OPTIONS ?= detect_leaks=1:halt_on_error=1
endif

.PHONY: clean memcheck

SOURCES = main.cpp LinearSystem.cpp solver.cpp block_ops.cpp matrix_io.cpp residuals.cpp
HEADERS = LinearSystem.h solver.h block_ops.h matrix_io.h residuals.h

a.out: $(SOURCES) $(HEADERS) Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FPFLAGS) $(SOURCES) $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) -r build a.out

memcheck: $(SOURCES) $(HEADERS) Makefile
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FPFLAGS) -O3 -fsanitize=address -fno-omit-frame-pointer $(SOURCES) $(LDFLAGS) $(LDLIBS) -o build/memcheck
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 40 2 4 3
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1340 130 10 1
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 20 10 10 1
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 21 10 10 1
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1000 7 10 1
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 201 40 10 1
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 2000 20 10 1
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 2000 1000 10 1
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1500 142 10 3
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1660 300 10 2
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 100 23 10 4
	ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 200 2 10 4

zip:
	zip Korenev_ID.zip *.cpp *.h Makefile
