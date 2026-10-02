CXX = g++
CXXFLAGS = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
MEMCHECK_ARGS ?= 40 2 4 3

.PHONY: clean memcheck

a.out: main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp LinearSystem.h solver.h matrix_io.h Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) -r build a.out

memcheck: main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp LinearSystem.h solver.h matrix_io.h Makefile
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -O1 -fsanitize=address -fno-omit-frame-pointer main.cpp LinearSystem.cpp solver.cpp matrix_io.cpp $(LDFLAGS) $(LDLIBS) -o build/memcheck
	ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 ./build/memcheck $(MEMCHECK_ARGS)
