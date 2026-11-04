UNAME_S := $(shell uname -s)

ifeq ($(origin CXX),default)
	ifeq ($(UNAME_S),Darwin)
		CXX = clang++
	else
		CXX = g++
	endif
endif

COMMON_FLAGS = -O3 -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -pthread

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

MEMCHECK_THREADS ?= 1 2 3 4 5 6 7 8

.PHONY: clean memcheck

SOURCES = *.cpp
HEADERS = *.h

a.out: $(SOURCES) $(HEADERS) Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FPFLAGS) $(SOURCES) $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) -r build a.out

memcheck: $(SOURCES) $(HEADERS) Makefile
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FPFLAGS) -O3 -fsanitize=address -fno-omit-frame-pointer $(SOURCES) $(LDFLAGS) $(LDLIBS) -o build/memcheck
	@set -e; \
	for p in $(MEMCHECK_THREADS); do \
		printf '\nmemcheck: p = %s\n' "$$p"; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 40 2 $$p 0 3; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1340 12 $$p 0 1; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 20 10 $$p 0 1; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 21 10 $$p 0 1; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1000 7 $$p 0 1; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 201 40 $$p 0 1; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 2000 20 $$p 0 1; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 2000 1000 $$p 0 1; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1500 142 $$p 0 3; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 1660 300 $$p 0 2; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 100 23 $$p 0 4; \
		ASAN_OPTIONS='$(ASAN_OPTIONS)' ./build/memcheck 200 2 $$p 0 4; \
	done

zip:
	zip Korenev_ID.zip *.cpp *.h Makefile
