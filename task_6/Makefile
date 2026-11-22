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

SOURCES = *.cpp
HEADERS = *.h

a.out: $(SOURCES) $(HEADERS) Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FPFLAGS) $(SOURCES) $(LDFLAGS) $(LDLIBS) -o $@

fpe_check: $(SOURCES) $(HEADERS) Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(FPFLAGS) $(SOURCES) $(LDFLAGS) $(LDLIBS) -DENABLE_FPE -o $@

clean:
	$(RM) -r build a.out fpe_check

zip:
	zip Korenev_ID.zip *.cpp *.h Makefile
