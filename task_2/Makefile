UNAME_S := $(shell uname -s)

CPPFLAGS =
COMMON_CXXFLAGS = -std=c++11 -O3 -fstack-protector-all -g \
           -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic \
           -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
           -Wmissing-format-attribute -Wwrite-strings -Wcast-align \
           -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual

ifeq ($(UNAME_S),Darwin)
CXX = clang++
# Apple Clang: -Wformat вместо -Wformat=1; без SSE для совместимости с ARM.
PLATFORM_CXXFLAGS = -Wformat
else ifeq ($(UNAME_S),Linux)
CXX = g++
# Исходный набор флагов для GCC на Linux x86/x86-64.
PLATFORM_CXXFLAGS = -mfpmath=sse -Wformat=1 -Wno-suggest-attribute=format
else
$(error Unsupported OS: $(UNAME_S))
endif

CXXFLAGS = $(COMMON_CXXFLAGS) $(PLATFORM_CXXFLAGS)
LDFLAGS =
LDLIBS =

BUILD_DIR = build
TARGET = $(BUILD_DIR)/a.out
SOURCES = main.cpp LinearSystem.cpp
OBJECTS = $(addprefix $(BUILD_DIR)/,$(SOURCES:.cpp=.o))
DEPS = $(OBJECTS:.o=.d)

.PHONY: all clean test
all: $(TARGET)

test: all $(BUILD_DIR)/tests/test_matrix
	./$(BUILD_DIR)/tests/test_matrix
	bash ./test.sh

$(BUILD_DIR)/tests/test_matrix: tests/test_matrix.cpp matrix.cpp matrix_io.cpp matrix_multiply.cpp residuals.cpp solver.cpp matrix.h solver.h Makefile
	mkdir -p $(BUILD_DIR)/tests
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) tests/test_matrix.cpp matrix.cpp matrix_io.cpp matrix_multiply.cpp residuals.cpp $(LDLIBS) -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: %.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPS)

clean:
	$(RM) -r -- "$(BUILD_DIR)"
