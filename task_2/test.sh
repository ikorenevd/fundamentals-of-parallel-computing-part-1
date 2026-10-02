#!/usr/bin/env bash
# сборка и проверка запуска программы
set -euo pipefail
export LC_ALL=C

task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
make -C "$task_dir"
mkdir -p "$task_dir/build"
"${CXX:-g++}" -std=c++11 -O3 -Wall -Wextra -Werror -pedantic \
    "$task_dir/tests/test_numerics.cpp" "$task_dir/solver.cpp" "$task_dir/matrix_io.cpp" \
    -o "$task_dir/build/test_numerics"
exec python3 "$task_dir/tests/test_cli.py"
