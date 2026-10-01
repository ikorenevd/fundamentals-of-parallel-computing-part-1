#!/usr/bin/env bash
# Check file input and block storage; the solver is not implemented yet.
set -euo pipefail
export LC_ALL=C

task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
program="$task_dir/build/a.out"
data_dir="$task_dir/test_data"
work_dir=$(mktemp -d "${TMPDIR:-/tmp}/task_2_tests.XXXXXX")
trap 'rm -rf -- "$work_dir"' EXIT

if ! make -C "$task_dir" all > "$work_dir/build.log" 2>&1; then
    cat "$work_dir/build.log" >&2
    exit 1
fi

failed=0
files=0
shopt -s nullglob

for file in "$data_dir"/*.txt; do
    files=$((files + 1))
    name=${file##*/}
    # Dimensions are encoded in filenames, except for the empty fixture.
    if [[ $name == invalid_empty.txt ]]; then
        n=2
    elif [[ $name =~ _([1-9][0-9]*)\.txt$ ]]; then
        n=$((10#${BASH_REMATCH[1]}))
    else
        printf 'FAIL %s: missing _N.txt dimension suffix\n' "$name"
        failed=$((failed + 1))
        continue
    fi

    invalid=false
    [[ $name != invalid_* ]] || invalid=true
    # These fixtures contain valid input despite their legacy names.
    case $name in
        invalid_inf_2.txt|invalid_separator_2.txt) invalid=false ;;
    esac
    # Overflow may be rejected during initialization or reported by the solver.
    computation_input=false
    case $name in
        invalid_overflow_*.txt|invalid_rhs_overflow_*.txt) computation_input=true ;;
    esac

    if [[ $invalid == false ]]; then
        # Reconstruct the printed matrix and b independently in row order.
        # b is the sum of columns 1, 3, 5, ... (one-based indexing).
        awk -v n="$n" '
            {
                rest = $0
                while (1) {
                    sub(/^[[:space:]]+/, "", rest)
                    if (rest == "") break
                    # Scan numbers, allowing adjacent signs such as 2-2.
                    # A sign inside an exponent belongs to the same number.
                    if (!match(rest, /^[+-]?(([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?|[iI][nN][fF]([iI][nN][iI][tT][yY])?)/)) {
                        bad = 1
                        exit
                    }
                    a[count++] = substr(rest, 1, RLENGTH)
                    rest = substr(rest, RLENGTH + 1)
                }
            }
            END {
                if (bad || count != n * n) exit 1
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++)
                        printf " %10.3e", a[i * n + j]
                    printf "\n"
                }
                for (i = 0; i < n; i++) {
                    sum = 0
                    for (j = 0; j < n; j += 2) sum += a[i * n + j]
                    printf " %10.3e", sum
                }
                printf "\n"
            }
        ' "$file" > "$work_dir/expected" || {
            printf 'FAIL %s: fixture does not contain n*n values\n' "$name"
            failed=$((failed + 1))
            continue
        }
    fi

    # Pass m unchanged to test that the program internally clamps m to n.
    for ((m = 1; m <= n + 1; m++)); do
        expected_m=$m
        if ((m > n)); then
            expected_m=$n
        fi
        status=0
        "$program" "$n" "$m" "$n" 0 "$file" \
            > "$work_dir/stdout" 2> "$work_dir/stderr" || status=$?

        unavailable=false
        minus_one='-1(\.0+)?([eE][+-]?0+)?'
        if ((status == 0 || status == 1)) &&
            grep -Eq " : Task = 11 Res1 = $minus_one Res2 = $minus_one T1 = .* S = 0 N = $n M = $expected_m$" "$work_dir/stdout"; then
            unavailable=true
        fi

        reason=
        if [[ $invalid == true ]]; then
            if [[ $computation_input == true && $unavailable == true ]]; then
                : # A complete -1/-1 report is an expected computation failure.
            elif ((status != 1)); then
                reason="expected rejection with exit 1, got $status"
            elif ! grep -q 'error:' "$work_dir/stderr"; then
                reason='missing input error diagnostic'
            elif grep -q 'Task =' "$work_dir/stdout"; then
                reason='invalid input reached the solver report'
            fi
        else
            head -n "$((n + 1))" "$work_dir/stdout" > "$work_dir/actual"
            if ((status != 0)) && [[ $unavailable == false ]]; then
                reason="expected success, got exit $status"
            elif [[ -s $work_dir/stderr && $unavailable == false ]]; then
                reason='unexpected error output'
            elif ! diff -u "$work_dir/expected" "$work_dir/actual" > "$work_dir/diff"; then
                reason='matrix or right-hand-side output differs'
            elif ! grep -Eq "Task = 11 .* S = 0 N = $n M = $expected_m$" "$work_dir/stdout"; then
                reason='missing or incorrect final report'
            fi
        fi

        if [[ -n $reason ]]; then
            printf 'FAIL %s (n=%d, m=%d): %s\n' "$name" "$n" "$m" "$reason"
            if [[ $reason == 'matrix or right-hand-side output differs' ]]; then
                cat "$work_dir/diff"
            fi
            cat "$work_dir/stdout" "$work_dir/stderr"
            failed=$((failed + 1))
        fi
    done
done

if ((files == 0)); then
    printf 'FAIL: no .txt test files found in %s\n' "$data_dir"
    failed=$((failed + 1))
fi

((failed == 0))
