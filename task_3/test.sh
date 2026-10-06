#!/usr/bin/env bash
# Каждый полный запуск добавляет одну строку T1/T2 для всех тестов в CSV.
set -euo pipefail
export LC_ALL=C

task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python=${PYTHON:-python3}
profile=full
results=
program="$task_dir/a.out"
build=1

usage() {
    cat <<'EOF'
Использование: bash task_3/test.sh [параметры]
  --profile quick|full  Набор тестов (по умолчанию full, n<=4000; n=4000 при s=1..3).
  --results FILE        CSV результатов (по умолчанию task_3/results_full_4000_m1_2000.csv).
  --no-build            Запустить уже собранную программу.
  --program FILE        Проверить другой исполняемый файл без сборки.
  --help                Показать справку.
В CSV: одна строка на запуск, сведения о системе, T1/T2, -1 при отказе.
Ограничения времени нет: скрипт ждёт завершения каждого теста.
Код 1 означает провал ожиданий теста; строка результатов всё равно сохраняется.
EOF
}

while (($#)); do
    case "$1" in
        --profile|--results|--program)
            if (($# < 2)); then
                printf 'Ошибка: для %s требуется значение.\n' "$1" >&2
                exit 2
            fi
            case "$1" in
                --profile) profile=$2 ;;
                --results) results=$2 ;;
                --program) program=$2; build=0 ;;
            esac
            shift 2 ;;
        --no-build) build=0; shift ;;
        --help|-h) usage; exit 0 ;;
        *) printf 'Неизвестный параметр: %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done
if [[ "$profile" != quick && "$profile" != full ]]; then
    printf 'Ошибка: --profile должен быть quick или full.\n' >&2
    exit 2
fi
if [[ -z "$results" ]]; then
    if [[ "$profile" == full ]]; then
        results="$task_dir/results_full_4000_m1_2000.csv"
    else
        results="$task_dir/results.csv"
    fi
fi
if [[ "$program" != /* ]]; then program="$PWD/$program"; fi
if [[ "$results" != /* ]]; then results="$PWD/$results"; fi

mkdir -p "$task_dir/build/test-runs"
run_dir=$(mktemp -d "$task_dir/build/test-runs/run.XXXXXXXX")
manifest="$run_dir/cases.tsv"
measurements="$run_dir/times.tsv"
system_info="$run_dir/system.json"
timestamp=$("$python" -c 'from datetime import datetime; print(datetime.now().astimezone().isoformat(timespec="microseconds"))')
"$python" "$task_dir/tests/benchmark_support.py" system --output "$system_info"
"$python" "$task_dir/tests/benchmark_cases.py" --profile "$profile" \
    --output-dir "$run_dir/matrices" > "$manifest"
"$python" "$task_dir/tests/benchmark_support.py" check \
    --manifest "$manifest" --results "$results"
if ((build)); then make -C "$task_dir" a.out; fi
if [[ ! -x "$program" ]]; then
    printf 'Ошибка: программа не найдена или не исполняемая: %s\n' "$program" >&2
    exit 2
fi

total=0
failed=0
unsolved=0
printf 'Тесты: %s. Без ограничения времени.\n' "$profile"
while IFS=$'\t' read -r name n m r s filename expectation; do
    command=("$program" "$n" "$m" "$r" "$s")
    if [[ "$filename" != '<none>' ]]; then command+=("$filename"); fi
    values=$("$python" "$task_dir/tests/benchmark_support.py" run \
        --log "$run_dir/$name.log" --expect "$expectation" \
        --n "$n" --m "$m" --s "$s" -- "${command[@]}")
    IFS=$'\t' read -r t1 t2 status passed <<< "$values"
    printf '%s\t%s\t%s\n' "$name" "$t1" "$t2" >> "$measurements"
    total=$((total + 1))
    if [[ "$t1" == -1 ]]; then unsolved=$((unsolved + 1)); fi
    if [[ "$passed" == 1 ]]; then mark=OK; else mark=FAIL; failed=$((failed + 1)); fi
    printf '[%s] %-42s T1=%s T2=%s (%s)\n' "$mark" "$name" "$t1" "$t2" "$status"
done < "$manifest"

"$python" "$task_dir/tests/benchmark_support.py" append \
    --manifest "$manifest" --measurements "$measurements" \
    --results "$results" --timestamp "$timestamp" --system-info "$system_info"
printf '\nТестов: %d; не выполнено ожиданий: %d; результатов -1: %d.\n' "$total" "$failed" "$unsolved"
printf 'Добавлена строка: %s\nПодробный вывод тестов: %s\n' "$results" "$run_dir"
if ((failed)); then exit 1; fi
