#!/usr/bin/env python3
"""Display task_3 benchmark timings from the wide CSV written by test.sh.

The CSV header is ``timestamp,system,architecture,cpu,hostname,<case>.T1,...``.
The older timestamp-only prefix is also supported. Values are seconds;
``-1`` means that the corresponding test did not produce a valid solution.
The ``--table`` mode uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import fnmatch
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence


SYSTEM_COLUMNS = ("system", "architecture", "cpu", "hostname")


class ResultsError(ValueError):
    """An unreadable or invalid benchmark result file."""


@dataclass(frozen=True)
class ResultRun:
    timestamp: str
    values: Mapping[str, Mapping[str, float]]
    system_info: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class BenchmarkResults:
    cases: tuple[str, ...]
    runs: tuple[ResultRun, ...]


def load_results(path: str | Path) -> BenchmarkResults:
    """Read and validate the benchmark CSV, without importing matplotlib.

    Every case must have one T1 column and one T2 column, with -1 in both
    columns for a failure. Incomplete rows, duplicate columns, nonfinite
    values and invalid negative values raise ResultsError with the file
    line and column where possible.
    """
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            fcntl.flock(stream, fcntl.LOCK_SH)
            reader = csv.reader(stream, strict=True)
            try:
                header = next(reader)
            except StopIteration:
                raise ResultsError(f"Файл {path} пуст.") from None

            if not header or header[0] != "timestamp":
                raise ResultsError("Первый столбец CSV должен называться timestamp.")
            if len(header) < 3:
                raise ResultsError("В CSV нет пары столбцов <тест>.T1 и <тест>.T2.")
            if len(set(header)) != len(header):
                raise ResultsError("В заголовке CSV повторяются имена столбцов.")

            timing_start = 1
            if any(name in SYSTEM_COLUMNS for name in header[1:]):
                if tuple(header[1:5]) != SYSTEM_COLUMNS:
                    raise ResultsError(
                        "Столбцы системы должны идти сразу после timestamp "
                        "в порядке system,architecture,cpu,hostname."
                    )
                timing_start = 5
            if len(header) < timing_start + 2:
                raise ResultsError("В CSV нет пары столбцов <тест>.T1 и <тест>.T2.")

            cases: list[str] = []
            columns: list[tuple[str, str]] = []
            metrics: dict[str, set[str]] = {}
            for name in header[timing_start:]:
                case, separator, metric = name.rpartition(".")
                if not separator or not case.strip() or metric not in {"T1", "T2"}:
                    raise ResultsError(f"Неверный столбец {name!r}: ожидается <тест>.T1 или <тест>.T2.")
                if case not in metrics:
                    metrics[case] = set()
                    cases.append(case)
                metrics[case].add(metric)
                columns.append((case, metric))
            incomplete = [case for case in cases if metrics[case] != {"T1", "T2"}]
            if incomplete:
                raise ResultsError("Нет полной пары T1/T2 для тестов: " + ", ".join(incomplete))

            runs: list[ResultRun] = []
            for row in reader:
                line = reader.line_num
                if len(row) != len(header):
                    raise ResultsError(
                        f"Строка {line}: найдено {len(row)} столбцов вместо {len(header)}."
                    )
                timestamp = row[0].strip()
                if not timestamp:
                    raise ResultsError(f"Строка {line}: пустая метка времени.")
                system_info: dict[str, str] = {}
                if timing_start == 5:
                    for column, raw_value in zip(SYSTEM_COLUMNS, row[1:5]):
                        value = raw_value.strip()
                        if not value:
                            raise ResultsError(f"Строка {line}, {column}: пустое описание системы.")
                        system_info[column] = value
                values: dict[str, dict[str, float]] = {case: {} for case in cases}
                for column, (case, metric), raw_value in zip(header[timing_start:], columns, row[timing_start:]):
                    try:
                        value = float(raw_value)
                    except ValueError:
                        raise ResultsError(
                            f"Строка {line}, {column}: {raw_value!r} не является числом."
                        ) from None
                    if not math.isfinite(value) or (value < 0 and value != -1):
                        raise ResultsError(
                            f"Строка {line}, {column}: ожидаются секунды >= 0 или -1, получено {raw_value!r}."
                        )
                    values[case][metric] = value
                for case in cases:
                    if (values[case]["T1"] == -1) != (values[case]["T2"] == -1):
                        raise ResultsError(
                            f"Строка {line}, {case}: для неудачного теста T1 и T2 должны оба равняться -1."
                        )
                runs.append(ResultRun(timestamp, values, system_info))
            if not runs:
                raise ResultsError(f"В файле {path} есть заголовок, но нет результатов запусков.")
    except ResultsError:
        raise
    except (OSError, UnicodeError, csv.Error) as error:
        raise ResultsError(f"Не удалось прочитать {path}: {error}") from None
    return BenchmarkResults(tuple(cases), tuple(runs))


def select_cases(cases: Sequence[str], patterns: Sequence[str]) -> tuple[str, ...]:
    """Select cases in CSV order by substring or shell wildcard patterns."""
    if not patterns:
        return tuple(cases)

    def matches(case: str, pattern: str) -> bool:
        if any(character in pattern for character in "*?["):
            return fnmatch.fnmatchcase(case, pattern)
        return pattern in case

    selected = tuple(case for case in cases if any(matches(case, pattern) for pattern in patterns))
    if not selected:
        raise ResultsError("Нет тестов, соответствующих --case: " + ", ".join(patterns))
    return selected


def format_timing(value: float) -> str:
    if value == -1:
        return "-1"
    return f"{value:.6g}"


def describe_system(run: ResultRun, include_details: bool = True) -> str:
    """Describe each run's machine, including unknown legacy metadata."""
    if not run.system_info:
        return "неизвестно"
    columns = SYSTEM_COLUMNS if include_details else SYSTEM_COLUMNS[:2]
    values = []
    for column in columns:
        value = run.system_info.get(column, "").strip()
        values.append(value if value and value.lower() != "unknown" else "неизвестно")
    return " / ".join(values) if any(value != "неизвестно" for value in values) else "неизвестно"


def print_table(runs: Sequence[ResultRun], cases: Sequence[str], metric: str) -> None:
    """Print one CSV run per table row, with seconds and failure counts."""
    header = ["Запуск", "Время запуска", "Система / архитектура / CPU / хост", *cases, "Ошибок"]
    rows = [
        [
            str(index),
            run.timestamp,
            describe_system(run),
            *(format_timing(run.values[case][metric]) for case in cases),
            str(sum(run.values[case][metric] == -1 for case in cases)),
        ]
        for index, run in enumerate(runs, start=1)
    ]
    widths = [max(len(row[column]) for row in [header, *rows]) for column in range(len(header))]
    print(f"{metric}, секунды; -1 = тест не решён / ошибка.")
    print(" | ".join(value.ljust(width) for value, width in zip(header, widths)))
    print("-+-".join("-" * width for width in widths))
    for row in rows:
        print(" | ".join(value.ljust(width) for value, width in zip(row, widths)))


def print_summary(runs: Sequence[ResultRun], cases: Sequence[str], metric: str) -> None:
    print(f"Запусков: {len(runs)}, тестов: {len(cases)}, показатель: {metric}.")
    first = max(0, len(runs) - 10)
    if first:
        print("Краткая сводка последних 10 запусков; график включает все выбранные запуски.")
    for index, run in enumerate(runs[first:], start=first + 1):
        values = [run.values[case][metric] for case in cases]
        solved = [value for value in values if value != -1]
        total = format_timing(sum(solved)) if solved else "—"
        print(
            f"{index}. {run.timestamp}: решено {len(solved)}/{len(cases)}, "
            f"ошибок {len(values) - len(solved)}, сумма {metric} = {total} с; "
            f"система / архитектура / CPU / хост: {describe_system(run)}"
        )


def plot_results(
    runs: Sequence[ResultRun], cases: Sequence[str], metric: str, save_path: Path | None
) -> None:
    """Plot one series per run; failed tests interrupt the series."""
    try:
        import matplotlib

        if save_path is not None:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError:
        raise ResultsError(
            "Для графика нужен matplotlib: python3 -m pip install matplotlib. "
            "Таблица доступна без зависимостей: python3 plot_results.py --table."
        ) from None

    figure, axis = plt.subplots(figsize=(min(24, max(10, len(cases) * 0.35)), 6.5))
    positions = list(range(len(cases)))
    maximum = max(
        (run.values[case][metric] for run in runs for case in cases if run.values[case][metric] >= 0),
        default=0.0,
    )
    failed_positions: set[int] = set()
    for index, run in enumerate(runs, start=1):
        values = [run.values[case][metric] for case in cases]
        successful = [math.nan if value == -1 else value for value in values]
        line, = axis.plot(
            positions,
            successful,
            marker="o" if len(cases) <= 80 else ".",
            markersize=4,
            linewidth=1,
            alpha=0.8,
            label=f"{index}. {run.timestamp} [{describe_system(run, include_details=False)}]",
        )
        zero_positions = [position for position, value in enumerate(values) if value == 0]
        if zero_positions:
            axis.scatter(
                zero_positions, [0] * len(zero_positions), marker="o", facecolors="white",
                edgecolors=line.get_color(), s=55, zorder=4,
            )
        failed_positions.update(position for position, value in enumerate(values) if value == -1)
    if failed_positions:
        axis.scatter(
            sorted(failed_positions), [0] * len(failed_positions), marker="x", color="crimson",
            s=65, linewidths=1.5, zorder=5,
        )

    stride = max(1, math.ceil(len(cases) / 35))
    ticks = list(range(0, len(cases), stride))
    if len(cases) - 1 not in ticks:
        if len(ticks) > 1 and len(cases) - 1 - ticks[-1] < stride / 2:
            ticks.pop()
        ticks.append(len(cases) - 1)
    axis.set_xticks(ticks)
    axis.set_xticklabels([cases[index] for index in ticks], rotation=60, ha="right", fontsize=8)
    axis.set_ylabel(f"{metric}, seconds")
    axis.set_xlabel("Test case")
    axis.set_title(f"task_3: {metric} — {len(cases)} cases, {len(runs)} runs")
    upper_limit = maximum * 1.08 if maximum > 0 else 1.0
    if not math.isfinite(upper_limit):
        upper_limit = maximum
    axis.set_ylim(-upper_limit * 0.035, upper_limit)
    axis.set_xlim(-0.6, len(cases) - 0.4)
    axis.grid(True, axis="y", alpha=0.25)

    handles, labels = axis.get_legend_handles_labels()
    # Keep a large history readable; --last narrows the history when needed.
    if len(handles) > 12:
        handles = handles[:3] + handles[-3:]
        labels = labels[:3] + labels[-3:]
    if failed_positions:
        handles.append(Line2D([], [], color="crimson", marker="x", linestyle="none"))
        labels.append("-1: unsolved / error (any run)")
    handles.append(Line2D([], [], color="gray", marker="o", markerfacecolor="white", linestyle="none"))
    labels.append("0: successful zero timing")
    axis.legend(handles, labels, loc="upper left", fontsize=8)
    if stride > 1 or len(runs) > 12:
        figure.text(0.01, 0.01, "Use --case PATTERN and --last N to inspect fewer cases or runs.", fontsize=8)
        figure.tight_layout(rect=(0, 0.04, 1, 1))
    else:
        figure.tight_layout()

    try:
        if save_path is not None:
            figure.savefig(save_path, dpi=160)
            print(f"График сохранён: {save_path.resolve()}")
        else:
            plt.show()
    except (OSError, ValueError, RuntimeError) as error:
        raise ResultsError(f"Не удалось отобразить или сохранить график: {error}") from None
    finally:
        plt.close(figure)


def positive_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("ожидается положительное целое число") from None
    if number <= 0:
        raise argparse.ArgumentTypeError("ожидается положительное целое число")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Таблица или график результатов task_3. Время в секундах, -1 означает неудачу.",
        epilog=(
            "Примеры: plot_results.py --table; "
            "plot_results.py results.csv --metric T2 --case 'random*' --last 3 --save timings.png"
        ),
    )
    parser.add_argument(
        "csv_file", nargs="?", type=Path,
        default=Path(__file__).resolve().with_name("results_full_4000_m1_2000.csv"),
        help="файл результатов (по умолчанию results_full_4000_m1_2000.csv рядом со скриптом)",
    )
    parser.add_argument("--metric", choices=("T1", "T2"), default="T1", help="показатель времени (по умолчанию T1)")
    parser.add_argument("--table", action="store_true", help="вывести таблицу без matplotlib; с --save также сохранить график")
    parser.add_argument("--case", action="append", default=[], metavar="PATTERN", help="подстрока имени или шаблон *, ?, []; можно повторять")
    parser.add_argument("--last", type=positive_integer, metavar="N", help="показать только последние N запусков")
    parser.add_argument("--save", type=Path, metavar="PATH", help="сохранить график в PNG/SVG/PDF вместо открытия окна")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        results = load_results(arguments.csv_file)
        cases = select_cases(results.cases, arguments.case)
        runs = results.runs[-arguments.last:] if arguments.last else results.runs
        if arguments.table:
            print_table(runs, cases, arguments.metric)
        else:
            print_summary(runs, cases, arguments.metric)
        if not arguments.table or arguments.save is not None:
            plot_results(runs, cases, arguments.metric, arguments.save)
    except ResultsError as error:
        print(f"Ошибка: {error}", file=sys.stderr)
        return 1
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Ошибка отображения результатов: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
