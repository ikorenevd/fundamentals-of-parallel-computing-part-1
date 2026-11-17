#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
import re
import shlex
import subprocess


TASK = Path(__file__).resolve().parent
DATA = TASK.parent / "test_data" / "task_6"

SPECTRA = [
    ("scalar_negative_1.txt", 1., (-7.5,)),
    ("zero_3.txt", 1., (0., 0., 0.)),
    ("diagonal_repeated_5.txt", 1., (-3., -3., 0., 7., 7.)),
    ("zero_pivot_3.txt", 1., (-math.sqrt(5.), 0., math.sqrt(5.))),
]

DENSE_SPECTRUM = (-11., -3., 5., 17.)
RANK_ONE_SPECTRUM = (0., 0., 9.)
SCALED_MATRICES = [
    ("dense_4.txt", 1., DENSE_SPECTRUM),
    ("dense_e20_4.txt", 1e20, DENSE_SPECTRUM),
    ("dense_em20_4.txt", 1e-20, DENSE_SPECTRUM),
    ("dense_e300_4.txt", 1e300, DENSE_SPECTRUM),
    ("dense_em300_4.txt", 1e-300, DENSE_SPECTRUM),
    ("dense_neg_e20_4.txt", 1e20, tuple(-value for value in reversed(DENSE_SPECTRUM))),
    ("rank_one_3.txt", 1., RANK_ONE_SPECTRUM),
    ("rank_one_e20_3.txt", 1e20, RANK_ONE_SPECTRUM),
    ("rank_one_em20_3.txt", 1e-20, RANK_ONE_SPECTRUM),
]

INVALID_FILES = [
    "empty_2.txt",
    "whitespace_2.txt",
    "short_2.txt",
    "extra_number_2.txt",
    "extra_text_2.txt",
    "text_2.txt",
    "comma_2.txt",
    "broken_exponent_2.txt",
    "nan_2.txt",
    "inf_2.txt",
    "overflow_2.txt",
    "broken_2.bin",
    "missing_2.txt",
]

INVALID_ARGUMENTS = [
    ("нечисловой размер", ["abc", "2", "1e-11", "1"], "Error:"),
    ("нулевой размер", ["0", "2", "1e-11", "1"], "Error:"),
    ("отрицательный размер печати", ["2", "-1", "1e-11", "1"], "Error:"),
    ("отрицательный eps", ["2", "2", "-1e-11", "1"], "Error:"),
    ("NaN вместо eps", ["2", "2", "nan", "1"], "Error:"),
    ("inf вместо eps", ["2", "2", "inf", "1"], "Error:"),
    ("формула вне диапазона", ["2", "2", "1e-11", "5"], "Error:"),
    ("размер с буквенным суффиксом", ["2x", "2", "1e-11", "1"], "Error:"),
    ("режим k=0 без имени файла", ["2", "2", "1e-11", "0"], "Usage:"),
    ("имя файла при k=1", ["2", "2", "1e-11", "1", str(DATA / "short_2.txt")], "Usage:"),
]


def check_spectrum(result, scale, expected):
    if result.returncode != 0 or result.stderr.strip():
        raise ValueError(f"ожидался успешный запуск, код возврата {result.returncode}")

    _, header, output = result.stdout.partition("Eigenvalues:\n")
    if not header or not output.strip():
        raise ValueError("нет строки с собственными значениями")
    values = [float(value) for value in output.splitlines()[0].split()]
    if len(values) != len(expected):
        raise ValueError(f"найдено {len(values)} собственных значений, ожидалось {len(expected)}")

    # Матрица и спектр выводятся с округлением в формате %.3e.
    for value, reference in zip(values, expected):
        if not math.isclose(value / scale, reference, rel_tol=5.1e-4, abs_tol=5e-11):
            raise ValueError(f"масштабированное собственное значение {value / scale:g}, ожидалось {reference:g}")

    residuals = re.search(r"\bResidual1 = (\S+)\s+Residual2 = (\S+)", output)
    if residuals is None:
        raise ValueError("нет итоговой строки с двумя невязками")
    if not all(0. <= float(value) <= 1e-9 for value in residuals.groups()):
        raise ValueError(f"ожидались две невязки не больше 1e-9, получено {residuals.groups()}")


def check_rejection(result, diagnostic):
    output = result.stdout + result.stderr
    if result.returncode <= 0 or diagnostic not in output:
        raise ValueError(f"ожидался отказ с сообщением {diagnostic!r}, код возврата {result.returncode}")
    if "Eigenvalues:" in output or "Residual1 =" in output:
        raise ValueError("программа вывела результат для некорректного ввода")


def main():
    parser = argparse.ArgumentParser(description="Проверка матриц и входных аргументов task_6")
    parser.add_argument("program", nargs="?", type=Path, default=TASK / "a.out",
                        help="путь к проверяемой программе (по умолчанию task_6/a.out)")
    args = parser.parse_args()
    program = args.program.resolve()
    if not program.is_file():
        parser.error(f"программа {program} не найдена, выполните make -C task_6")

    cases = []
    for name, scale, expected in SPECTRA + SCALED_MATRICES:
        n = str(len(expected))
        cases.append((name, [n, n, f"{1e-11 * scale:.16e}", "0", str(DATA / name)],
                      check_spectrum, (scale, expected)))

    for name in INVALID_FILES + ["non_symmetric_2.txt"]:
        diagnostic = "Error:" if name == "non_symmetric_2.txt" else "Error: couldn't init matrix"
        cases.append((name, ["2", "2", "1e-11", "0", str(DATA / name)],
                      check_rejection, (diagnostic,)))

    for name, arguments, diagnostic in INVALID_ARGUMENTS:
        cases.append((name, arguments, check_rejection, (diagnostic,)))

    failures = 0
    for name, arguments, check, expected in cases:
        command = [str(program), *arguments]
        result = None
        try:
            result = subprocess.run(command, capture_output=True, text=True, errors="replace", timeout=5)
            check(result, *expected)
        except (OSError, subprocess.TimeoutExpired, ValueError) as error:
            failures += 1
            print(f"FAIL {name} — {error}\n  {shlex.join(command)}")
            if result is not None:
                print((result.stdout + result.stderr).rstrip())
        else:
            print(f"OK   {name}")

    print(f"Пройдено {len(cases) - failures} из {len(cases)} проверок")
    return int(failures != 0)


if __name__ == "__main__":
    raise SystemExit(main())
