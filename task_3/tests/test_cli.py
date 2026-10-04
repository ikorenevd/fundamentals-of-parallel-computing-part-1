#!/usr/bin/env python3
"""Проверки запуска с заглушкой решателя."""
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

TASK = Path(__file__).resolve().parents[1]
PROGRAM = TASK / "a.out"
DATA = TASK / "test_data"
NUMBER = r"[+-]?(?:\d+\.\d+e[+-]\d+|inf|nan)"
REPORT = re.compile(
    re.escape(str(PROGRAM)) + r" : Task = 11 Res1 = (" + NUMBER
    + r") Res2 = (" + NUMBER + r") T1 = (-?\d+\.\d{2})"
    + r" T2 = (-?\d+\.\d{2}) S = (\d+) N = (\d+) M = (\d+)\n"
)


def run(*args):
    return subprocess.run([str(PROGRAM), *map(str, args)], capture_output=True,
                          text=True, timeout=10)


def line(values):
    return "".join(f" {value:10.3e}" for value in values) + "\n"


class CommandLineTests(unittest.TestCase):
    def test_norm_solver_and_rectangular_print(self):
        result = subprocess.run([str(TASK / "build/test_numerics")],
                                capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        expected = ""
        for rows in (1, 3, 5):
            for cols in (1, 3, 5):
                for _m in (1, 2, 4, 7):
                    for r in (0, 1, 3, 7):
                        expected += "".join(
                            line([10. * i + j + 1. for j in range(min(r, cols))])
                            for i in range(min(r, rows)))
        self.assertEqual(result.stdout, expected)

    def check_output(self, a, m, r, s, filename=None, finite=True):
        n = len(a)
        args = [n, m, r, s]
        if filename is not None:
            args.append(filename)
        result = run(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        limit = min(r, n)
        rhs = [sum(row[::2]) for row in a]
        expected = "".join(line(row[:limit]) for row in a[:limit])
        if limit:
            expected += line(rhs[:limit])
        if finite and limit:
            expected += line([1. if i % 2 == 0 else 0. for i in range(limit)])
        self.assertTrue(result.stdout.startswith(expected),
                        f"expected prefix: {expected!r}\nactual: {result.stdout!r}")
        report = REPORT.fullmatch(result.stdout[len(expected):])
        self.assertIsNotNone(report, result.stdout)
        res1, res2, t1, t2, actual_s, actual_n, actual_m = report.groups()
        self.assertEqual((int(actual_s), int(actual_n), int(actual_m)),
                         (s, n, m))
        for elapsed in (t1, t2):
            self.assertGreaterEqual(float(elapsed), 0.)
            self.assertFalse(elapsed.startswith("-"), elapsed)
        if finite:
            self.assertTrue(math.isfinite(float(res1)), res1)
            self.assertAlmostEqual(float(res1), 0., delta=1e-12)
            self.assertAlmostEqual(float(res2), 0., delta=1e-12)
        else:
            # при бесконечной норме или NaN метод неприменим
            self.assertEqual(float(res1), -1.)
            self.assertEqual(float(res2), -1.)
            self.assertEqual(float(t2), 0.)

    def test_formulas_and_print_limits(self):
        formulas = {
            1: lambda n, i, j: n - max(i, j),
            2: lambda n, i, j: max(i, j) + 1,
            3: lambda n, i, j: abs(i - j),
            4: lambda n, i, j: 1. / (i + j + 1),
        }
        for n in (1, 2, 5, 6):
            for s, formula in formulas.items():
                a = [[formula(n, i, j) for j in range(n)] for i in range(n)]
                for m in range(1, n + 1):
                    for r in sorted({0, 1, n - 1, n, n + 2}):
                        with self.subTest(n=n, m=m, r=r, s=s):
                            self.check_output(a, m, r, s)

    def test_valid_files(self):
        names = ("scalar_1.txt", "identity_3.txt", "dense_4.txt", "mixed_3.txt",
                 "singular_3.txt", "zero_3.txt", "matrix with spaces_2.txt",
                 "invalid_separator_2.txt", "invalid_inf_2.txt",
                 "invalid_overflow_2.txt")
        for name in names:
            path = DATA / name
            n = int(path.stem.rsplit("_", 1)[1])
            if name == "invalid_separator_2.txt":
                a = [[2., -2.], [2., 1.]]
            else:
                values = list(map(float, path.read_text().split()))
                self.assertEqual(len(values), n * n)
                a = [values[i:i + n] for i in range(0, len(values), n)]
            finite = all(math.isfinite(value) for row in a for value in row)
            for m in range(1, n + 1):
                for r in sorted({0, 1, n, n + 1}):
                    with self.subTest(file=name, m=m, r=r):
                        self.check_output(a, m, r, 0, path, finite)

    def test_generated_matrices(self):
        # несимметричные матрицы выявляют ошибки порядка элементов
        with tempfile.TemporaryDirectory(prefix="linear-system-tests-") as directory:
            path = Path(directory) / "matrix.txt"
            for n in range(1, 10):
                matrices = (
                    [[float((i + 1) * 10 - 3 * j) for j in range(n)]
                     for i in range(n)],
                    [[0. if j % 2 == 0 else float(i + j + 1)
                      for j in range(n)] for i in range(n)],
                )
                for variant, a in enumerate(matrices):
                    path.write_text("\n".join(
                        " ".join(map(str, row)) for row in a) + "\n")
                    for m in range(1, n + 2):
                        with self.subTest(n=n, m=m, variant=variant):
                            self.check_output(a, m, n, 0, path)

    def test_report_preserves_m_argument(self):
        # в отчёте сохраняются исходные S, N, M (PDF, стр. 13)
        self.check_output([[2., 1.], [1., 1.]], 3, 2, 1)

    def test_submission_executable(self):
        # make должен создавать ./a.out (mail.txt, п. 7)
        self.assertTrue((TASK / "a.out").is_file(),
                        "make must produce task_2/a.out, not only build/a.out")

    def test_invalid_files(self):
        names = ("invalid_empty.txt", "invalid_short_3.txt", "invalid_extra_2.txt",
                 "invalid_text_2.txt", "invalid_suffix_2.txt", "invalid_nan_2.txt",
                 "invalid_rhs_overflow_3.txt", "does_not_exist_2.txt")
        for name in names:
            n = 2 if name == "invalid_empty.txt" else int(Path(name).stem.rsplit("_", 1)[1])
            for m in range(1, n + 1):
                with self.subTest(file=name, m=m):
                    result = run(n, m, n, 0, DATA / name)
                    self.assertEqual(result.returncode, 1)
                    self.assertIn("error:", result.stderr)
                    self.assertEqual(result.stdout, "")

    def test_invalid_arguments(self):
        cases = [(), (3,), (3, 2, 3), (3, 2, 3, 0),
                 (0, 2, 3, 1), (-1, 2, 3, 1), (3, 0, 3, 1),
                 (3, -1, 3, 1), (3, 2, -1, 1), (3, 2, 3, -1),
                 (3, 2, 3, 5), (3, 2, 3, 1, DATA / "identity_3.txt"),
                 (3, 2, 3, 0, DATA / "identity_3.txt", "extra")]
        for index in range(4):
            args = [3, 2, 3, 1]
            args[index] = "text"
            cases.append(tuple(args))
        for args in cases:
            with self.subTest(args=args):
                result = run(*args)
                self.assertEqual(result.returncode, 1)
                self.assertTrue((result.stdout + result.stderr).strip())
                self.assertNotIn("Task =", result.stdout)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]], verbosity=2)
