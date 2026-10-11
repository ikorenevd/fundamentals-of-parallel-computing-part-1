"""Regressions for scaled systems and rejection of singular matrices."""

import math
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


TASK = Path(__file__).resolve().parents[1]


class SolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["make", "-C", str(TASK), "a.out"],
                       check=True, capture_output=True, text=True)

    def check_matrix(self, filename, n, m, singular=False):
        result = subprocess.run(
            [str(TASK / "a.out"), str(n), str(m), "0", "0", str(filename)],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = re.search(r"Res1 = (\S+) Res2 = (\S+)", result.stdout)
        self.assertIsNotNone(report, result.stdout)
        residuals = tuple(map(float, report.groups()))
        if singular:
            self.assertEqual(residuals, (-1., -1.), result.stdout)
        else:
            self.assertTrue(all(math.isfinite(value) and 0 <= value <= limit
                                for value, limit in zip(residuals, (1e-8, 1e-6))),
                            result.stdout)

    def test_new_input_files_for_every_block_size(self):
        for name, n in (("a", 4), ("a20", 4), ("a40", 4), ("b", 4),
                        ("c", 6), ("d", 6), ("e", 6), ("f", 4)):
            for m in range(1, n + 2):
                with self.subTest(matrix=name, m=m):
                    self.check_matrix(TASK / "test_data" / (name + ".txt"),
                                      n, m, singular=name in {"b", "f"})

    def test_common_scale_preserves_solution_and_singularity(self):
        with tempfile.TemporaryDirectory(prefix="task3 scaling ") as directory:
            filename = Path(directory) / "scaled.txt"
            for name in ("a", "b"):
                values = list(map(float, (TASK / "test_data" / (name + ".txt"))
                                  .read_text().split()))
                for scale in (1e-20, 1e40):
                    filename.write_text(" ".join(format(value * scale, ".17g")
                                                 for value in values))
                    for m in range(1, 6):
                        with self.subTest(matrix=name, scale=scale, m=m):
                            self.check_matrix(filename, 4, m, singular=name == "b")

    def test_hilbert_12_all_blocks_and_input_sources(self):
        with tempfile.TemporaryDirectory(prefix="task3 hilbert inputs ") as directory:
            filename = Path(directory) / "matrix.txt"
            # .17g сохраняет те же значения double, которые создаёт s=4.
            filename.write_text(" ".join(format(1. / (i + j + 1), ".17g")
                                         for i in range(12) for j in range(12)))
            for m in range(1, 13):
                with self.subTest(m=m):
                    residuals = []
                    for mode in (("4",), ("0", str(filename))):
                        result = subprocess.run(
                            [str(TASK / "a.out"), "12", str(m), "0", *mode],
                            capture_output=True, text=True)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        report = re.search(r"Res1 = (\S+) Res2 = (\S+)", result.stdout)
                        self.assertIsNotNone(report, result.stdout)
                        res1, res2 = map(float, report.groups())
                        # При n=12 допустим численный отказ. При успехе обе
                        # невязки диагностические: double не гарантирует точность.
                        self.assertTrue((res1, res2) == (-1., -1.) or
                                        all(math.isfinite(value) and value >= 0
                                            for value in (res1, res2)), result.stdout)
                        residuals.append(report.groups())
                    self.assertEqual(residuals[0], residuals[1])

    def test_hilbert_numerical_rejection_preserves_scale_and_input_source(self):
        with tempfile.TemporaryDirectory(prefix="task3 hilbert rank ") as directory:
            filename = Path(directory) / "matrix.txt"
            for n in (5, 20):
                values = [1. / (i + j + 1) for i in range(n) for j in range(n)]
                for scale in (1e-20, 1., 1e40):
                    filename.write_text(" ".join(format(value * scale, ".17g")
                                                 for value in values))
                    for m in range(1, n + 2):
                        with self.subTest(n=n, m=m, scale=scale):
                            self.check_matrix(filename, n, m, singular=n == 20)
                            if scale == 1.:
                                formula = subprocess.run(
                                    [str(TASK / "a.out"), str(n), str(m), "0", "4"],
                                    capture_output=True, text=True)
                                self.assertEqual(formula.returncode, 0, formula.stderr)
                                file_result = subprocess.run(
                                    [str(TASK / "a.out"), str(n), str(m), "0", "0", str(filename)],
                                    capture_output=True, text=True)
                                pattern = r"Res1 = (\S+) Res2 = (\S+)"
                                formula_report = re.search(pattern, formula.stdout)
                                file_report = re.search(pattern, file_result.stdout)
                                self.assertIsNotNone(formula_report, formula.stdout)
                                self.assertIsNotNone(file_report, file_result.stdout)
                                self.assertEqual(formula_report.groups(), file_report.groups())

    def test_small_pivot_is_rejected_including_one_element_tail(self):
        with tempfile.TemporaryDirectory(prefix="task3 small pivot ") as directory:
            filename = Path(directory) / "diagonal.txt"
            for last_diagonal in (1., 1e-20):
                for scale in (1e-20, 1., 1e40):
                    values = (scale, 0., 0., 0., scale, 0.,
                              0., 0., last_diagonal * scale)
                    filename.write_text(" ".join(format(value, ".17g") for value in values))
                    for m in range(1, 5):
                        with self.subTest(last_diagonal=last_diagonal, scale=scale, m=m):
                            self.check_matrix(filename, 3, m, singular=last_diagonal < 1.)

    def test_underflow_with_floating_point_traps_enabled(self):
        # Эти размеры раньше давали цепочки очень малых остатков округления.
        for n, m, s in ((1000, 18, 1), (1000, 7, 1), (1340, 130, 1),
                        (1500, 142, 3), (1660, 300, 2)):
            with self.subTest(n=n, m=m, s=s):
                result = subprocess.run(
                    [str(TASK / "a.out"), str(n), str(m), "0", str(s)],
                    capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                report = re.search(r"Res1 = (\S+) Res2 = (\S+)", result.stdout)
                self.assertIsNotNone(report, result.stdout)
                residuals = tuple(map(float, report.groups()))
                self.assertTrue(all(math.isfinite(value) and 0 <= value <= limit
                                    for value, limit in zip(residuals, (1e-8, 1e-6))),
                                result.stdout)


if __name__ == "__main__":
    unittest.main()
