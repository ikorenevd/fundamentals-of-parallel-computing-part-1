#!/usr/bin/env python3
"""Проверка обычного запуска полного набора без больших вычислений."""

import csv
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


TASK = Path(__file__).resolve().parents[1]


class RunnerTests(unittest.TestCase):
    def test_default_profile_runs_full_catalog_and_keeps_makefile(self):
        makefile_before = (TASK / "Makefile").read_bytes()
        with tempfile.TemporaryDirectory(prefix="task3 runner with spaces ") as name:
            directory = Path(name)
            program = directory / "fake solver.sh"
            results = directory / "results with spaces.csv"
            output = directory / "runner output.txt"
            program.write_text(
                "#!/usr/bin/env bash\n"
                "printf '%s : Task = 11 Res1 = -1 Res2 = -1 T1 = 0.00 "
                "T2 = 0.00 S = %s N = %s M = %s\\n' \"$0\" \"$4\" \"$1\" \"$2\"\n",
                encoding="utf-8")
            program.chmod(0o755)
            environment = dict(os.environ, PYTHON=sys.executable,
                               PYTHONDONTWRITEBYTECODE="1")
            with output.open("w", encoding="utf-8") as target:
                completed = subprocess.run(
                    ["bash", str(TASK / "test.sh"), "--program", str(program),
                     "--results", str(results)],
                    cwd=directory, env=environment, stdout=target,
                    stderr=subprocess.STDOUT, timeout=60, check=False)
            self.assertEqual((TASK / "Makefile").read_bytes(), makefile_before)
            contents = output.read_text(encoding="utf-8")
            self.assertEqual(completed.returncode, 1, contents[-4000:])
            self.assertIn("Тесты: full.", contents)
            with results.open(encoding="utf-8", newline="") as source:
                rows = list(csv.reader(source))
            self.assertEqual(len(rows), 2, "One invocation must append one data row")
            header, measurements = rows
            self.assertEqual(header[:5],
                             ["timestamp", "system", "architecture", "cpu", "hostname"])
            self.assertEqual(len(header), 5 + 535 * 2)
            self.assertEqual(len(measurements), len(header))
            self.assertTrue(all(value.strip() for value in measurements[:5]))
            self.assertEqual(measurements[5:], ["-1"] * (535 * 2))
            t1_columns, t2_columns = header[5::2], header[6::2]
            self.assertTrue(all(column.endswith(".T1") for column in t1_columns))
            self.assertEqual(t2_columns,
                             [column.removesuffix(".T1") + ".T2"
                              for column in t1_columns])
            case_names = [column.removesuffix(".T1") for column in t1_columns]
            self.assertEqual(len(set(case_names)), 535)
            for s in range(1, 4):
                with self.subTest(s=s):
                    self.assertNotIn(f"benchmark_formula_s{s}_n4000_m1", case_names)
                    self.assertIn(f"benchmark_formula_s{s}_n2000_m1", case_names)
                    for m in (8, 16, 31, 32, 64, 4000, 4001):
                        self.assertIn(f"benchmark_formula_s{s}_n4000_m{m}", case_names)
            formula_sizes = [int(match.group(1)) for case in case_names
                             if (match := re.search(r"(?:^|_)formula_s\d+_n(\d+)_m", case))]
            self.assertTrue(all(n <= 4000 for n in formula_sizes),
                            "The default catalog must exclude n > 4000")
            self.assertFalse(any("_n8000_" in case for case in case_names))
            hilbert_sizes = [int(match.group(1)) for case in case_names
                             if (match := re.search(r"(?:^|_)formula_s4_n(\d+)_m", case))]
            self.assertTrue(hilbert_sizes, "Small s=4 cases must remain in the catalog")
            self.assertTrue(all(n <= 15 for n in hilbert_sizes),
                            "The default catalog must exclude s=4 with n > 15")
            logs_line = next(line for line in contents.splitlines()
                             if line.startswith("Подробный вывод тестов: "))
            run_directory = Path(logs_line.removeprefix("Подробный вывод тестов: "))
            with (run_directory / "system.json").open(encoding="utf-8") as source:
                system_info = json.load(source)
            self.assertEqual([system_info[key] for key in header[1:5]],
                             measurements[1:5])


if __name__ == "__main__":
    unittest.main(verbosity=2)
