#!/usr/bin/env python3
"""Регрессии сбора результатов без зависимости от конкретного решателя."""

import csv
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import benchmark_support as support


SYSTEM_COLUMNS = ("system", "architecture", "cpu", "hostname")
METADATA_HEADER = ["timestamp", *SYSTEM_COLUMNS]


def report(res1="1.000000e-12", res2="2.000000e-12", t1="12.34",
           t2="56.78", n=3, m=2, s=1):
    return (f"program with spaces : Task = 11 Res1 = {res1} Res2 = {res2} "
            f"T1 = {t1} T2 = {t2} S = {s} N = {n} M = {m}\n")


class ReportTests(unittest.TestCase):
    def test_reads_last_nonempty_line_and_keeps_reported_times(self):
        stdout = ("Matrix A:\n  1.000e+00  0.000e+00\nSolution x:\n"
                  + report(t1="0.01", t2="0.02")
                  + report() + "\n \t\n")
        self.assertEqual(support.parse_report(stdout, 3, 2, 1),
                         ("12.34", "56.78", "solved"))

    def test_zero_seconds_is_a_valid_solution(self):
        self.assertEqual(support.parse_report(report(t1="0.00", t2="0.00"),
                                              3, 2, 1),
                         ("0.00", "0.00", "solved"))

    def test_original_block_size_is_preserved(self):
        self.assertEqual(support.parse_report(report(n=3, m=8), 3, 8, 1)[2],
                         "solved")

    def test_unsolved_report_is_detected_even_with_zero_exit_status(self):
        self.assertEqual(support.parse_report(report(res1="-1", res2="-1"),
                                              3, 2, 1),
                         ("-1", "-1", "unsolved"))

    def test_invalid_or_stale_reports_are_rejected(self):
        samples = [
            "", "Matrix A:\n", report() + "unrelated final message\n",
            report().replace(" T2 = 56.78", ""),
            report().replace("Task = 11", "Task = 12"),
            report(n=4), report(m=3), report(s=2),
            report(t1="-0.01"), report(t2="-0.01"),
            report(res1="-1", res2="0"), report(res1="0", res2="-1"),
            report(res1="-2", res2="-2"),
        ]
        for field in ("res1", "res2", "t1", "t2"):
            for value in ("nan", "inf", "-inf", "1e999"):
                samples.append(report(**{field: value}))
        for stdout in samples:
            with self.subTest(stdout=stdout):
                with self.assertRaises(ValueError):
                    support.parse_report(stdout, 3, 2, 1)


class ExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="benchmark with spaces ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.log = self.directory / "complete log.txt"

    def execute(self, stdout, *, stderr="", exit_code=0, expectation="any"):
        program = (f"import sys; sys.stdout.write({stdout!r}); "
                   f"sys.stderr.write({stderr!r}); sys.exit({exit_code})")
        return support.execute_case([sys.executable, "-c", program],
                                    self.log, expectation, 3, 2, 1)

    def test_timing_comes_from_report_and_log_contains_both_streams(self):
        stdout = "Matrix A:\nSolution x:\n" + report()
        result = self.execute(stdout, stderr="diagnostic text\n", expectation="solved")
        self.assertEqual(result, ("12.34", "56.78", "solved", True))
        contents = self.log.read_text()
        self.assertIn(stdout, contents)
        self.assertIn("diagnostic text\n", contents)

    def test_failure_outcomes_write_minus_one(self):
        samples = (
            (report(res1="-1", res2="-1"), 0, "unsolved", "unsolved", True),
            (report(res1="-1", res2="-1"), 0, "solved", "unsolved", False),
            ("Wrong input\n", 1, "error", "error", True),
            (report(), 1, "any", "error", False),
            ("Matrix A:\n", 0, "any", "invalid_report", False),
        )
        for stdout, code, expectation, status, passed in samples:
            with self.subTest(status=status, expectation=expectation, code=code):
                self.assertEqual(self.execute(stdout, exit_code=code,
                                              expectation=expectation),
                                 ("-1", "-1", status, passed))

    def test_solved_expectation_checks_accuracy(self):
        for stdout in (report(res1="1e-3"), report(res2="1e-3")):
            with self.subTest(stdout=stdout):
                result = self.execute(stdout, expectation="solved")
                self.assertEqual(result[2:], ("solved", False))
                self.assertTrue(self.execute(stdout, expectation="any")[3])

    def test_process_crash_is_not_an_expected_input_error(self):
        program = "import os, signal; os.kill(os.getpid(), signal.SIGTERM)"
        result = support.execute_case([sys.executable, "-c", program],
                                      self.log, "error", 3, 2, 1)
        self.assertEqual(result, ("-1", "-1", "crash", False))

    def test_waits_until_completion_and_uses_reported_times(self):
        program = ("import sys, time; print('before delay', flush=True); "
                   f"time.sleep(0.2); sys.stdout.write({report()!r})")
        started = time.monotonic()
        result = support.execute_case([sys.executable, "-c", program],
                                      self.log, "solved", 3, 2, 1)
        elapsed = time.monotonic() - started
        self.assertEqual(result, ("12.34", "56.78", "solved", True))
        self.assertGreaterEqual(elapsed, 0.15)
        self.assertIn("before delay", self.log.read_text())
        self.assertIn(report(), self.log.read_text())

    def test_executable_and_arguments_with_spaces(self):
        program = self.directory / "fake solver.py"
        argument = self.directory / "matrix with spaces.txt"
        program.write_text("import sys\n"
                           f"assert sys.argv[1] == {str(argument)!r}\n"
                           f"sys.stdout.write({report()!r})\n")
        result = support.execute_case([sys.executable, str(program), str(argument)],
                                      self.log, "solved", 3, 2, 1)
        self.assertEqual(result, ("12.34", "56.78", "solved", True))


class ResultFileTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="results with spaces ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.results = self.directory / "results.csv"
        self.measurements = self.directory / "measurements.tsv"
        self.cases = [{"name": "alpha"}, {"name": "beta"}]
        self.measurements.write_text("alpha\t1.25\t3.50\nbeta\t-1\t-1\n")
        self.system_info = dict(system="Darwin 26.0", architecture="arm64",
                                cpu="Sample, CPU", hostname="sample-host")
        collector = patch.object(support, "collect_system_info",
                                 return_value=self.system_info)
        collector.start()
        self.addCleanup(collector.stop)

    def test_appends_one_full_row_per_run_and_one_header(self):
        support.append_results(self.results, self.cases, self.measurements,
                               "2026-10-06T01:02:03+03:00")
        support.append_results(self.results, self.cases, self.measurements,
                               "2026-10-06T01:02:04+03:00")
        with self.results.open(newline="") as source:
            rows = list(csv.reader(source))
        metadata = [self.system_info[key] for key in SYSTEM_COLUMNS]
        self.assertEqual(rows, [
            METADATA_HEADER + ["alpha.T1", "alpha.T2", "beta.T1", "beta.T2"],
            ["2026-10-06T01:02:03+03:00", *metadata, "1.25", "3.50", "-1", "-1"],
            ["2026-10-06T01:02:04+03:00", *metadata, "1.25", "3.50", "-1", "-1"],
        ])
        support.check_results(self.results, self.cases)

    def test_legacy_history_is_migrated_without_changing_reported_times(self):
        original = ("timestamp,alpha.T1,alpha.T2,beta.T1,beta.T2\n"
                    "old1,0.00,3.500,-1,-1\n"
                    "old2,4.20,8.00,12.34,56.78\n").encode()
        self.results.write_bytes(original)
        support.check_results(self.results, self.cases)
        self.assertEqual(self.results.read_bytes(), original)
        support.append_results(self.results, self.cases, self.measurements, "new",
                               system_info=self.system_info)
        with self.results.open(newline="") as source:
            rows = list(csv.reader(source))
        metadata = [self.system_info[key] for key in SYSTEM_COLUMNS]
        self.assertEqual(rows, [
            METADATA_HEADER + ["alpha.T1", "alpha.T2", "beta.T1", "beta.T2"],
            ["old1", *(["unknown"] * 4), "0.00", "3.500", "-1", "-1"],
            ["old2", *(["unknown"] * 4), "4.20", "8.00", "12.34", "56.78"],
            ["new", *metadata, "1.25", "3.50", "-1", "-1"],
        ])

    def test_incompatible_header_does_not_overwrite_results(self):
        original = b"timestamp,other.T1,other.T2\nold,1,2\n"
        self.results.write_bytes(original)
        with self.assertRaises(ValueError):
            support.check_results(self.results, self.cases)
        with self.assertRaises(ValueError):
            support.append_results(self.results, self.cases, self.measurements, "new")
        self.assertEqual(self.results.read_bytes(), original)

    def test_incomplete_measurements_do_not_append_partial_row(self):
        support.append_results(self.results, self.cases, self.measurements, "first")
        original = self.results.read_bytes()
        self.measurements.write_text("alpha\t1.25\t3.50\n")
        with self.assertRaises(ValueError):
            support.append_results(self.results, self.cases, self.measurements, "second")
        self.assertEqual(self.results.read_bytes(), original)

    def test_corrupted_existing_rows_are_rejected_without_changes(self):
        header = ("timestamp,system,architecture,cpu,hostname,"
                  "alpha.T1,alpha.T2,beta.T1,beta.T2\n")
        metadata = "old,Darwin,arm64,Sample CPU,sample-host,"
        for row in (metadata + "1,2\n", metadata + "nan,2,-1,-1\n",
                    metadata + "-2,2,-1,-1\n", metadata + "1,2,-1,-1",
                    metadata.removeprefix("old") + "1,2,-1,-1\n",
                    metadata + "-1,2,-1,-1\n", metadata + "1,-1,-1,-1\n",
                    "old,Darwin,arm64,,sample-host,1,2,-1,-1\n",
                    'old,Darwin,arm64,"Sample\nCPU",sample-host,1,2,-1,-1\n'):
            with self.subTest(row=row):
                original = (header + row).encode()
                self.results.write_bytes(original)
                with self.assertRaises(ValueError):
                    support.append_results(self.results, self.cases,
                                           self.measurements, "new")
                self.assertEqual(self.results.read_bytes(), original)

    def test_invalid_system_metadata_does_not_damage_history(self):
        support.append_results(self.results, self.cases, self.measurements, "first")
        original = self.results.read_bytes()
        samples = ({}, dict(self.system_info, cpu=""),
                   dict(self.system_info, cpu="bad\ncpu"),
                   dict(self.system_info, system="bad\rrelease"),
                   dict(self.system_info, hostname=123))
        for system_info in samples:
            with self.subTest(system_info=system_info):
                with self.assertRaises(ValueError):
                    support.append_results(self.results, self.cases,
                                           self.measurements, "new",
                                           system_info=system_info)
                self.assertEqual(self.results.read_bytes(), original)

    def test_invalid_measurements_do_not_create_a_result_row(self):
        for measurements in (
                "beta\t1\t2\nalpha\t3\t4\n",
                "alpha\t-1\t2\nbeta\t3\t4\n",
                "alpha\tnan\t2\nbeta\t3\t4\n",
                "alpha\t1\t2\textra\nbeta\t3\t4\n"):
            with self.subTest(measurements=measurements):
                self.measurements.write_text(measurements)
                with self.assertRaises(ValueError):
                    support.append_results(self.results, self.cases,
                                           self.measurements, "new")
                self.assertFalse(self.results.exists())


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="manifest with spaces ")
        self.addCleanup(self.temporary.cleanup)
        self.manifest = Path(self.temporary.name) / "cases.tsv"

    def test_tab_delimiters_preserve_filenames_with_spaces(self):
        self.manifest.write_text("formula\t3\t2\t0\t1\t<none>\tany\n"
                                 "from_file\t3\t2\t3\t0\tmatrix with spaces.txt\tsolved\n")
        cases = support.read_manifest(self.manifest)
        self.assertEqual(cases[0], dict(name="formula", n=3, m=2, r=0, s=1,
                                      filename="<none>", expectation="any"))
        self.assertEqual(cases[1]["filename"], "matrix with spaces.txt")
        self.assertEqual(support.result_header(cases),
                         METADATA_HEADER + ["formula.T1", "formula.T2",
                                            "from_file.T1", "from_file.T2"])

    def test_invalid_or_duplicate_cases_are_rejected(self):
        valid = "formula\t3\t2\t0\t1\t<none>\tany\n"
        samples = ("", valid + valid, valid.replace("formula", "bad name"),
                   valid.replace("\t3\t", "\tthree\t"),
                   valid.replace("<none>", ""), valid.replace("any", "unknown"),
                   valid.rstrip() + "\textra\n")
        for text in samples:
            with self.subTest(text=text):
                self.manifest.write_text(text)
                with self.assertRaises(ValueError):
                    support.read_manifest(self.manifest)


class SystemInfoTests(unittest.TestCase):
    def test_current_system_has_four_nonempty_metadata_fields(self):
        self.assertEqual(tuple(support.SYSTEM_COLUMNS), SYSTEM_COLUMNS)
        system_info = support.collect_system_info()
        self.assertEqual(set(system_info), set(SYSTEM_COLUMNS))
        for key, value in system_info.items():
            with self.subTest(key=key):
                self.assertIsInstance(value, str)
                self.assertTrue(value.strip())
                self.assertNotIn("\n", value)
                self.assertNotIn("\r", value)


if __name__ == "__main__":
    unittest.main(verbosity=2)
