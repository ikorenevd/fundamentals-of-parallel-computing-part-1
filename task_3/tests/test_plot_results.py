#!/usr/bin/env python3
"""Regression checks for reading and displaying benchmark result history."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "plot_results.py"
sys.path.insert(0, str(SCRIPT.parent))
import plot_results


class ResultFileTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="plot results ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.csv_file = self.directory / "results.csv"

    def write(self, content):
        self.csv_file.write_text(content, encoding="utf-8")
        return self.csv_file

    def sample(self):
        return self.write(
            "timestamp,regular.T1,regular.T2,singular.T1,singular.T2\n"
            "2026-10-06T02:00:00.123456+03:00,0,0,-1,-1\n"
            "2026-10-06T02:00:00.654321+03:00,1.2e-5,2.3e-6,-1,-1\n"
        )

    def system_sample(self):
        return self.write(
            "timestamp,system,architecture,cpu,hostname,regular.T1,regular.T2,singular.T1,singular.T2\n"
            '2026-10-06T02:00:00.123456+03:00,"Test OS, version 4",x86_64,"Example CPU, 4 cores",host.example,0,0,-1,-1\n'
        )

    def test_zero_is_success_and_failure_remains_minus_one(self):
        results = plot_results.load_results(self.sample())
        self.assertEqual(results.cases, ("regular", "singular"))
        self.assertEqual(len(results.runs), 2)
        self.assertEqual(results.runs[0].timestamp, "2026-10-06T02:00:00.123456+03:00")
        self.assertEqual(results.runs[0].values["regular"], {"T1": 0, "T2": 0})
        self.assertEqual(results.runs[0].values["singular"], {"T1": -1, "T2": -1})
        self.assertEqual(results.runs[1].values["regular"]["T1"], 1.2e-5)
        self.assertEqual(results.runs[0].system_info, {})
        self.assertEqual(plot_results.ResultRun("t", {}).system_info, {})

    def test_loads_system_metadata_with_csv_quoted_commas(self):
        results = plot_results.load_results(self.system_sample())
        self.assertEqual(results.cases, ("regular", "singular"))
        self.assertEqual(results.runs[0].system_info, {
            "system": "Test OS, version 4", "architecture": "x86_64",
            "cpu": "Example CPU, 4 cores", "hostname": "host.example",
        })
        self.assertEqual(results.runs[0].values["regular"], {"T1": 0, "T2": 0})

    def test_legacy_or_unknown_system_is_displayed_as_unknown(self):
        self.assertEqual(plot_results.describe_system(plot_results.ResultRun("t", {})), "неизвестно")
        migrated = plot_results.ResultRun("t", {}, {column: "unknown" for column in plot_results.SYSTEM_COLUMNS})
        self.assertEqual(plot_results.describe_system(migrated), "неизвестно")

    def test_rejects_partial_misordered_unknown_or_empty_system_metadata(self):
        contents = (
            "timestamp,system,a.T1,a.T2\nt,OS,0,0\n",
            "timestamp,architecture,system,cpu,hostname,a.T1,a.T2\nt,x86_64,OS,CPU,host,0,0\n",
            "timestamp,os,a.T1,a.T2\nt,OS,0,0\n",
            "timestamp,system,architecture,cpu,hostname,a.T1,a.T2\nt,OS,x86_64,,host,0,0\n",
        )
        for content in contents:
            with self.subTest(content=content):
                with self.assertRaises(plot_results.ResultsError):
                    plot_results.load_results(self.write(content))

    def test_rejects_partial_failures_and_incomplete_case_pairs(self):
        contents = (
            "timestamp,a.T1,a.T2\nt,-1,0\n",
            "timestamp,a.T1,a.T2\nt,0,-1\n",
            "timestamp,a.T1,b.T2\nt,0,0\n",
            "timestamp,a.T1,a.T1\nt,0,0\n",
        )
        for content in contents:
            with self.subTest(content=content):
                with self.assertRaises(plot_results.ResultsError):
                    plot_results.load_results(self.write(content))

    def test_rejects_damaged_rows_and_invalid_times(self):
        rows = ("t,0", "t,0,0,0", ",0,0", "t,nan,0", "t,-0.1,0", "t,abc,0")
        for row in rows:
            with self.subTest(row=row):
                with self.assertRaises(plot_results.ResultsError):
                    plot_results.load_results(self.write("timestamp,a.T1,a.T2\n" + row + "\n"))
        with self.assertRaises(plot_results.ResultsError):
            plot_results.load_results(self.write("timestamp,a.T1,a.T2\n"))

    def test_filters_preserve_csv_order_without_duplicates(self):
        cases = ("random_1", "singular_2", "random_3", "formula_4")
        self.assertEqual(plot_results.select_cases(cases, []), cases)
        self.assertEqual(
            plot_results.select_cases(cases, ["random*", "singular", "random_1"]),
            ("random_1", "singular_2", "random_3"),
        )
        with self.assertRaises(plot_results.ResultsError):
            plot_results.select_cases(cases, ["missing"])

    def test_table_cli_needs_only_standard_library_and_supports_last_filter(self):
        result = subprocess.run(
            [sys.executable, "-S", str(SCRIPT), str(self.sample()), "--table",
             "--metric", "T2", "--case", "regular", "--last", "1"],
            text=True, capture_output=True, timeout=10,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("2026-10-06T02:00:00.654321+03:00", result.stdout)
        self.assertNotIn("2026-10-06T02:00:00.123456+03:00", result.stdout)
        self.assertIn("2.3e-06", result.stdout)
        self.assertNotIn("singular", result.stdout)
        self.assertIn("Ошибок", result.stdout)

    def test_cli_reports_invalid_input_without_traceback(self):
        result = subprocess.run(
            [sys.executable, "-S", str(SCRIPT), str(self.write("bad header\n")), "--table"],
            text=True, capture_output=True, timeout=10,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("Ошибка:", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_table_cli_displays_system_metadata_without_dependencies(self):
        result = subprocess.run(
            [sys.executable, "-S", str(SCRIPT), str(self.system_sample()), "--table"],
            text=True, capture_output=True, timeout=10,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for value in ("Test OS, version 4", "x86_64", "Example CPU, 4 cores", "host.example"):
            self.assertIn(value, result.stdout)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "optional matplotlib is unavailable")
    def test_saves_headless_graph_with_successful_zero_and_failure(self):
        output = self.directory / "timings.png"
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(self.directory / "matplotlib-cache")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), str(self.system_sample()), "--save", str(output)],
            env=environment, text=True, capture_output=True, timeout=45,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertIn("Test OS, version 4", result.stdout)
        self.assertIn("host.example", result.stdout)


if __name__ == "__main__":
    unittest.main()
