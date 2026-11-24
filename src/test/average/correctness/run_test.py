#!/usr/bin/env python3
"""Focused average crash and percentile regressions."""
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

EXECUTABLE = str(pathlib.Path(sys.argv[1]).resolve())

class AverageTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = pathlib.Path(self.tmp.name)

    def write(self, name, content):
        (self.directory / name).write_text(content)
        return name

    def run_tool(self, *args, success=True):
        result = subprocess.run([EXECUTABLE, *args], cwd=self.directory,
                                capture_output=True, text=True, timeout=20)
        self.assertGreaterEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("AddressSanitizer", result.stderr)
        self.assertNotIn("runtime error:", result.stderr)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result.stdout

    def test_single_value(self):
        name = self.write("one", "7\n")
        output = self.run_tool("-c", "1", "-e", "0,50,100", name)
        self.assertIn("ave 7 std 0", output)
        for pct in (0, 50, 100):
            self.assertIn(f"{pct}% 7", output)

    def test_small_sample_percentiles(self):
        for content, median in [("2\n1\n", 2), ("30\n10\n20\n", 20)]:
            with self.subTest(content=content):
                name = self.write("sample", content)
                output = self.run_tool("-c", "1", "-e", "0,50,99,100", name)
                values = [int(x) for x in content.split()]
                for pct, expected in [(0, min(values)), (50, median),
                                      (99, max(values)), (100, max(values))]:
                    self.assertRegex(output, rf"(?:^| ){pct}% {expected}(?: |$)")

    def test_bounded_sample_default_and_explicit_columns(self):
        name = self.write("many", "7\n" * 20000)
        for columns in [[], ["-c", "1"]]:
            with self.subTest(columns=columns):
                output = self.run_tool(*columns, "-r", "10", "-e", "0,50,100", name)
                self.assertIn("20000 values", output)
                self.assertIn("0% 7 50% 7 100% 7", output)

    def test_sample_resets_between_files(self):
        self.write("first", "1\n" * 1000)
        self.write("second", "9\n" * 1000)
        output = self.run_tool("-c", "1", "-r", "10", "-e", "50", "first", "second")
        self.assertEqual(re.findall(r"50% (\S+)", output), ["1", "9"])

    def test_descriptor_columns_are_resolved_for_each_file(self):
        self.write("first", "ID a x\na 1 9\nb 2 7\n")
        for header, data in [("ID x a", "a 3 30\nb 5 50\n"),
                             ("ID x", "a 3\nb 5\n"),
                             ("ID a b x", "a 30 300 3\nb 50 500 5\n")]:
            with self.subTest(header=header):
                self.write("second", header + "\n" + data)
                output = self.run_tool("-d", "x", "first", "second")
                lines = output.splitlines()
                self.assertEqual(len(lines), 2, output)
                self.assertTrue(all(line.startswith("x 2 values") for line in lines), output)
                self.assertIn("ave 8", lines[0])
                self.assertIn("ave 4", lines[1])

    def test_numeric_column_names_do_not_survive_previous_file(self):
        self.write("first", "ID first_name\na 1\nb 3\n")
        self.write("second", "ID second_name\na 5\nb 7\n")
        output = self.run_tool("-c", "2", "-j", "first", "second")
        self.assertTrue(output.splitlines()[1].startswith("second_name "), output)

    def test_scaled_prevalence_overflow_is_rejected(self):
        self.write("counts", "1 1e300\n")
        self.run_tool("-c", "1", "-h", "2", "-h", "mult=10", "counts", success=False)
        self.write("counts", "1 1e308\n")
        self.run_tool("-c", "1", "-h", "2", "-h", "mult=10", "counts", success=False)

    def test_numeric_prevalence_applies_to_every_file(self):
        self.write("counts", "1 1\n9 3\n")
        output = self.run_tool("-c", "1", "-h", "2", "counts", "counts")
        self.assertEqual(output.count("4 values"), 2, output)
        self.assertEqual(output.count("ave 7"), 2, output)

    def test_prevalence_multiplier_applies_to_every_file(self):
        self.write("counts", "1 0.1\n9 0.3\n")
        output = self.run_tool("-c", "1", "-h", "2", "-h", "mult=10", "counts", "counts")
        self.assertEqual(output.count("4 values"), 2, output)
        self.assertEqual(output.count("ave 7"), 2, output)

    def test_named_prevalence_is_resolved_for_each_header(self):
        self.write("first", "ID value count\na 1 1\nb 9 3\n")
        self.write("second", "ID count value\na 3 2\nb 1 6\n")
        output = self.run_tool("-d", "value", "-h", "count", "first", "second")
        lines = output.splitlines()
        self.assertEqual(len(lines), 2, output)
        self.assertIn("4 values", lines[0])
        self.assertIn("ave 7", lines[0])
        self.assertIn("4 values", lines[1])
        self.assertIn("ave 3", lines[1])

    def test_named_prevalence_with_numeric_selection(self):
        self.write("first", "value count other\n1 1 100\n9 3 100\n")
        self.write("second", "value other count\n2 100 3\n6 100 1\n")
        for header_option in [["-j"], ["-s", "1"]]:
            with self.subTest(header_option=header_option):
                output = self.run_tool("-c", "1", *header_option, "-h", "count", "first", "second")
                self.assertEqual(output.count("4 values"), 2, output)
                self.assertIn("ave 7", output.splitlines()[0])
                self.assertIn("ave 3", output.splitlines()[1])

    def test_named_prevalence_multiplier_survives_moving_columns(self):
        self.write("first", "ID value count\na 1 0.1\nb 9 0.3\n")
        self.write("second", "ID count value\na 0.3 2\nb 0.1 6\n")
        output = self.run_tool("-d", "value", "-h", "count", "-h", "mult=10", "first", "second")
        self.assertEqual(output.count("4 values"), 2, output)
        self.assertIn("ave 7", output.splitlines()[0])
        self.assertIn("ave 3", output.splitlines()[1])

    def test_missing_numeric_prevalence_in_later_file_is_rejected(self):
        self.write("first", "1 1\n9 3\n")
        self.write("second", "2\n6\n")
        self.run_tool("-c", "1", "-h", "2", "first", "second", success=False)

    def test_missing_prevalence_name_in_second_file_is_rejected(self):
        self.write("first", "ID value count\na 1 1\nb 9 3\n")
        self.write("second", "ID value other\na 2 3\nb 6 1\n")
        self.run_tool("-d", "value", "-h", "count", "first", "second", success=False)
        self.write("duplicate", "ID value count count\na 1 1 1\nb 9 3 3\n")
        self.run_tool("-d", "value", "-h", "count", "duplicate", success=False)

    def test_prevalence_configuration_requires_a_column(self):
        self.write("counts", "1 1\n9 3\n")
        for args in [["-h", "mult=10"], ["-h", "count"], ["-h", "0"],
                     ["-h", "-1"], ["-h", "1", "-h", "2"]]:
            with self.subTest(args=args):
                self.run_tool("-c", "1", *args, "counts", success=False)

    def test_zero_prevalence_is_independent_of_row_order(self):
        baseline = "2 1\n4 1\n"
        args = ["-c", "1", "-h", "2", "-p", "99", "-e", "0,50,100"]
        self.write("baseline", baseline)
        expected = self.run_tool(*args, "baseline")
        for content in ["99 0\n" + baseline, baseline + "99 0\n",
                        "2 1\n99 0\n4 1\n"]:
            with self.subTest(content=content):
                self.write("counts", content)
                self.assertEqual(self.run_tool(*args, "counts"), expected)

    def test_zero_prevalence_has_no_output_in_any_selection_mode(self):
        self.write("counts", "99 0\n1 0\n")
        self.assertEqual(self.run_tool("-c", "1", "-h", "2", "-e", "50", "counts"), "")
        self.write("labelled", "value 99 0\nvalue 1 0\n")
        self.assertEqual(self.run_tool("-R", "^value$", "-h", "3", "labelled"), "")
        self.write("header", "value count\n99 0\n1 0\n")
        self.assertEqual(self.run_tool("-d", "value", "-h", "count", "header"), "")

    def test_zero_prevalence_after_scaling(self):
        self.write("counts", "99 0.01\n2 0.1\n4 0.1\n")
        output = self.run_tool("-c", "1", "-h", "2", "-h", "mult=10", "-p", "99",
                               "-e", "0,50,100", "counts")
        self.assertIn("2 values between 2 and 4", output)
        self.assertIn("0 instances of 99 fraction 0", output)
        self.assertIn("0% 2 50% 4 100% 4", output)

    def test_zero_only_file_does_not_affect_next_file(self):
        self.write("zero", "99 0\n1 0\n")
        self.write("positive", "2 1\n4 1\n")
        args = ["-c", "1", "-h", "2", "-e", "0,50,100"]
        expected = self.run_tool(*args, "positive")
        self.assertEqual(self.run_tool(*args, "zero", "positive"), expected)

    def test_weighted_special_value_count(self):
        self.write("counts", "1 1\n9 3\n")
        output = self.run_tool("-c", "1", "-h", "2", "-p", "9", "counts")
        self.assertIn("3 instances of 9 fraction 0.75", output)

    def test_weighted_percentiles_match_expanded_population(self):
        cases = [[(1, 1), (9, 3)], [(8, 4), (-3, 2), (0, 5)],
                 [(7, 12)], [(10, 2), (1, 1), (10, 3), (0, 0)]]
        for entries in cases:
            with self.subTest(entries=entries):
                self.write("weighted", "".join(f"{v} {count}\n" for v, count in entries))
                self.write("expanded", "".join(f"{v}\n" * count for v, count in entries))
                args = ["-c", "1", "-b", "-p", "0", "-e", "0,1,25,50,75,99,100"]
                expected = self.run_tool(*args, "expanded")
                self.assertEqual(self.run_tool(*args, "-h", "2", "weighted"), expected)

    def test_large_prevalences_do_not_expand_storage_or_truncate_counts(self):
        for counts in [(2000000000, 2000000000, 2000000000),
                       (2147483647, 2147483647, 2)]:
            with self.subTest(counts=counts):
                self.write("weighted", "".join(f"9 {n}\n" for n in counts))
                output = self.run_tool("-c", "1", "-h", "2", "-p", "9",
                                       "-e", "0,50,100", "weighted")
                total = sum(counts)
                self.assertIn(f"{total} values", output)
                self.assertIn(f"{total} instances of 9 fraction 1", output)
                self.assertIn("0% 9 50% 9 100% 9", output)

    def test_weighted_percentiles_and_special_counts_reset_across_files(self):
        self.write("first", "ID value count\na 1 1\nb 9 3\n")
        self.write("second", "ID count value\na 3 2\nb 1 9\n")
        output = self.run_tool("-d", "value", "-h", "count", "-p", "9",
                               "-e", "0,50,100", "first", "second")
        first, second = output.splitlines()
        self.assertIn("3 instances of 9 fraction 0.75", first)
        self.assertIn("0% 1 50% 9 100% 9", first)
        self.assertIn("1 instances of 9 fraction 0.25", second)
        self.assertIn("0% 2 50% 2 100% 9", second)

    def test_weighted_percentiles_apply_prevalence_multiplier(self):
        self.write("counts", "1 0.1\n9 0.3\n")
        output = self.run_tool("-c", "1", "-h", "2", "-h", "mult=10", "-p", "9",
                               "-e", "0,50,100", "counts")
        self.assertIn("3 instances of 9 fraction 0.75", output)
        self.assertIn("0% 1 50% 9 100% 9", output)

    def test_bounded_sampling_with_prevalence_is_rejected(self):
        self.write("counts", "1 1\n9 3\n")
        self.run_tool("-c", "1", "-h", "2", "-r", "10", "-e", "50", "counts", success=False)

    def test_missing_header_is_rejected(self):
        self.write("empty", "")
        self.run_tool("-d", "x", "empty", success=False)
        self.run_tool("-c", "1", "-s", "1", "empty", success=False)

if __name__ == "__main__":
    suite = unittest.main(argv=[sys.argv[0]], exit=False)
    if not suite.result.wasSuccessful():
        sys.exit(1)
    print("average correctness tests passed")
