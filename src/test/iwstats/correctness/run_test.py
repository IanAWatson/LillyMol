#!/usr/bin/env python3
"""Numerical references and input-validity regressions for iwstats."""
import math
import pathlib
import re
import statistics
import subprocess
import sys
import tempfile
import unittest

EXE = str(pathlib.Path(sys.argv[1]).resolve())

class StatsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = pathlib.Path(self.tmp.name)

    def write(self, name, text):
        (self.directory / name).write_text(text)
        return name

    def run_stats(self, text, *options, success=True):
        self.write("data", text)
        args = [EXE, "-e", "1", "-p", "2", "-c", "3", "-Y", "allequals", *options, "data"]
        result = subprocess.run(args, cwd=self.directory, text=True, capture_output=True, timeout=20)
        self.assertGreaterEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("AddressSanitizer", result.stderr)
        self.assertNotIn("runtime error:", result.stderr)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result.stdout

    def metric(self, output, key):
        match = re.search(r"(?:^| )" + re.escape(key) + r"[ =]+([^ \n]+)", output)
        self.assertIsNotNone(match, (key, output))
        value = float(match[1])
        self.assertTrue(math.isfinite(value), (key, output))
        return value

    def check_reference(self, observed, predicted, *options):
        text = "".join(f"{x:.17g} {y:.17g}\n" for x, y in zip(observed, predicted))
        output = self.run_stats(text, "-w", "-h", *options)
        ox, py = statistics.mean(observed), statistics.mean(predicted)
        sx = sum((x-ox)**2 for x in observed)
        sy = sum((y-py)**2 for y in predicted)
        cov = sum((x-ox)*(y-py) for x,y in zip(observed,predicted))
        rho = cov / math.sqrt(sx*sy) if sx and sy else 0
        squared_errors = sum((x-y)**2 for x,y in zip(observed,predicted))
        expected = {"R2": math.copysign(rho*rho, rho),
                    "trQ2": 1-squared_errors/sx if sx else (1 if not squared_errors else 0),
                    "Bias": ox-py,
                    "AAE": statistics.mean(abs(x-y) for x,y in zip(observed,predicted)),
                    "RMS": math.sqrt(squared_errors/len(observed))}
        for key, value in expected.items():
            self.assertAlmostEqual(self.metric(output,key), value,
                                   delta=max(1e-5,abs(value)*1e-5), msg=(key,output))
        return output

    def test_double_precision_and_centered_correlation(self):
        values = [1e9 + i/10 for i in range(25)]
        self.check_reference(values, values)
        self.check_reference(values, [v + .025 for v in values])

    def test_general_reference_and_negative_correlation(self):
        self.check_reference([1,2,4,7,11], [2,1,5,6,10])
        self.check_reference([1,2,4,7,11], [10,8,5,3,1])

    def test_missing_measured_value_is_not_zero(self):
        text = "".join(f"{i} {i+.2}\n" for i in range(1,25))
        expected = self.run_stats(text, "-w", "-d")
        self.assertEqual(self.run_stats(". 1\n" + text, "-w", "-d", "-M", "."), expected)

    def test_missing_predictions_do_not_enter_metrics(self):
        text = "".join(f"{i} {i+.2}\n" for i in range(1,25))
        expected = self.run_stats(text, "-w", "-d", "-G", "1", "-D", "1")
        self.assertEqual(self.run_stats("99 .\n" + text, "-w", "-d", "-M", ".",
                                        "-G", "1", "-D", "1"), expected)

    def test_filtering_matches_physical_subset(self):
        text = "".join(f"{i} {i+.2}\n" for i in range(25))
        subset = "".join(f"{i} {i+.2}\n" for i in range(20,25))
        expected = self.run_stats(subset, "-w", "-d", "-G", "1", "-D", "1")
        self.assertEqual(self.run_stats(text, "-w", "-d", "-u", "20", "-G", "1", "-D", "1"), expected)
        self.run_stats(text, "-u", "100", success=False)

    def test_truncation_applies_to_all_error_metrics(self):
        observed = list(range(1,11))
        predicted = [-10] + list(range(2,10)) + [20]
        text = "".join(f"{x} {y}\n" for x,y in zip(observed,predicted))
        output = self.run_stats(text,"-T", "-w", "-h")
        for key, expected in [("R2",1),("trQ2",1),("Bias",0),("AAE",0),("RMS",0)]:
            self.assertEqual(self.metric(output,key), expected)

    def test_odd_residual_median_against_reference(self):
        residuals = [5,15,5,-20,6,-15,1]
        text = "".join(f"{100*i} {100*i-d}\n" for i,d in enumerate(residuals))
        output = self.run_stats(text,"-d","-w")
        bias = abs(statistics.mean(residuals))
        mad = statistics.median(abs(d-statistics.median(residuals)) for d in residuals)*1.4826
        k = 1.6449+.3325*math.exp(-4.1564*bias/mad)
        self.assertAlmostEqual(self.metric(output,"cMSD95"),bias+k*mad,places=4)

    def test_constant_vectors_keep_error_metrics(self):
        self.check_reference(list(range(25)), [100]*25)
        self.check_reference([7]*25, list(range(25)))
        self.check_reference([7]*25, [7]*25)

    def test_absolute_error_percentile_100(self):
        output = self.run_stats("1 2\n2 4\n3 6\n", "-B", "100", "-w")
        self.assertEqual(self.metric(output,"AE100"),3)

    def test_within_option_without_lowercase_w(self):
        output = self.run_stats("1 1.1\n2 2.2\n3 3.3\n4 5\n", "-W", "0.25")
        self.assertEqual(self.metric(output,"WITHIN:0.25"),.5)

    def test_ranges_exclude_missing_predictions(self):
        text = "".join(f"{i} {i+.2}\n" for i in range(1,25))
        expected = self.run_stats(text,"-w","-d","-r","1","-T")
        output = self.run_stats("999 .\n"+text,"-M",".","-w","-d","-r","1","-T")
        self.assertEqual(output,expected)

    def test_correction_formula_uses_retained_pairs(self):
        text = "999 .\n" + "".join(f"{i} {2*i+3}\n" for i in range(10))
        self.run_stats(text,"-M",".","-u","5","-C","correction","-w")
        proto = (self.directory/"correction.textproto").read_text()
        self.assertIn("slope: 2", proto)
        self.assertIn("intercept: 3", proto)

    def test_top_fraction_tiny_population(self):
        self.write("activities", "ID measured\na 1\nb 2\nc 3\n")
        self.write("predictions", "ID prediction\na 1\nb 2\nc 3\n")
        result = subprocess.run([EXE,"-E","activities","-e","2","-p","2","-j","-c","3",
                                 "-U","0.001","predictions"],cwd=self.directory,
                                 text=True,capture_output=True)
        self.assertGreaterEqual(result.returncode,0,result.stderr)
        self.assertNotIn("AddressSanitizer", result.stderr)
        self.assertNotIn("runtime error:",result.stderr)

    def test_optional_outputs_cannot_overwrite_multiple_analyses(self):
        text = "1 2 3\n2 3 4\n3 4 5\n"
        for option in ["-H","-C"]:
            self.run_stats(text,option,"out","-p","3",success=False)
            self.assertFalse((self.directory/"out.txt").exists())
            self.assertFalse((self.directory/"out.textproto").exists())

    def test_input_read_errors(self):
        (self.directory/"input_directory").mkdir()
        result = subprocess.run([EXE,"-e","1","-p","2","input_directory"],
                                cwd=self.directory,text=True,capture_output=True)
        self.assertGreater(result.returncode,0,result.stderr)

    @unittest.skipUnless(pathlib.Path("/dev/full").exists(), "requires /dev/full")
    def test_output_write_errors(self):
        text = "".join(f"{i} {i+((i%3)+1)*.1}\n" for i in range(150))
        self.run_stats(text,"-L","/dev/full",success=False)
        (self.directory/"coverage.txt").symlink_to("/dev/full")
        self.run_stats(text,"-H","coverage",success=False)
        self.write("data",text)
        with open("/dev/full","w") as full:
            result = subprocess.run([EXE,"-e","1","-p","2","data"],cwd=self.directory,
                                    stdout=full,stderr=subprocess.PIPE,text=True)
        self.assertGreater(result.returncode,0,result.stderr)

    def test_optional_distributions_use_valid_pairs(self):
        text = "".join(f"{i} {i+.2}\n" for i in range(1,25))
        expected = self.run_stats(text,"-w","-b","5")
        self.assertEqual(self.run_stats("999 .\n"+text,"-w","-M",".","-b","5"),expected)

    def test_coverage_percentages_are_complete_and_monotone(self):
        text = "".join(f"{i} {i+(i%7)*.1}\n" for i in range(150))
        self.run_stats(text,"-w","-H","coverage")
        values = [line.split() for line in (self.directory/"coverage.txt").read_text().splitlines()[1:]]
        self.assertEqual([int(v[0]) for v in values],list(range(101)))
        widths = [float(v[1]) for v in values]
        self.assertTrue(all(math.isfinite(x) for x in widths))
        self.assertEqual(widths,sorted(widths))
        self.assertAlmostEqual(widths[-1],1.2)

    def test_full_percentage_does_not_exceed_large_input(self):
        text = "".join(f"{i} {i}\n" for i in range(100000))
        output = self.run_stats(text,"-w","-P","100")
        self.assertIn("first 100% observations (100000 records)", output)

    def test_top_fraction_requires_external_activities(self):
        self.run_stats("1 2\n2 4\n3 6\n", "-U", "0.1", success=False)

    def test_missing_experimental_value_from_external_file(self):
        self.write("experimental", "ID measured\na .\nb 1\nc 2\nd 3\n")
        # Missing activity values follow the selected missing-value convention.
        result = subprocess.run([EXE,"-E","experimental","-e","2","-p","2","-c","3",
                                 "-M",".","-j",self.write("predictions","ID prediction\na 99\nb 1\nc 2\nd 3\n")],
                                 cwd=self.directory,text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn("N = 3", result.stdout)

if __name__ == "__main__":
    suite = unittest.main(argv=[sys.argv[0]], exit=False)
    if not suite.result.wasSuccessful(): sys.exit(1)
    print("iwstats correctness tests passed")
