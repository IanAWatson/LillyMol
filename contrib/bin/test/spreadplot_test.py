"""Run with python3 contrib/bin/test/spreadplot_test.py."""
import argparse
import csv
import importlib.util
import io
import pathlib
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "spreadplot.py"
SPEC = importlib.util.spec_from_file_location("spreadplot", SCRIPT)
spreadplot = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = spreadplot
SPEC.loader.exec_module(spreadplot)


class SpreadplotTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.directory = pathlib.Path(self._tmp.name)
        self.input = self.write(
            "trajectory.csv", "sel,distance\n0,1\n10,0.8\n20,0.5\n30,0.2\n40,0.1\n")

    def write(self, name, contents):
        path = self.directory / name
        path.write_text(contents)
        return path

    def run_tool(self, *args, input_text=None):
        return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                              input=input_text, capture_output=True, text=True, timeout=30)

    def trajectory(self):
        return spreadplot.read_trajectory(str(self.input))

    def test_parse_thresholds(self):
        self.assertEqual(spreadplot.parse_thresholds("0.1, 0.25,1"), [0.1, 0.25, 1])
        for value in ("", "0.1,", ",0.1", "abc", "-0.1", "1.1", "NaN", "Inf"):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    spreadplot.parse_thresholds(value)

    def test_crossings_and_regular_points(self):
        values = [1.0, 0.8, 0.5, 0.2, 0.0]
        self.assertEqual(spreadplot.crossing_point(values, 1.0), 0)
        self.assertEqual(spreadplot.crossing_point(values, 0.5), 2)
        self.assertIsNone(spreadplot.crossing_point(values, -0.1))
        self.assertEqual(spreadplot.crossing_points(values, 0.1), [1, 2, 3, 4])
        self.assertEqual(spreadplot.crossing_points(values, 0.25), [2, 3])
        self.assertEqual(spreadplot.crossing_points(values, 1.0), [])
        self.assertEqual(spreadplot.crossing_points([1, 0.1, 0], 0.1), [1])
        self.assertEqual(spreadplot.crossing_points(values, 1e-8), [1, 2, 3, 4])
        self.assertEqual(spreadplot.crossing_points(values, 1e-300), [1, 2, 3, 4])
        for delta in (0, -0.1, 1.1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                spreadplot.crossing_points(values, delta)

    def test_read_trajectory_preserves_positions(self):
        trajectory = self.trajectory()
        self.assertEqual(trajectory.selected, [0, 10, 20, 30, 40])
        self.assertEqual(trajectory.distance, [1, 0.8, 0.5, 0.2, 0.1])
        self.assertAlmostEqual(trajectory.fraction(1), 11 / 41)
        self.assertEqual(trajectory.fraction(4), 1)

    def test_bad_trajectories(self):
        cases = [
            "", "count,distance\n0,1\n1,0.5\n", "sel,value\n0,1\n1,0.5\n",
            "sel,distance\n0,1\n", "sel,distance\n0,1\n1,\n",
            "sel,distance\n0,1\n1.5,0.5\n", "sel,distance\n-1,1\n0,0.5\n",
            "sel,distance\n0,1\n0,0.5\n", "sel,distance\n1,1\n0,0.5\n",
            "sel,distance\n0,1\n1,1.1\n", "sel,distance\n0,1\n1,-0.1\n",
            "sel,distance\n0,0.5\n1,0.6\n", "sel,distance\n0,1\n1,nan\n",
            "sel,distance\n0,1\n1,inf\n", "sel,sel,distance\n0,0,1\n1,1,0.5\n",
            "sel,distance\n0,1,extra\n1,0.5\n",
            'sel,distance\n0,1\n1,"unterminated\n']
        for contents in cases:
            with self.subTest(contents=contents):
                path = self.write("bad.csv", contents)
                with self.assertRaises(ValueError):
                    spreadplot.read_trajectory(str(path))

    def test_extra_columns_and_stdin(self):
        path = self.write("extra.csv", "sel,distance,note\n0,1,first\n10,0.5,last\n")
        self.assertEqual(spreadplot.read_trajectory(str(path)).selected, [0, 10])
        with patch.object(sys, "stdin", io.StringIO(path.read_text())):
            self.assertEqual(spreadplot.read_trajectory("-").selected, [0, 10])

    def test_summary_status_and_positions(self):
        rows = spreadplot.threshold_rows([self.trajectory()], [1, 0.5, 0.05])
        self.assertEqual([r["status"] for r in rows], ["at_start", "crossed", "not_reached"])
        self.assertEqual([r["sel"] for r in rows], [0, 20, ""])
        self.assertAlmostEqual(rows[1]["fraction_selected"], 21 / 41)
        self.assertEqual(rows[2]["fraction_selected"], "")

    def test_reduced_trajectory_crossings_match_full(self):
        full = self.write("full.csv",
                          "sel,distance\n0,1\n1,0.8\n2,0.8\n3,0.8\n"
                          "4,0.5\n5,0.5\n6,0.1\n")
        reduced = self.write("reduced.csv",
                             "sel,distance\n0,1\n1,0.8\n3,0.8\n4,0.5\n5,0.5\n6,0.1\n")
        a, b = [spreadplot.read_trajectory(str(path)) for path in (full, reduced)]
        for threshold in (i / 100 for i in range(101)):
            ra, rb = [spreadplot.threshold_rows([t], [threshold])[0] for t in (a, b)]
            self.assertEqual(ra["sel"], rb["sel"])
            self.assertEqual(ra["fraction_selected"], rb["fraction_selected"])

    def test_second_distance_controls_y_limit(self):
        figure = spreadplot.make_plot([self.trajectory()], [])
        self.assertEqual(figure.axes[0].get_ylim(), (0, 0.8))
        self.assertEqual(list(figure.axes[0].lines[0].get_ydata())[0], 1)
        figure = spreadplot.make_plot([self.trajectory()], [0.9])
        self.assertEqual(figure.axes[0].get_ylim(), (0, 0.9))

    def test_default_connects_distinct_levels_without_changing_crossings(self):
        path = self.write("plateaus.csv",
                          "sel,distance\n0,1\n1,0.8\n4,0.8\n5,0.7\n"
                          "9,0.7\n10,0.5\n14,0.5\n")
        trajectory = spreadplot.read_trajectory(str(path))
        before = spreadplot.threshold_rows([trajectory], [0.75, 0.6])
        figure = spreadplot.make_plot([trajectory], [])
        line = figure.axes[0].lines[0]
        self.assertEqual(list(line.get_xdata()), [0, 1, 5, 10, 14])
        self.assertEqual(list(line.get_ydata()), [1, 0.8, 0.7, 0.5, 0.5])
        self.assertEqual(figure.axes[0].get_ylim(), (0, 0.8))
        self.assertEqual(spreadplot.threshold_rows([trajectory], [0.75, 0.6]), before)
        self.assertEqual([row["sel"] for row in before], [5, 10])
        exact = spreadplot.make_plot([trajectory], [], plateaus=True)
        self.assertEqual(list(exact.axes[0].lines[0].get_xdata()), trajectory.selected)
        self.assertEqual(list(exact.axes[0].lines[0].get_ydata()), trajectory.distance)

    def test_plateaus_option_renders_exact_curve(self):
        result = self.run_tool("-S", self.directory / "plateaus", "--plateaus",
                               "--format", "svg", self.input)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.directory / "plateaus.svg").exists())

    def test_all_zero_distances(self):
        path = self.write("zero.csv", "sel,distance\n0,0\n100,0\n")
        figure = spreadplot.make_plot([spreadplot.read_trajectory(str(path))], [])
        self.assertGreater(figure.axes[0].get_ylim()[1], 0)

    def test_normalization_and_overlays(self):
        a = self.trajectory()
        path = self.write("other.csv", "sel,distance\n0,1\n99,0.6\n199,0.2\n")
        b = spreadplot.read_trajectory(str(path), "other")
        figure = spreadplot.make_plot([a, b], [], normalize=True)
        self.assertEqual(figure.axes[0].get_ylim(), (0, 0.8))
        for line, total in zip(figure.axes[0].lines, [41, 200]):
            self.assertAlmostEqual(line.get_xdata()[0], 1 / total)
            self.assertEqual(line.get_xdata()[-1], 1)

    def test_log_axis_omits_zero_but_keeps_summary(self):
        trajectory = self.trajectory()
        figure = spreadplot.make_plot([trajectory], [], logx=True)
        self.assertEqual(figure.axes[0].get_xscale(), "log")
        self.assertGreater(min(figure.axes[0].lines[0].get_xdata()), 0)
        self.assertEqual(spreadplot.threshold_rows([trajectory], [1])[0]["sel"], 0)

    def test_monotonic_spline_and_original_markers(self):
        trajectory = self.trajectory()
        figure = spreadplot.make_plot([trajectory], [], spline=50, point=0.25)
        values = figure.axes[0].lines[0].get_ydata()
        self.assertEqual(len(values), 50)
        self.assertTrue(all(a >= b for a, b in zip(values, values[1:])))
        offsets = figure.axes[0].collections[0].get_offsets()
        self.assertEqual(list(offsets[:, 0]), [20, 30])

    def test_spline_warns_below_input_size_but_still_writes_plot(self):
        for samples, warning in ((4, True), (5, False), (6, False)):
            with self.subTest(samples=samples):
                stem = self.directory / f"spline_{samples}"
                result = self.run_tool("-S", stem, "--spline", samples, self.input)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(stem.with_suffix(".png").exists())
                self.assertEqual("Warning:" in result.stderr, warning)
                if warning:
                    self.assertIn(str(self.input), result.stderr)
                    self.assertIn("5 input datapoints", result.stderr)
                    self.assertIn("at least --spline 5", result.stderr)

    def test_png_svg_pdf_rendering(self):
        signatures = {"png": b"\x89PNG", "svg": b"<?xml", "pdf": b"%PDF"}
        for extension, signature in signatures.items():
            with self.subTest(extension=extension):
                stem = self.directory / extension
                result = self.run_tool("-S", stem, "--format", extension, "--dist",
                                       "0.9,0.5,0.05", "--point", "0.1", self.input)
                self.assertEqual(result.returncode, 0, result.stderr)
                output = stem.with_suffix("." + extension)
                self.assertGreater(output.stat().st_size, 100)
                self.assertTrue(output.read_bytes().startswith(signature))

    def test_cli_multiple_curves_spline_and_summary(self):
        other = self.write("other.csv", "sel,distance\n0,1\n99,0.6\n199,0.2\n")
        summary = self.directory / "thresholds.csv"
        stem = self.directory / "combined"
        result = self.run_tool(
            "-S", stem, "--format", "svg", "--normalize", "--logx",
            "--label", "Set A", "--label", "Set B", "--color", "green",
            "--color", "blue", "--dist", "1,0.5,0.05", "--summary", summary,
            "--spline", "60", "--width", "8", "--height", "5", "--dpi", "100",
            "--title", "Collections", self.input, other)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(stem.with_suffix(".svg").exists())
        with summary.open() as input_file:
            rows = list(csv.DictReader(input_file))
        self.assertEqual(len(rows), 6)
        self.assertEqual([row["label"] for row in rows], ["Set A"] * 3 + ["Set B"] * 3)
        self.assertEqual(rows[1]["sel"], "20")
        self.assertEqual(rows[4]["sel"], "199")
        self.assertEqual(rows[-1]["status"], "not_reached")
        self.assertEqual(rows[-1]["sel"], "")

    def test_cli_stdin_and_summary_stdout(self):
        result = self.run_tool("-S", self.directory / "stdin", "--dist", "0.5",
                               "--summary", "-", "-", input_text=self.input.read_text())
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = list(csv.DictReader(io.StringIO(result.stdout)))
        self.assertEqual(rows[0]["sel"], "20")
        self.assertEqual(rows[0]["file"], "-")

    def test_invalid_arguments_do_not_create_output(self):
        stem = self.directory / "bad"
        for extra in (["--summary", "x.csv"], ["--point", "0"], ["--point", "2"],
                      ["--spline", "1"], ["--spline", "0"], ["--width", "-1"],
                      ["--height", "nan"], ["--dpi", "0"], ["--lwd", "inf"],
                      ["--dist", "0.1,"], ["--label", "a", "--label", "b"],
                      ["--color", "red", "--color", "blue"], ["--color", "not_a_colour"]):
            with self.subTest(extra=extra):
                result = self.run_tool("-S", stem, *extra, self.input)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("Traceback", result.stderr)
                self.assertFalse(stem.with_suffix(".png").exists())

    def test_bad_csv_reports_error(self):
        path = self.write("bad.csv", "sel,distance\n0,1\n5,NaN\n")
        result = self.run_tool("-S", self.directory / "bad", path)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CSV line 3", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_does_not_overwrite_input_or_alias_outputs(self):
        original = self.input.read_text()
        result = self.run_tool("-S", self.directory / "result", "--dist", "0.5",
                               "--summary", self.input, self.input)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.input.read_text(), original)
        disguised = self.write("trajectory.png", original)
        result = self.run_tool("-S", self.directory / "trajectory", disguised)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(disguised.read_text(), original)
        result = self.run_tool("-S", self.directory / "result", "--dist", "0.5",
                               "--summary", self.directory / "result.png", self.input)
        self.assertNotEqual(result.returncode, 0)
        alias = self.directory / "alias.csv"
        alias.symlink_to(self.input)
        result = self.run_tool("-S", self.directory / "result", "--dist", "0.5",
                               "--summary", alias, self.input)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.input.read_text(), original)

    def test_help_does_not_require_plotting_packages(self):
        result = subprocess.run([sys.executable, "-S", str(SCRIPT), "--help"],
                                capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("nplotnn", result.stdout)


if __name__ == "__main__":
    unittest.main()
