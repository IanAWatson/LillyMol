"""Reduction of nplotnn -S trajectories without changing the plotted line."""
import csv
import io
import pathlib
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class SpreadReduction(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._input = self._directory / "input.nn"

    def run_trajectory(self, values, reduce=False, precision=None, verbose=False):
        with self._input.open("w") as output:
            for i, distance in enumerate(values):
                output.write(f"$SMI<C>\nPCN<m{i}>\n$SMI<N>\nPCN<n{i}>\n"
                             f"DIST<{distance}>\n|\n")
        profile = self._directory / ("reduced.csv" if reduce else "full.csv")
        command = [_EXECUTABLE, "-S", str(profile)]
        if reduce:
            command.extend(["-X", "spreadreduce"])
        if precision is not None:
            command.extend(["-c", str(precision)])
        if verbose:
            command.append("-v")
        result = subprocess.run([*command, str(self._input)],
                                capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        text = profile.read_text()
        self.assertTrue(text.startswith("sel,distance\n"))
        rows = [(int(row["sel"]), row["distance"])
                for row in csv.DictReader(io.StringIO(text))]
        return result, rows

    @staticmethod
    def expected_endpoints(rows):
        result = []
        first = 0
        while first < len(rows):
            end = first + 1
            while end < len(rows) and rows[end][1] == rows[first][1]:
                end += 1
            result.append(rows[first])
            if end - first > 1:
                result.append(rows[end - 1])
            first = end
        return result

    def test_plateaus_and_original_selection_positions(self):
        values = [1.0] + [0.81] * 5 + [0.70] * 4 + [0.20] * 6 + [0.0] * 4
        full_result, full = self.run_trajectory(values)
        reduced_result, reduced = self.run_trajectory(values, reduce=True, verbose=True)
        self.assertEqual([i for i, _ in full], list(range(len(values))))
        self.assertEqual(reduced, self.expected_endpoints(full))
        self.assertEqual([i for i, _ in reduced], [0, 1, 5, 6, 9, 10, 15, 16, 19])
        self.assertEqual(full_result.stdout, reduced_result.stdout)
        self.assertIn("reduced from 20 to 9 rows", reduced_result.stderr)
        for threshold in (i / 100 for i in range(101)):
            first_full = next((i for i, value in full if float(value) <= threshold), None)
            first_reduced = next((i for i, value in reduced
                                  if float(value) <= threshold), None)
            self.assertEqual(first_full, first_reduced, threshold)

    def test_small_and_non_monotonic_sequences(self):
        for values in ([], [1.0], [0.5, 0.5], [0.5] * 3, [0.5] * 8,
                       [1.0, 0.8, 0.6, 0.4, 0.0],
                       [1.0, 0.8, 0.8, 0.9, 0.9, 0.5]):
            with self.subTest(values=values):
                _, full = self.run_trajectory(values)
                _, reduced = self.run_trajectory(values, reduce=True)
                self.assertEqual(reduced, self.expected_endpoints(full))
                if full:
                    self.assertEqual(reduced[0], full[0])
                    self.assertEqual(reduced[-1], full[-1])

    def test_precision_controls_plateaus(self):
        values = [1.0, 0.804, 0.803, 0.802, 0.801, 0.5]
        _, full = self.run_trajectory(values)
        _, reduced = self.run_trajectory(values, reduce=True)
        self.assertEqual(reduced, self.expected_endpoints(full))
        self.assertEqual([i for i, _ in reduced], [0, 1, 4, 5])
        _, fine = self.run_trajectory(values, precision=3)
        _, fine_reduced = self.run_trajectory(values, reduce=True, precision=3)
        self.assertEqual(fine_reduced, fine)
        self.assertEqual(len(fine_reduced), len(values))

    def test_two_decimal_trajectory_has_at_most_202_points(self):
        values = [level / 100 for level in range(100, -1, -1) for _ in range(40)]
        _, full = self.run_trajectory(values)
        _, reduced = self.run_trajectory(values, reduce=True)
        self.assertEqual(len(full), 4040)
        self.assertEqual(reduced, self.expected_endpoints(full))
        self.assertLessEqual(len(reduced), 202)

    def test_preserves_polyline_at_all_original_points(self):
        values = [1.0, 0.9, 0.9, 0.9, 0.8, 0.8, 0.3, 0.3, 0.3, 0.0]
        _, full = self.run_trajectory(values)
        _, reduced = self.run_trajectory(values, reduce=True)
        segment = 0
        for x, y in full:
            while segment + 1 < len(reduced) - 1 and x > reduced[segment + 1][0]:
                segment += 1
            x0, y0 = reduced[segment]
            x1, y1 = reduced[segment + 1]
            interpolated = float(y0) + (float(y1) - float(y0)) * (x - x0) / (x1 - x0)
            self.assertAlmostEqual(float(y), interpolated, places=12)

    def test_reduction_requires_summary_output(self):
        self._input.write_text("")
        result = subprocess.run(
            [_EXECUTABLE, "-X", "spreadreduce", str(self._input)],
            capture_output=True, text=True, timeout=10)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires", result.stderr)
        self.assertIn("-S", result.stderr)


if __name__ == "__main__":
    unittest.main()
