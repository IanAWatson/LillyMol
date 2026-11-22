"""Regression checks for integral/spatial distance accumulation and histogram cleanup."""

import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class DbfDistancesTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        # C-C-C: bond distances 1, 1, 2; spatial distances 1.25, 1.75, 3.
        record = """propane
  LillyMol

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    0.0000    1.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8000    0.0000    2.4000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""
        (self._directory / "input.sdf").write_text(record)

    def run_tool(self, *options, expected_records=1, inputs=("input.sdf",)):
        result = subprocess.run(
            [_EXECUTABLE, "-s", "[#6]", "-d", *options, *inputs],
            cwd=self._directory, capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = result.stdout.splitlines()
        self.assertEqual(len(rows), expected_records + 1, result.stdout)
        header, values = (row.split() for row in rows[:2])
        self.assertEqual(len(header), len(values), result.stdout)
        self.assertEqual(values[0], "propane")
        return dict(zip(header[1:], map(float, values[1:]))), result

    # The existing output table has a 0.003-wide distance grid.
    def test_topological_distances_keep_fractional_mean(self):
        values, _ = self.run_tool("-p", "2d")
        self.assertEqual(values["dbf_minaa"], 1)
        self.assertEqual(values["dbf_maxaa"], 2)
        self.assertAlmostEqual(values["dbf_aveaa"], 4 / 3, delta=0.002)

    def test_spatial_distances_keep_fractional_values(self):
        values, _ = self.run_tool("-p", "3d")
        self.assertAlmostEqual(values["dbf_3minaa"], 1.25, delta=0.002)
        self.assertAlmostEqual(values["dbf_3maxaa"], 3, delta=0.002)
        self.assertAlmostEqual(values["dbf_3aveaa"], 2, delta=0.002)

    def test_spatial_topological_ratios(self):
        values, _ = self.run_tool("-p", "ratio")
        self.assertAlmostEqual(values["dbf_rminaa"], 1.25, delta=0.002)
        self.assertAlmostEqual(values["dbf_rmaxaa"], 1.75, delta=0.002)
        self.assertAlmostEqual(values["dbf_raveaa"], 1.5, delta=0.002)

    def test_histogram_setup_and_destruction_preserve_distances(self):
        for mode in ("2d", "3d", "ratio"):
            with self.subTest(mode=mode):
                baseline, _ = self.run_tool("-p", mode)
                values, _ = self.run_tool("-p", mode, "-T", "0,10,1")
                self.assertEqual(values, baseline)

    def histograms(self, stderr):
        """Extract labeled histogram totals and bucket populations."""
        reports = {}
        label = None
        for line in stderr.splitlines():
            match = re.fullmatch(r"(.+) i = 0 j = 0", line)
            if match:
                label = match[1]
            match = re.fullmatch(r"Histogram from .* with (\d+) samples", line)
            if match and label is not None:
                reports[label] = {"samples": int(match[1]), "buckets": {}}
            match = re.fullmatch(r"([\d.eE+-]+) (\d+) values", line)
            if match and label in reports:
                reports[label]["buckets"][float(match[1])] = int(match[2])
        return reports

    def test_histograms_include_single_final_molecule(self):
        _, result = self.run_tool("-p", "ratio", "-T", "0,10,1")
        reports = self.histograms(result.stderr)
        self.assertEqual(set(reports), {"Topological", "Spatial", "Spatial/topological ratio"})
        expected = {
            "Topological": {1: 2, 2: 1},
            "Spatial": {1: 1, 2: 1, 3: 1},
            "Spatial/topological ratio": {1: 1, 2: 2},
        }
        for label, counts in expected.items():
            with self.subTest(label=label):
                self.assertEqual(reports[label]["samples"], 3)
                nonempty = {k: v for k, v in reports[label]["buckets"].items() if v}
                self.assertEqual(nonempty, counts)

    def test_histograms_count_each_molecule_once(self):
        path = self._directory / "input.sdf"
        record = path.read_text()
        path.write_text(record * 3)
        _, result = self.run_tool("-p", "ratio", "-T", "0,10,1", expected_records=3)
        reports = self.histograms(result.stderr)
        self.assertEqual(len(reports), 3)
        for report in reports.values():
            self.assertEqual(report["samples"], 9)
            self.assertEqual(sum(report["buckets"].values()), 9)

    def test_histograms_accumulate_across_input_files(self):
        (self._directory / "second.sdf").write_text(
            (self._directory / "input.sdf").read_text())
        _, result = self.run_tool("-p", "3d", "-T", "0,10,1",
                                  expected_records=2, inputs=("input.sdf", "second.sdf"))
        report = self.histograms(result.stderr)["Spatial"]
        self.assertEqual(report["samples"], 6)
        self.assertEqual(sum(report["buckets"].values()), 6)

    def test_histograms_report_only_enabled_distance_types(self):
        for mode, label in (("2d", "Topological"), ("3d", "Spatial")):
            with self.subTest(mode=mode):
                _, result = self.run_tool("-p", mode, "-T", "0,10,1")
                reports = self.histograms(result.stderr)
                self.assertEqual(set(reports), {label})
                self.assertEqual(reports[label]["samples"], 3)

    def test_histogram_excludes_out_of_range_distances(self):
        _, result = self.run_tool("-p", "3d", "-T", "0,2,1")
        report = self.histograms(result.stderr)["Spatial"]
        self.assertEqual(report["samples"], 2)
        self.assertEqual(sum(report["buckets"].values()), 2)


if __name__ == "__main__":
    unittest.main()
