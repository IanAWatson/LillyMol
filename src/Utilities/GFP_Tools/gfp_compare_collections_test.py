"""Input, sampling, and output regressions for gfp_compare_collections."""

import csv
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest


_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())
_FIXTURE = pathlib.Path(sys.argv.pop(1)).resolve()
_A, _B = [record + "|\n" for record in _FIXTURE.read_text().split("|\n")[:2]]


def remove_tag(record, tag):
    return "".join(line + "\n" for line in record.splitlines()
                   if not line.startswith(tag + "<"))


def zero_fingerprint(tag, nbits):
    nbytes = nbits // 8
    encoded = "." * (((nbytes + 2) // 3) * 4) + str(nbytes % 3 or 3)
    return f"{tag}<{encoded};{nbits};0;{nbits};0;1>\n"


class CompareCollections(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)

    def run_tool(self, reference=_A, queries=_A, *options, stem="out", stdin=False, threads=1):
        (self._directory / "reference.gfp").write_text(reference)
        (self._directory / "query.gfp").write_text(queries)
        return subprocess.run(
            [_EXECUTABLE, "-h", str(threads), "-v", "-p", "reference.gfp", "-S", stem,
             *options, "-" if stdin else "query.gfp"],
            cwd=self._directory, input=queries if stdin else None,
            capture_output=True, text=True, timeout=20)

    def distribution(self, name="out0.csv"):
        with (self._directory / name).open() as inp:
            return {float(row["Dist"]): float(row["Fraction"])
                    for row in csv.DictReader(inp)}

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stderr)

    def assert_failure(self, result):
        # An ordinary failure, not SIGABRT or SIGSEGV.
        self.assertGreater(result.returncode, 0, result.stderr)

    def test_short_stream_writes_normalised_final_distribution(self):
        result = self.run_tool(_A, _A + _B)
        self.assert_success(result)
        dist = self.distribution()
        self.assertAlmostEqual(dist[0], 0.5)
        self.assertAlmostEqual(sum(dist.values()), 1.0)
        self.assertEqual(len([value for value in dist.values() if value]), 2)
        mean = float(re.search(r"Average (\S+)", result.stderr).group(1))
        self.assertAlmostEqual(mean, sum(d * f for d, f in dist.items()), places=5)

    def test_final_partial_batch_is_included(self):
        result = self.run_tool(_A, _A * 10 + _B, "-f", "10", "-n", "100", "-b", "7")
        self.assert_success(result)
        self.assertEqual(self.distribution()[0], 1)
        self.assertAlmostEqual(self.distribution("out1.csv")[0], 10 / 11, places=5)
        self.assertIn("without convergence after 11", result.stderr)

    def test_exact_sample_boundary_is_not_written_twice(self):
        result = self.run_tool(_A, _A * 10, "-f", "10")
        self.assert_success(result)
        self.assertTrue((self._directory / "out0.csv").exists())
        self.assertFalse((self._directory / "out1.csv").exists())

    def test_sampling_counts_do_not_depend_on_batch_size(self):
        for batch in (1, 7, 100):
            with self.subTest(batch=batch):
                result = self.run_tool(_A, _A * 30, "-f", "10", "-n", "3",
                                       "-b", str(batch), stem=f"batch{batch}_")
                self.assert_success(result)
                self.assertEqual(re.findall(r"Sample after (\d+)", result.stderr), ["10", "13"])
                self.assertIn("Converged after 13", result.stderr)

    def test_relative_tolerance_rejects_newly_populated_bin(self):
        result = self.run_tool(_A, _A * 10 + _B * 10, "-f", "10", "-n", "10", "-r", "0.01")
        self.assert_success(result)
        self.assertIn("without convergence after 20", result.stderr)
        self.assertEqual(self.distribution("out1.csv")[0], 0.5)

    def test_absolute_tolerance_converges(self):
        result = self.run_tool(_A, _A * 10 + _B * 10 + _A * 10,
                               "-f", "10", "-n", "10", "-a", "0.6")
        self.assert_success(result)
        self.assertIn("Converged after 20", result.stderr)
        self.assertEqual(self.distribution("out1.csv")[0], 0.5)

    def test_reference_limit_equal_to_or_larger_than_file(self):
        automatic = self.run_tool(_A + _B, _A, stem="automatic")
        self.assert_success(automatic)
        for count in (2, 3):
            with self.subTest(count=count):
                result = self.run_tool(_A + _B, _A, "-s", str(count), stem="explicit")
                self.assert_success(result)
                self.assertIn("Read 2 fingerprints", result.stderr)
                self.assertEqual(self.distribution("automatic0.csv"),
                                 self.distribution("explicit0.csv"))

    def test_reference_limit_uses_only_prefix(self):
        result = self.run_tool(_A + _B, _A, "-s", "1")
        self.assert_success(result)
        self.assertIn("Read 1 fingerprints", result.stderr)
        self.assertEqual(self.distribution()[0], 1)
        self.assertEqual(sum(self.distribution().values()), 1)

    def test_reference_tail_is_not_read_after_limit(self):
        result = self.run_tool(_A + "PCN<malformed tail>\n", _A, "-s", "1")
        self.assert_success(result)
        self.assertEqual(self.distribution()[0], 1)

    def test_empty_inputs_are_errors(self):
        for reference, queries, options in (("", _A, ()), ("", _A, ("-s", "1")),
                                            (_A, "", ())):
            with self.subTest(reference=bool(reference), queries=bool(queries), options=options):
                self.assert_failure(self.run_tool(reference, queries, *options))
        self.assertFalse((self._directory / "out0.csv").exists())

    def test_no_retained_distances_is_an_error(self):
        result = self.run_tool(_A, _B * 10, "-f", "10", "-T", "0.0001")
        self.assert_failure(result)
        self.assertIn("No distances retained", result.stderr)
        self.assertFalse((self._directory / "out0.csv").exists())

    def test_later_batches_can_supply_first_retained_distances(self):
        result = self.run_tool(_A, _B * 10 + _A * 10, "-f", "10", "-n", "10", "-T", "0.0001")
        self.assert_success(result)
        self.assertEqual(self.distribution(), {0: 1})
        self.assertIn("without convergence after 20", result.stderr)

    def test_open_failure_returns_nonzero(self):
        for queries in (_A, _A * 20):
            with self.subTest(count=queries.count("|")):
                result = self.run_tool(_A, queries, "-f", "10", "-n", "10", stem="missing/out")
                self.assert_failure(result)
                self.assertIn("cannot open", result.stderr)

    @unittest.skipUnless(pathlib.Path("/dev/full").exists(), "requires /dev/full")
    def test_buffered_write_failure_returns_nonzero(self):
        (self._directory / "out0.csv").symlink_to("/dev/full")
        result = self.run_tool()
        self.assert_failure(result)
        self.assertIn("error writing", result.stderr)

    def test_missing_components_are_errors_without_abort(self):
        for tag in ("FPIW", "FPMK", "FPMK2", "MPR"):
            with self.subTest(tag=tag):
                self.assert_failure(self.run_tool(remove_tag(_A, tag), _A))

    def test_wrong_fingerprint_width_is_rejected(self):
        for tag, width in (("FPIW", 192), ("FPMK", 2048), ("FPMK2", 256), ("MPR", 32)):
            with self.subTest(tag=tag):
                bad = remove_tag(_A, tag).replace("|\n", zero_fingerprint(tag, width) + "|\n")
                self.assert_failure(self.run_tool(bad, bad))

    def test_reordered_tags_preserve_distances(self):
        normal = self.run_tool(_A, _B, stem="normal")
        lines = _A.splitlines()
        reordered = "\n".join(reversed(lines[:-1])) + "\n|\n"
        result = self.run_tool(reordered, _B, stem="reordered")
        self.assert_success(normal)
        self.assert_success(result)
        self.assertEqual(self.distribution("normal0.csv"), self.distribution("reordered0.csv"))

    def test_query_stdin(self):
        result = self.run_tool(_A, _A + _B, stdin=True)
        self.assert_success(result)
        self.assertEqual(self.distribution()[0], 0.5)

    def test_truncated_tdt_is_not_clean_eof(self):
        for reference, queries in ((_A + _B.removesuffix("|\n"), _A),
                                   (_A, _A + _B.removesuffix("|\n"))):
            with self.subTest(reference=reference.count("|"), queries=queries.count("|")):
                result = self.run_tool(reference, queries)
                self.assert_failure(result)
                self.assertIn("Malformed", result.stderr)

    def test_invalid_numeric_options(self):
        for option, value in (("-a", "nan"), ("-r", "inf"), ("-T", "nan"),
                              ("-T", "0"), ("-b", "0"), ("-n", "0"), ("-f", "9")):
            with self.subTest(option=option, value=value):
                self.assert_failure(self.run_tool(_A, _A, option, value))

    def test_threads_give_identical_distributions(self):
        for threads in (1, 4):
            result = self.run_tool(_A + _B, _A * 11 + _B * 6,
                                   stem=f"threads{threads}_", threads=threads)
            self.assert_success(result)
        self.assertEqual(self.distribution("threads1_0.csv"), self.distribution("threads4_0.csv"))

    def test_output_cannot_overwrite_query(self):
        query = self._directory / "out0.csv"
        query.write_text(_A)
        result = subprocess.run(
            [_EXECUTABLE, "-h", "1", "-p", "out0.csv", "-S", "out", "out0.csv"],
            cwd=self._directory, capture_output=True, text=True, timeout=20)
        self.assert_failure(result)
        self.assertIn("would overwrite input", result.stderr)
        self.assertEqual(query.read_text(), _A)


if __name__ == "__main__":
    unittest.main()
