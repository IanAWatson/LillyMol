"""Output handling in the Ruby tsubstructure parallel wrapper."""

import pathlib
import subprocess
import tempfile
import unittest


_WRAPPER = pathlib.Path(__file__).resolve().parents[1] / "tsubstructure_parallel.rb"


class ParallelOutput(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._input = self._directory / "input.smi"
        self._smiles = "".join(f"CN molecule_{i}\n" for i in range(30))
        self._input.write_text(self._smiles)
        # A minimal worker that honours the wrapper's byte ranges. No LillyMol
        # installation is needed to test splitting and joining named output.
        self._worker = self._directory / "worker"
        self._worker.write_text('''#!/usr/bin/env python3
import pathlib
import sys
args = sys.argv[1:]
start, stop = 0, None
stem = None
for i, arg in enumerate(args[:-1]):
    if arg == '-i':
        key, value = args[i + 1].split('=')
        if key == 'seek':
            start = int(value)
        elif key == 'stop':
            stop = int(value)
    elif arg == '-m' and args[i + 1] != 'QDT':
        stem = args[i + 1]
data = pathlib.Path(args[-1]).read_bytes()[start:stop]
if stem:
    pathlib.Path(stem + '.smi').write_bytes(data)
''')
        self._worker.chmod(0o755)

    def run_tool(self, *args):
        return subprocess.run(
            ["ruby", str(_WRAPPER), "-thr", "3", "-tmpdir", str(self._directory),
             "-tsub", str(self._worker), *args, str(self._input)],
            cwd=self._directory, capture_output=True, text=True, timeout=15)

    def assert_stdout_rejected(self, *args):
        result = self.run_tool(*args)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("stdout, which is not supported", result.stderr)
        self.assertIn("Use a named output stem", result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertFalse((self._directory / "-.smi").exists())
        self.assertFalse((self._directory / "results.smi").exists())

    def test_matching_stdout(self):
        self.assert_stdout_rejected("-m", "-")

    def test_nonmatching_stdout(self):
        self.assert_stdout_rejected("-n", "-")

    def test_all_occurrences_checked(self):
        for values in (("-", "QDT"), ("QDT", "-"), ("-", "results"),
                       ("results", "-")):
            with self.subTest(values=values):
                self.assert_stdout_rejected("-m", values[0], "-m", values[1])

    def test_combined_stdout(self):
        for modifier in ("QDT", "QDTVB", "SEPA", "NONM", "NONMX"):
            with self.subTest(modifier=modifier):
                self.assert_stdout_rejected("-m", modifier + ":-")

    def test_no_join_also_rejects_stdout(self):
        self.assert_stdout_rejected("-nj", "-m", "-")

    def test_stem_that_normalises_to_stdout(self):
        self.assert_stdout_rejected("-m", "-.smi")

    def test_named_output_joins_all_workers(self):
        for stem in ("results", "results.smi"):
            with self.subTest(stem=stem):
                result = self.run_tool("-m", stem, "-m", "QDT")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual((self._directory / "results.smi").read_text(),
                                 self._smiles)
                self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
