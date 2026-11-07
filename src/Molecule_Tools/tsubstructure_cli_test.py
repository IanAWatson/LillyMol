"""Output modifier and filename handling in tsubstructure."""

import pathlib
import subprocess
import sys
import tempfile
import unittest


_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class OutputModifiers(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        (self._directory / "input.smi").write_text("CN methylamine\nCC ethane\n")

    def run_tool(self, *args):
        return subprocess.run(
            [_EXECUTABLE, "-s", "N nitrogen", *args, "input.smi"],
            cwd=self._directory, capture_output=True, text=True, timeout=15)

    def outputs(self, stem):
        return {p.name.replace(stem, "output", 1): p.read_text()
                for p in self._directory.glob(stem + "*.smi")}

    def test_combined_matches_both_separate_orders(self):
        for modifier in ("QDT", "QDTVB", "SEPA", "NONM", "NONMX"):
            with self.subTest(modifier=modifier):
                for args in (("-m", "before", "-m", modifier),
                             ("-m", modifier, "-m", "after"),
                             ("-m", modifier + ":combined")):
                    result = self.run_tool(*args)
                    self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(self.outputs("combined"))
                self.assertEqual(self.outputs("before"), self.outputs("combined"))
                self.assertEqual(self.outputs("after"), self.outputs("combined"))
                for stem in ("before", "after", "combined"):
                    for path in self._directory.glob(stem + "*.smi"):
                        path.unlink()

    def test_absolute_filename_and_existing_suffix(self):
        path = self._directory / "hits.smi"
        result = self.run_tool("-m", "QDT:" + str(path))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("methylamine", path.read_text())
        self.assertFalse(path.with_suffix(".smi.smi").exists())

    def test_stdout_and_additional_modifier(self):
        combined = self.run_tool("-m", "QDT:-", "-m", "NONMX")
        separate = self.run_tool("-m", "NONMX", "-m", "-", "-m", "QDT")
        self.assertEqual(combined.returncode, 0, combined.stderr)
        self.assertEqual(separate.returncode, 0, separate.stderr)
        self.assertEqual(combined.stdout, separate.stdout)
        self.assertIn("ethane", combined.stdout)
        self.assertIn("matches to", combined.stdout)

    def test_modifier_requires_filename(self):
        for modifier in ("QDT", "QDTVB", "SEPA", "NONM", "NONMX"):
            with self.subTest(modifier=modifier):
                result = self.run_tool("-m", modifier)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("require an output filename", result.stderr)

    def test_empty_combined_filename(self):
        for modifier in ("QDT", "QDTVB", "SEPA", "NONM", "NONMX"):
            with self.subTest(modifier=modifier):
                result = self.run_tool("-m", modifier + ":")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Missing output filename", result.stderr)

    def test_conflicting_filenames(self):
        for args in (("-m", "first", "-m", "QDT:second"),
                     ("-m", "QDT:first", "-m", "second"),
                     ("-m", "QDT:first", "-m", "NONMX:second")):
            with self.subTest(args=args):
                result = self.run_tool(*args)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("first", result.stderr)
                self.assertIn("second", result.stderr)
                self.assertFalse(self.outputs("first"))

    def test_csr_keeps_separator_meaning(self):
        result = self.run_tool("-m", "hits", "-m", "CSR:,")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(",nitrogen:1", (self._directory / "hits.smi").read_text())

    def test_reserved_word_as_literal_filename(self):
        result = self.run_tool("-m", "./QDT")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self._directory / "QDT.smi").exists())

    def test_help_explains_combined_form(self):
        result = self.run_tool("-m", "help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("MODIFIER:<filename>", result.stderr)
        self.assertIn("separator, not a filename", result.stderr)


if __name__ == "__main__":
    unittest.main()
