"""Command-line regression tests for additive chemotype isotope labels."""
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class ChemotypeIsotopes(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._input = self._directory / "input.smi"
        self._input.write_text("CC1(C)CCNCC1 example\n")

    def run_tool(self, *args):
        return subprocess.run(
            [_EXECUTABLE, "-s", "[N]", "-i", "smi", *args, str(self._input)],
            cwd=self._directory, text=True, capture_output=True, check=False)

    def labels(self, text):
        return sorted(int(x) for x in re.findall(r"\[(\d+)", text))

    def test_additive_labels_complement_and_qualifier_order(self):
        for qualifiers in [("-I", "ring=10,20", "-I", "2"),
                           ("-I", "2", "-I", "ring=10,20")]:
            result = self.run_tool(*qualifiers, "-C", "complement")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(self.labels(result.stdout), [10] * 5 + [12])
            self.assertEqual(self.labels(
                (self._directory / "complement.smi").read_text()), [2, 2])

    def test_repeat_last_base_and_zero_base(self):
        self._input.write_text("N1CCC(c2ccccc2)C(C2CC2)C1 example\n")
        result = self.run_tool("-n", "2", "-I", "ring=10,20")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.labels(result.stdout), [10] * 6 + [20] * 9)
        result = self.run_tool("-n", "2", "-I", "ring=0,20")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.labels(result.stdout), [20] * 9)

    def test_invalid_lists_conflicts_and_overflow(self):
        for value in ["ring=", "ring=,10", "ring=10,", "ring=10,,20",
                      "ring=-1", "ring=abc", "ring=4294967296", "0"]:
            with self.subTest(value=value):
                result = self.run_tool("-I", value)
                self.assertNotEqual(result.returncode, 0)
        for extra in [("-P", "Z"), ("-J", "FPCTYPE"),
                      ("-I", "ring=20"), ("-I", "4294967295")]:
            with self.subTest(extra=extra):
                result = self.run_tool("-I", "ring=10", *extra)
                self.assertNotEqual(result.returncode, 0)

    def test_whole_molecule_labels_and_zero(self):
        for inside, outside in [(1, 0), (10, 20), (0, 2), (0, 0)]:
            result = self.run_tool("-I", f"whole={inside},{outside}")
            self.assertEqual(result.returncode, 0, result.stderr)
            expected = ([inside] * 6 if inside else []) + ([outside] * 2 if outside else [])
            self.assertEqual(self.labels(result.stdout), sorted(expected))
            self.assertEqual(len(result.stdout.splitlines()), 1)
        result = self.run_tool("-I", "whole=10,20", "-C", "complement")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.labels(
            (self._directory / "complement.smi").read_text()), [20, 20])

    def test_invalid_whole_labels_and_conflicts(self):
        for value in ["whole=", "whole=1", "whole=1,", "whole=,2",
                      "whole=1,2,3", "whole=1,,2", "whole=-1,0",
                      "whole=1,4294967296", "whole=x,0"]:
            result = self.run_tool("-I", value)
            self.assertNotEqual(result.returncode, 0, value)
        for extra in [("-I", "ring=10"), ("-I", "2"), ("-I", "whole=2,0"),
                      ("-P", "Z"), ("-J", "FPCTYPE"), ("-F", "summary"),
                      ("-H", "def"), ("-p", ".")]:
            result = self.run_tool("-I", "whole=1,0", *extra)
            self.assertNotEqual(result.returncode, 0, extra)


if __name__ == "__main__":
    unittest.main()
