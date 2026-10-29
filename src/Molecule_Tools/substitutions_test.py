"""Regression tests for substitutions matching, properties and output safety."""
import csv
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

_EXE = str(Path(sys.argv.pop(1)).resolve())


class SubstitutionsTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._dir = Path(self._tmp.name)
        self._input = self._dir / "input.smi"

    def run_tool(self, data, *args):
        self._input.write_text(data)
        return subprocess.run([_EXE, *args, str(self._input)], cwd=self._dir,
                              text=True, capture_output=True, check=False)

    def table(self):
        text = (self._dir / "groups.tsv").read_text()
        rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
        self.assertTrue(rows)
        self.assertTrue(all(len(row) == len(rows[0]) for row in rows), text)
        return [dict(zip(rows[0], row)) for row in rows[1:]]

    def test_library_tpsa(self):
        data = ("ClCO alcohol\nClCSC thioether\nClc1ccsc1 thiophene\n"
                "ClCN([H])[H] explicit_amine\nClCP(C)C phosphine\n")
        result = self.run_tool(data, "-s", "Cl", "-x", "-w", "-d", "groups.tsv")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.table()
        self.assertEqual(len(rows), 5)
        for row, expected in zip(rows, (20.23, 25.30, 28.24, 26.02, 13.59)):
            self.assertAlmostEqual(float(row["R1_nvrtspsa"]), expected, places=2)

    def test_missing_fragments_keep_their_columns(self):
        result = self.run_tool("CN short\nCCN long\n", "-s", "N", "-x", "-C", "2",
                               "-w", "-d", "groups.tsv")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.table()
        self.assertNotEqual(rows[0]["R1_SMILES"], ".")
        self.assertEqual(rows[1]["R1_SMILES"], ".")
        self.assertEqual(rows[1]["R1_nvrtspsa"], ".")

    def test_multiple_queries_and_different_column_counts(self):
        result = self.run_tool("CCNC amine\nCCO alcohol\n", "-s", "CN", "-s", "O",
                               "-j", "0,1", "-j", "0", "-z", "first", "-z", "ignore",
                               "-x", "-w", "-d", "groups.tsv")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.table()
        self.assertEqual(len(rows), 2)
        self.assertNotEqual(rows[0]["R2_SMILES"], ".")
        self.assertEqual(rows[1]["R2_SMILES"], ".")
        result = self.run_tool("CCN sample\n", "-s", "C", "-s", "N", "-j", "0",
                               "-j", "0", "-z", "first", "-x", "-d", "groups.tsv")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.table()

    def test_distinct_matches_not_silently_collapsed(self):
        data = "NC.CN(C)C two_sites\n"
        result = self.run_tool(data, "-s", "N", "-z", "ignmmatch", "-z", "oknomatch")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("two_sites 0\n", result.stdout)
        result = self.run_tool(data, "-s", "N", "-z", "first")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("two_sites 3\n", result.stdout)

    def test_ring_closure_labels_and_parent_isotopes(self):
        result = self.run_tool("C1CCCC1 ring\n", "-s", "CCC", "-j", "0", "-r", "-x",
                               "-c", "99", "-z", "first", "-S", "fragments", "-y")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("99", (self._dir / "fragments.smi").read_text())
        result = self.run_tool("[13CH3]N labelled\n", "-s", "N", "-x", "-J", "joined.txt",
                               "-S", "fragments", "-v")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("13", (self._dir / "joined.txt").read_text().split()[0])
        parent = (self._dir / "fragments.smi").read_text().splitlines()[0]
        self.assertIn("13", parent.split()[0])
        self.assertIn("1 molecules had 1 fragments", result.stderr)
        self.assertIn("1 molecules matched 1 queries", result.stderr)

    def test_no_substituents_and_parent_suppression(self):
        result = self.run_tool("N ammonia\n", "-s", "N", "-x", "-w", "-d", "groups.tsv")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.table()
        self.assertEqual(len(rows), 1)
        self.assertNotIn("R1_SMILES", rows[0])
        result = self.run_tool("CN amine\n", "-s", "N", "-S", "fragments", "-y")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len((self._dir / "fragments.smi").read_text().splitlines()), 1)

    def test_capping_reference_for_unsubstituted_groups(self):
        cap = self._dir / "cap.smi"
        cap.write_text("C cap\n")
        result = self.run_tool("CN substituted\nN unsubstituted\n", "-s", "N", "-x",
                               "-P", str(cap), "-w", "-d", "groups.tsv")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = self.table()
        self.assertEqual(rows[1]["R1_SMILES"], "C")
        self.assertEqual(float(rows[1]["R1_natoms"]), 1)
        self.assertEqual(float(rows[1]["R1_nvrtspsa"]), 0)

    def test_output_protection_and_write_failure(self):
        data = "CN input\n"
        for flag in ("-J", "-T", "-d"):
            result = self.run_tool(data, "-s", "N", flag, str(self._input))
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(self._input.read_text(), data)
        for link in ("symlink", "hardlink"):
            alias = self._dir / (link + ".smi")
            if link == "symlink":
                alias.symlink_to(self._input)
            else:
                os.link(self._input, alias)
            result = self.run_tool(data, "-s", "N", "-J", str(alias))
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(self._input.read_text(), data)
        result = self.run_tool(data, "-s", "N", "-d", "missing/groups.tsv")
        self.assertNotEqual(result.returncode, 0)

    def test_pair_attachment_keeps_single_atom_linker(self):
        result = self.run_tool("NCN linker\n", "-s", "N.N", "-j", "0-1", "-z", "first",
                               "-x", "-S", "fragments", "-y")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len((self._dir / "fragments.smi").read_text().splitlines()), 1)
        self.assertIn("CH4", (self._dir / "fragments.smi").read_text())
        result = self.run_tool("NN direct\n", "-s", "N.N", "-j", "0-1", "-z", "first",
                               "-x", "-J", "joined.txt")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len((self._dir / "joined.txt").read_text().splitlines()), 1)

    def test_invalid_isotope_query_combinations(self):
        for args in [("-s", "N", "-m", "missing.smi"),
                     ("-m", "missing.smi", "-m", "another.smi")]:
            result = self.run_tool("CN input\n", *args)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("one -m query file", result.stderr)

    def test_malformed_input_returns_failure(self):
        result = self.run_tool("CN good\n[ bad\nCCN later\n", "-s", "N",
                               "-d", "groups.tsv")
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self._dir / "groups.tsv").exists())

    def test_invalid_attachment_index(self):
        result = self.run_tool("CN input\n", "-s", "N", "-j", "3", "-d", "groups.tsv")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Attachment index", result.stderr)


if __name__ == "__main__":
    unittest.main()
