"""Executable regression tests for section-size filtering and aromatic counts."""
import pathlib
import subprocess
import sys
import tempfile
import unittest

_EXE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class HydrophobicSections(unittest.TestCase):
    def run_tool(self, smiles, *args):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "input.smi"
            path.write_text(smiles)
            result = subprocess.run([_EXE, *args, str(path)], cwd=tmp,
                                    capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = [line.split() for line in result.stdout.splitlines()]
            self.assertTrue(lines)
            self.assertTrue(all(len(row) == len(lines[0]) for row in lines[1:]))
            rows = {row[0]: dict(zip(lines[0][1:], map(float, row[1:]))) for row in lines[1:]}
            labelled = (pathlib.Path(tmp) / "labelled.smi")
            return rows, labelled.read_text() if labelled.exists() else ""

    def test_simple_rules_exclude_multiple_bonded_heteroatom_carbons(self):
        rows, labelled = self.run_tool(
            "CC=O aldehyde\nC=O formaldehyde\nCC=N imine\nCC(=O)C ketone\n"
            "CC#N nitrile\nC=C alkene\nc1ccccc1 aromatic\n",
            "-G", "simple", "-G", "iso=7", "-S", "labelled")
        for name, expected in (("aldehyde", 1), ("formaldehyde", 0), ("imine", 1),
                               ("ketone", 2), ("nitrile", 1), ("alkene", 2),
                               ("aromatic", 6)):
            self.assertEqual(rows[name]["hpo_HPOnumbr"], expected, name)
        aldehyde = next(line for line in labelled.splitlines() if line.endswith(" aldehyde"))
        self.assertTrue(aldehyde.startswith("[7CH3]C=O"), aldehyde)
        # Explicit queries can deliberately override the built-in selection.
        rows, _ = self.run_tool("CC=O aldehyde\n", "-G", "simple", "-G", "[#6]")
        self.assertEqual(rows["aldehyde"]["hpo_HPOnumbr"], 2)

    def test_aromatic_halogens_remain_connected_to_ring(self):
        rows, labelled = self.run_tool(
            "Fc1ccccc1 fluoro\nClc1ccccc1 chloro\nBrc1ccccc1 bromo\n"
            "Ic1ccccc1 iodo\nOc1ccccc1 phenol\nNc1ccccc1 aniline\n"
            "C1=CC(O)=C(C=C1Br)CNC1=CC=C(C=C1)OC CHEMBL1808838\n",
            "-G", "simple", "-G", "minsize=2", "-G", "iso=1", "-S", "labelled")
        for name in ("fluoro", "chloro", "bromo", "iodo"):
            self.assertEqual(rows[name]["hpo_HPOnumbr"], 7, name)
            self.assertEqual(rows[name]["hpo_HPOnsect"], 1, name)
            self.assertEqual(rows[name]["hpo_HPOmaxsz"], 7, name)
        for name in ("phenol", "aniline"):
            self.assertEqual(rows[name]["hpo_HPOnumbr"], 5, name)
        for element in ("F", "Cl", "Br", "I"):
            self.assertIn("[1" + element + "]", labelled)
        example = next(line for line in labelled.splitlines()
                       if line.endswith(" CHEMBL1808838"))
        self.assertIn("[1Br]", example)
        self.assertTrue(example.endswith("OC CHEMBL1808838"), example)

    def test_substituted_hydrocarbon_unsaturation(self):
        rows, labelled = self.run_tool(
            "C=C(C)C isobutene\nCC(C)=C(C)C substituted_alkene\n"
            "CC#CC alkyne\nC=C1CCCCC1 exocyclic\n"
            "C1=CC(=O)C2=C(C1=O)C[C@@]1(C)[C@H](C(=C)CCC1)C2 CHEMBL401659\n",
            "-G", "simple", "-G", "minsize=2", "-G", "iso=1", "-S", "labelled")
        for name, expected in (("isobutene", 4), ("substituted_alkene", 6),
                               ("alkyne", 4), ("exocyclic", 7)):
            self.assertEqual(rows[name]["hpo_HPOnumbr"], expected, name)
            self.assertEqual(rows[name]["hpo_HPOnsect"], 1, name)
        example = next(line for line in labelled.splitlines()
                       if line.endswith(" CHEMBL401659"))
        self.assertIn("[1C](=[1CH2])", example)

    def test_aromatic_carbon_next_to_polar_carbon(self):
        rows, labelled = self.run_tool(
            "N=Cc1ccccc1 imine\nO=Cc1ccccc1 aldehyde\nN#Cc1ccccc1 nitrile\n"
            "NC(=O)c1ccccc1 amide\n"
            "O=C(NN=CC1=CC=C(N(C)C)C=C1)C1=CC=CC=N1 CHEMBL1987655\n",
            "-G", "simple", "-G", "minsize=2", "-G", "iso=1", "-S", "labelled")
        for name in ("imine", "aldehyde", "nitrile", "amide"):
            self.assertEqual(rows[name]["hpo_HPOnumbr"], 6, name)
            self.assertEqual(rows[name]["hpo_HPOmxarm"], 6, name)
            line = next(line for line in labelled.splitlines() if line.endswith(" " + name))
            self.assertEqual(line.count("[1"), 6, line)
        self.assertEqual(rows["CHEMBL1987655"]["hpo_HPOnumbr"], 11)
        example = next(line for line in labelled.splitlines()
                       if line.endswith(" CHEMBL1987655"))
        self.assertIn("NN=C[1C]1", example)
        self.assertTrue(example.startswith("O=C("), example)

    def test_neutral_aromatic_sulfur(self):
        rows, labelled = self.run_tool(
            "c1ccsc1 thiophene\nc1ccc2sccc2c1 benzothiophene\n"
            "c1nccs1 thiazole\n[s+]1ccccc1 charged\nCS(=O)C sulfoxide\n"
            "CS(=O)(=O)C sulfone\n",
            "-G", "simple", "-G", "minsize=2", "-G", "iso=1", "-L", "simple",
            "-S", "labelled")
        for name, expected in (("thiophene", 5), ("benzothiophene", 9), ("thiazole", 4)):
            self.assertEqual(rows[name]["hpo_HPOnumbr"], expected, name)
            self.assertEqual(rows[name]["hpo_HPOnsect"], 1, name)
            self.assertEqual(rows[name]["hpo_HPInumbr"], name == "thiazole", name)
            line = next(line for line in labelled.splitlines() if line.endswith(" " + name))
            self.assertIn("[1S]", line)
        self.assertEqual(rows["charged"]["hpo_HPInumbr"], 1)
        self.assertEqual(rows["sulfoxide"]["hpo_HPInumbr"], 2)
        self.assertEqual(rows["sulfone"]["hpo_HPInumbr"], 3)
        for name in ("charged", "sulfoxide", "sulfone"):
            line = next(line for line in labelled.splitlines() if line.endswith(" " + name))
            self.assertNotIn("[1S", line)

    def test_isotope_mask_splits_sections_and_clears_labels(self):
        rows, labelled = self.run_tool(
            "[13CH3][13CH2]C[13CH2][13CH3] sample\n", "-I", "iso",
            "-G", "simple", "-G", "minsize=2", "-G", "iso=7", "-S", "labelled")
        row = rows["sample"]
        self.assertEqual(row["hpo_HPOnumbr"], 4)
        self.assertEqual(row["hpo_HPOratio"], 1)
        self.assertEqual(row["hpo_HPOnsect"], 2)
        self.assertEqual(row["hpo_HPOmaxsz"], 2)
        self.assertEqual(row["hpo_HPOmxntr"], 0.5)
        self.assertEqual(labelled.count("[7"), 4)
        self.assertNotIn("[13", labelled)

    def test_unlabelled_mask_and_neither_counts(self):
        rows, labelled = self.run_tool(
            "[13CH3]CCO sample\n", "-I", "noiso", "-G", "simple", "-L", "[N]",
            "-G", "iso=7", "-S", "labelled")
        row = rows["sample"]
        self.assertEqual(row["hpo_HPOnumbr"], 2)
        self.assertAlmostEqual(row["hpo_HPOratio"], 0.667, places=3)
        self.assertEqual(row["hpo_NTnhphpi"], 1)
        self.assertAlmostEqual(row["hpo_NTfnhopi"], 0.333, places=3)
        self.assertEqual(labelled.count("[7"), 2)
        self.assertNotIn("[13", labelled)

    def test_isotope_mask_constrains_complement_and_queries(self):
        for args in (("-x",), ("-L", "[O]")):
            rows, labelled = self.run_tool(
                "[13CH3][18OH] sample\nCO ignored\n", "-I", "iso",
                "-G", "[#6]", *args, "-G", "iso=7", "-L", "iso=9",
                "-S", "labelled")
            self.assertEqual(rows["sample"]["hpo_HPOnumbr"], 1)
            self.assertEqual(rows["sample"]["hpo_HPInumbr"], 1)
            self.assertEqual(rows["sample"]["hpo_HPOratio"], 0.5)
            self.assertEqual(rows["sample"]["hpo_HPIratio"], 0.5)
            self.assertTrue(all(v == 0 for v in rows["ignored"].values()))
            self.assertNotIn("[13", labelled)
            self.assertNotIn("[18", labelled)
            self.assertIn("[7CH3][9OH]", labelled)
        rows, _ = self.run_tool("[13CH3][13CH3] ignored\n", "-I", "noiso", "-x")
        self.assertTrue(all(v == 0 for v in rows["ignored"].values()))

    def test_input_isotopes_cleared_without_output_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "input.smi"
            path.write_text("[13CH3][13CH3] sample\n")
            for args, expected in ((("-I", "iso"), "CC"), ((), "[13CH3][13CH3]")):
                result = subprocess.run([_EXE, *args, "-Y", "wpipe", str(path)],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.splitlines()[1].split()[0], expected)

    def test_invalid_isotope_selection(self):
        for args in (("-I", "iso", "-I", "noiso"), ("-I", "bad")):
            result = subprocess.run([_EXE, *args], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("-I", result.stderr)

    def test_stable_descriptor_formatting(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "input.smi"
            path.write_text("CCCOCCOCCO sample\n")
            result = subprocess.run([_EXE, "-G", "[#6]", "-L", "[O]", str(path)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            header, values = (line.split() for line in result.stdout.splitlines())
            row = dict(zip(header, values))
            self.assertEqual(row["hpo_HPOnumbr"], "7")
            self.assertEqual(row["hpo_HPOratio"], "0.7")
            self.assertEqual(row["hpo_HPOavesz"], "2.333")
            self.assertEqual(row["hpo_HPOmaxsz"], "3")
            self.assertEqual(row["hpo_HPIratio"], "0.3")
            for value in values[1:]:
                self.assertRegex(value, r"^-?\d+(?:\.\d{1,3})?$")
            # The separately written hydrophobic/hydrophilic ratio uses the
            # same formatter as the descriptor arrays.
            self.assertEqual(row["hpo_OIratio"], "2.333")

    def test_aromatic_atom_counts(self):
        rows, _ = self.run_tool("c1ccccc1 aromatic\nCCCCCC aliphatic\nCCCCc1ccccc1 mixed\n",
                                "-G", "[#6]")
        for name, expected in (("aromatic", 6), ("aliphatic", 0), ("mixed", 6)):
            self.assertEqual(rows[name]["hpo_HPOmxarm"], expected)
            self.assertEqual(rows[name]["hpo_HPOavarm"], expected)

    def test_minimum_affects_counts_components_and_isotope_labels(self):
        rows, labelled = self.run_tool("COCCC sample\n", "-G", "[#6]", "-G", "minsize=2",
                                       "-G", "iso=7", "-S", "labelled")
        row = rows["sample"]
        self.assertEqual(row["hpo_HPOnumbr"], 3)
        self.assertAlmostEqual(row["hpo_HPOratio"], 0.6)
        self.assertEqual(row["hpo_HPOnsect"], 1)
        for descriptor in ("minsz", "avesz", "maxsz", "lgtrm"):
            self.assertEqual(row["hpo_HPO" + descriptor], 3)
        self.assertEqual(labelled.count("[7"), 3)
        self.assertTrue(labelled.startswith("CO"), labelled)

    def test_all_sections_below_threshold(self):
        rows, _ = self.run_tool("COCC sample\n", "-G", "[#6]", "-G", "minsize=3")
        self.assertTrue(all(value == 0 for value in rows["sample"].values()))

    def test_hydrophilic_minimum_and_complement(self):
        rows, _ = self.run_tool("COCCC sample\n", "-G", "[#6]", "-G", "minsize=2", "-x")
        self.assertEqual(rows["sample"]["hpo_HPInumbr"], 2)
        self.assertEqual(rows["sample"]["hpo_HPImaxsz"], 2)
        rows, _ = self.run_tool("OCOCC sample\n", "-G", "[#6]", "-L", "[O]",
                                "-L", "minsize=2")
        self.assertEqual(rows["sample"]["hpo_HPOnumbr"], 3)
        for name, value in rows["sample"].items():
            if name.startswith("hpo_HPI"):
                self.assertEqual(value, 0, name)
        self.assertEqual(rows["sample"]["hpo_OInumsep"], 0)


if __name__ == "__main__":
    unittest.main()
