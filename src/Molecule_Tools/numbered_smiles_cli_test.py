"""Error reporting and buffered output checks for numbered_smiles."""

import pathlib
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class NumberedSmilesErrors(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._input = self._directory / "input.smi"
        self._input.write_text("CCO example\n")

    def run_tool(self, *options, inputs=("input.smi",)):
        return subprocess.run([_EXECUTABLE, *options, *inputs],
                              cwd=self._directory, capture_output=True,
                              text=True, timeout=15)

    def test_valid_input_succeeds(self):
        result = self.run_tool("-S", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("example", result.stdout)

    def test_missing_input_does_not_create_output(self):
        result = self.run_tool("-S", "numbered", inputs=("missing.smi",))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Cannot open", result.stderr)
        self.assertFalse((self._directory / "numbered.smi").exists())

    def test_missing_later_input_is_not_success(self):
        result = self.run_tool("-S", "-", inputs=("input.smi", "missing.smi"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("example", result.stdout)

    def test_missing_first_input_stops_processing(self):
        result = self.run_tool("-S", "-", inputs=("missing.smi", "input.smi"))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_malformed_input_is_not_success(self):
        self._input.write_text("CCO valid\ninvalid_smiles bad\nCC later\n")
        result = self.run_tool("-S", "-")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("valid", result.stdout)
        self.assertNotIn("later", result.stdout)

    def test_requested_input_error_tolerance_still_works(self):
        self._input.write_text("invalid_smiles bad\nCC valid\n")
        result = self.run_tool("-i", "ICTE", "-S", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("valid", result.stdout)

    def test_output_open_failure_is_not_success(self):
        result = self.run_tool("-S", "missing_directory/output")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Cannot set stem", result.stderr)

    def test_default_output_open_failure_is_not_success(self):
        (self._directory / "input_num.smi").mkdir()
        result = self.run_tool()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Cannot set new stem", result.stderr)

    def test_invalid_options_are_errors(self):
        for options in (("-K", "not_a_smiles_option"), ("-Y", "not_an_option")):
            with self.subTest(options=options):
                result = self.run_tool(*options, "-S", "numbered")
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((self._directory / "numbered.smi").exists())

    def test_help_still_succeeds(self):
        result = self.run_tool("-Y", "help", inputs=())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("spiro", result.stderr)

    def test_removed_options_are_rejected(self):
        for option in ("-t", "-p", "-q", "-s"):
            self.assertNotEqual(self.run_tool(option).returncode, 0)

    @unittest.skipUnless(pathlib.Path("/dev/full").exists(), "requires /dev/full")
    def test_buffered_file_write_failure_is_not_success(self):
        (self._directory / "full.smi").symlink_to("/dev/full")
        result = self.run_tool("-S", "full")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("error writing output", result.stderr)

    @unittest.skipUnless(pathlib.Path("/dev/full").exists(), "requires /dev/full")
    def test_output_error_is_checked_before_switching_files(self):
        (self._directory / "input_num.smi").symlink_to("/dev/full")
        (self._directory / "second.smi").write_text("CC second\n")
        result = self.run_tool(inputs=("input.smi", "second.smi"))
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self._directory / "second_num.smi").exists())

    @unittest.skipUnless(pathlib.Path("/dev/full").exists(), "requires /dev/full")
    def test_buffered_stdout_failure_is_not_success(self):
        with open("/dev/full", "w") as output:
            result = subprocess.run([_EXECUTABLE, "-S", "-", "input.smi"],
                                    cwd=self._directory, stdout=output,
                                    stderr=subprocess.PIPE, text=True, timeout=15)
        self.assertNotEqual(result.returncode, 0)


class NumberedSmilesNumbering(unittest.TestCase):
    setUp = NumberedSmilesErrors.setUp
    run_tool = NumberedSmilesErrors.run_tool

    def convert(self, smiles, *options):
        self._input.write_text(smiles + " example\n")
        result = self.run_tool(*options, "-S", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout

    def sdf_atoms(self, text):
        lines = text.splitlines()
        natoms = int(lines[3][:3])
        nbonds = int(lines[3][3:6])
        elements = [line.split()[3] for line in lines[4:4 + natoms]]
        isotopes = [0] * natoms
        degrees = [0] * natoms
        for line in lines[4 + natoms:4 + natoms + nbonds]:
            a, b = int(line[:3]) - 1, int(line[3:6]) - 1
            degrees[a] += 1
            degrees[b] += 1
        for line in lines:
            if line.startswith("M  ISO"):
                fields = line.split()
                for i in range(int(fields[2])):
                    isotopes[int(fields[3 + 2 * i]) - 1] = int(fields[4 + 2 * i])
        return elements, isotopes, degrees

    def test_atom_number_step_is_an_increment(self):
        self.assertIn("[1CH3][4CH2][7OH]", self.convert("CCO", "-n", "1", "-e", "3"))
        self.assertIn("C[3CH2][6OH]", self.convert("CCO", "-n", "0", "-e", "3"))

    def test_numeric_start_selects_atom_number_mode(self):
        self.assertIn("[1CH3][4OH]", self.convert("CO", "-n", "symm", "-n", "1", "-e", "3"))

    def test_map_numbering_is_independent_of_isotopes(self):
        for extra in ((), ("-I",)):
            with self.subTest(extra=extra):
                text = self.convert("[13CH3:99][18OH:88]", "-m", "-n", "1", "-e", "3", *extra)
                self.assertIn("[13CH3:1][18OH:4]", text)

    def test_existing_isotopes_are_preserved_unless_overwritten(self):
        self.assertIn("[13CH3][4OH]", self.convert("[13CH3]O", "-n", "1", "-e", "3"))
        self.assertIn("[1CH3][4OH]", self.convert("[13CH3]O", "-n", "1", "-e", "3", "-I"))

    def test_other_numbering_modes_use_step_as_multiplier(self):
        self.assertIn("[18CH3][24OH]", self.convert("CO", "-n", "z", "-e", "3"))
        self.assertIn("[CH3:18][OH:24]", self.convert("CO", "-n", "z", "-e", "3", "-m"))

    def test_negative_starts_are_rejected(self):
        result = self.run_tool("-n", "-1", "-S", "-")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("non-negative", result.stderr)
        self.assertEqual(result.stdout, "")

    def test_excessive_labels_are_rejected_without_wraparound(self):
        for options in (("-n", "2147483647", "-e", "2147483647"),
                        ("-m", "-n", "2147483647"),
                        ("-n", "z", "-e", "2147483647")):
            with self.subTest(options=options):
                result = self.run_tool(*options, "-S", "-")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("supported label range", result.stderr)
                self.assertEqual(result.stdout, "")

    def test_sorting_with_nonzero_start_and_step(self):
        text = self.convert("OCN", "-n", "1", "-e", "3", "-r", "-o", "sdf")
        self.assertEqual(self.sdf_atoms(text)[1], [1, 4, 7])

    def test_atomic_number_sort_is_stable_for_ties(self):
        text = self.convert("OCC", "-n", "z", "-r", "-o", "sdf")
        elements, isotopes, degrees = self.sdf_atoms(text)
        self.assertEqual(elements, ["C", "C", "O"])
        self.assertEqual(isotopes, [6, 6, 8])
        self.assertEqual(degrees, [2, 1, 1])

    def test_map_sort_preserves_isotopes_and_uses_maps(self):
        text = self.convert("[11OH][12CH2][13NH2]", "-m", "-n", "z", "-r", "-o", "sdf")
        elements, isotopes, _ = self.sdf_atoms(text)
        self.assertEqual(elements, ["C", "N", "O"])
        self.assertEqual(isotopes, [12, 13, 11])

    def test_sorting_uses_preserved_isotope_labels(self):
        text = self.convert("[13CH3]O", "-n", "1", "-r", "-o", "sdf")
        self.assertEqual(self.sdf_atoms(text)[:2], (["O", "C"], [2, 13]))

    def test_sorting_handles_every_numbering_scheme(self):
        cases = (("-n", "symm"), ("-n", "canon"), ("-n", "ring"),
                 ("-n", "z"), ("-P", "UST:Z"))
        for options in cases:
            with self.subTest(options=options):
                text = self.convert("CC1CC1O.CN", *options, "-r", "-o", "sdf")
                _, isotopes, _ = self.sdf_atoms(text)
                self.assertEqual(isotopes, sorted(isotopes))

    def test_sorting_acyclic_ring_labels_keeps_input_order(self):
        text = self.convert("OCN", "-n", "ring", "-r", "-o", "sdf")
        self.assertEqual(self.sdf_atoms(text)[:2], (["O", "C", "N"], [0, 0, 0]))


class NumberedSmilesOutputNames(unittest.TestCase):
    setUp = NumberedSmilesErrors.setUp
    run_tool = NumberedSmilesErrors.run_tool

    def test_filename_retains_earlier_periods(self):
        (self._directory / "sample.v2.smi").write_text("CC example\n")
        result = self.run_tool(inputs=("sample.v2.smi",))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self._directory / "sample.v2_num.smi").exists())
        self.assertFalse((self._directory / "sample_num.smi").exists())

    def test_directory_and_filename_periods_are_preserved(self):
        directory = self._directory / "run.v2"
        directory.mkdir()
        (directory / "sample.v3.smi").write_text("CC example\n")
        result = self.run_tool(inputs=("run.v2/sample.v3.smi",))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((directory / "sample.v3_num.smi").exists())
        self.assertFalse((self._directory / "run_num.smi").exists())

    def test_dotted_directories_produce_distinct_outputs(self):
        for name in ("run.a", "run.b"):
            directory = self._directory / name
            directory.mkdir()
            (directory / "sample.smi").write_text(f"CC {name}\n")
        result = self.run_tool(inputs=("run.a/sample.smi", "run.b/sample.smi"))
        self.assertEqual(result.returncode, 0, result.stderr)
        for name in ("run.a", "run.b"):
            self.assertIn(name, (self._directory / name / "sample_num.smi").read_text())

    def test_input_without_extension(self):
        (self._directory / "plain").write_text("CC example\n")
        result = self.run_tool("-i", "smi", inputs=("plain",))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self._directory / "plain_num.smi").exists())

    def test_multiple_output_formats_keep_the_same_dotted_stem(self):
        (self._directory / "sample.v2.smi").write_text("CC example\n")
        result = self.run_tool("-o", "smi", "-o", "sdf", inputs=("sample.v2.smi",))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self._directory / "sample.v2_num.smi").exists())
        self.assertTrue((self._directory / "sample.v2_num.sdf").exists())

    def test_explicit_stem_combines_input_files(self):
        (self._directory / "second.smi").write_text("CN second\n")
        result = self.run_tool("-S", "combined", inputs=("input.smi", "second.smi"))
        self.assertEqual(result.returncode, 0, result.stderr)
        text = (self._directory / "combined.smi").read_text()
        self.assertIn("example", text)
        self.assertIn("second", text)
        self.assertFalse((self._directory / "input_num.smi").exists())

    def test_default_output_must_not_replace_later_input(self):
        later = self._directory / "input_num.smi"
        original = "N original\n"
        later.write_text(original)
        result = self.run_tool(inputs=("input.smi", "input_num.smi"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("would overwrite input", result.stderr)
        self.assertEqual(later.read_text(), original)
        self.assertFalse((self._directory / "input_num_num.smi").exists())

    def test_two_default_outputs_must_not_coincide(self):
        for suffix in ("smi", "sdf"):
            (self._directory / f"pair.{suffix}").write_text("CC example\n")
        result = self.run_tool("-i", "smi", inputs=("pair.smi", "pair.sdf"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("coincide", result.stderr)
        self.assertFalse((self._directory / "pair_num.smi").exists())

    def test_default_output_symlink_must_not_replace_input(self):
        (self._directory / "input_num.smi").symlink_to(self._input)
        original = self._input.read_text()
        result = self.run_tool()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self._input.read_text(), original)

    def test_explicit_output_must_not_replace_input(self):
        original = self._input.read_text()
        result = self.run_tool("-S", "./input")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self._input.read_text(), original)

    def test_output_comparison_error_does_not_open_output(self):
        (self._directory / "loop").symlink_to("loop")
        result = self.run_tool("-S", "loop/output")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Cannot compare file paths", result.stderr)
        self.assertEqual(self._input.read_text(), "CCO example\n")

    def test_bond_symmetry_uses_stderr(self):
        baseline = self.run_tool("-S", "-")
        self.assertEqual(baseline.returncode, 0, baseline.stderr)
        for mode in ("small", "large"):
            with self.subTest(mode=mode):
                result = self.run_tool("-b", mode, "-S", "-")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout, baseline.stdout)
                lines = result.stderr.splitlines()
                self.assertEqual(lines[0], "example")
                self.assertEqual(lines[3:], baseline.stderr.splitlines())
                for line in lines[1:3]:
                    self.assertRegex(line, r"^\d+ \(\d+\) \d+ \(\d+\) \d+$")


if __name__ == "__main__":
    unittest.main()
