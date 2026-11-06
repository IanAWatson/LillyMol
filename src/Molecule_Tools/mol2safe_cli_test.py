"""Optional ring-system complexity and numbering guards for mol2SAFE."""
import pathlib
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class SafeRingGuards(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._input = self._directory / "input.smi"
        self._library = self._directory / "library.textproto"

    def run_tool(self, smiles, *args):
        self._input.write_text(smiles)
        result = subprocess.run(
            [_EXECUTABLE, "-I", "7", "-v", "-S", str(self._library),
             *args, str(self._input)],
            capture_output=True, text=True, check=False, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result, self._library.read_text()

    def test_small_molecule_default_is_unchanged(self):
        molecule = "COc1ccccc1 anisole\n"
        default, default_library = self.run_tool(molecule)
        checked, checked_library = self.run_tool(molecule, "-M", "maxringsys=7")
        self.assertTrue(default.stdout)
        self.assertEqual(default.stdout, checked.stdout)
        self.assertEqual(sorted(default_library.splitlines()),
                         sorted(checked_library.splitlines()))
        self.assertNotIn("skipped for ring systems", default.stderr)
        self.assertIn("0 molecules skipped for ring systems", checked.stderr)
        self.assertIn("0 molecules skipped for unsafe ring numbering", checked.stderr)

    def test_fused_ring_limit_is_optional_and_inclusive(self):
        molecule = "COc1ccc2ccccc2c1 naphthyl\n"
        default, library = self.run_tool(molecule)
        self.assertTrue(default.stdout)
        self.assertTrue(library)
        rejected, library = self.run_tool(molecule, "-M", "maxringsys=1", "-p")
        self.assertEqual(rejected.stdout, "")
        self.assertEqual(library, "")
        self.assertIn("1 molecules skipped for ring systems", rejected.stderr)
        accepted, library = self.run_tool(molecule, "-M", "maxringsys=2")
        self.assertTrue(accepted.stdout)
        self.assertTrue(library)

    def test_spiro_rings_are_one_system(self):
        molecule = "COC1CC12CC23CC3 spiro\n"
        rejected, library = self.run_tool(molecule, "-M", "maxringsys=2")
        self.assertEqual(rejected.stdout, "")
        self.assertEqual(library, "")
        self.assertIn("1 molecules skipped for ring systems", rejected.stderr)
        accepted, library = self.run_tool(molecule, "-M", "maxringsys=3")
        self.assertTrue(accepted.stdout)
        self.assertTrue(library)

    def test_seven_ring_cutoff(self):
        molecule = "COC1CC12CC23CC34CC45CC56CC67CC78CC8 eight_spiro_rings\n"
        rejected, library = self.run_tool(molecule, "-M", "maxringsys=7")
        self.assertEqual(rejected.stdout, "")
        self.assertEqual(library, "")
        self.assertIn("1 molecules skipped for ring systems", rejected.stderr)
        accepted, library = self.run_tool(molecule, "-M", "maxringsys=8")
        self.assertTrue(accepted.stdout)
        self.assertTrue(library)

    def test_limit_is_per_ring_system(self):
        accepted, library = self.run_tool(
            "COc1ccc(-c2ccccc2)cc1 biphenyl\n", "-M", "maxringsys=1")
        self.assertTrue(accepted.stdout)
        self.assertTrue(library)
        self.assertIn("0 molecules skipped for ring systems", accepted.stderr)

    def test_processing_continues_after_rejection(self):
        result, library = self.run_tool(
            "COC1CC12CC23CC3 rejected\nCOc1ccccc1 accepted\n",
            "-M", "maxringsys=2")
        self.assertNotIn("rejected", result.stdout)
        self.assertIn("accepted", result.stdout)
        self.assertNotIn('par: "rejected"', library)
        self.assertIn('par: "accepted"', library)

    def test_number_reuse_and_other_M_qualifiers(self):
        molecule = "COc1ccccc1 anisole\n"
        for extra in ((), ("-y",)):
            result, library = self.run_tool(
                molecule, "-M", "16", "-M", "maxnr=9", "-M", "maxringsys=7", *extra)
            self.assertTrue(result.stdout)
            self.assertTrue(library)
            self.assertIn("0 molecules skipped for unsafe ring numbering", result.stderr)

    def test_excessive_connection_numbers_reject_before_library_recording(self):
        molecule = "CO" * 70 + "C many_cuts\n"
        default, _ = self.run_tool(molecule)
        self.assertTrue(default.stdout)
        checked, library = self.run_tool(molecule, "-M", "maxringsys=7")
        self.assertEqual(checked.stdout, "")
        self.assertEqual(library, "")
        self.assertIn("1 molecules skipped for unsafe ring numbering", checked.stderr)

    def test_invalid_limits(self):
        self._input.write_text("COc1ccccc1 anisole\n")
        for value in ("", "0", "-1", "abc", "2147483648"):
            with self.subTest(value=value):
                result = subprocess.run(
                    [_EXECUTABLE, "-M", "maxringsys=" + value, str(self._input)],
                    capture_output=True, text=True, check=False, timeout=10)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Invalid maxringsys=", result.stderr)


if __name__ == "__main__":
    unittest.main()
