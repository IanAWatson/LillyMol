"""Regression tests for SAFE attachment markers alongside other isotopes."""
import json
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class AttachmentIsotopes(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._input = self._directory / "input.smi"
        self._library = self._directory / "library.textproto"
        self._input.write_text(
            "[13CH3][1CH2]%10.[1CH2]%10[2CH2][13CH3] parent\n")
        self._library.write_text(
            'smi: "[1NH][2CH2][15CH3]" par: "replacement"\n')

    def run_tool(self, *args, library=True, limit=100, mode="e"):
        command = [_EXECUTABLE, "-" + mode, str(limit)]
        if library:
            command.extend(["-L", str(self._library)])
        return subprocess.run(
            [*command, *args, str(self._input)],
            cwd=self._directory, capture_output=True, text=True, check=False, timeout=5)

    def labels(self, result):
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 1, result.stderr)
        smiles = lines[0].split()[0]
        return sorted(map(int, re.findall(r"\[(\d+)", smiles)))

    def test_preserves_unrelated_isotopes(self):
        self.assertEqual(self.labels(self.run_tool("-j", "1")),
                         [1, 1, 2, 13, 15])

    def test_remove_only_attachment_isotopes(self):
        self.assertEqual(self.labels(self.run_tool("-j", "1", "-I")),
                         [2, 13, 15])

    def test_safe_output_retains_labels_and_connections(self):
        result = self.run_tool("-j", "1", "-I", "-s")
        self.assertEqual(self.labels(result), [1, 1, 2, 13, 15])
        self.assertEqual(result.stdout.count("%10"), 2)

    def test_configuration_library_uses_marker(self):
        config = self._directory / "config.textproto"
        config.write_text("library: " + json.dumps(str(self._library)) + "\n")
        self.assertEqual(self.labels(self.run_tool(
            "-C", str(config), "-j", "1", "-I", library=False)), [2, 13, 15])

    def test_default_behavior_is_preserved(self):
        self._input.write_text("[1CH3]%10.[1CH2]%10CC parent\n")
        self._library.write_text('smi: "[1NH]CC" par: "replacement"\n')
        default = self.run_tool()
        explicit = self.run_tool("-j", "1")
        self.assertEqual(self.labels(default), [1, 1])
        self.assertEqual(default.stdout, explicit.stdout)
        self.assertEqual(self.labels(self.run_tool("-I")), [])

    def test_every_entry_in_one_bucket_is_visited(self):
        self._library.write_text(
            'smi: "[1NH]CC" par: "amine"\n'
            'smi: "[1O]CC" par: "ether"\n')
        result = self.run_tool("-j", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 2, result.stdout)
        self.assertIn(" %% amine.1.3", result.stdout)
        self.assertIn(" %% ether.1.3", result.stdout)

    def test_every_bucket_is_visited(self):
        self._library.write_text(
            'smi: "[1NH]CC" par: "amine3"\n'
            'smi: "[1O]CC" par: "ether3"\n'
            'smi: "[1NH]CCC" par: "amine4"\n'
            'smi: "[1O]CCC" par: "ether4"\n')
        config = self._directory / "config.textproto"
        config.write_text("extra_atoms: 1\nfewer_atoms: 0\n")
        result = self.run_tool("-j", "1", "-C", str(config))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 6, result.stdout)
        for name, size in (("amine3", 3), ("ether3", 3),
                           ("amine4", 4), ("ether4", 4)):
            self.assertIn(f" %% {name}.1.{size}", result.stdout)

    def test_exhaustive_limits(self):
        self._library.write_text(
            'smi: "[1NH]CC" par: "amine"\n'
            'smi: "[1NH]CC" par: "duplicate"\n'
            'smi: "[1O]CC" par: "ether"\n')
        for limit, expected in ((0, 2), (1, 1), (2, 2), (10, 2)):
            with self.subTest(limit=limit):
                result = self.run_tool("-j", "1", limit=limit)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(len(result.stdout.splitlines()), expected, result.stdout)

    def test_exhaustive_limit_applies_to_each_library(self):
        self._library.write_text(
            'smi: "[1NH]CC" par: "first_amine"\n'
            'smi: "[1O]CC" par: "first_ether"\n')
        second = self._directory / "second.textproto"
        second.write_text(
            'smi: "[1CH2]CO" par: "second_alcohol"\n'
            'smi: "[1CH2]CN" par: "second_amine"\n')
        result = self.run_tool("-j", "1", "-L", str(second), limit=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 2, result.stdout)
        self.assertTrue(any(" %% first_" in line for line in lines), result.stdout)
        self.assertTrue(any(" %% second_" in line for line in lines), result.stdout)

    def config(self, text):
        config = self._directory / "config.textproto"
        config.write_text(text)
        return str(config)

    def test_invalid_generation_inputs(self):
        result = self.run_tool("-j", "1", library=False, mode="n", limit=1)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("nonempty fragment library", result.stderr)
        result = self.run_tool("-j", "1", library=False, mode="b", limit=1)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("at least two", result.stderr)

    def test_breeding_with_two_parents_needs_no_library(self):
        self._input.write_text(
            "C[1CH2]%10.[1NH]%10C amine\n"
            "C[1CH2]%10.[1O]%10C ether\n")
        result = self.run_tool("-j", "1", library=False, mode="b", limit=1)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1, result.stderr)

    def test_query_no_match_is_fatal(self):
        self._input.write_text(
            "[1CH3]%10.[1CH2]%10CC first\n"
            "C[1CH2]%10.[1CH2]%10N second\n")
        result = self.run_tool("-j", "1", "-Y", "SMARTS:[N]")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_query_ignore_continues_with_later_molecules(self):
        self._input.write_text(
            "[1CH3]%10.[1CH2]%10CC first\n"
            "C[1CH2]%10.[1CH2]%10N second\n")
        self._library.write_text('smi: "[1O]C" par: "ether"\n')
        result = self.run_tool("-j", "1", "-Y", "SMARTS:[N]", "-z", "i")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1)
        self.assertIn(" second ", result.stdout)
        self.assertNotIn(" first ", result.stdout)
        result = self.run_tool("-j", "1", "-Y", "SMARTS:[P]", "-z", "i")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")

    def test_identifier_with_period(self):
        self._input.write_text(
            "[13CH3][1CH2]%10.[1CH2]%10[2CH2][13CH3] parent.v1\n")
        result = self.run_tool("-j", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(" parent.v1 ", result.stdout)

    def test_exhaustive_uses_small_and_multiple_source_fragments(self):
        self._input.write_text("[1CH3]%10.[1CH2]%10CC parent\n")
        self._library.write_text('smi: "[1NH2]" par: "amine"\n')
        result = self.run_tool("-j", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1)
        self._input.write_text("C[1CH2]%10.[1CH2]%10N parent\n")
        self._library.write_text('smi: "[1O]C" par: "ether"\n')
        result = self.run_tool("-j", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 2, result.stdout)
        result = self.run_tool("-j", "1", "-N", "SMARTS:[N]")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1)

    def test_fragment_size_bounds_apply_to_exhaustive(self):
        result = self.run_tool("-j", "1", "-C", self.config(
            "min_atoms_in_fragment: 1\nmax_atoms_in_fragment: 2\nextra_atoms: 1\n"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        for text in ("min_atoms_in_fragment: 4\nmax_atoms_in_fragment: 2\n",
                     "min_atoms_in_fragment: 0\n"):
            result = self.run_tool("-C", self.config(text))
            self.assertNotEqual(result.returncode, 0)

    def test_atom_delta_defaults_and_cli(self):
        self._input.write_text("[1CH3]%10.[1CH2]%10CC parent\n")
        self._library.write_text('smi: "[1NH]CCC" par: "larger"\n')
        result = self.run_tool("-j", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        result = self.run_tool("-j", "1", "-x", "extra=1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1, result.stdout)
        self._library.write_text('smi: "[1NH]C" par: "smaller"\n')
        result = self.run_tool("-j", "1", "-C", self.config("fewer_atoms: 1\n"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1, result.stdout)
        result = self.run_tool("-j", "1", "-x", "extra=0")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")

    def test_ring_count_constraints(self):
        self._input.write_text("[1CH3]%10.[1CH2]%10CC parent\n")
        self._library.write_text('smi: "[1CH]1CC1" par: "cyclopropyl"\n')
        for allowed, expected in ((0, 0), (1, 1)):
            result = self.run_tool("-j", "1", "-C", self.config(
                f"extra_rings: {allowed}\n"))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(result.stdout.splitlines()), expected, result.stdout)
        self._input.write_text("[1CH3]%10.[1CH]%10C1CC1 parent\n")
        self._library.write_text('smi: "[1NH]CCC" par: "chain"\n')
        for allowed, expected in ((0, 0), (1, 1), (2, 1)):
            result = self.run_tool("-j", "1", "-C", self.config(
                f"fewer_rings: {allowed}\n"))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(result.stdout.splitlines()), expected, result.stdout)

    def test_ring_size_constraints(self):
        self._input.write_text("[1CH3]%10.[1CH]%101CCC1 parent\n")
        self._library.write_text('smi: "[1CH]1CC1" par: "smaller_ring"\n')
        for allowed, expected in ((0, 0), (1, 1)):
            result = self.run_tool("-j", "1", "-C", self.config(
                f"fewer_atoms: 1\nmax_decrease_ring_size: {allowed}\n"))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(result.stdout.splitlines()), expected, result.stdout)

    def test_futility_counts_consecutive_failures(self):
        self._input.write_text("C[1CH2]%10.[1NH]%10C parent\n")
        self._library.write_text('smi: "[1O]C" par: "ether"\n')
        for cutoff, attempts in ((2, 3), (0, 100)):
            result = self.run_tool("-j", "1", "-Y", "SMARTS:[N]",
                                   "-k", str(cutoff), "-v", mode="n", limit=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(result.stdout.splitlines()), 1)
            self.assertIn(f"generated {attempts} molecules", result.stderr)

    def test_deduplicate_after_attachment_isotope_removal(self):
        self._input.write_text(
            "[1CH3]%10.[1CH2]%10CC first\n"
            "C[1CH2]%10.[1CH2]%10C second\n")
        self._library.write_text(
            'smi: "[1CH3]" par: "methyl"\n'
            'smi: "[1CH2]C" par: "ethyl"\n')
        result = self.run_tool("-j", "1", "-I")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1, result.stdout)

    def test_unsupported_settings_and_negative_limits_fail(self):
        for args in (("-l",), ("-c",), ("-g", "all"), ("-T", "Br=Cl"),
                     ("-x", "extra=-1"), ("-x", "fewer=-1"), ("-k", "-1"),
                     ("-z", "unknown")):
            with self.subTest(args=args):
                result = self.run_tool(*args)
                self.assertNotEqual(result.returncode, 0)
        for text in ('atom_type: "UST:Z"\n',
                     'etrans { etrans { efrom: "Br" eto: "Cl" } }\n'):
            result = self.run_tool("-C", self.config(text))
            self.assertNotEqual(result.returncode, 0)

    def test_adjacent_heteroatoms_are_actually_rejected(self):
        self._input.write_text("[1CH3]%10.[1CH2]%10CNOC parent\n")
        self._library.write_text('smi: "[1NH2]" par: "amine"\n')
        result = self.run_tool("-j", "1", "-v")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        self.assertIn("1 rejected for adjacent heteroatoms", result.stderr)

    def test_parent_output_also_works_for_exhaustive(self):
        result = self.run_tool("-j", "1", "-p")
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].endswith(" parent parent"))
        self.assertIn(" %% replacement.", lines[1])

    def test_missing_generation_mode_or_library(self):
        for options in ([], ["-e", "0"], ["-n", "1"]):
            result = subprocess.run(
                [_EXECUTABLE, *options, str(self._input)],
                capture_output=True, text=True, check=False, timeout=5)
            self.assertNotEqual(result.returncode, 0)

    def test_internal_two_digit_library_ring_is_rejected(self):
        self._input.write_text(
            "[1CH3]%10.[1CH2]%10%11.[1CH2]%11CC parent\n")
        self._library.write_text('smi: "[1CH]%10CC%10" par: "ring"\n')
        result = self.run_tool("-j", "1")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("internal two-digit ring closures", result.stderr)
        self.assertEqual(result.stdout, "")

    def test_output_failure(self):
        if not pathlib.Path("/dev/full").exists():
            self.skipTest("/dev/full is unavailable")
        with open("/dev/full", "w") as output:
            result = subprocess.run(
                [_EXECUTABLE, "-j", "1", "-e", "1", "-L", str(self._library),
                 str(self._input)], stdout=output, stderr=subprocess.PIPE,
                text=True, check=False, timeout=5)
        self.assertNotEqual(result.returncode, 0)

    def test_empty_library(self):
        self._library.write_text("")
        result = self.run_tool("-j", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")

    def test_invalid_markers(self):
        for value in ("0", "-1", "abc", "4294967296"):
            with self.subTest(value=value):
                result = self.run_tool("-j", value)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Invalid attachment isotope", result.stderr)
        result = self.run_tool("-j", "1", "-j", "7")
        self.assertNotEqual(result.returncode, 0)

    def test_missing_library_marker(self):
        result = self.run_tool("-j", "7")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no attachment isotopes", result.stderr)

    def test_missing_configuration_library_marker(self):
        config = self._directory / "config.textproto"
        config.write_text("library: " + json.dumps(str(self._library)) + "\n")
        self._input.write_text("[7CH3]%10.[7CH2]%10CC parent\n")
        result = self.run_tool("-C", str(config), "-j", "7", library=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no attachment isotopes", result.stderr)

    def test_missing_input_marker(self):
        self._library.write_text('smi: "[7NH][2CH2][15CH3]" par: "replacement"\n')
        result = self.run_tool("-j", "7")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no attachment isotopes", result.stderr)


if __name__ == "__main__":
    unittest.main()
