"""Focused CLI regressions for fetch_smiles_quick."""

import os
import pathlib
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class FetchSmilesQuickTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._haystack = self._directory / "haystack.smi"
        # Keep the match last so these tests do not depend on the separate
        # early-termination issue affecting unmatched output after a match.
        self._haystack.write_text("N other\nC wanted\n")
        (self._directory / "ids.txt").write_text("wanted\n")

    def run_tool(self, *options, inputs=("haystack.smi",), success=True, env=None):
        result = subprocess.run([_EXECUTABLE, *options, *inputs], cwd=self._directory,
                                capture_output=True, text=True, timeout=20, env=env)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertGreater(result.returncode, 0, result.stderr)
        return result

    def test_explicit_identifier_with_unselected_output(self):
        result = self.run_tool("-K", "wanted", "-X", "unselected.smi")
        self.assertEqual(result.stdout, "C wanted\n")
        self.assertEqual((self._directory / "unselected.smi").read_text(), "N other\n")

    def test_explicit_identifier_with_missing_output(self):
        result = self.run_tool("-K", "wanted", "-Y", "missing.txt")
        self.assertEqual(result.stdout, "C wanted\n")
        self.assertEqual((self._directory / "missing.txt").read_text(), "")

    def test_explicit_identifiers_with_both_outputs(self):
        result = self.run_tool("-K", "wanted", "-K", "missing", "-X", "unselected.smi",
                               "-Y", "missing.txt")
        self.assertEqual(result.stdout, "C wanted\n")
        self.assertEqual((self._directory / "unselected.smi").read_text(), "N other\n")
        self.assertEqual((self._directory / "missing.txt").read_text(), "missing\n")

    def test_external_identifier_file_with_one_positional_input(self):
        for flag in ("-X", "-Y"):
            with self.subTest(flag=flag):
                result = self.run_tool("-F", "ids.txt", flag, "extra.txt")
                self.assertEqual(result.stdout, "C wanted\n")
                expected = "N other\n" if flag == "-X" else ""
                self.assertEqual((self._directory / "extra.txt").read_text(), expected)

    def test_explicit_identifier_output_cannot_replace_input(self):
        original = self._haystack.read_text()
        for flag in ("-X", "-Y"):
            with self.subTest(flag=flag):
                result = self.run_tool("-K", "wanted", flag, "haystack.smi", success=False)
                self.assertIn("Cannot overwrite input", result.stderr)
                self.assertEqual(self._haystack.read_text(), original)

    def test_output_cannot_replace_later_positional_input(self):
        later = self._directory / "later.smi"
        later.write_text("O later\n")
        for flag in ("-X", "-Y"):
            with self.subTest(flag=flag):
                self.run_tool("-K", "wanted", flag, "later.smi",
                              inputs=("haystack.smi", "later.smi"), success=False)
                self.assertEqual(later.read_text(), "O later\n")

    def test_output_aliases_cannot_replace_inputs(self):
        original = self._haystack.read_text()
        alias = self._directory / "alias.smi"
        for hardlink in (False, True):
            if hardlink:
                alias.hardlink_to(self._haystack)
            else:
                alias.symlink_to(self._haystack)
            for flag in ("-X", "-Y"):
                with self.subTest(hardlink=hardlink, flag=flag):
                    self.run_tool("-K", "wanted", flag, "alias.smi", success=False)
                    self.assertEqual(self._haystack.read_text(), original)
            alias.unlink()

    def test_external_identifier_input_is_protected(self):
        for flag in ("-X", "-Y"):
            with self.subTest(flag=flag):
                self.run_tool("-F", "ids.txt", flag, "ids.txt", success=False)
                self.assertEqual((self._directory / "ids.txt").read_text(), "wanted\n")

    def test_side_outputs_must_differ_before_either_is_opened(self):
        self.run_tool("-K", "wanted", "-X", "extra.txt", "-Y", "./extra.txt", success=False)
        self.assertFalse((self._directory / "extra.txt").exists())

    def test_side_outputs_cannot_be_hard_links_to_one_another(self):
        first = self._directory / "first.txt"
        first.write_text("original\n")
        (self._directory / "second.txt").hardlink_to(first)
        self.run_tool("-K", "wanted", "-X", "first.txt", "-Y", "second.txt", success=False)
        self.assertEqual(first.read_text(), "original\n")

    def test_w_suffix_applies_only_to_missing_identifier_output(self):
        self.run_tool("-K", "wanted", "-K", "absent", "-w", "-X", "extra", "-Y", "missing")
        self.assertEqual((self._directory / "extra").read_text(), "N other\n")
        self.assertFalse((self._directory / "extra.smi").exists())
        self.assertEqual((self._directory / "missing.smi").read_text(), "absent\n")

    def test_output_suffix_is_included_in_conflict_checks(self):
        self.run_tool("-K", "wanted", "-w", "-X", "extra.smi", "-Y", "extra", success=False)
        self.assertFalse((self._directory / "extra.smi").exists())

    def test_output_comparison_error_is_not_ignored(self):
        (self._directory / "loop").symlink_to("loop")
        self.run_tool("-K", "wanted", "-X", "loop/output", success=False)
        self.assertEqual(self._haystack.read_text(), "N other\nC wanted\n")

    def test_empty_output_names_are_rejected(self):
        for flag in ("-X", "-Y"):
            with self.subTest(flag=flag):
                self.run_tool("-K", "wanted", flag, "", success=False)
        self.run_tool("-S", "", inputs=("ids.txt", "haystack.smi"), success=False)

    def test_subset_output_cannot_replace_identifier_input(self):
        (self._directory / "subset0").write_text("wanted\n")
        self.run_tool("-S", "subset", inputs=("subset0", "haystack.smi"), success=False)
        self.assertEqual((self._directory / "subset0").read_text(), "wanted\n")

    def test_subset_append_prefix_cannot_bypass_input_check(self):
        path = self._directory / "subset0"
        path.write_text("wanted\n")
        self.run_tool("-S", ">>subset", inputs=("subset0", "haystack.smi"), success=False)
        self.assertEqual(path.read_text(), "wanted\n")

    def test_subset_expanded_output_cannot_replace_input(self):
        path = self._directory / "subset0"
        path.write_text("wanted\n")
        env = dict(os.environ, FETCH_QUICK_TEST_DIR=str(self._directory))
        self.run_tool("-S", "${FETCH_QUICK_TEST_DIR}/subset",
                      inputs=("subset0", "haystack.smi"), success=False, env=env)
        self.assertEqual(path.read_text(), "wanted\n")

    def test_gzip_output_uses_literal_append_prefix(self):
        (self._directory / ">>subset0.gz").symlink_to(self._haystack)
        original = self._haystack.read_text()
        self.run_tool("-S", ">>subset", "-u", ".gz",
                      inputs=("ids.txt", "haystack.smi"), success=False)
        self.assertEqual(self._haystack.read_text(), original)

    def test_expanded_input_cannot_be_overwritten(self):
        env = dict(os.environ, FETCH_QUICK_TEST_DIR=str(self._directory))
        original = self._haystack.read_text()
        self.run_tool("-K", "wanted", "-X", "haystack.smi",
                      inputs=("${FETCH_QUICK_TEST_DIR}/haystack.smi",),
                      success=False, env=env)
        self.assertEqual(self._haystack.read_text(), original)

    def test_all_subset_outputs_are_checked_before_opening(self):
        (self._directory / "subset1").write_text("wanted\n")
        self.run_tool("-S", "subset", inputs=("ids.txt", "subset1", "haystack.smi"), success=False)
        self.assertFalse((self._directory / "subset0").exists())
        self.assertEqual((self._directory / "subset1").read_text(), "wanted\n")

    def test_subset_output_alias_cannot_replace_haystack(self):
        (self._directory / "subset0.smi").symlink_to(self._haystack)
        original = self._haystack.read_text()
        self.run_tool("-S", "subset", "-u", ".smi", inputs=("ids.txt", "haystack.smi"), success=False)
        self.assertEqual(self._haystack.read_text(), original)

    def test_subset_output_roles_must_use_distinct_files(self):
        self.run_tool("-S", "subset", "-X", "subset", inputs=("ids.txt", "haystack.smi"), success=False)
        self.assertFalse((self._directory / "subset0").exists())

    def test_subset_output_collision_across_different_indices(self):
        (self._directory / "second.txt").write_text("wanted\n")
        self.run_tool("-S", "out", "-u", "1", "-X", "out0",
                      inputs=("ids.txt", "second.txt", "haystack.smi"), success=False)
        self.assertFalse((self._directory / "out01").exists())
        self.assertFalse((self._directory / "out00").exists())

    def test_subset_sequence_index_is_validated(self):
        for index in ("abc", "-1"):
            with self.subTest(index=index):
                self.run_tool("-S", "subset", "-g", index,
                              inputs=("ids.txt", "haystack.smi"), success=False)
        (self._directory / "second.txt").write_text("wanted\n")
        self.run_tool("-S", "subset", "-g", "2147483647",
                      inputs=("ids.txt", "second.txt", "haystack.smi"), success=False)

    def test_haystack_read_failure_returns_nonzero(self):
        (self._directory / "directory").mkdir()
        self.run_tool("-K", "wanted", inputs=("directory",), success=False)

    def test_subset_identifier_read_failure_returns_nonzero(self):
        (self._directory / "directory").mkdir()
        self.run_tool("-S", "subset", inputs=("directory", "haystack.smi"), success=False)

    def test_normal_output_write_failure_returns_nonzero(self):
        if not pathlib.Path("/dev/full").exists():
            self.skipTest("Requires /dev/full")
        for options in (("-X", "/dev/full"), ("-Y", "/dev/full")):
            with self.subTest(options=options):
                self.run_tool("-K", "wanted", "-K", "absent", *options, success=False)

    def test_stdout_write_failure_returns_nonzero(self):
        if not pathlib.Path("/dev/full").exists():
            self.skipTest("Requires /dev/full")
        for block_size, length in (("1", 10000), ("8192", 1)):
            with self.subTest(block_size=block_size):
                self._haystack.write_text("C" * length + " wanted\n")
                with open("/dev/full", "w") as sink:
                    result = subprocess.run([_EXECUTABLE, "-K", "wanted", "-B", block_size,
                                             "haystack.smi"], cwd=self._directory, stdout=sink,
                                            stderr=subprocess.PIPE, text=True, timeout=20)
                self.assertGreater(result.returncode, 0, result.stderr)

    def test_subset_write_failure_returns_nonzero(self):
        if not pathlib.Path("/dev/full").exists():
            self.skipTest("Requires /dev/full")
        (self._directory / "ids.txt").write_text("wanted\nabsent\n")
        for target, options in (("subset0", ()), ("extra0", ("-X", "extra")),
                                ("missing0", ("-Y", "missing"))):
            with self.subTest(target=target):
                path = self._directory / target
                path.symlink_to("/dev/full")
                self.run_tool("-S", "subset", *options,
                              inputs=("ids.txt", "haystack.smi"), success=False)
                path.unlink()

    def test_removed_flush_option_is_rejected(self):
        result = self.run_tool("-h", "-K", "wanted", success=False)
        self.assertIn("Unrecognised options", result.stderr)


if __name__ == "__main__":
    unittest.main()
