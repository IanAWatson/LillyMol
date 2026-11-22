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


class ProtoFileQueries(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        query_dir = self._directory / "query files"
        query_dir.mkdir()
        self._list = query_dir / "queries"
        (query_dir / "nitrogen#query.textproto").write_text(
            'name: "nitrogen"\nquery { smarts: "N" }\n')
        (query_dir / "oxygen query.textproto").write_text(
            'name: "oxygen"\nquery { smarts: "O" }\n')
        (self._directory / "input.smi").write_text("CN amine\nCO alcohol\nCC ethane\n")

    def run_list(self, text):
        self._list.write_text(text)
        return subprocess.run(
            [_EXECUTABLE, "-q", "PROTOFILE:" + str(self._list), "-a", "input.smi"],
            cwd=self._directory, capture_output=True, text=True, timeout=15)

    def test_comments_and_blank_lines_preserve_queries_and_order(self):
        clean = self.run_list("nitrogen#query.textproto\noxygen query.textproto\n")
        decorated = self.run_list(
            "# Heading\n\n  # Indented comment ${UNDEFINED_QUERY_VARIABLE}\n \t\n"
            "  nitrogen#query.textproto  \n\t# More comments\n"
            "oxygen query.textproto\n\n# Trailing comment\n")
        self.assertEqual(clean.returncode, 0, clean.stderr)
        self.assertEqual(decorated.returncode, 0, decorated.stderr)
        self.assertEqual(decorated.stdout, clean.stdout)
        self.assertIn("Name nitrogen oxygen", decorated.stdout)
        self.assertIn("amine 1 0", decorated.stdout)
        self.assertIn("alcohol 0 1", decorated.stdout)

    def test_noncomment_missing_file_is_still_an_error(self):
        result = self.run_list("# Skip this\n\nmissing.textproto\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no file 'missing.textproto'", result.stderr)

    def test_comment_only_list_has_no_queries(self):
        result = self.run_list("# No queries\n\n \t\n  # Still no queries\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("no file", result.stderr)


class Correctness(unittest.TestCase):
    setUp = OutputModifiers.setUp
    run_tool = OutputModifiers.run_tool
    def test_malformed_input_returns_failure(self):
        (self._directory / "input.smi").write_text("CN valid\ninvalid_smiles broken\n")
        self.assertNotEqual(self.run_tool("-m", "-").returncode, 0)

    def test_timing_works_for_short_runs(self):
        result = self.run_tool("-M", "time", "-m", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("molecules per second", result.stderr)

    def test_break_clears_unsearched_query_results(self):
        (self._directory / "input.smi").write_text("CCCC first\nN second\n")
        result = self.run_tool("-s", "CC carbon", "-b", "-G", "atoms.txt")
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = (self._directory / "atoms.txt").read_text().splitlines()
        self.assertEqual(lines[1], "second (0),")

    def test_break_at_first_nonmatch_clears_old_results(self):
        (self._directory / "input.smi").write_text("NCC first\nC second\n")
        result = self.run_tool("-s", "CC carbon", "-B", "-G", "atoms.txt")
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = (self._directory / "atoms.txt").read_text().splitlines()
        self.assertEqual(lines[1], "second ,")

    def test_self_match_skip_clears_old_results(self):
        (self._directory / "input.smi").write_text("CCCC first\nN carbon\n")
        result = self.run_tool("-s", "CC carbon", "-M", "nosm", "-G", "atoms.txt")
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = (self._directory / "atoms.txt").read_text().splitlines()
        self.assertEqual(lines[1], "carbon (0),")

    def test_negative_isotope_offset(self):
        (self._directory / "input.smi").write_text("[5NH3] example\n")
        result = self.run_tool("-j", "-2", "-m", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("[3NH3]", result.stdout)

    def test_sequential_labels_start_at_one(self):
        (self._directory / "input.smi").write_text("N example\n")
        result = self.run_tool("-j", "sequential", "-m", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("[1NH3]", result.stdout)

    def test_qmatch_counts_queries_not_overlapping_embeddings(self):
        (self._directory / "input.smi").write_text("CCC example\n")
        result = self.run_tool("-s", "CC", "-u", "-j", "qmatch", "-m", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("[1CH3][1CH2][1CH3]", result.stdout)

    def test_qmatch_counts_distinct_queries(self):
        (self._directory / "input.smi").write_text("CCC example\n")
        result = self.run_tool("-s", "CC first", "-s", "CC second", "-u", "-j", "qmatch", "-m", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("[2CH3][2CH2][2CH3]", result.stdout)

    def test_qmatch_respects_atom_limit(self):
        (self._directory / "input.smi").write_text("CO example\n")
        for limit, expected in ((1, "[1CH3]O"), (2, "[1CH3][1OH]")):
            with self.subTest(limit=limit):
                result = self.run_tool("-s", "CO", "-j", "qmatch", "-j", f"n={limit}", "-m", "-")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(expected, result.stdout)

    def test_qmatch_atom_limit_with_numbering_gaps(self):
        (self._directory / "input.smi").write_text("CO example\n")
        query = self._directory / "limited.textproto"
        query.write_text('query { respect_initial_atom_numbering: true '
                         'query_atom { id: 0 atom_properties { atomic_number: 6 } } '
                         'query_atom { id: 2 atom_properties { atomic_number: 8 } single_bond: 0 } }')
        for limit, expected in ((2, "[1CH3]O"), (3, "[1CH3][1OH]")):
            with self.subTest(limit=limit):
                result = self.run_tool("-q", "PROTO:" + str(query), "-j", "qmatch", "-j", f"n={limit}", "-m", "-")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(expected, result.stdout)

    def test_element_labels_skip_gaps_in_query_atom_numbers(self):
        query = self._directory / "gaps.textproto"
        query.write_text('query {\n respect_initial_atom_numbering: true\n query_atom { id: 0 atom_properties { atomic_number: 6 } }\n query_atom { id: 2 atom_properties { atomic_number: 7 } single_bond: 0 }\n}\n')
        result = self.run_tool("-q", "PROTO:" + str(query), "-j", "O", "-m", "-")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("AddressSanitizer", result.stderr)
        self.assertIn("methylamine", result.stdout)

    def test_invalid_label_configuration_is_an_error(self):
        for label in ("NOT_AN_ELEMENT", "mqoffset=2"):
            result = self.run_tool("-j", label, "-m", "-")
            self.assertNotEqual(result.returncode, 0, result.stderr)

    def test_tdt_counts_each_molecule_once(self):
        (self._directory / "input.smi").write_text("$SMI<N>\nPCN<one>\n|\n")
        result = self.run_tool("-M", "TDT", "-J", "FPTEST")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("1 molecules read, 1 molecules match fraction 1", result.stderr)
        self.assertEqual(result.stdout.count("FPTEST<"), 1)

    def test_output_filenames_must_not_replace_input(self):
        original = (self._directory / "input.smi").read_text()
        for options in (("-R", "input.smi"), ("-G", "input.smi"),
                        ("-M", "DCF=input.smi")):
            result = self.run_tool(*options)
            self.assertNotEqual(result.returncode, 0, result.stderr)
            self.assertEqual((self._directory / "input.smi").read_text(), original)
        # SEPA's query zero output is input0.smi.
        (self._directory / "input0.smi").write_text(original)
        result = subprocess.run([_EXECUTABLE, "-s", "N", "-m", "input", "-m", "SEPA", "input0.smi"],
                                cwd=self._directory, capture_output=True, text=True, timeout=15)
        self.assertNotEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self._directory / "input0.smi").read_text(), original)

    def test_match_and_nonmatch_destinations_must_differ(self):
        result = self.run_tool("-m", "hits", "-n", "./hits.smi")
        self.assertNotEqual(result.returncode, 0, result.stderr)

    @unittest.skipUnless(pathlib.Path("/dev/full").exists(), "requires /dev/full")
    def test_buffered_output_failures_are_errors(self):
        (self._directory / "full.smi").symlink_to("/dev/full")
        for options in (("-m", "full"), ("-n", "full"), ("-R", "/dev/full"),
                        ("-G", "/dev/full"), ("-j", "1", "-j", "writeach=full"),
                        ("-m", "separate", "-m", "SEPA")):
            if "SEPA" in options:
                (self._directory / "separate0.smi").symlink_to("/dev/full")
            result = self.run_tool(*options)
            self.assertNotEqual(result.returncode, 0, (options, result.stderr))


if __name__ == "__main__":
    unittest.main()
