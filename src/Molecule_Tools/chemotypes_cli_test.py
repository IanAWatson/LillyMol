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

    def rejection_input(self):
        self._input.write_text('C1CCCCC1 no_match\nNCC no_ring\n'
                               'N1CCCCC1CCN2CCCCC2 ambiguous\nN1CCCCC1 accepted\n')

    def rejected(self, name='rejected.smi'):
        return (self._directory / name).read_text().splitlines()

    def test_rejections_continue_and_append_smi(self):
        self.rejection_input()
        for filename, actual in [('rejected', 'rejected.smi'),
                                 ('rejected.smi', 'rejected.smi'),
                                 ('rejected.txt', 'rejected.txt.smi')]:
            result = self.run_tool('-B', filename, '-v')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(result.stdout.splitlines()), 1)
            self.assertIn(' accepted', result.stdout)
            lines = self.rejected(actual)
            self.assertEqual(len(lines), 3)
            for line, identifier, reason in zip(lines,
                    ['no_match', 'no_ring', 'ambiguous'],
                    ['no_query_match', 'no_ring_atom', 'ambiguous_ring_systems']):
                self.assertIn(f' {identifier} chemotype_reject={reason}', line)
            self.assertIn('Read 4 molecules', result.stderr)
            self.assertIn('Wrote 3 rejected molecules', result.stderr)
        self.assertNotEqual(self.run_tool().returncode, 0)
        self.assertNotEqual(self.run_tool('-z', 'i').returncode, 0)
        result = self.run_tool('-B', 'rejected', '-z', 'i', '-z', 'first')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 2)
        self.assertEqual(len(self.rejected()), 2)

    def test_rejections_min_rings_and_output_format(self):
        self._input.write_text('N1CCCCC1 one_ring\nN1CCCCC1CC2CCCCC2 two_rings\n')
        result = self.run_tool('-B', 'rejected', '-r', '2', '-o', 'sdf', '-S', 'accepted')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.rejected()), 1)
        self.assertIn(' one_ring chemotype_reject=below_min_rings', self.rejected()[0])
        self.assertTrue((self._directory / 'accepted.sdf').is_file())
        self.assertNotIn('$$$$', self.rejected()[0])

    def test_rejections_preserve_preprocessed_parent_and_input_isotopes(self):
        self._input.write_text('[13CH3]C1CCCCC1.[Na+] original name\n')
        result = self.run_tool('-B', 'rejected', '-l', '-P', 'UST:Z')
        self.assertEqual(result.returncode, 0, result.stderr)
        line = self.rejected()[0]
        self.assertIn(' original name chemotype_reject=no_query_match', line)
        self.assertEqual(self.labels(line), [13])
        self.assertNotIn('Na', line)
        self._input.write_text('N1CCCCC1 accepted\nC1CCCCC1 rejected\n')
        result = self.run_tool('-B', 'rejected', '-p', 'def')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 2)  # Parent plus chemotype only.

    def test_rejections_keep_fingerprint_records_aligned(self):
        self.rejection_input()
        for tag in ['FPTEST', 'NCTEST']:
            result = self.run_tool('-B', 'rejected', '-J', tag)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.splitlines().count('|'), 4)
            self.assertEqual(result.stdout.count(f'{tag}<>\n'), 3)
            self.assertEqual(len(self.rejected()), 3)
            self.assertNotIn('chemotype_reject=', result.stdout)

    def test_rejections_tdt_identifiers_metadata_and_empty_fingerprints(self):
        records = [('C1CCCCC1', 'no_match idA'), ('NCC', 'no_ring idB'),
                   ('N1CCCCC1CCN2CCCCC2', 'ambiguous idC'), ('N1CCCCC1', 'accepted idD')]
        self._input.write_text(''.join(f'$SMI<{smi}>\nPCN<{name}>\nEXTRA<keep>\n|\n'
                                       for smi, name in records))
        result = self.run_tool('-B', 'rejected', '-J', 'FPTEST', '-f')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.count('FPTEST<>\n'), 3)
        self.assertEqual(result.stdout.count('EXTRA<keep>'), 4)
        self.assertEqual(result.stdout.splitlines().count('|'), 4)
        for _, name in records:
            self.assertIn(f'PCN<{name}>', result.stdout)
        for line, (_, name) in zip(self.rejected(), records):
            self.assertIn(f' {name} chemotype_reject=', line)
        self._input.write_text('$SMI<N1CCCCC1>\nPCN<incomplete>\n')
        self.assertNotEqual(self.run_tool('-B', 'rejected', '-J', 'FPTEST', '-f').returncode, 0)
        self._input.write_text('$SMI<N1CCCCC1X\nPCN<unterminated>\n|\n')
        self.assertNotEqual(self.run_tool('-B', 'rejected', '-J', 'FPTEST', '-f').returncode, 0)

    def test_rejections_supply_identifier_for_unnamed_input(self):
        self._input.write_text('C1CCCCC1\n')
        result = self.run_tool('-B', 'rejected')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(' molecule_1 chemotype_reject=no_query_match', self.rejected()[0])

    def test_rejections_do_not_overwrite_input_or_other_outputs(self):
        original = self._input.read_text()
        for filename in [str(self._input), str(self._input.with_suffix(''))]:
            self.assertNotEqual(self.run_tool('-B', filename).returncode, 0)
            self.assertEqual(self._input.read_text(), original)
        alias = self._directory / 'alias.smi'
        alias.symlink_to(self._input)
        self.assertNotEqual(self.run_tool('-B', str(alias)).returncode, 0)
        self.assertEqual(self._input.read_text(), original)
        sentinel = self._directory / 'shared.smi'
        sentinel.write_text('preserve me\n')
        for args in [('-S', 'shared'), ('-C', 'shared'), ('-F', 'shared.smi'),
                     ('-H', 'xref=shared.smi')]:
            result = self.run_tool('-B', 'shared', *args)
            self.assertNotEqual(result.returncode, 0, args)
            self.assertEqual(sentinel.read_text(), 'preserve me\n')
        self.assertNotEqual(self.run_tool('-B', '-').returncode, 0)
        self.assertNotEqual(self.run_tool('-B', 'a', '-B', 'b').returncode, 0)

    def test_rejections_invalid_input_options_and_output_errors_stay_fatal(self):
        result = self.run_tool('-B', 'invalid_options', '-I', 'ring=1', '-P', 'Z')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self._directory / 'invalid_options.smi').exists())
        self._input.write_text('not_a_smiles broken\n')
        for extra in [[], ['-J', 'FPTEST']]:
            self.assertNotEqual(self.run_tool('-B', 'rejected', *extra).returncode, 0)
        self.rejection_input()
        self.assertNotEqual(self.run_tool('-B', 'absent/destination').returncode, 0)
        if pathlib.Path('/dev/full').exists():
            (self._directory / 'full.smi').symlink_to('/dev/full')
            self.assertNotEqual(self.run_tool('-B', 'full').returncode, 0)


if __name__ == "__main__":
    unittest.main()
