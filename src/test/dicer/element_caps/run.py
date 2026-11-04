#!/usr/bin/env python3
"""Element-cap regression checks; no Python chemistry bindings required."""
import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(pathlib.Path(sys.argv[1]).resolve())
sys.argv = sys.argv[:1]


class ElementCaps(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input = pathlib.Path(self.tmp.name) / 'input.smi'

    def run_dicer(self, smiles, extra=(), success=True, depth=3):
        self.input.write_text(smiles + ' example\n')
        args = [EXE, '-B', 'bscb', '-B', 'proto', '-k', str(depth),
                '-m', '1', '-M', '30', *extra, str(self.input)]
        result = subprocess.run(args, text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        self.assertNotRegex(result.stderr, r'AddressSanitizer|runtime error:')
        return [(smi, int(nat)) for smi, nat in re.findall(
            r'fragment \{ smi: "([^"]+)" nat: (\d+)', result.stdout)]

    def provenance(self, smi):
        return sorted(int(x) % 100 for x in re.findall(r'\[(\d+)', smi))

    def test_endpoint_elements_and_original_numbering(self):
        for smiles in ['[1NH2][2CH2][3C]([4F])([5F])[6F]',
                       '[1NH2][2CH2][3CH]=[4CH2]',
                       '[1NH2][2CH2][3C]#[4CH]',
                       '[1NH2][2CH2][3c]1[4cH][5cH][6cH][7cH][8cH]1']:
            with self.subTest(smiles=smiles):
                fragments = self.run_dicer(smiles, ['-I', 'cap', '-I', 'inc=100'])
                ncc = [s for s, n in fragments if n == 3 and self.provenance(s) == [1, 2, 3]]
                self.assertTrue(ncc, fragments)
                for smi in ncc:
                    self.assertEqual(len(re.findall(r'\[\d+CH2\]', smi)), 1)
                    self.assertNotRegex(smi, r'\[\d+c')
                # Caps must never create duplicate parent atoms within a fragment.
                for smi, nat in fragments:
                    labels = self.provenance(smi)
                    self.assertEqual(len(labels), nat)
                    self.assertEqual(len(labels), len(set(labels)), smi)

    def test_recursive_fragments_keep_previous_caps(self):
        smiles = '[1NH2][2CH2][3CH2][4O][5CH3]'
        one = self.run_dicer(smiles, ['-I', 'cap'], depth=1)
        recursive = self.run_dicer(smiles, ['-I', 'cap'], depth=3)
        # Atom 3 gets a carbon cap at one cut and an oxygen cap at another.
        target = [2, 3, 4]
        self.assertFalse(any(self.provenance(s) == target for s, _ in one))
        self.assertTrue(any(self.provenance(s) == target for s, _ in recursive))
        # Two-atom cores must not collapse to one atom by cutting a cap bond.
        self.assertTrue(all(n >= 2 for _, n in recursive))

    def test_default_is_unchanged_and_cap_needs_no_other_label(self):
        smiles = '[1NH2][2CH2][3C]([4F])([5F])[6F]'
        plain = self.run_dicer(smiles)
        capped = self.run_dicer(smiles, ['-I', 'cap'])
        target = [1, 2, 3]
        self.assertFalse(any(self.provenance(s) == target for s, _ in plain))
        self.assertTrue(any(self.provenance(s) == target for s, _ in capped))

    def test_charge_and_conflicting_modes(self):
        fragments = self.run_dicer('[1NH3+][2CH2][3CH3]', ['-I', 'cap'])
        self.assertTrue(any('+' in s and n == 2 and self.provenance(s) == [1, 2]
                            for s, n in fragments))
        for extra in [['-I', 'enva'], ['-B', 'chiraldummy'], ['-B', 'BRB'],
                      ['-B', 'BRF'], ['-B', 'BRS'], ['-B', 'dearom'],
                      ['-B', 'recap'], ['-C', 'def'], ['-I', 'fragatype']]:
            with self.subTest(extra=extra):
                self.run_dicer('NCC', ['-I', 'cap', *extra], success=False)

    def test_both_modes_preserve_each_candidate_family(self):
        for smiles in ['[1NH2][2CH2][3C]([4F])([5F])[6F]',
                       '[1NH2][2CH2][3CH2][4O][5CH3]',
                       '[1NH2][2CH]([3CH3])[4CH2][5CH3]']:
            for labels in [[], ['-I', 'inc=100']]:
                with self.subTest(smiles=smiles, labels=labels):
                    plain = set(self.run_dicer(smiles, labels))
                    capped = set(self.run_dicer(smiles, ['-I', 'cap', *labels]))
                    both = set(self.run_dicer(smiles, ['-I', 'capboth', *labels]))
                    self.assertTrue(plain <= both, plain - both)
                    self.assertTrue(capped <= both, capped - both)

    def test_both_modes_filter_each_size_independently(self):
        # The capped amine-side fragment has three atoms; its bare form has two.
        fragments = self.run_dicer('[1NH2][2CH2][3C]([4F])([5F])[6F]',
                                   ['-I', 'capboth', '-M', '2'])
        self.assertTrue(any(self.provenance(s) == [1, 2] for s, _ in fragments))
        self.assertTrue(all(n <= 2 for _, n in fragments))
        self.run_dicer('NCC', ['-I', 'cap', '-I', 'capboth'], success=False)

    def test_marker_identifies_caps_and_survives_recursion(self):
        smiles = '[1NH2][2CH2][3CH2][4O][5CH3]'
        capped = self.run_dicer(smiles, ['-I', 'cap', '-I', 'inc=100'])
        self.assertTrue(all(any(int(i) >= 10000 for i in re.findall(r'\[(\d+)', s))
                            for s, _ in capped))
        # Both previous carbon and new oxygen caps survive in the central fragment.
        middle = next(s for s, _ in capped if self.provenance(s) == [2, 3, 4])
        self.assertEqual(sorted(int(i) % 100 for i in re.findall(r'\[(\d+)', middle)
                                if int(i) >= 10000), [2, 4])
        both = self.run_dicer(smiles, ['-I', 'capboth', '-I', 'inc=100'])
        self.assertTrue(any(all(int(i) < 10000 for i in re.findall(r'\[(\d+)', s))
                            for s, _ in both))

    def test_marker_range_is_reserved_only_in_cap_modes(self):
        self.run_dicer('[10000NH2]CC', ['-I', 'cap'], success=False)
        self.run_dicer('NCC', ['-I', 'cap', '-I', '10000'], success=False)
        self.run_dicer('NCC', ['-I', 'cap', '-I', 'inc=10000'], success=False)
        self.run_dicer('[9999NH2]CC', ['-I', 'cap', '-I', 'INC=100'], success=False)
        self.run_dicer('NCC', ['-I', 'cap', '-I', 'atype'], success=False)
        self.run_dicer('NCC', ['-I', 'cap', '-I', 'atype', '-P', 'UST:Z'])
        self.run_dicer('[10000NH2]CC')


if __name__ == '__main__':
    result = unittest.main(exit=False)
    if not result.result.wasSuccessful():
        sys.exit(1)
    print('element cap checks passed')
