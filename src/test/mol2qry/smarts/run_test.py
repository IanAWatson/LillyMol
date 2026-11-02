#!/usr/bin/env python3
"""Validate generated atom SMARTS by searching molecules with LillyMol."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = Path(sys.argv[1]).resolve()
SEARCH = EXE.parent / 'tsubstructure'


class Smarts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def generate(self, smiles, *args):
        source = self.directory / 'source.smi'
        source.write_text(smiles + ' source\n')
        output = self.directory / 'query.sma'
        result = subprocess.run([str(EXE), *map(str, args), '-B', str(output), str(source)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        text = output.read_text().strip()
        smarts, name = text.split(maxsplit=1)
        self.assertEqual(name, 'source')
        return smarts

    def matches(self, smarts, molecules):
        source = self.directory / 'targets.smi'
        source.write_text(''.join(f'{smi} {name}\n' for name, smi in molecules.items()))
        result = subprocess.run([str(SEARCH), '-s', smarts, '-m', '-', str(source)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return {line.split()[1] for line in result.stdout.splitlines() if line.strip()}

    def test_exact_and_minimum_degree(self):
        molecules = {'same': 'CCO', 'extended': 'CCCO', 'branched': 'CC(C)O', 'ether': 'CCOC'}
        smarts = self.generate('CCO')
        self.assertEqual(smarts.count(';D1'), 2)
        self.assertIn(';D2', smarts)
        self.assertEqual(self.matches(smarts, molecules), {'same'})
        minimum = self.generate('CCO', '-m')
        self.assertEqual(self.matches(minimum, molecules), set(molecules))
        self.assertEqual(self.matches(self.generate('C'), {'same': 'C', 'larger': 'CC'}), {'same'})
        self.assertEqual(self.matches(self.generate('C', '-m'), {'same': 'C', 'larger': 'CC'}), {'same', 'larger'})

    def test_isotopic_substitution(self):
        smarts = self.generate('[1CH3]CO', '-s')
        self.assertIn('D>1', smarts)
        self.assertNotIn('[1', smarts)
        molecules = {'original': 'CCO', 'substituted': 'CCCO', 'other_site': 'CC(C)O'}
        self.assertEqual(self.matches(smarts, molecules), {'substituted'})
        self.assertEqual(self.matches(self.generate('[1CH3]CO', '-t'), molecules), {'substituted'})
        # Even a labelled atom with no implicit H must obey the degree condition.
        smarts = self.generate('C[1C](C)(C)C', '-s')
        self.assertIn('D>4', smarts)
        self.assertEqual(self.matches(smarts, {'same': 'CC(C)(C)C'}), set())

    def test_extra_connections_and_wildcard_are_independent(self):
        smarts = self.generate('[2CH3]CO', '-c')
        self.assertIn(';D3', smarts)
        self.assertNotIn('*', smarts)
        self.assertEqual(self.matches(smarts, {'two': 'CC(C)CO', 'one': 'CCCO'}), {'two'})
        smarts = self.generate('[9CH3]O', '-x', 9)
        self.assertIn('[*;D1]', smarts)
        self.assertEqual(self.matches(smarts, {'carbon': 'CO', 'nitrogen': 'NO'}), {'carbon', 'nitrogen'})

    def test_aromaticity_and_ring_constraints(self):
        targets = {'chain': 'CC', 'ring': 'C1CCCCC1'}
        for args in (('-n',), ('-j',)):
            smarts = self.generate('CC', '-m', *args)
            self.assertEqual(self.matches(smarts, targets), {'chain'})
        self.assertEqual(self.matches(self.generate('CC', '-m'), targets), set(targets))
        targets = {'aliphatic': 'C', 'aromatic': 'c1ccccc1'}
        self.assertEqual(self.matches(self.generate('C', '-m'), targets), set(targets))
        self.assertEqual(self.matches(self.generate('C', '-m', '-a'), targets), {'aliphatic'})
        smarts = self.generate('c1ccccc1', '-j')
        self.assertEqual(smarts.count(';x2'), 6)
        self.assertIn(':', smarts)
        self.assertEqual(self.matches(smarts, {'aromatic': 'c1ccccc1', 'aliphatic': 'C1CCCCC1'}), {'aromatic'})

    def test_saturation(self):
        targets = {'saturated': 'CCC', 'unsaturated': 'CC=C'}
        self.assertEqual(self.matches(self.generate('CC', '-m'), targets), set(targets))
        self.assertEqual(self.matches(self.generate('CC', '-m', '-d'), targets), {'saturated'})

    def test_isotopes_charges_and_preprocessing(self):
        targets = {'labelled': '[13CH4]', 'plain': 'C'}
        self.assertEqual(self.matches(self.generate('[13CH4]'), targets), {'labelled'})
        self.assertEqual(self.matches(self.generate('[13CH4]', '-Y', 'rmiso'), targets), set(targets))
        self.assertEqual(self.matches(self.generate('[O-]'), {'charged': '[O-]', 'neutral': 'O'}), {'charged'})
        smarts = self.generate('CC(=O)[O-]', '-g', 'all')
        self.assertNotIn(';-' , smarts)
        self.assertEqual(self.matches(smarts, {'acid': 'CC(=O)O'}), {'acid'})

    def test_element_override(self):
        smarts = self.generate('[SiH3]O', '-f', 'Si=[N]')
        self.assertTrue(smarts.startswith('[N]'), smarts)
        self.assertEqual(self.matches(smarts, {'match': 'NO', 'different': 'CO'}), {'match'})

    def test_relational_ranges(self):
        targets = {'original': 'CCO', 'substituted': 'CCCO', 'other_site': 'CC(C)O'}
        for mode in (['-s'], ['-m']):
            ordinary = self.generate('[1CH3]CO', *mode)
            ranged = self.generate('[1CH3]CO', *mode, '-Y', 'smtrange')
            self.assertNotIn('D>', ranged)
            self.assertRegex(ranged, r'D\{\d+-\}')
            self.assertEqual(self.matches(ordinary, targets), self.matches(ranged, targets))

    def test_failed_conversion_returns_failure(self):
        source = self.directory / 'bad.smi'
        source.write_text('invalid_smiles bad\n')
        result = subprocess.run([str(EXE), '-B', str(self.directory / 'bad.sma'), str(source)],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Smarts))
    if not result.wasSuccessful():
        sys.exit(1)
    print('mol2qry SMARTS checks passed')
