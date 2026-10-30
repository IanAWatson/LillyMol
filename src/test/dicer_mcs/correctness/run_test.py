#!/usr/bin/env python3
"""Focused correctness checks, isolated in a temporary directory per test."""
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
INDIR = Path(sys.argv[2]).resolve()


def record(name, parent, *fragments):
    return (f'name: {json.dumps(name)} smiles: {json.dumps(parent)} ' +
            ' '.join(f'fragment {{ smi: {json.dumps(smi)} }}' for smi in fragments) + '\n')


def numbered_chain(n):
    return 'C' + ''.join(f'[{i}CH2]' for i in range(1, n - 1)) + f'[{n - 1}CH3]'


class Correctness(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def run_tool(self, data, *options, success=True):
        path = self.root / 'input.textproto'
        path.write_text(data)
        result = subprocess.run([EXE, *options, str(path)], cwd=self.root,
                                capture_output=True, text=True, timeout=30)
        self.assertGreaterEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result

    def test_basic_atom_zero_mapping_and_multiple_increments(self):
        basic = (INDIR / 'basic.textproto').read_text()
        expected = self.run_tool(basic).stdout
        self.assertEqual([line.split()[1] for line in expected.splitlines()], ['one', 'two'])
        self.assertTrue(all(line.startswith('[1CH3]') for line in expected.splitlines()))
        for isotope in (200, 300):
            self.assertEqual(self.run_tool(basic.replace('[100CH4]', f'[{isotope}CH4]')).stdout,
                             expected)

    def test_invalid_parent_labels(self):
        for parent in ('C[99CH3]', 'CC', '[1CH3][2CH3]'):
            result = self.run_tool(record('bad', parent, '[100CH4]'), success=False)
            self.assertIn('isotopic', result.stderr)

    def test_no_surviving_fragments(self):
        for data, options in ((record('none', 'C[1CH3]'), ()),
                              (record('filtered', 'C[1CH3]', '[100CH4]'), ('-m', '2'))):
            result = self.run_tool(data, *options, success=False)
            self.assertIn('No fragments survive', result.stderr)
            self.assertEqual(result.stdout, '')

    def test_parent_without_fragments_does_not_prevent_other_matches(self):
        data = record('none', 'C[1CH3]') + (INDIR / 'basic.textproto').read_text()
        result = self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
        self.assertNotIn(' none', result.stdout)
        pairs = (self.root / 'pairs.smi').read_text()
        self.assertNotIn(' none', pairs)
        self.assertIn('one LHS', pairs)

    def test_distinct_parent_support(self):
        data = (record('one', 'C[1CH2][2CH3]', 'C[101CH3]', '[101CH3][2CH3]', '[100CH4]') +
                record('two', 'C[1CH2][2NH2]', '[100CH4]'))
        result = self.run_tool(data)
        self.assertIn('[1CH3]C in 1 molecules', result.stderr)
        self.assertIn('Most common fragment [1CH4] found 2 times', result.stderr)
        self.assertEqual([line.split()[1] for line in result.stdout.splitlines()], ['one', 'two'])

    def test_fragment_size_boundary(self):
        parent = numbered_chain(99)
        result = self.run_tool(record('ninety_nine', parent, parent))
        self.assertEqual(len(re.findall(r'\[\d+', result.stdout)), 99)
        parent = numbered_chain(100)
        result = self.run_tool(record('hundred', parent, parent), success=False)
        self.assertIn('fragments must contain fewer than 100 atoms', result.stderr)
        # A 100-atom parent with a small fragment remains supported.
        self.run_tool(record('hundred', parent, '[100CH4]'))

    def test_ambiguous_parent_encoding(self):
        result = self.run_tool(record('too_big', numbered_chain(101), '[100CH4]'), success=False)
        self.assertIn('attachment encoding supports at most 100 atoms', result.stderr)

    def test_invalid_minimum(self):
        data = (INDIR / 'basic.textproto').read_text()
        for value in ('0', '-1', 'invalid'):
            result = self.run_tool(data, '-m', value, success=False)
            self.assertIn('min fragment size', result.stderr)
            self.assertEqual(result.stdout, '')

    def test_pairwise_output_open_failure(self):
        result = self.run_tool((INDIR / 'basic.textproto').read_text(),
                               '-Y', 'missing_directory/pairs.smi', success=False)
        self.assertIn('cannot open', result.stderr)
        self.assertEqual(result.stdout, '')

    def test_expanded_pairwise_counts(self):
        data = (record('one', 'C[1CH2][2CH3]', '[101CH4]') +
                record('two', 'C[1CH2][2CH3]', '[101CH4]'))
        self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
        lines = (self.root / 'pairs.smi').read_text().splitlines()
        self.assertEqual(len(lines), 4)
        for line, count in zip(lines, (1, 1, 3, 3)):
            fields = line.split()
            self.assertEqual(int(fields[-2]), count)
            self.assertAlmostEqual(float(fields[-1]), count / 3, places=6)
            self.assertEqual(len(re.findall(r'\[\d+', fields[0])), count)

    def test_tied_selection_is_deterministic(self):
        fragments = ('C[101CH3]', '[101CH3][2NH2]')
        outputs = set()
        for order in (fragments, fragments[::-1]):
            data = (record('one', 'C[1CH2][2NH2]', *order) +
                    record('two', 'C[1CH2][2NH2]', *order))
            for _ in range(4):
                result = self.run_tool(data, '-Y', 'pairs.smi')
                chosen = next(line for line in result.stderr.splitlines()
                              if line.startswith('Most common fragment '))
                outputs.add((result.stdout, chosen, (self.root / 'pairs.smi').read_text()))
        self.assertEqual(len(outputs), 1)

    def support_example(self):
        return (record('one', 'C[1CH2][2CH2][3NH2]', '[100CH4]',
                       'C[101CH3]', 'C[1CH2][102CH3]') +
                record('two', 'C[1CH2][2OH]', '[100CH4]', 'C[101CH3]') +
                record('three', 'C[1CH2][2NH2]', '[100CH4]'))

    def test_minimum_support_selects_largest_qualifying_fragment(self):
        data = self.support_example()
        default = self.run_tool(data)
        self.assertIn('Most common fragment [1CH4] found 3 times', default.stderr)
        for threshold, atoms, parents in (('1', 3, ['one']),
                                          ('2', 2, ['one', 'two']),
                                          ('3', 1, ['one', 'two', 'three']),
                                          ('50%', 2, ['one', 'two']),
                                          ('67%', 1, ['one', 'two', 'three']),
                                          ('100%', 1, ['one', 'two', 'three'])):
            with self.subTest(threshold=threshold):
                result = self.run_tool(data, '-s', threshold)
                lines = result.stdout.splitlines()
                self.assertEqual([line.split()[1] for line in lines], parents)
                self.assertTrue(all(len(re.findall(r'\[\d+', line.split()[0])) == atoms
                                    for line in lines))
        result = self.run_tool(data, '-s', '4', success=False)
        self.assertIn('No fragment meets minimum support 4 of 3', result.stderr)
        self.assertEqual(result.stdout, '')

    def test_minimum_support_ties_use_support_before_smiles(self):
        data = (record('one', 'C[1CH2][2NH2]', 'C[101CH3]', '[101CH3][2NH2]') +
                record('two', 'C[1CH2][2NH2]', 'C[101CH3]', '[101CH3][2NH2]') +
                record('three', 'C[1CH2][2NH2]', '[101CH3][2NH2]'))
        result = self.run_tool(data, '-s', '2')
        self.assertIn('Selected fragment N[1CH3] found 3 times', result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 3)

    def test_percentage_support_uses_all_parents_and_exact_boundaries(self):
        # 14% of 100 is exactly 14, not 15 due to floating point rounding.
        data = ''.join(record(str(i), 'C[1CH2][2NH2]', '[100CH4]', 'C[101CH3]')
                       if i < 14 else record(str(i), 'C[1CH2][2NH2]', '[100CH4]')
                       for i in range(100))
        result = self.run_tool(data, '-s', '14%')
        self.assertEqual(len(result.stdout.splitlines()), 14)
        self.assertIn('Selected fragment [1CH3]C found 14 times', result.stderr)
        data = record('one', 'C[1CH3]', '[100CH4]') + record('empty', 'C[1CH3]')
        result = self.run_tool(data, '-s', '100%', success=False)
        self.assertIn('No fragment meets minimum support 2 of 2', result.stderr)
        self.run_tool(data, '-s', '50%')

    def test_invalid_support_thresholds(self):
        for threshold in ('0', '-1', '0%', '101%', 'invalid', '1.5', '0.5%', '%'):
            result = self.run_tool(self.support_example(), '-s', threshold, success=False)
            self.assertIn('support', result.stderr)
            self.assertEqual(result.stdout, '')
        self.run_tool(self.support_example(), '-s', '1', '-s', '2', success=False)

    def test_no_pairwise_match_is_not_an_error(self):
        data = record('carbon', 'C[1CH3]', '[100CH4]') + record('nitrogen', 'N[1NH2]', '[100NH3]')
        self.run_tool(data, '-Y', 'pairs.smi')
        self.assertEqual((self.root / 'pairs.smi').read_text(), '')


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Correctness)
    result = unittest.TextTestRunner().run(suite)
    if not result.wasSuccessful():
        sys.exit(1)
    print('dicer_mcs correctness checks passed')
