#!/usr/bin/env python3
"""Focused correctness checks, isolated in a temporary directory per test."""
import csv
import json
import os
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

    def assert_smarts_matches(self, result, molecules, expected):
        smarts = next(line[:-len(' MCS smarts')]
                      for line in result.stderr.splitlines() if line.endswith(' MCS smarts'))
        self.assertNotIn('D>0', smarts)
        query_tool = Path(EXE).with_name('tsubstructure')
        self.assertTrue(query_tool.is_file(), 'Build/install tsubstructure alongside dicer_mcs')
        path = self.root / 'query_targets.smi'
        path.write_text(molecules)
        matched = self.root / 'matched.smi'
        search = subprocess.run([str(query_tool), '-s', smarts, '-m', str(matched), str(path)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(search.returncode, 0, search.stderr)
        names = {line.split()[1] for line in matched.read_text().splitlines()}
        self.assertEqual(names, set(expected), (smarts, search.stderr))
        return smarts

    def assert_highlight_count(self, line, count):
        self.assertEqual(re.findall(r'\[(\d+)', line.split()[0]), ['1'] * count)

    def test_select_fragment_identity_whole_parent_and_no_parent_search(self):
        data = (record('cut example', '[3OH][2CH2][1CH2]N', 'N[1CH2][102CH3]') +
                record('whole', 'N[1CH2][2CH3]') +
                record('contains_but_not_emitted', 'N[1CH2][2CH2][3CH3]'))
        result = self.run_tool(data, '-F', 'CCN')
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 2)
        self.assertIn(' cut example nmatch=1', lines[0])
        self.assertIn(' whole nmatch=1', lines[1])
        for line in lines:
            self.assert_highlight_count(line, 3)
        self.assertIn('found in 2 of 3 parents; 2 distinct occurrences', result.stderr)

    def test_select_union_separate_and_duplicate_atom_sets(self):
        left, right = 'N[1CH2][102CH3]', '[6NH2][5CH2][104CH3]'
        cap_left = 'N[1CH2][10102CH3]'
        parent = 'N[1CH2][2CH2][3CH2][4CH2][5CH2][6NH2]'
        outputs = []
        for fragments in [(left, right, left, cap_left), (cap_left, left, right, left)]:
            data = record('two', parent, *fragments)
            union = self.run_tool(data, '-F', 'NCC')
            self.assertEqual(len(union.stdout.splitlines()), 1)
            self.assertIn('nmatch=2', union.stdout)
            self.assert_highlight_count(union.stdout, 6)
            separate = self.run_tool(data, '-F', 'NCC', '-F', 'separate')
            lines = separate.stdout.splitlines()
            self.assertEqual(len(lines), 2)
            for i, line in enumerate(lines, 1):
                self.assertIn(f'nmatch=2 match={i}/2', line)
                self.assert_highlight_count(line, 3)
            outputs.append((union.stdout, separate.stdout))
        self.assertEqual(outputs[0], outputs[1])

    def test_select_overlapping_occurrences(self):
        data = record('overlap', 'C[1CH2][2CH2][3CH3]',
                      'C[1CH2][102CH3]', '[1CH3][2CH2][3CH3]')
        result = self.run_tool(data, '-F', 'CCC')
        self.assertIn('nmatch=2', result.stdout)
        self.assert_highlight_count(result.stdout, 4)
        result = self.run_tool(data, '-F', 'CCC', '-F', 'separate')
        for line in result.stdout.splitlines():
            self.assert_highlight_count(line, 3)

    def test_select_attachment_positions_and_unlabelled_atoms_unconstrained(self):
        parent = 'N[1CH2][2CH2][3OH]'
        middle = '[100NH2][101CH2][2CH3]'
        terminal = 'N[1CH2][102CH3]'
        data = (record('middle', parent, middle) + record('terminal', parent, terminal) +
                record('plain', parent, 'N[1CH2][2CH3]'))
        for smiles, expected in [('N[1CH2]C', 'middle'), ('NC[1CH3]', 'terminal'),
                                  ('N[7CH2]C', 'middle')]:
            result = self.run_tool(data, '-F', smiles)
            self.assertEqual(len(result.stdout.splitlines()), 1)
            self.assertIn(f' {expected} nmatch=1', result.stdout)
        self.assertEqual(len(self.run_tool(data, '-F', 'NCC').stdout.splitlines()), 3)
        self.run_tool(data, '-F', 'N[1CH2][1CH3]', success=False)
        # Every occurrence must be considered, not just the first representative.
        for fragments in [(middle, terminal), (terminal, middle)]:
            result = self.run_tool(record('both', parent, *fragments), '-F', 'NC[1CH3]')
            self.assertIn(' both nmatch=1', result.stdout)

    def test_select_attachment_symmetry(self):
        data = record('symmetric', 'C[1CH2][2CH2][3CH3]',
                      'C[101CH3]', '[102CH3][3CH3]')
        outputs = []
        for pattern in ['[1CH3]C', 'C[1CH3]']:
            result = self.run_tool(data, '-F', pattern)
            self.assertIn('nmatch=2', result.stdout)
            self.assert_highlight_count(result.stdout, 4)
            outputs.append(result.stdout)
        self.assertEqual(outputs[0], outputs[1])
        self.run_tool(data, '-F', '[1CH3][1CH3]', success=False)

    def test_select_caps_highlight_original_parent_atoms(self):
        data = record('cap', 'N[1CH2][2C]([3F])([4F])[5F]',
                      'N[101CH2][10102CH3]')
        result = self.run_tool(data, '-F', 'NC[1CH3]')
        self.assertIn('[1C]', result.stdout)
        self.assertNotIn('[1F]', result.stdout)
        self.assertEqual(result.stdout.count('F'), 3)
        self.assert_highlight_count(result.stdout, 3)

    def test_select_no_match_invalid_input_and_conflicting_options(self):
        data = record('parent', 'N[1CH2][2CH3]')
        result = self.run_tool(data, '-F', 'CCC', success=False)
        self.assertEqual(result.stdout, '')
        self.assertIn('No parents contributed', result.stderr)
        for args in [('-F', 'separate'), ('-F', ''), ('-F', 'invalid'), ('-F', 'C.C'),
                     ('-F', 'NC', '-F', 'CC')]:
            self.run_tool(data, *args, success=False)
        for args in [('-m', '2'), ('-s', '50%'), ('-R', 'report'), ('-X', 'missing'),
                     ('-Y', 'pairs'), ('-P', 'UST:Z')]:
            self.run_tool(data, '-F', 'NCC', *args, success=False)
        self.run_tool('not a proto\n', '-F', 'NC', success=False)
        self.run_tool(record('bad', 'N[2CH3]'), '-F', 'NC', success=False)
        self.run_tool('', '-F', 'NC', success=False)
        result = subprocess.run([EXE, '-F', 'NC', str(self.root / 'absent')],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)

    def test_select_stdin_and_multiple_input_files(self):
        data = record('stdin', 'N[1CH3]')
        result = subprocess.run([EXE, '-F', 'NC', '-'], input=data,
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(' stdin nmatch=1', result.stdout)
        first, second = self.root / 'first.textproto', self.root / 'second.textproto'
        first.write_text(record('absent', 'N[1NH2]'))
        second.write_text(data)
        result = subprocess.run([EXE, '-F', 'NC', str(first), str(second)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 1)
        self.assertIn('found in 1 of 2 parents', result.stderr)

    def test_attachment_can_match_hydrogen(self):
        plain = record('plain', 'C[1CH3]', 'C[1CH3]')
        cut = record('cut', 'C[1CH2][2CH3]', 'C[101CH3]')
        results = []
        for data in (plain + cut, cut + plain):
            result = self.run_tool(data, '-s', '100%', '-Y', 'pairs.smi')
            self.assertIn('Selected fragment CC found 2 times', result.stderr)
            self.assertEqual(len(result.stdout.splitlines()), 2)
            pairs = (self.root / 'pairs.smi').read_text().splitlines()
            self.assertEqual(len(pairs), 2)
            self.assertTrue(all(int(line.split()[-2]) == 2 for line in pairs))
            results.append(self.assert_smarts_matches(
                result, 'CC ethane\nCCC propane\nC(C)(C)C branched\nCO wrong_element\nC=C wrong_bond\n',
                ('ethane', 'propane', 'branched')))
        self.assertEqual(results[0], results[1])

    def test_different_cut_positions_preserve_atom_correspondence(self):
        data = (record('plain', 'C[1CH2][2OH]', 'C[1CH2][2OH]') +
                record('carbon_cut', 'C[1CH]([3CH3])[2OH]', 'C[101CH2][2OH]') +
                record('oxygen_cut', 'C[1CH2][2O][3CH3]', 'C[1CH2][102OH]'))
        result = self.run_tool(data, '-s', '100%')
        self.assertIn('Selected fragment OCC found 3 times', result.stderr)
        lines = result.stdout.splitlines()
        self.assertEqual(len(lines), 3)
        mappings = [dict(re.findall(r'\[(\d+)([A-Z][a-z]?)', line.split()[0])) for line in lines]
        self.assertTrue(all(len(mapping) == 3 for mapping in mappings))
        self.assertTrue(all(mapping == mappings[0] for mapping in mappings))
        self.assert_smarts_matches(result,
                                  'CCO plain\nCC(C)O carbon_cut\nCCOC oxygen_cut\n'
                                  'CCN wrong_element\nCC=O wrong_bond\n',
                                  ('plain', 'carbon_cut', 'oxygen_cut'))

    def test_aromatic_fragment_with_and_without_cut(self):
        data = (record('benzene', 'c1[1cH][2cH][3cH][4cH][5cH]1',
                       'c1[1cH][2cH][3cH][4cH][5cH]1') +
                record('toluene', 'C[1c]1[2cH][3cH][4cH][5cH][6cH]1',
                       '[101cH]1[2cH][3cH][4cH][5cH][6cH]1'))
        result = self.run_tool(data, '-s', '100%')
        self.assertEqual(len(result.stdout.splitlines()), 2)
        self.assertTrue(all(len(re.findall(r'\[\d+', line)) == 6
                            for line in result.stdout.splitlines()))
        self.assert_smarts_matches(result,
                                  'c1ccccc1 benzene\nCc1ccccc1 toluene\n'
                                  'Cc1ccc(C)cc1 xylene\nC1CCCCC1 nonaromatic\n',
                                  ('benzene', 'toluene', 'xylene'))

    def test_single_atom_smarts_allows_zero_connections(self):
        result = self.run_tool((INDIR / 'basic.textproto').read_text())
        self.assert_smarts_matches(result, 'C methane\nCC ethane\nCO methanol\nN nitrogen\n',
                                  ('methane', 'ethane', 'methanol'))

    def test_genuine_terminal_atoms_retain_their_symbols(self):
        result = self.run_tool(record('methylamine', 'N[1CH3]'))
        smarts = self.assert_smarts_matches(
            result,
            'NC methylamine\nNc1ccccc1 aniline\nCn1cccc1 methylpyrrole\n'
            'OC wrong_element\nN=C wrong_bond\n',
            ('methylamine',))
        self.assertIn('[C]', smarts)
        self.assertIn('[N]', smarts)
        self.assertNotIn('#', smarts)

    def test_only_marked_caps_match_aromatic_and_aliphatic_elements(self):
        examples = [
            ('N[1c]1[2cH][3cH][4cH][5cH][6cH]1', '[100NH2][10101CH3]',
             '[#6]', '[N]', ('methylamine', 'aniline')),
            ('C[1n]1[2cH][3cH][4cH][5cH]1', '[100CH3][10101NH2]',
             '[#7]', '[C]', ('methylamine', 'methylpyrrole'))]
        for parent, fragment, cap, genuine, expected in examples:
            with self.subTest(cap=cap):
                data = record('capped', parent, fragment) + record('plain', 'N[1CH3]')
                result = self.run_tool(data, '-s', '100%', '-R', 'fragments.out')
                smarts = self.assert_smarts_matches(
                    result, 'NC methylamine\nNc1ccccc1 aniline\nCn1cccc1 methylpyrrole\n'
                    'OC wrong_element\nN=C wrong_bond\n', expected)
                self.assertIn(cap, smarts)
                self.assertIn(genuine, smarts)
                row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'NC')
                self.assertEqual(row['SMARTS'], smarts)

    def test_internal_atom_aromaticity_is_not_relaxed(self):
        data = (record('benzylamine', 'N[1CH2][2c]1[3cH][4cH][5cH][6cH][7cH]1',
                       'N[101CH2][10102CH3]') + record('ethylamine', 'N[1CH2][2CH3]'))
        result = self.run_tool(data, '-s', '100%')
        smarts = self.assert_smarts_matches(
            result,
            'NCC ethylamine\nNCc1ccccc1 benzylamine\nNc1ccccc1 aniline\n'
            'NCO wrong_element\n', ('ethylamine', 'benzylamine'))
        self.assertIn('[#6]', smarts)
        self.assertIn('[CD>1]', smarts)

    def test_same_parent_subset_prefers_fewer_caps_for_smarts(self):
        two_caps = '[10100NH2][101CH2][10102CH3]'
        one_cap = '[100NH2][101CH2][10102CH3]'
        for fragments in [(two_caps, one_cap), (one_cap, two_caps)]:
            data = (record('fluorinated', 'N[1CH2][2C]([3F])([4F])[5F]', *fragments) +
                    record('ethylamine', 'N[1CH2][2CH3]'))
            result = self.run_tool(data, '-s', '100%', '-R', 'fragments.out')
            self.assertIn('[N]', result.stderr)
            self.assertNotIn('[#7]', result.stderr)
            row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'NCC')
            self.assertIn('[N]', row['SMARTS'])
            self.assertIn('[#6]', row['SMARTS'])

    def test_whole_parent_has_genuine_atoms_even_if_capped_identity_exists(self):
        data = record('propylamine', 'N[1CH2][2CH2][3CH3]',
                      '[10100NH2][101CH2][2CH2][3CH3]')
        result = self.run_tool(data, '-R', 'fragments.out')
        self.assertNotIn('[#7]', result.stderr)
        row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'NCCC')
        self.assertIn('[N]', row['SMARTS'])
        self.assertNotIn('#', row['SMARTS'])
        self.assertEqual(row['NWHOLE'], '1')
        self.assertEqual(row['NCUT'], '1')

    def assert_induced_pair(self, lines):
        self.assertEqual(len(lines), 2)
        fields = [line.split() for line in lines]
        isotopes = [set(map(int, re.findall(r'\[(\d+)', row[0]))) for row in fields]
        self.assertEqual(isotopes[0], isotopes[1])
        labels = sorted(isotopes[0])
        self.assertTrue(all(len(labels) == int(row[-2]) for row in fields))
        if len(labels) < 2:
            return
        # Ask LillyMol about each edge and bond order between labelled atoms;
        # do not duplicate a SMILES parser in the test. Aromaticity takes
        # precedence over a particular Kekule single/double assignment.
        args = [str(Path(EXE).with_name('tsubstructure')), '-a']
        for i, iso1 in enumerate(labels):
            for iso2 in labels[:i]:
                for bond in ('~', '-', '=', '#', ':'):
                    args.extend(('-s', f'[{iso1}*]{bond}[{iso2}*]'))
        path = self.root / 'pair_targets.smi'
        path.write_text(''.join(f'{row[0]} {row[1]}\n' for row in fields))
        result = subprocess.run([*args, str(path)], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [line.split() for line in result.stdout.splitlines()]
        self.assertEqual(len(rows), 3, result.stdout)
        values = [list(map(int, row[1:])) for row in rows[1:]]
        for offset in range(0, len(values[0]), 5):
            lhs = values[0][offset:offset + 5]
            rhs = values[1][offset:offset + 5]
            if lhs[-1] or rhs[-1]:
                self.assertEqual((lhs[0], lhs[-1]), (rhs[0], rhs[-1]), lines)
            else:
                self.assertEqual(lhs, rhs, lines)

    def test_induced_ring_chain_in_both_directions(self):
        cases = (('C1[1CH2][2CH2]1', 'C[1CH2][2CH3]', 2),
                 ('C1[1CH2][2CH2][3CH2]1', 'C[1CH2][2CH2][3CH3]', 3))
        for ring, chain, expected in cases:
            ring_record = record('ring', ring, '[100CH4]')
            chain_record = record('chain', chain, '[101CH4]')
            for data in (ring_record + chain_record, chain_record + ring_record):
                self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
                lines = (self.root / 'pairs.smi').read_text().splitlines()
                self.assertEqual([int(line.split()[-2]) for line in lines],
                                 [1, 1, expected, expected])
                self.assert_induced_pair(lines[2:])

    def test_induced_ring_closures_and_bond_orders(self):
        triangle = 'C1[1CH2][2CH2]1'
        double_triangle = 'C1[1CH]=[2CH]1'
        for lhs, rhs, expected in ((triangle, 'C1[1CH2][2CH]1[3CH3]', 3),
                                   (triangle, double_triangle, 2),
                                   (double_triangle, triangle, 2),
                                   ('C=[1CH][2CH3]', 'C=[1CH][2CH2][3OH]', 3),
                                   ('C#[1C][2CH3]', 'C#[1C][2CH2][3OH]', 3)):
            self.run_tool(record('lhs', lhs, '[100CH4]') + record('rhs', rhs, '[100CH4]'),
                          '-P', 'UST:Z', '-Y', 'pairs.smi')
            lines = (self.root / 'pairs.smi').read_text().splitlines()
            self.assertEqual([int(line.split()[-2]) for line in lines],
                             [1, 1, expected, expected])
            self.assert_induced_pair(lines[2:])

    def test_noninduced_seed_is_skipped_and_smaller_valid_seed_is_tried(self):
        for extra in ((), ('[100CH4]',)):
            data = (record('ring', 'C1[1CH2][2CH2]1', '[100CH3][1CH2][102CH3]', *extra) +
                    record('chain', 'C[1CH2][2CH3]', 'C[1CH2][2CH3]', *extra))
            self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
            lines = (self.root / 'pairs.smi').read_text().splitlines()
            if not extra:
                self.assertEqual(lines, [])
            else:
                self.assertEqual([int(line.split()[-2]) for line in lines], [1, 1, 2, 2])
                self.assert_induced_pair(lines[2:])

    def test_expansion_tries_next_rhs_neighbor_after_incompatible_closure(self):
        # At rhs atom 1, the ring-closure neighbour precedes the methyl branch.
        # The first would add an extra edge; the methyl is a valid continuation.
        data = (record('chain', 'C[1CH2][2CH3]', 'C[101CH3]') +
                record('ring', 'C1[1CH]([3CH2]1)[2CH3]', '[100CH3][101CH3]'))
        self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
        lines = (self.root / 'pairs.smi').read_text().splitlines()
        self.assertEqual([int(line.split()[-2]) for line in lines], [2, 2, 3, 3])
        self.assert_induced_pair(lines[2:])

    def test_induced_aromatic_bonds_ignore_kekule_assignment(self):
        lhs = 'C1=[1CH][2CH]=[3CH][4CH]=[5CH]1'
        rhs = 'C1[1CH]=[2CH][3CH]=[4CH][5CH]=1'
        self.run_tool(record('lhs', lhs, lhs) + record('rhs', rhs, rhs),
                      '-P', 'UST:Z', '-Y', 'pairs.smi')
        lines = (self.root / 'pairs.smi').read_text().splitlines()
        self.assertEqual([int(line.split()[-2]) for line in lines], [6, 6])
        self.assert_induced_pair(lines)

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
        for data, options in ((record('none', 'C[1CH3]'), ('-m', '3')),
                              (record('filtered', 'C[1CH3]', '[100CH4]'), ('-m', '3'))):
            result = self.run_tool(data, *options, success=False)
            self.assertIn('No fragments survive', result.stderr)
            self.assertEqual(result.stdout, '')

    def test_parent_without_emitted_fragments_contributes_zero_cut_candidate(self):
        data = (record('whole', 'N[1CH2][2CH3]') +
                record('cut', 'N[1CH2][2CH2][3CH3]', 'N[1CH2][102CH3]'))
        result = self.run_tool(data, '-m', '2', '-s', '100%', '-P', 'UST:Z', '-Y', 'pairs.smi')
        self.assertEqual([line.split()[1] for line in result.stdout.splitlines()], ['whole', 'cut'])
        self.assertTrue(all(len(re.findall(r'\[\d+', line)) == 3
                            for line in result.stdout.splitlines()))
        pairs = (self.root / 'pairs.smi').read_text().splitlines()
        self.assertEqual([int(line.split()[-2]) for line in pairs], [3, 3])
        self.assert_induced_pair(pairs)

    def test_zero_cut_candidate_does_not_double_count_explicit_whole_fragment(self):
        data = record('one', 'C[1CH3]', 'C[1CH3]') + record('two', 'C[1CH3]')
        result = self.run_tool(data, '-s', '100%', '-Y', 'pairs.smi')
        self.assertIn('Selected fragment CC found 2 times', result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 2)
        self.assertEqual(len((self.root / 'pairs.smi').read_text().splitlines()), 2)

    def test_ethylamine_and_tertbutylamine_share_whole_smaller_parent(self):
        data = (INDIR / 'ethyl_tertbutyl.textproto').read_text()
        for options in (('-m', '2'), ('-m', '2', '-s', '100%')):
            result = self.run_tool(data, *options, '-Y', 'pairs.smi')
            self.assertEqual([line.split()[1] for line in result.stdout.splitlines()], ['p1a', 'p1b'])
            self.assertTrue(all(len(re.findall(r'\[\d+', line)) == 3
                                for line in result.stdout.splitlines()))
            self.assertIn('fragment NCC found 2 times', result.stderr)
            self.assert_smarts_matches(result, 'NCC ethyl\nNC(C)(C)C tertbutyl\nCCC wrong_element\n',
                                      ('ethyl', 'tertbutyl'))
            pairs = (self.root / 'pairs.smi').read_text().splitlines()
            self.assertEqual([int(line.split()[-2]) for line in pairs], [3, 3])
            self.assert_induced_pair(pairs)

    def test_distinct_parent_support(self):
        data = (record('one', 'C[1CH2][2CH3]', 'C[101CH3]', '[101CH3][2CH3]', '[100CH4]') +
                record('two', 'C[1CH2][2NH2]', '[100CH4]'))
        result = self.run_tool(data, '-R', 'fragments.tsv')
        with (self.root / 'fragments.tsv').open() as report:
            rows = {row['SMILES']: row for row in csv.DictReader(report, delimiter='\t')}
        self.assertEqual(rows['CC']['NPARENTS'], '1')
        self.assertEqual(rows['CC']['NCUT'], '1')
        self.assertNotIn(' in 1 molecules', result.stderr)
        self.assertIn('Most common fragment C found 2 times', result.stderr)
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
        result = self.run_tool(record('hundred', parent), success=False)
        self.assertIn('No fragments survive', result.stderr)
        self.run_tool(record('whole_99', numbered_chain(99)))
        self.run_tool(record('filtered_99', numbered_chain(99)), '-m', '100', success=False)

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
                record('two', 'C[1CH2][2CH2][3OH]', '[101CH4]'))
        self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
        lines = (self.root / 'pairs.smi').read_text().splitlines()
        self.assertEqual(len(lines), 4)
        for line, count, denominator in zip(lines, (1, 1, 3, 3), (3, 4, 3, 4)):
            fields = line.split()
            self.assertEqual(int(fields[-2]), count)
            self.assertAlmostEqual(float(fields[-1]), count / denominator, places=6)
            self.assertEqual(len(re.findall(r'\[\d+', fields[0])), count)

    def test_best_occurrence_is_found_on_either_side(self):
        reference = {}
        for fragments in (('[103CH4]', '[101CH4]'), ('[101CH4]', '[103CH4]')):
            repeated = record('repeated', 'O[1CH2][2CH2][3CH2][4NH2]', *fragments)
            single = record('single', 'O[1CH2][2CH3]', '[101CH4]')
            for repeated_first in (True, False):
                data = repeated + single if repeated_first else single + repeated
                self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
                pairs = (self.root / 'pairs.smi').read_text()
                lines = pairs.splitlines()
                self.assertEqual(len(lines), 4)
                self.assertEqual([int(line.split()[-2]) for line in lines], [1, 1, 3, 3])
                # The better occurrence includes the OH end, not the amine end.
                for line in lines[2:]:
                    smiles = line.split()[0]
                    self.assertRegex(smiles, r'^\[\d+OH\]')
                    self.assertNotRegex(smiles, r'\[\d+N')
                    self.assertEqual(len(re.findall(r'\[\d+', smiles)), 3)
                if repeated_first in reference:
                    self.assertEqual(pairs, reference[repeated_first])
                reference[repeated_first] = pairs

    def test_all_occurrence_pairings_and_trial_reset(self):
        reference = None
        lhs_occurrences = ('[103CH4]', '[101CH4]')
        rhs_occurrences = ('[103CH4]', '[101CH4]')
        for lhs_order in (lhs_occurrences, lhs_occurrences[::-1]):
            for rhs_order in (rhs_occurrences, rhs_occurrences[::-1]):
                data = (record('lhs', 'O[1CH2][2CH2][3CH2][4NH2]', *lhs_order) +
                        record('rhs', 'O[1CH2][2CH3].[3CH3][4F]', *rhs_order))
                result = self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
                pairs = (self.root / 'pairs.smi').read_text()
                lines = pairs.splitlines()
                self.assertEqual(len(lines), 4)
                self.assertEqual([int(line.split()[-2]) for line in lines], [1, 1, 3, 3])
                for line in lines[2:]:
                    self.assertEqual(len(re.findall(r'\[\d+', line.split()[0])), 3)
                    self.assertNotRegex(line.split()[0], r'\[\d+[NF]')
                self.assertIn('Most common fragment C found 2 times', result.stderr)
                if reference is not None:
                    self.assertEqual(pairs, reference)
                reference = pairs

    def test_occurrence_selection_without_expansion_is_deterministic(self):
        outputs = set()
        for order in (('[103CH4]', '[101CH4]'), ('[101CH4]', '[103CH4]')):
            data = (record('repeated', 'O[1CH2][2CH2][3CH2][4NH2]', *order) +
                    record('single', 'O[1CH2][2CH3]', '[101CH4]'))
            result = self.run_tool(data, '-Y', 'pairs.smi')
            pairs = (self.root / 'pairs.smi').read_text()
            self.assertEqual(len(pairs.splitlines()), 2)
            self.assertTrue(all(int(line.split()[-2]) == 1 for line in pairs.splitlines()))
            outputs.add((result.stdout, pairs))
        self.assertEqual(len(outputs), 1)

    def test_largest_shared_seed_takes_priority_over_smaller_seed_expansion(self):
        data = (record('lhs', 'C[1CH3].[2OH][3CH2][4CH3]',
                       'C[1CH3]', '[103CH4]') +
                record('rhs', 'C[1CH2][2OH]', 'C[101CH3]', '[101CH4]'))
        self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
        lines = (self.root / 'pairs.smi').read_text().splitlines()
        # The ethane seed cannot expand on lhs. The smaller C seed could match
        # ethanol, but is intentionally outside the largest-seed search.
        self.assertEqual(len(lines), 2)
        self.assertTrue(all(int(line.split()[-2]) == 2 for line in lines))

    def test_equal_size_fragment_identities_compete_by_expansion(self):
        data = (record('lhs', 'C[1CH2][2OH].[3NH2][4CH3]', '[100CH4]', '[103NH3]') +
                record('rhs', 'C[1F].[2NH2][3CH2][4OH]', '[100CH4]', '[102NH3]'))
        self.run_tool(data, '-P', 'UST:Z', '-Y', 'pairs.smi')
        lines = (self.root / 'pairs.smi').read_text().splitlines()
        self.assertEqual(len(lines), 4)
        self.assertEqual([int(line.split()[-2]) for line in lines], [1, 1, 2, 2])
        self.assertTrue(all('[1NH2]' in line.split()[0] for line in lines[:2]))

    def test_tied_selection_is_deterministic(self):
        fragments = ('C[101CH3]', '[101CH3][2NH2]')
        outputs = set()
        for order in (fragments, fragments[::-1]):
            data = (record('one', 'C[1CH2][2NH2]', *order) +
                    record('two', 'C[1CH2][2NH][3CH3]', *order))
            for _ in range(4):
                result = self.run_tool(data, '-Y', 'pairs.smi')
                chosen = next(line for line in result.stderr.splitlines()
                              if line.startswith('Most common fragment '))
                outputs.add((result.stdout, chosen, (self.root / 'pairs.smi').read_text()))
        self.assertEqual(len(outputs), 1)

    def test_nonmatching_output_uses_candidate_membership(self):
        data = (record('hit', 'C[1CH2][2OH]', 'C[101CH3]') +
                record('second_hit', 'C[1CH2][2F]', 'C[101CH3]') +
                record('contains_core_but_not_emitted', 'C[1CH2][2NH2]', '[100CH4]') +
                record('no_fragments', 'N[1NH2]'))
        result = self.run_tool(data, '-s', '2', '-X', 'missing.smi', '-Y', 'pairs.smi')
        self.assertEqual([line.split()[1] for line in result.stdout.splitlines()], ['hit', 'second_hit'])
        self.assertEqual((self.root / 'missing.smi').read_text(),
                         'CCN contains_core_but_not_emitted\nNN no_fragments\n')
        self.assertEqual(len((self.root / 'pairs.smi').read_text().splitlines()), 2)

    def read_fragment_report(self, delimiter='\t', filename='fragments.out'):
        with (self.root / filename).open(newline='') as report:
            return list(csv.DictReader(report, delimiter=delimiter))

    def test_fragment_report_fields_sorting_and_provenance(self):
        data = (INDIR / 'basic.textproto').read_text()
        baseline = self.run_tool(data)
        result = self.run_tool(data, '-R', 'fragments.out')
        self.assertEqual(result.stdout, baseline.stdout)
        rows = self.read_fragment_report()
        self.assertEqual([row['SMILES'] for row in rows], ['C', 'CC', 'NC'])
        carbon, ethane, methylamine = rows
        self.assertEqual(carbon, {'SMILES': 'C', 'NATOMS': '1', 'NPARENTS': '2',
                                  'FPARENTS': '1', 'FIRST_PARENT': 'one', 'SMARTS': '[C]',
                                  'NWHOLE': '0', 'NCUT': '2', 'SELECTED': '1',
                                  'ATTACHMENT_SMILES': '[1CH4]',
                                  'EXAMPLE_SMILES': '[1CH3]C'})
        self.assertEqual(ethane['NWHOLE'], '1')
        self.assertEqual(ethane['NCUT'], '0')
        self.assertEqual(ethane['NPARENTS'], '1')
        self.assertEqual(ethane['FPARENTS'], '0.5')
        self.assertEqual(ethane['FIRST_PARENT'], 'one')
        self.assertEqual(methylamine['FIRST_PARENT'], 'two')
        self.assertEqual(sum(int(row['SELECTED']) for row in rows), 1)
        self.assertIn('Found 3 unique candidate fragments', result.stderr)
        self.assertNotIn(' in 1 molecules', result.stderr)

    def test_fragment_report_whole_and_cut_support(self):
        data = (INDIR / 'ethyl_tertbutyl.textproto').read_text()
        self.run_tool(data, '-m', '2', '-s', '100%', '-R', 'fragments.out')
        rows = {row['SMILES']: row for row in self.read_fragment_report()}
        core = rows['NCC']
        self.assertEqual(core['NPARENTS'], '2')
        self.assertEqual(core['NWHOLE'], '1')
        self.assertEqual(core['NCUT'], '1')  # Multiple occurrences in p1b count once.
        self.assertEqual(core['FIRST_PARENT'], 'p1a')
        self.assertEqual(core['SELECTED'], '1')
        self.assertTrue(all(int(row['NATOMS']) >= 2 for row in rows.values()))
        self.assertTrue(any(int(row['NPARENTS']) == 1 for row in rows.values()))
        data = record('first', 'C[1CH3]', 'C[1CH3]') + record('second', 'C[1CH3]')
        self.run_tool(data, '-R', 'fragments.out')
        row = self.read_fragment_report()[0]
        self.assertEqual((row['NPARENTS'], row['NWHOLE'], row['NCUT']), ('2', '2', '0'))

    def test_fragment_report_example_maps_caps_to_original_parent(self):
        # The endpoint cap is carbon 2, which remains a CF3 carbon in the example.
        data = record('fluorinated', 'N[1CH2][2C]([3F])([4F])[5F]',
                      'N[101CH2][102CH3]', 'N[101CH3]')
        baseline = self.run_tool(data).stdout
        result = self.run_tool(data, '-R', 'fragments.out')
        self.assertEqual(result.stdout, baseline)
        rows = {row['SMILES']: row for row in self.read_fragment_report()}
        example = rows['NCC']['EXAMPLE_SMILES']
        self.assertIn('[1NH2]', example)
        self.assertIn('[1CH2]', example)
        self.assertIn('[1C]', example)
        self.assertNotIn('[1F]', example)
        self.assertEqual(re.findall(r'\[(\d+)', example), ['1'] * 3)
        self.assertEqual(example.count('F'), 3)
        # Each row starts from a clean parent copy, not the preceding row's labels.
        for row in rows.values():
            self.assertEqual(re.findall(r'\[(\d+)', row['EXAMPLE_SMILES']),
                             ['1'] * int(row['NATOMS']))

    def test_fragment_report_example_uses_first_contributing_parent(self):
        big = record('big', 'N[1CH2][2CH3]', 'N[101CH3]')
        small = record('small', 'N[1CH3]')
        for data, name, carbons in [(big + small, 'big', 2),
                                    (small + big, 'small', 1)]:
            self.run_tool(data, '-R', 'fragments.out')
            row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'NC')
            self.assertEqual(row['FIRST_PARENT'], name)
            self.assertEqual(row['EXAMPLE_SMILES'].count('C'), carbons)
            self.assertEqual(re.findall(r'\[(\d+)', row['EXAMPLE_SMILES']), ['1'] * 2)

    def test_cap_marker_preserves_identity_mapping_and_attachment_payload(self):
        parent = 'N[1CH2][2CH2][3OH]'
        baseline = self.run_tool(record('parent', parent, '[100NH2][101CH2][302CH3]'),
                                 '-R', 'fragments.out')
        rows = self.read_fragment_report()
        for caps in ['[10100NH2][101CH2][10302CH3]',
                     '[100NH2][101CH2][10302CH3]']:
            result = self.run_tool(record('parent', parent, caps), '-R', 'fragments.out')
            self.assertEqual(result.stdout, baseline.stdout)
            current = self.read_fragment_report()
            self.assertEqual([{k: v for k, v in r.items() if k != 'SMARTS'} for r in current],
                             [{k: v for k, v in r.items() if k != 'SMARTS'} for r in rows])
        # A marker can have a zero payload, mapping a cap to parent atom zero.
        self.run_tool(record('parent', parent, '[10000NH2][101CH2][102CH3]'),
                      '-R', 'fragments.out')
        row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'NCC')
        self.assertEqual(row['NPARENTS'], '1')
        self.assertIn('[1NH2]', row['EXAMPLE_SMILES'])

    def test_cap_markers_do_not_split_support(self):
        data = (record('capped', 'N[1CH2][2CH2][3OH]', '[100NH2][101CH2][10102CH3]') +
                record('plain', 'N[1CH2][2CH3]'))
        result = self.run_tool(data, '-m', '2', '-s', '100%', '-R', 'fragments.out')
        self.assertIn('Selected fragment NCC found 2 times', result.stderr)
        row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'NCC')
        self.assertEqual(row['NPARENTS'], '2')

    def test_invalid_cap_marker_payloads_and_parent_markers(self):
        for fragment in ['N[10101CH2][2CH3]', 'N[1CH2][20102CH3]']:
            self.run_tool(record('bad', 'N[1CH2][2CH3]', fragment), success=False)
        self.run_tool(record('bad_parent', '[10000NH2][1CH3]'), success=False)

    def test_sparse_support_gaps_duplicates_and_parent_flags(self):
        contributors = {1, 7, 22, 38}
        data = ''
        for i in range(40):
            if i == 1:
                data += record(str(i), 'C[1CH3]', '[100CH3][101CH3]',
                               '[100CH3][101CH3]')
            elif i == 7:
                data += record(str(i), 'C[1CH2][2NH2]', '[100CH3][101CH3]',
                               'C[101CH3]')
            elif i == 22:
                data += record(str(i), 'C[1CH3]')
            elif i == 38:
                data += record(str(i), 'C[1CH2][2NH2]', 'C[1CH3]')
            else:
                data += record(str(i), 'N[1NH2]')
        result = self.run_tool(data, '-R', 'fragments.out', '-X', 'missing')
        row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'CC')
        self.assertEqual((row['NPARENTS'], row['NWHOLE'], row['NCUT']), ('4', '2', '2'))
        self.assertEqual(row['FIRST_PARENT'], '1')
        self.assertEqual({int(s.split()[1]) for s in
                          (self.root / 'missing.smi').read_text().splitlines()}, contributors)
        self.assertEqual(len(result.stdout.splitlines()), 36)

    def test_inverse_numbering_handles_parent_order_and_zero(self):
        data = record('permuted', '[2OH][1CH2]N', 'N[101CH3]')
        self.run_tool(data, '-R', 'fragments.out')
        row = next(r for r in self.read_fragment_report() if r['SMILES'] == 'NC')
        self.assertIn('[1NH2]', row['EXAMPLE_SMILES'])
        self.assertIn('[1CH2]', row['EXAMPLE_SMILES'])
        self.assertNotIn('[1OH]', row['EXAMPLE_SMILES'])

    def test_minimum_size_filter_still_validates_discarded_fragments(self):
        for fragment in ['[99CH4]', '[10000CH4]', '']:
            with self.subTest(fragment=fragment):
                self.run_tool(record('invalid', 'C[1NH2]', fragment), '-m', '2', success=False)
        self.run_tool(record('valid', 'C[1NH2]', '[100CH4]'), '-m', '2',
                      '-R', 'fragments.out')
        self.assertEqual([r['SMILES'] for r in self.read_fragment_report()], ['NC'])

    def test_fragment_report_separators_and_quoted_fields(self):
        first_name = 'zzz,| "quoted"\tparent'
        data = (record(first_name, 'C[1CH3]', '[100CH4]') +
                record('a_second', 'C[1NH2]', '[100CH4]'))
        baseline = self.run_tool(data).stdout
        for name, separator in (('tab', '\t'), ('comma', ','), ('vbar', '|'),
                                ('space', ' '), ('semic', ';')):
            result = self.run_tool(data, '-R', f'sep={name}', '-R', 'fragments.out')
            self.assertEqual(result.stdout, baseline)
            rows = self.read_fragment_report(separator)
            self.assertTrue(all(len(row) == 11 and None not in row for row in rows))
            self.assertEqual(rows[0]['FIRST_PARENT'], first_name)
            self.assertEqual([row['SMILES'] for row in rows], ['C', 'CC', 'NC'])

    def test_fragment_report_when_no_candidate_qualifies(self):
        data = (INDIR / 'basic.textproto').read_text()
        self.run_tool(data, '-s', '3', '-R', 'fragments.out', success=False)
        rows = self.read_fragment_report()
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row['SELECTED'] == '0' for row in rows))
        self.run_tool(data, '-m', '3', '-R', 'fragments.out', success=False)
        self.assertEqual(self.read_fragment_report(), [])
        self.assertEqual(len((self.root / 'fragments.out').read_text().splitlines()), 1)

    def test_fragment_report_flushes_multiple_buffers(self):
        fragments = ['[100CH4]']
        for length in range(2, 99):
            fragments.append('C' + ''.join(f'[{i}CH2]' for i in range(1, length - 1)) +
                             f'[{100 + length - 1}CH3]')
        self.run_tool(record('long', numbered_chain(99), *fragments), '-R', 'fragments.out')
        rows = self.read_fragment_report()
        self.assertGreater((self.root / 'fragments.out').stat().st_size, 4096)
        self.assertEqual([int(row['NATOMS']) for row in rows], list(range(99, 0, -1)))
        self.assertTrue(all(row['NPARENTS'] == '1' for row in rows))
        self.assertEqual(sum(int(row['NWHOLE']) for row in rows), 1)
        self.assertEqual(sum(int(row['NCUT']) for row in rows), 98)
        self.assertEqual(rows[0]['SELECTED'], '1')

    def test_fragment_report_options_and_path_errors(self):
        data = (INDIR / 'basic.textproto').read_text()
        for options in (('-R', 'sep=tab'), ('-R', 'a', '-R', 'b'),
                        ('-R', 'fragments.out', '-R', 'sep=unknown'),
                        ('-R', 'fragments.out', '-R', 'sep='),
                        ('-R', 'fragments.out', '-R', 'sep=dquote'),
                        ('-R', 'fragments.out', '-R', 'sep=\n'), ('-R', '-')):
            self.run_tool(data, *options, success=False)
        result = self.run_tool(data, '-R', 'missing_directory/report.out', success=False)
        self.assertIn('Cannot open fragment report', result.stderr)
        self.run_tool(data, '-R', 'input.textproto', success=False)
        self.assertEqual((self.root / 'input.textproto').read_text(), data)
        os.link(self.root / 'input.textproto', self.root / 'report_alias')
        self.run_tool(data, '-R', 'report_alias', success=False)
        self.assertEqual((self.root / 'input.textproto').read_text(), data)
        self.run_tool(data, '-R', './missing.smi', '-X', 'missing', success=False)
        self.run_tool(data, '-R', './pairs.out', '-Y', 'pairs.out', success=False)
        self.assertFalse((self.root / 'missing.smi').exists())
        self.assertFalse((self.root / 'pairs.out').exists())

    def test_nonmatching_filename_suffix(self):
        data = record('hit', 'C[1CH2][2OH]') + record('miss', 'N[1NH2]')
        for name, expected in (('missing', 'missing.smi'),
                               ('missing.smi', 'missing.smi'),
                               ('missing.other', 'missing.other.smi')):
            self.run_tool(data, '-s', '1', '-X', name)
            self.assertEqual((self.root / expected).read_text(), 'NN miss\n')
        self.assertFalse((self.root / 'missing').exists())
        self.assertFalse((self.root / 'missing.smi.smi').exists())

    def test_nonmatching_output_empty_when_all_parents_contribute(self):
        self.run_tool((INDIR / 'basic.textproto').read_text(), '-X', 'missing.smi')
        self.assertEqual((self.root / 'missing.smi').read_text(), '')

    def test_nonmatching_output_contains_all_parents_when_no_core_selected(self):
        for data, options in ((record('none', 'C[1CH3]'), ('-m', '3')),
                              (record('filtered', 'C[1CH3]', '[100CH4]'), ('-m', '3')),
                              ((INDIR / 'basic.textproto').read_text(), ('-s', '3'))):
            result = self.run_tool(data, *options, '-X', 'missing.smi', success=False)
            self.assertEqual(result.stdout, '')
            lines = (self.root / 'missing.smi').read_text().splitlines()
            self.assertEqual(len(lines), data.count('name:'))
            self.assertTrue(all('[' not in line.split()[0] for line in lines))

    def test_nonmatching_output_flushes_multiple_buffers(self):
        data = record('hit', 'C[1CH2][2CH3]', 'C[1CH2][2CH3]')
        names = [f'missing_parent_{i:04d}' for i in range(300)]
        data += ''.join(record(name, 'N[1NH2]', '[100NH3]') for name in names)
        self.run_tool(data, '-s', '1', '-X', 'missing.smi')
        self.assertEqual((self.root / 'missing.smi').read_text(),
                         ''.join(f'NN {name}\n' for name in names))

    def test_nonmatching_output_errors_and_path_collisions(self):
        data = (INDIR / 'basic.textproto').read_text()
        result = self.run_tool(data, '-X', 'nonexistent/missing.smi', success=False)
        self.assertIn('Cannot open non-matching output', result.stderr)
        # Check collisions using the actual filename after suffix addition.
        (self.root / 'input.textproto.smi').symlink_to('input.textproto')
        for path in ('input.textproto', './input.textproto'):
            self.run_tool(data, '-X', path, success=False)
            self.assertEqual((self.root / 'input.textproto').read_text(), data)
        os.link(self.root / 'input.textproto', self.root / 'hardlink.smi')
        (self.root / 'symlink.smi').symlink_to('input.textproto')
        for path in ('hardlink', 'symlink'):
            self.run_tool(data, '-X', path, success=False)
            self.assertEqual((self.root / 'input.textproto').read_text(), data)
        result = self.run_tool(data, '-X', 'missing', '-Y', './missing.smi', success=False)
        self.assertIn('different output files', result.stderr)
        self.assertFalse((self.root / 'missing.smi').exists())
        self.run_tool(data, '-X', 'missing.smi', '-Y', 'input.textproto', success=False)
        self.assertEqual((self.root / 'input.textproto').read_text(), data)
        self.assertFalse((self.root / 'missing.smi').exists())
        self.run_tool(data, '-X', '-', success=False)

    def support_example(self):
        return (record('one', 'C[1CH2][2CH2][3NH2]', '[100CH4]',
                       'C[101CH3]', 'C[1CH2][102CH3]') +
                record('two', 'C[1CH2][2OH]', '[100CH4]', 'C[101CH3]') +
                record('three', 'C[1CH2][2NH2]', '[100CH4]'))

    def test_minimum_support_selects_largest_qualifying_fragment(self):
        data = self.support_example()
        default = self.run_tool(data)
        self.assertIn('Most common fragment C found 3 times', default.stderr)
        for threshold, atoms, parents in (('1', 4, ['one']),
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
                record('two', 'C[1CH2][2NH][3CH3]', 'C[101CH3]', '[101CH3][102NH2]') +
                record('three', 'C[1CH2][2NH][3OH]', '[101CH3][102NH2]'))
        result = self.run_tool(data, '-s', '2')
        self.assertIn('Selected fragment NC found 3 times', result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 3)

    def test_percentage_support_uses_all_parents_and_exact_boundaries(self):
        # 14% of 100 is exactly 14, not 15 due to floating point rounding.
        # Spread implicit parent identities over eight structures so that
        # each whole occurs fewer than 14 times and cannot win this test.
        termini = ('NH2', 'OH', 'F', 'Cl', 'Br', 'I', 'SH', 'PH2')
        data = ''.join(record(str(i), f'C[1CH2][2{termini[i % 8]}]',
                              *(['[100CH4]', 'C[101CH3]'] if i < 14 else ['[100CH4]']))
                       for i in range(100))
        result = self.run_tool(data, '-s', '14%')
        self.assertEqual(len(result.stdout.splitlines()), 14)
        self.assertIn('Selected fragment CC found 14 times', result.stderr)
        data = record('one', 'C[1CH3]', '[100CH4]') + record('empty', 'C[1NH2]')
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
