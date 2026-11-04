#!/usr/bin/env python3
"""Executable regression tests using explicit synthetic SDF coordinates."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
SQUARE = [(0, 1), (1, 2), (2, 3), (3, 0)]
CLEAN = [(1, 1), (3, 1), (3, 3), (1, 3)]
BOWTIE = [(1, 1), (3, 3), (1, 3), (3, 1)]


def record(name, points, edges, z=0):
    text = [name, '  LillyMol regression 2D', '',
            f'{len(points):3d}{len(edges):3d}  0  0  0  0            999 V2000']
    text.extend(f'{x:10.4f}{y:10.4f}{z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0'
                for x, y in points)
    text.extend(f'{a + 1:3d}{b + 1:3d}  1  0  0  0  0' for a, b in edges)
    return '\n'.join([*text, 'M  END', '$$$$', ''])


class Crossings(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def run_tool(self, sdf, *args, success=True, stdin=False):
        source = self.directory / 'input.sdf'
        source.write_text(sdf)
        result = subprocess.run([EXE, *map(str, args), '-' if stdin else str(source)],
                                input=sdf if stdin else None, capture_output=True,
                                text=True, timeout=30)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def rows(self, text):
        lines = [line.split() for line in text.splitlines()]
        self.assertTrue(lines)
        self.assertEqual(lines[0][0], 'Name')
        rows = {}
        for line in lines[1:]:
            self.assertEqual(len(line), len(lines[0]))
            for value in line[1:]:
                self.assertRegex(value, r'^\d+(?:\.\d{1,4})?$')
            rows[line[0]] = dict(zip(lines[0][1:], map(float, line[1:])))
        return rows

    def test_counts_distances_and_stdin(self):
        sdf = record('clean', CLEAN, SQUARE) + record('crossed', BOWTIE, SQUARE)
        result = self.run_tool(sdf)
        rows = self.rows(result.stdout)
        self.assertEqual(rows['clean']['bc_score'], 0)
        crossed = rows['crossed']
        self.assertEqual(crossed['bc_crossings'], 1)
        self.assertEqual(crossed['bc_score'], 1)
        self.assertEqual(crossed['bc_score_per_atom'], 0.25)
        self.assertEqual(crossed['bc_min_2d'], 2)
        self.assertEqual(crossed['bc_max_2d'], 2.8284)
        self.assertEqual(crossed['bc_topological_diameter'], 2)
        self.assertEqual(crossed['bc_max_2d_per_topology'], 1.4142)
        self.assertEqual(self.run_tool(sdf, stdin=True).stdout, result.stdout)

    def test_local_congestion_and_spiro_separation(self):
        points = [(0, 0), (4, 0), (4, 4), (0, 4), (2, 6), (1, 6), (3, 6)]
        sdf = ''
        for count in (1, 2, 3):
            edges = SQUARE + [(0, 4)]
            if count >= 2: edges.append((1, 5))
            if count >= 3: edges.append((0, 6))
            sdf += record(f'local{count}', points[:4 + count], edges)
        sdf += record('spiro', BOWTIE + [(-1, -1), (1, -1), (-1, 1)],
                      SQUARE + [(0, 4), (4, 5), (5, 6), (6, 0)])
        rows = self.rows(self.run_tool(sdf).stdout)
        for count in (1, 2, 3):
            self.assertEqual(rows[f'local{count}']['bc_ring_chain'], count)
            self.assertEqual(rows[f'local{count}']['bc_score'], count ** 3)
        self.assertEqual(rows['spiro']['bc_score'], 2)
        self.assertEqual(rows['spiro']['bc_nrs'], 2)
        self.assertEqual(rows['spiro']['bc_nratoms'], 7)

    def test_chain_chain_crossings_are_not_scored(self):
        rows = self.rows(self.run_tool(record('chain', BOWTIE, SQUARE[:3]), '-d').stdout)
        self.assertEqual(rows['chain']['bc_crossings'], 0)
        self.assertEqual(rows['chain']['bc_score'], 0)
        self.assertEqual(rows['chain']['bc_near_misses'], 0)

    def test_dimensionality_policy_and_statistics(self):
        sdf = record('three_d', CLEAN, SQUARE, z=1)
        sdf += record('flat', [(0, 0), (1, 0)], [(0, 1)])
        sdf += record('clean', CLEAN, SQUARE)
        failure = self.run_tool(sdf, success=False)
        self.assertIn('dimensionality', failure.stderr)
        skipped = self.run_tool(sdf, '-z', 'i')
        self.assertEqual(set(self.rows(skipped.stdout)), {'clean'})
        self.assertIn('skipped 2', skipped.stderr)

    def test_disconnected_requires_largest_fragment(self):
        sdf = record('salt', CLEAN + [(100, 100)], SQUARE)
        failure = self.run_tool(sdf, success=False)
        self.assertIn('use -l', failure.stderr)
        self.assertIn('use -l', self.run_tool(sdf, '-z', 'i', success=False).stderr)
        rows = self.rows(self.run_tool(sdf, '-l').stdout)
        self.assertEqual(rows['salt']['bc_natoms'], 4)
        self.assertEqual(rows['salt']['bc_max_2d'], 2.8284)

    def test_optional_diagnostics(self):
        points = [(0, 0), (4, 0), (4, 4), (0, 4)]
        sdf = record('near', points + [(-1, -0.01), (5, -0.01)], SQUARE + [(0, 4), (4, 5)])
        sdf += record('overlap', points + [(-1, 0), (5, 0)], SQUARE + [(0, 4), (4, 5)])
        sdf += record('degenerate', points + [(0, 0)], SQUARE + [(0, 4)])
        basic = self.run_tool(sdf)
        self.assertNotIn('bc_overlaps', basic.stdout)
        detailed = self.rows(self.run_tool(sdf, '-d').stdout)
        self.assertGreater(detailed['near']['bc_near_misses'], 0)
        self.assertGreater(detailed['overlap']['bc_overlaps'], 0)
        self.assertEqual(detailed['degenerate']['bc_degenerate_bonds'], 1)
        self.assertEqual(detailed['degenerate']['bc_coincident_atoms'], 1)
        for row in detailed.values(): self.assertEqual(row['bc_score'], 0)
        disabled_near = self.rows(self.run_tool(sdf, '-d', '-t', 0).stdout)
        self.assertEqual(disabled_near['near']['bc_near_misses'], 0)

    def test_bad_options_names_and_missing_files(self):
        sdf = record('name with spaces', CLEAN, SQUARE)
        self.assertEqual(set(self.rows(self.run_tool(sdf).stdout)), {'name'})
        for args in (['-t', '0.1'], ['-d', '-t', '-1'], ['-z', 'invalid']):
            self.run_tool(sdf, *args, success=False)
        result = subprocess.run([EXE, str(self.directory / 'missing.sdf')],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Cannot open', result.stderr)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Crossings))
    if not result.wasSuccessful(): sys.exit(1)
    print('bond_crossings checks passed')
