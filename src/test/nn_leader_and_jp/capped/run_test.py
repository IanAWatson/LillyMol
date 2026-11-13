#!/usr/bin/env python3
"""Check capped cluster membership, headers and verbose size statistics."""
from pathlib import Path
from collections import Counter
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())


class Capped(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input = Path(self.tmp.name) / 'input.nn'

    def graph(self, limits, self_neighbours=False, distances=None):
        records = []
        n = len(limits)
        for i, limit in enumerate(limits):
            lines = ['$SMI<C>', f'PCN<m{i} {limit}>']
            neighbours = []
            for j in range(n):
                if i == j and not self_neighbours:
                    continue
                d = distances[i][j] if distances else (0 if i == j else 0.1)
                neighbours.append((d, j))
            for d, j in sorted(neighbours):
                lines.extend([f'PCN<m{j}>', f'DIST<{d}>'])
            records.append('\n'.join([*lines, '|', '']))
        self.input.write_text(''.join(records))

    def run_tool(self, *args):
        result = subprocess.run([EXE, '-v', *map(str, args), str(self.input)],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def clusters(self, result):
        clusters = []
        for record in result.stdout.split('|'):
            if not record.strip():
                continue
            ids = re.findall(r'^PCN<(m\d+)', record, re.MULTILINE)
            size = int(re.search(r'^CSIZE<(\d+)>', record, re.MULTILINE).group(1))
            self.assertEqual(size, len(ids), record)
            clusters.append(ids)
        # These statistics use the count returned by cluster formation.
        reported = {int(size): int(count) for count, size in re.findall(
            r'^(\d+) clusters of size (\d+) members$', result.stderr, re.MULTILINE)}
        self.assertEqual(reported, dict(Counter(map(len, clusters))))
        return clusters

    def test_global_limits(self):
        for mode in ([], ['-T']):
            for self_neighbours in (False, True):
                for limit in (2, 3, 4, 8):
                    with self.subTest(mode=mode, self_neighbours=self_neighbours, limit=limit):
                        self.graph([limit] * 5, self_neighbours)
                        result = self.run_tool(*mode, '-h', '-t', '0.2', '-m', limit)
                        clusters = self.clusters(result)
                        self.assertEqual(list(map(len, clusters)),
                                         [min(limit, 5 - i) for i in range(0, 5, limit)])
                        self.assertEqual(sorted(m for c in clusters for m in c),
                                         [f'm{i}' for i in range(5)])
                        self.assertNotIn('In clusters ', result.stderr)

    def test_column_limits_include_one(self):
        for mode in ([], ['-T']):
            for self_neighbours in (False, True):
                for limits in ([1] * 5, [1, 2, 3, 1, 2], [3] * 5):
                    with self.subTest(mode=mode, self_neighbours=self_neighbours, limits=limits):
                        self.graph(limits, self_neighbours)
                        result = self.run_tool(*mode, '-h', '-t', '0.2', '-M', 'col=2')
                        clusters = self.clusters(result)
                        for cluster in clusters:
                            leader = int(cluster[0][1:])
                            self.assertLessEqual(len(cluster), limits[leader])
                        self.assertEqual(sorted(m for c in clusters for m in c),
                                         [f'm{i}' for i in range(5)])
                        self.assertNotIn('In clusters ', result.stderr)

    def test_per_molecule_limit_overrides_global(self):
        for mode in ([], ['-T']):
            self.graph([1] * 3, True)
            clusters = self.clusters(self.run_tool(*mode, '-h', '-t', '0.2',
                                                  '-m', '3', '-M', 'col=2'))
            self.assertEqual(list(map(len, clusters)), [1, 1, 1])

    def test_self_neighbour_output_matches_uncapped(self):
        # Without -h, retain the established convention of writing a listed
        # self neighbour. Large caps should agree with the uncapped path.
        for mode in ([], ['-T']):
            self.graph([8] * 3, True)
            uncapped = self.run_tool(*mode, '-t', '0.2')
            capped = self.run_tool(*mode, '-t', '0.2', '-m', '8')
            self.assertEqual(capped.stdout, uncapped.stdout)
            self.assertEqual(len(self.clusters(capped)[0]), 4)
            for limit in (2, 3):
                clusters = self.clusters(self.run_tool(*mode, '-t', '0.2', '-m', limit))
                self.assertTrue(all(len(c) <= limit for c in clusters))
                self.assertEqual(set(m for c in clusters for m in c), {'m0', 'm1', 'm2'})

    def test_threshold_inclusive_and_isolated_molecule(self):
        distances = [[0, 0.2, 0.8], [0.2, 0, 0.8], [0.8, 0.8, 0]]
        for mode in ([], ['-T']):
            for self_neighbours in (False, True):
                self.graph([8] * 3, self_neighbours, distances)
                clusters = self.clusters(self.run_tool(*mode, '-h', '-t', '0.2', '-m', '8'))
                self.assertEqual(clusters, [['m0', 'm1'], ['m2']])
                clusters = self.clusters(self.run_tool(*mode, '-h', '-t', '0.1', '-m', '8'))
                self.assertEqual(clusters, [['m0'], ['m1'], ['m2']])

    def test_smiles_headers(self):
        for mode in ([], ['-T']):
            self.graph([3] * 5, True)
            result = self.run_tool(*mode, '-h', '-a', '-t', '0.2', '-m', '3')
            sizes, actual = [], []
            for line in result.stdout.splitlines():
                if ' CLUSTER ' in line:
                    sizes.append(int(re.search(r'\((\d+) members\)', line).group(1)))
                    actual.append(1)
                else:
                    actual[-1] += 1
            self.assertEqual(sizes, [3, 2])
            self.assertEqual(sizes, actual)
            self.assertNotIn('In clusters ', result.stderr)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Capped))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp capped cluster checks passed')
