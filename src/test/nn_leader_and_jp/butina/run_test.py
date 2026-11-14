#!/usr/bin/env python3
"""Compare dynamic Butina clustering with an independent full recount."""
from collections import Counter
from pathlib import Path
import random
import re
import struct
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tfrecord_utils import tfrecord

EXE = str(Path(sys.argv[1]).resolve())


def f32(value):
    return struct.unpack('<f', struct.pack('<f', value))[0]


def recount(rows, threshold, caps, suppress_self=True, global_cap=0,
            column=False, retained=None, max_clusters=None, static_order=False):
    n = len(rows) if retained is None else retained
    threshold = f32(threshold)
    neighbours = [sorted([(j, f32(d)) for j, d in row if j < n and
                           (not suppress_self or j != i)], key=lambda e: e[1])
                  for i, row in enumerate(rows[:n])]
    initial_counts = [sum(d <= threshold for _, d in row) for row in neighbours]
    selected = set()
    clusters = []
    while len(selected) < n:
        # Deliberately recompute every count after each cluster. No heap or
        # reverse adjacency is used by this reference implementation.
        centre = min((i for i in range(n) if i not in selected),
                     key=lambda i: (-(initial_counts[i] if static_order else
                                      sum(j not in selected and d <= threshold
                                          for j, d in neighbours[i])), i))
        members = [centre]
        cap = caps[centre] if column else global_cap
        for j, d in neighbours[centre]:
            if d > threshold:
                break
            if j in selected:
                continue
            if cap and len(members) >= cap:
                break
            members.append(j)
            selected.add(j)
        selected.add(centre)
        clusters.append(members)
        if max_clusters and len(clusters) >= max_clusters:
            break
    return clusters


class Butina(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def write(self, rows, caps, serialized):
        ids = [f'm{len(rows) - i}' for i in range(len(rows))]
        path = self.directory / 'input'
        if serialized:
            path.write_bytes(b''.join(tfrecord(f'{ids[i]} {caps[i]}',
                                      [(ids[j], d) for j, d in sorted(row, key=lambda e: e[1])])
                                      for i, row in enumerate(rows)))
        else:
            lines = []
            for i, row in enumerate(rows):
                lines.extend(['$SMI<C>', f'PCN<{ids[i]} {caps[i]}>'])
                for j, d in sorted(row, key=lambda e: e[1]):
                    lines.extend([f'PCN<{ids[j]}>', f'DIST<{d}>'])
                lines.append('|')
            path.write_text('\n'.join(lines) + '\n')
        return path, ids

    def run_tool(self, path, args, serialized):
        result = subprocess.run([EXE, '-T', '-v', *args, *(['-r'] if serialized else []), str(path)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def assert_clusters(self, text, expected, smiles=False):
        if smiles:
            actual, sizes = [], []
            for line in text.splitlines():
                name = re.search(r'\b(m\d+)\b', line).group(1)
                if ' CLUSTER ' in line:
                    number, size = map(int, re.search(r'CLUSTER (\d+) \((\d+) members\)', line).groups())
                    self.assertEqual(number, len(actual))
                    actual.append([name])
                    sizes.append(size)
                else:
                    actual[-1].append(name)
            self.assertEqual(sizes, list(map(len, actual)))
        else:
            actual = []
            for record in text.split('|'):
                if not record.strip():
                    continue
                names = re.findall(r'^PCN<(m\d+)', record, re.MULTILINE)
                self.assertEqual(int(re.search(r'CLUSTER<(\d+)>', record).group(1)), len(actual))
                self.assertEqual(int(re.search(r'CSIZE<(\d+)>', record).group(1)), len(names))
                actual.append(names)
        self.assertEqual(actual, expected)

    def check(self, rows, threshold, serialized, caps=None, suppress_self=True,
              global_cap=0, column=False, retained=None, max_clusters=None, smiles=False,
              static_order=False):
        caps = caps if caps is not None else [4] * len(rows)
        path, ids = self.write(rows, caps, serialized)
        groups = recount(rows, threshold, caps, suppress_self, global_cap, column,
                         retained, max_clusters, static_order)
        expected = [[ids[i] for i in group] for group in groups]
        args = ['-t', str(threshold)]
        if suppress_self:
            args.append('-h')
        if global_cap:
            args.extend(['-m', str(global_cap)])
        if column:
            args.extend(['-M', 'col=2'])
        if retained is not None:
            args.extend(['-s', str(retained)])
        if max_clusters:
            args.extend(['-C', str(max_clusters)])
        if smiles:
            args.append('-a')
        result = self.run_tool(path, args, serialized)
        self.assert_clusters(result.stdout, expected, smiles)
        histogram = {int(size): int(count) for count, size in re.findall(
            r'^(\d+) clusters of size (\d+) members$', result.stderr, re.MULTILINE)}
        self.assertEqual(histogram, dict(Counter(map(len, groups))))
        return groups

    def test_counts_change_after_cluster_and_directed_edges(self):
        rows = [[(1, .125), (2, .125)], [(0, .125), (2, .125)],
                [(0, .125), (1, .125)], [(0, .125), (1, .125)], [(5, .125)], []]
        self.assertEqual(recount(rows, .25, [4] * 6), [[0, 1, 2], [4, 5], [3]])
        for serialized in (False, True):
            self.check(rows, .25, serialized)
            self.check(rows, .25, serialized, global_cap=2)
            self.check(rows, .25, serialized, caps=[1, 2, 3, 4, 1, 2], column=True, global_cap=5)

    def test_empty_self_only_and_dense_lists(self):
        cases = [[[]], [[(0, 0)]], [[] for _ in range(12)],
                 [[(i, 0)] for i in range(12)],
                 [[(j, 0) for j in range(12)] for _ in range(12)]]
        for serialized in (False, True):
            for rows in cases:
                for suppress in (False, True):
                    self.check(rows, 0, serialized, suppress_self=suppress)
                    self.check(rows, 1, serialized, caps=[1] * len(rows), column=True,
                               suppress_self=suppress)

    def random_rows(self, seed, n):
        rng = random.Random(seed)
        rows = [[(j, rng.choice((0, .125, .25, .5, .75, 1)))
                 for j in rng.sample(range(n), rng.randrange(14))] for _ in range(n)]
        caps = [rng.randrange(1, 7) for _ in range(n)]
        return rows, caps

    def test_random_directed_lists_and_caps(self):
        for seed in range(6):
            rows, caps = self.random_rows(seed, 48)
            for serialized in (False, True):
                for threshold in (0, .25, .75, 1):
                    for suppress in (False, True):
                        for global_cap, column in ((0, False), (3, False), (5, True)):
                            with self.subTest(seed=seed, serialized=serialized, threshold=threshold,
                                              suppress=suppress, global_cap=global_cap, column=column):
                                self.check(rows, threshold, serialized, caps=caps,
                                           suppress_self=suppress, global_cap=global_cap,
                                           column=column, smiles=bool(seed % 2))

    def test_subset_and_early_stop(self):
        rows, caps = self.random_rows(71, 64)
        for serialized in (False, True):
            for threshold in (0, .5, 1):
                self.check(rows, threshold, serialized, caps=caps, column=True, retained=37)
                self.check(rows, threshold, serialized, caps=caps, column=True, max_clusters=7)

    def test_parallel_thresholds_with_shared_reverse_index(self):
        rows, caps = self.random_rows(92, 64)
        for serialized in (False, True):
            path, ids = self.write(rows, caps, serialized)
            stem = self.directory / 'thresholds'
            # The largest threshold is first, not last; smaller thresholds must
            # not count edges solely because they occur in the shared index.
            thresholds = (1, 0, .25, 0)
            args = ['-h', '-M', 'col=2', '-B', str(stem)]
            for threshold in thresholds:
                args.extend(['-t', str(threshold)])
            result = self.run_tool(path, args, serialized)
            self.assertEqual(result.stdout, '')
            for i, threshold in enumerate(thresholds):
                groups = recount(rows, threshold, caps, column=True)
                expected = [[ids[j] for j in group] for group in groups]
                self.assert_clusters(Path(f'{stem}{i}.ldr').read_text(), expected)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Butina))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp Butina checks passed')
