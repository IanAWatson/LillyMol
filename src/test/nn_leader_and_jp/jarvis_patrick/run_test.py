#!/usr/bin/env python3
"""Check JP components and ordering against an independent graph traversal."""
from collections import Counter, deque
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tfrecord_utils import tfrecord

EXE = str(Path(sys.argv[1]).resolve())


def components(lists, spec, retained):
    length, common = spec.split(',')
    limit = None if length == '.' else int(length)
    neighbours = [set([j for j in row if j < retained][:limit]) for row in lists[:retained]]
    # Model JP as edges followed by breadth-first traversal, independently of
    # the tool's union-find implementation or its internal representatives.
    edges = [[] for _ in range(retained)]
    for i in range(retained):
        for j in range(i + 1, retained):
            if len(neighbours[i] & neighbours[j]) >= int(common):
                edges[i].append(j)
                edges[j].append(i)
    seen = set()
    result = []
    for i in range(retained):
        if i in seen:
            continue
        seen.add(i)
        queue = deque([i])
        group = []
        while queue:
            item = queue.popleft()
            group.append(item)
            for j in edges[item]:
                if j not in seen:
                    seen.add(j)
                    queue.append(j)
        result.append(sorted(group))
    return result


class JarvisPatrick(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def write(self, lists, serialized):
        # Reverse lexical/numeric ID order to distinguish it from input order.
        names = [f'm{len(lists) - i}' for i in range(len(lists))]
        path = self.directory / 'input'
        if serialized:
            path.write_bytes(b''.join(tfrecord(names[i], [(names[j], 0) for j in row])
                                      for i, row in enumerate(lists)))
        else:
            text = []
            for i, row in enumerate(lists):
                text.extend(['$SMI<C>', f'PCN<{names[i]}>'])
                for j in row:
                    text.extend([f'PCN<{names[j]}>', 'DIST<0>'])
                text.append('|')
            path.write_text('\n'.join(text) + '\n')
        return path, names

    def run_tool(self, path, args, serialized):
        result = subprocess.run([EXE, '-v', *args, *(['-r'] if serialized else []), str(path)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def assert_output(self, text, expected, smiles=False):
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
                self.assertEqual(int(re.search(r'CLUSTER<(\d+)>', record).group(1)), len(actual))
                names = re.findall(r'^PCN<([^>]+)>', record, re.MULTILINE)
                self.assertEqual(int(re.search(r'CSIZE<(\d+)>', record).group(1)), len(names))
                actual.append(names)
            self.assertNotIn('DIST<', text)
        self.assertEqual(actual, expected)

    def check(self, lists, spec, serialized, smiles=False, retained=None):
        path, names = self.write(lists, serialized)
        retained = len(lists) if retained is None else retained
        groups = components(lists, spec, retained)
        expected = [[names[i] for i in group] for group in groups]
        args = ['-J', spec]
        if smiles:
            args.append('-a')
        if retained != len(lists):
            args.extend(['-s', str(retained)])
        result = self.run_tool(path, args, serialized)
        self.assert_output(result.stdout, expected, smiles)
        self.assertIn(f'Clustered {retained} items into {len(groups)} clusters', result.stderr)
        histogram = {int(size): int(count) for count, size in re.findall(
            r'^(\d+) clusters of size (\d+)$', result.stderr, re.MULTILINE)}
        self.assertEqual(histogram, dict(Counter(map(len, groups))))
        self.assertIn(f'{retained} items in clusters', result.stderr)

    def test_singletons_and_repeated_merges(self):
        for serialized in (False, True):
            for smiles in (False, True):
                for lists in ([[]], [[] for _ in range(12)], [[19] for _ in range(20)]):
                    self.check(lists, '.,1', serialized, smiles)

    def test_late_bridge_and_nonminimal_root(self):
        # The two groups grow separately, then merge. Union-by-size can choose
        # a representative other than the first molecule of the final group.
        lists = [[0], [1, 2], [1], [2, 0], [], [3]]
        self.assertEqual(components(lists, '.,1', 6), [[0, 1, 2, 3], [4], [5]])
        for serialized in (False, True):
            for smiles in (False, True):
                self.check(lists, '.,1', serialized, smiles)

    def test_interleaved_components_and_list_limits(self):
        cases = [[[6], [7], [6], [7], [6], [7], [], []],
                 [[2, 3], [3, 2], [0], [0], [], []]]
        for serialized in (False, True):
            for lists in cases:
                for spec in ('.,1', '.,2', '1,1', '2,1', '2,2'):
                    self.check(lists, spec, serialized)

    def test_random_directed_lists(self):
        for seed in range(6):
            rng = random.Random(seed)
            n = 36
            lists = [rng.sample(range(n), rng.randrange(11)) for _ in range(n)]
            for serialized in (False, True):
                for spec in ('.,1', '.,3', '1,1', '2,2', '4,3'):
                    self.check(lists, spec, serialized, smiles=bool(seed % 2))

    def test_subset(self):
        lists = [[4, 1, 3], [4, 0, 3], [5, 3], [5, 2], [], []]
        for serialized in (False, True):
            for spec in ('.,1', '2,1', '.,2'):
                self.check(lists, spec, serialized, retained=4)

    def test_multiple_parameter_files(self):
        lists = [[2, 3], [3, 2], [0], [0], [], []]
        for serialized in (False, True):
            path, names = self.write(lists, serialized)
            specs = ('.,1', '.,2', '1,1', '.,1')
            stem = self.directory / 'jp'
            args = [value for spec in specs for value in ('-J', spec)]
            result = self.run_tool(path, [*args, '-B', str(stem)], serialized)
            self.assertEqual(result.stdout, '')
            for i, spec in enumerate(specs):
                expected = [[names[j] for j in group] for group in components(lists, spec, len(lists))]
                self.assert_output(Path(f'{stem}{i}.jpc').read_text(), expected)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(JarvisPatrick))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp Jarvis Patrick checks passed')
