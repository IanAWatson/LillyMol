#!/usr/bin/env python3
"""Normalise neighbour lists once, before clustering starts."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tfrecord_utils import tfrecord

EXE = str(Path(sys.argv[1]).resolve())
MODES = (['-h', '-t', '0.2'], ['-T', '-h', '-t', '0.2'], ['-J', '.,1'])


class Ingestion(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / 'input'

    def write(self, records, serialized):
        if serialized:
            self.path.write_bytes(b''.join(tfrecord(name, nbrs) for name, nbrs in records))
        else:
            text = []
            for name, neighbours in records:
                text.extend(['$SMI<C>', f'PCN<{name}>'])
                for neighbour, distance in neighbours:
                    text.extend([f'PCN<{neighbour}>', f'DIST<{distance}>'])
                text.append('|')
            self.path.write_text('\n'.join(text) + '\n')

    def run_tool(self, args, serialized, success=True):
        result = subprocess.run([EXE, *args, *(['-r'] if serialized else []), str(self.path)],
                                capture_output=True, text=True, timeout=15)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertGreater(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, '')
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def pairs(self, self_ids=(), self_last=False):
        records = []
        for i in range(4):
            # Different IDs may have distance zero; remove self by index,
            # never by distance. Equal distances permit self anywhere in the list.
            neighbours = [(f'm{i ^ 1}', 0)]
            if i in self_ids:
                neighbours.insert(len(neighbours) if self_last else 0, (f'm{i}', 0))
            records.append((f'm{i} 2', neighbours))
        return records

    def test_mixed_self_neighbours_do_not_change_clustering(self):
        for serialized in (False, True):
            for mode in (['-h', '-t', '0'], ['-T', '-h', '-t', '0']):
                for cap in ([], ['-m', '2'], ['-M', 'col=2']):
                    args = [*mode, *cap]
                    self.write(self.pairs(), serialized)
                    expected = self.run_tool(args, serialized).stdout
                    for self_ids in ((1,), (1, 3), (0, 1, 2, 3)):
                        for self_last in (False, True):
                            with self.subTest(serialized=serialized, args=args,
                                              self_ids=self_ids, self_last=self_last):
                                self.write(self.pairs(self_ids, self_last), serialized)
                                result = self.run_tool(args, serialized)
                                self.assertEqual(result.stdout, expected)
                    self.assertEqual(re.findall(r'^PCN<(m\d+)', expected, re.MULTILINE),
                                     ['m0', 'm1', 'm2', 'm3'])

    def test_self_only_lists_become_singletons(self):
        for serialized in (False, True):
            self.write([(f'm{i}', [(f'm{i}', 0)]) for i in range(3)], serialized)
            for mode in MODES[:2]:
                result = self.run_tool(mode, serialized)
                self.assertEqual(re.findall(r'^CSIZE<(\d+)>', result.stdout, re.MULTILINE),
                                 ['1', '1', '1'])
                self.assertNotIn('DIST<', result.stdout)

    def test_self_pruning_combines_with_subset_filtering(self):
        records = [('m0', [('m2', 0), ('m0', 0), ('m1', 0)]),
                   ('m1', [('m0', 0), ('m2', 0), ('m1', 0)]), ('m2', [])]
        for serialized in (False, True):
            self.write(records, serialized)
            for mode in MODES[:2]:
                result = self.run_tool([*mode, '-s', '2'], serialized)
                self.assertEqual(re.findall(r'^PCN<(m\d+)', result.stdout, re.MULTILINE),
                                 ['m0', 'm1'])
                self.assertIn('CSIZE<2>', result.stdout)

    def test_h_remains_specific_to_leader_and_butina(self):
        records = [('m0', [('m0', 0), ('m1', 0)]),
                   ('m1', [('m0', 0), ('m1', 0)])]
        for serialized in (False, True):
            self.write(records, serialized)
            expected = self.run_tool(['-J', '2,2'], serialized).stdout
            self.assertIn('CSIZE<2>', expected)
            self.assertEqual(self.run_tool(['-J', '2,2', '-h'], serialized).stdout, expected)

    def test_nonfinite_distances_rejected_at_load_time(self):
        for serialized in (False, True):
            for distance in (float('nan'), float('inf'), float('-inf')):
                for neighbour in ('m0', 'm1'):
                    for mode in MODES:
                        with self.subTest(serialized=serialized, distance=distance,
                                          neighbour=neighbour, mode=mode):
                            self.write([('m0', [(neighbour, distance)]), ('m1', [])], serialized)
                            result = self.run_tool(mode, serialized, success=False)
                            self.assertIn('invalid distance', result.stderr)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Ingestion))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp ingestion checks passed')
