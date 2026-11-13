#!/usr/bin/env python3
"""Check global pool limits for TDT and TFDataRecord neighbour inputs."""
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




class PoolLimit(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def write(self, filename, names, serialized):
        path = self.directory / filename
        if serialized:
            path.write_bytes(b''.join(tfrecord(name) for name in names))
        else:
            path.write_text(''.join(f'$SMI<C>\nPCN<{name}>\n|\n' for name in names))
        return path

    def run_tool(self, files, mode, serialized, limit=None, success=True):
        args = [EXE, *mode]
        if serialized:
            args.append('-r')
        if limit is not None:
            args.extend(['-s', str(limit)])
        result = subprocess.run([*args, *map(str, files)], capture_output=True,
                                text=True, timeout=15)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
            # An oversized -s is not an error or a warning.
            self.assertEqual(''.join(line for line in result.stderr.splitlines(True)
                                        if not line.startswith("Warning: distance coverage in '")), '')
        else:
            self.assertGreater(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, '')
        return result

    def ids(self, result):
        return re.findall(r'^PCN<([^>]+)>', result.stdout, re.MULTILINE)

    def test_single_file_without_limit(self):
        for serialized in (False, True):
            path = self.write('all', ['m0', 'm1', 'm2'], serialized)
            for mode in MODES:
                self.assertEqual(self.ids(self.run_tool([path], mode, serialized)),
                                 ['m0', 'm1', 'm2'])

    def test_multiple_files_require_pool_size(self):
        for serialized in (False, True):
            first = self.write('first', ['m0'], serialized)
            second = self.write('second', ['m1'], serialized)
            for mode in MODES:
                result = self.run_tool([first, second], mode, serialized, success=False)
                self.assertIn('Multiple input files require', result.stderr)
                self.assertIn('-s', result.stderr)

    def test_stop_within_one_file(self):
        for serialized in (False, True):
            path = self.write('all', [f'm{i}' for i in range(5)], serialized)
            for mode in MODES:
                self.assertEqual(self.ids(self.run_tool([path], mode, serialized, 2)),
                                 ['m0', 'm1'])

    def test_stop_before_opening_another_file(self):
        for serialized in (False, True):
            first = self.write('first', ['m0', 'm1'], serialized)
            missing = self.directory / 'must_not_be_opened'
            for mode in MODES:
                self.assertEqual(self.ids(self.run_tool([first, missing], mode, serialized, 2)),
                                 ['m0', 'm1'])

    def test_stop_midway_through_later_file(self):
        for serialized in (False, True):
            first = self.write('first', ['m0', 'm1'], serialized)
            second = self.write('second', ['m2', 'm3'], serialized)
            missing = self.directory / 'must_not_be_opened'
            for mode in MODES:
                result = self.run_tool([first, second, missing], mode, serialized, 3)
                self.assertEqual(self.ids(result), ['m0', 'm1', 'm2'])

    def test_oversized_limit(self):
        for serialized in (False, True):
            first = self.write('first', ['m0', 'm1'], serialized)
            second = self.write('second', ['m2', 'm3'], serialized)
            for mode in MODES:
                for files, expected in (([first], ['m0', 'm1']),
                                        ([first, second], ['m0', 'm1', 'm2', 'm3'])):
                    self.assertEqual(self.ids(self.run_tool(files, mode, serialized, 999)),
                                     expected)

    def write_graph(self, filename, records, serialized):
        path = self.directory / filename
        if serialized:
            path.write_bytes(b''.join(tfrecord(name, nbrs) for name, nbrs in records))
        else:
            text = []
            for name, nbrs in records:
                text.extend(['$SMI<C>', f'PCN<{name}>'])
                for nbr, distance in nbrs:
                    text.extend([f'PCN<{nbr}>', f'DIST<{distance}>'])
                text.append('|')
            path.write_text('\n'.join(text) + '\n')
        return path

    def test_neighbours_outside_subset(self):
        # Outside neighbours precede retained neighbours, so indices and
        # distances must be compacted together without changing their order.
        records = [('m0', [('m2', 0.05), ('m1', 0.15)]),
                   ('m1', [('m2', 0.05), ('m0', 0.15)]), ('m2', [])]
        for serialized in (False, True):
            path = self.write_graph('graph', records, serialized)
            for mode in MODES:
                result = self.run_tool([path], mode, serialized, 1)
                self.assertEqual(self.ids(result), ['m0'])
            for mode in (['-h'], ['-T', '-h']):
                for threshold, sizes in (('0.1', ['1', '1']), ('0.2', ['2'])):
                    result = self.run_tool([path], [*mode, '-t', threshold], serialized, 2)
                    self.assertEqual(self.ids(result), ['m0', 'm1'])
                    self.assertEqual(re.findall(r'^CSIZE<(\d+)>', result.stdout, re.MULTILINE),
                                     sizes)
                    distances = re.findall(r'^DIST<([^>]+)>', result.stdout, re.MULTILINE)
                    self.assertEqual(distances, [] if threshold == '0.1' else ['0.15'])
            # In JP, the omitted common neighbour no longer links m0 and m1.
            result = self.run_tool([path], ['-J', '.,1'], serialized, 2)
            self.assertEqual(re.findall(r'^CSIZE<(\d+)>', result.stdout, re.MULTILINE),
                             ['1', '1'])

    def test_neighbours_across_files(self):
        for serialized in (False, True):
            first = self.write_graph('first', [('m0', [('m2', 0.05), ('m1', 0.15)])], serialized)
            second = self.write_graph('second', [('m1', [('m2', 0.05), ('m0', 0.15)]),
                                                 ('m2', [])], serialized)
            for mode in MODES:
                result = self.run_tool([first, second], mode, serialized, 2)
                self.assertEqual(self.ids(result), ['m0', 'm1'])

    def test_unknown_neighbours_still_rejected_without_truncation(self):
        for serialized in (False, True):
            path = self.write_graph('graph', [('m0', [('missing', 0.1)])], serialized)
            for mode in MODES:
                for limit in (None, 999):
                    result = self.run_tool([path], mode, serialized, limit, success=False)
                    self.assertIn('missing', result.stderr)

    def test_unread_tail_is_not_parsed(self):
        for serialized in (False, True):
            # The duplicate ID after the cutoff would otherwise prevent loading.
            path = self.write('all', ['m0', 'm1', 'm0'], serialized)
            for mode in MODES:
                self.assertEqual(self.ids(self.run_tool([path], mode, serialized, 2)),
                                 ['m0', 'm1'])


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(PoolLimit))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp pool limit checks passed')
