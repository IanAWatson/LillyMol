#!/usr/bin/env python3
"""Compare threaded threshold output with independent serial runs."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tfrecord_utils import tfrecord

EXE = str(Path(sys.argv[1]).resolve())
THRESHOLDS = ('0', '0.1234567', '0.5', '1', '0.1234567')


class ParallelThresholds(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def records(self, n):
        result = []
        for i in range(n):
            distances = [(f'm{i}', 0.0)] if i % 2 else []
            for j in range((i // 8) * 8, min(n, (i // 8 + 1) * 8)):
                if j != i:
                    d = (0.0, 0.05, 0.1234567, 0.5, 1.0)[(i + j) % 5]
                    distances.append((f'm{j}', d))
            # Directed links between blocks also exercise shared selections.
            if n > 32:
                distances.append((f'm{(i + 17) % n}', 0.75))
            distances.sort(key=lambda item: item[1])
            result.append((f'm{i} {1 if i % 5 == 0 else 4}', distances))
        return result

    def write(self, name, records, serialized):
        path = self.directory / name
        if serialized:
            path.write_bytes(b''.join(tfrecord(name, nbrs) for name, nbrs in records))
        else:
            lines = []
            for name, nbrs in records:
                lines.extend(['$SMI<C>', f'PCN<{name}>'])
                for neighbour, distance in nbrs:
                    lines.extend([f'PCN<{neighbour}>', f'DIST<{distance}>'])
                lines.append('|')
            path.write_text('\n'.join(lines) + '\n')
        return path

    def run_tool(self, files, args, serialized, success=True):
        result = subprocess.run([EXE, *args, *(['-r'] if serialized else []), *map(str, files)],
                                capture_output=True, text=True, timeout=30)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertGreater(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def compare(self, files, options, serialized):
        # No -B and a single threshold executes the serial path to stdout.
        expected = {t: self.run_tool(files, [*options, '-t', t], serialized).stdout
                    for t in set(THRESHOLDS)}
        for repeat in range(2):
            stem = self.directory / f'parallel{repeat}'
            if repeat == 0:
                thresholds = ['-t', ','.join(THRESHOLDS)]
            else:
                thresholds = [value for t in THRESHOLDS for value in ('-t', t)]
            result = self.run_tool(files, [*options, *thresholds, '-B', str(stem)], serialized)
            self.assertEqual(result.stdout, '')
            for i, t in enumerate(THRESHOLDS):
                self.assertEqual(Path(f'{stem}{i}.ldr').read_text(), expected[t],
                                 f'threshold={t}, index={i}, options={options}')
        return expected

    def test_outputs_match_serial_runs(self):
        # Large enough to exercise periodic output-buffer writes, with variable
        # neighbour counts, caps, thresholds, duplicate thresholds and self IDs.
        for serialized in (False, True):
            path = self.write('input', self.records(1024), serialized)
            for options in (['-h'], ['-T', '-h'], ['-h', '-m', '3'],
                            ['-T', '-h', '-M', 'col=2'], ['-h', '-C', '7'], []):
                with self.subTest(serialized=serialized, options=options):
                    self.compare([path], options, serialized)

    def test_subset_and_multiple_input_files(self):
        for serialized in (False, True):
            records = self.records(64)
            first = self.write('first', records[:32], serialized)
            second = self.write('second', records[32:], serialized)
            for mode in ([], ['-T']):
                self.compare([first, second], [*mode, '-h', '-s', '48'], serialized)

    def test_verbose_reports_and_files(self):
        path = self.write('input', self.records(32), False)
        self.compare([path], ['-T', '-h', '-v', '-v'], False)
        stem = self.directory / 'verbose'
        result = self.run_tool([path], ['-h', '-v', '-t', '0,0.5,1', '-B', str(stem)], False)
        self.assertEqual(len(re.findall(r'^Threshold [^:]+:$', result.stderr, re.MULTILINE)), 3)
        self.assertEqual(len(re.findall(r'^Clustered 32 fingerprints into \d+ clusters$',
                                        result.stderr, re.MULTILINE)), 3)

    def test_one_failed_output_does_not_abandon_other_workers(self):
        path = self.write('input', self.records(32), False)
        expected0 = self.run_tool([path], ['-h', '-t', '0'], False).stdout
        expected1 = self.run_tool([path], ['-h', '-t', '1'], False).stdout
        stem = self.directory / 'failure'
        Path(f'{stem}1.ldr').mkdir()  # This worker cannot open its output.
        result = self.run_tool([path], ['-h', '-t', '0,0.5,1', '-B', str(stem)],
                               False, success=False)
        self.assertEqual(result.stdout, '')
        self.assertIn('Clustering at 0.5 failed', result.stderr)
        self.assertEqual(Path(f'{stem}0.ldr').read_text(), expected0)
        self.assertEqual(Path(f'{stem}2.ldr').read_text(), expected1)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(ParallelThresholds))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp parallel threshold checks passed')
