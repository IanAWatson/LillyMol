#!/usr/bin/env python3
"""Check member distance formatting for all clustering and output modes."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())


def neighbours():
    smiles = ['C', 'CC', 'CCC']
    distances = [[0, 0.1234567, 0.8765432],
                 [0.1234567, 0, 0.7654321],
                 [0.8765432, 0.7654321, 0]]
    records = []
    for i in range(3):
        lines = [f'$SMI<{smiles[i]}>', f'PCN<m{i} 3>']
        for j in sorted(range(3), key=lambda j: distances[i][j]):
            lines.extend([f'$SMI<{smiles[j]}>', f'PCN<m{j}>', f'DIST<{distances[i][j]}>'])
        records.append('\n'.join([*lines, '|', '']))
    return ''.join(records)


class Formatting(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.input = Path(self.tmp.name) / 'input.nn'
        self.input.write_text(neighbours())

    def run_tool(self, *args):
        result = subprocess.run([EXE, *map(str, args), str(self.input)],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result.stdout

    def distances(self, text):
        values = re.findall(r'^DIST<([^>]+)>$', text, re.MULTILINE)
        for value in values:
            self.assertRegex(value, r'^\d+(?:\.\d{1,4})?$')
        return values

    def test_leader_and_butina(self):
        for mode in ([], ['-T']):
            text = self.run_tool('-h', '-t', '0.9', *mode)
            values = self.distances(text)
            self.assertEqual(len(values), 2)
            self.assertIn('0.1235', values)
            self.assertTrue(set(values) <= {'0.1235', '0.7654', '0.8765'})
            self.assertEqual(text.count('PCN<'), 3)
            self.assertEqual(text.count('|\n'), 1)
            self.assertIn('CSIZE<3>', text)

    def test_threshold_excludes_distant_members(self):
        for mode in ([], ['-T']):
            for threshold in ('0.2', '0.1234567'):
                with self.subTest(mode=mode, threshold=threshold):
                    text = self.run_tool('-h', '-t', threshold, *mode)
                    clusters = [record for record in text.split('|') if record.strip()]
                    self.assertEqual(len(clusters), 2)
                    membership = []
                    for cluster in clusters:
                        ids = re.findall(r'^PCN<(m[0-9]+)', cluster, re.MULTILINE)
                        size = int(re.search(r'CSIZE<([0-9]+)>', cluster).group(1))
                        self.assertEqual(size, len(ids))
                        membership.append(set(ids))
                    self.assertIn({'m0', 'm1'}, membership)
                    self.assertIn({'m2'}, membership)
                    self.assertEqual(self.distances(text), ['0.1235'])

    def test_smiles_output(self):
        for mode in ([], ['-T']):
            text = self.run_tool('-h', '-a', '-t', '0.9', *mode)
            self.assertNotIn('DIST<', text)
            members = [line for line in text.splitlines() if 'CLUSTER' not in line and line.strip()]
            self.assertEqual(len(members), 2)
            values = [line.split()[-1] for line in members]
            self.assertIn('0.1235', values)
            self.assertTrue(set(values) <= {'0.1235', '0.7654', '0.8765'})

    def test_self_neighbours(self):
        text = self.run_tool('-t', '0.2')
        values = self.distances(text)
        self.assertIn('0', values)
        self.assertIn('0.1235', values)

    def test_multiple_threshold_files(self):
        stem = Path(self.tmp.name) / 'clusters'
        self.run_tool('-h', '-t', '0.2,0.9,0.2', '-B', stem)
        first = Path(str(stem) + '0.ldr').read_text()
        second = Path(str(stem) + '1.ldr').read_text()
        self.assertEqual(first, self.run_tool('-h', '-t', '0.2'))
        self.assertEqual(second, self.run_tool('-h', '-t', '0.9'))
        self.assertEqual(Path(str(stem) + '2.ldr').read_text(), first)
        self.assertEqual(self.distances(first), ['0.1235'])
        self.assertEqual(first.count('|\n'), 2)
        self.assertEqual(second.count('PCN<'), 3)
        self.assertEqual(second.count('|\n'), 1)

    def test_jarvis_patrick_files(self):
        stem = Path(self.tmp.name) / 'jp'
        expected = self.run_tool('-J', '.,1')
        self.run_tool('-J', '.,1', '-B', stem)
        self.assertEqual(Path(str(stem) + '.jpc').read_text(), expected)
        self.run_tool('-J', '.,1', '-J', '.,1', '-B', stem)
        for i in range(2):
            self.assertEqual(Path(str(stem) + str(i) + '.jpc').read_text(), expected)

    def test_output_larger_than_buffer(self):
        # Disjoint pairs exercise periodic writes and the final partial buffer.
        records = []
        for i in range(1000):
            partner = i ^ 1
            records.append(f'$SMI<C>\nPCN<m{i} 2>\n$SMI<C>\nPCN<m{i}>\n'
                           f'DIST<0>\n$SMI<C>\nPCN<m{partner}>\nDIST<0.1234567>\n|\n')
        self.input.write_text(''.join(records))
        text = self.run_tool('-h', '-t', '0.9')
        self.assertGreater(len(text), 32768)
        self.assertEqual(sorted(re.findall(r'^PCN<(m[0-9]+)', text, re.MULTILINE)),
                         sorted(f'm{i}' for i in range(1000)))
        self.assertEqual(text.count('|\n'), 500)
        self.assertEqual(self.distances(text), ['0.1235'] * 500)
        self.assertTrue(text.endswith('|\n'))

    def test_jarvis_patrick_omits_distances(self):
        text = self.run_tool('-J', '.,1')
        self.assertEqual(self.distances(text), [])
        self.assertEqual(text.count('PCN<'), 3)
        self.assertEqual(text.count('|\n'), 1)
        text = self.run_tool('-J', '.,1', '-a')
        self.assertNotIn('DIST<', text)
        self.assertIn('CLUSTER', text)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Formatting))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp formatting checks passed')
