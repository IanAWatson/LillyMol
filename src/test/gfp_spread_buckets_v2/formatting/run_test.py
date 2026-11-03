#!/usr/bin/env python3
"""Check stable distances without changing bucket selection or weighting."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
POOL = Path(__file__).resolve().parents[2] / 'gfp_pairwise_distances/mpr_only/in/mpr.gfp'


class Formatting(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.pool = Path(self.tmp.name) / 'pool.gfp'
        records = POOL.read_text().split('|')
        for i in range(5):
            records[i] = records[i].rstrip() + f'\nBUCKET<{i}>\n'
        self.pool.write_text('|'.join(records))

    def run_tool(self, *args):
        result = subprocess.run([EXE, '-B', 'BUCKET', *map(str, args), str(self.pool)],
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

    def test_normal_selection(self):
        text = self.run_tool()
        self.assertEqual(self.distances(text), ['0.4038', '0.3214', '0.1923', '0.113'])
        records = [r for r in text.split('|') if r.strip()]
        ids = [re.search(r'^PCN<([^>]+)>', r, re.MULTILINE).group(1) for r in records]
        self.assertEqual(ids, ['aspirin', 'decane', 'caffeine', 'benzene', 'paracetamol'])
        self.assertNotIn('DIST<', records[0])
        self.assertEqual(sorted(re.findall(r'^BUCKET<([0-9]+)>', text, re.MULTILINE)),
                         ['0', '1', '2', '3', '4'])
        limited = self.run_tool('-n', 3)
        self.assertEqual(limited, '|'.join(text.split('|')[:3]) + '|\n')

    def test_previously_selected(self):
        selected = Path(self.tmp.name) / 'selected.gfp'
        selected.write_text(self.pool.read_text().split('|')[0] + '|\n')
        for options in ([], ['-a']):
            text = self.run_tool('-A', selected, *options)
            values = self.distances(text)
            self.assertEqual(len(values), 5)
            self.assertIn('0.4038', values)
            self.assertIn('0', values)
            self.assertEqual(text.count('|\n'), 5)

    def test_scaled_distances_above_one(self):
        self.pool.write_text(self.pool.read_text().replace('\n|', '\nWEIGHT<3.123456>\n|'))
        text = self.run_tool('-p', 'WEIGHT')
        values = self.distances(text)
        self.assertEqual(len(values), 4)
        self.assertGreater(float(values[0]), 1)
        self.assertAlmostEqual(float(values[0]), 0.403846145 * 3.123456, delta=0.000051)
        selected = Path(self.tmp.name) / 'selected.gfp'
        selected.write_text(self.pool.read_text().split('|')[0] + '|\n')
        text = self.run_tool('-p', 'WEIGHT', '-A', selected)
        values = self.distances(text)
        self.assertEqual(len(values), 5)
        self.assertGreater(max(map(float, values)), 1)

    def test_bucket_rotation_and_penalties(self):
        records = self.pool.read_text().split('|')
        for i in range(5):
            records[i] = records[i].replace(f'BUCKET<{i}>', f'BUCKET<{i % 2}>')
        self.pool.write_text('|'.join(records))
        for options in (['-o'], ['-b', '0.5']):
            text = self.run_tool(*options)
            self.assertEqual(len(self.distances(text)), 4)
            self.assertEqual(text.count('|\n'), 5)
            self.assertEqual(text.count('BUCKET<0>\n'), 3)
            self.assertEqual(text.count('BUCKET<1>\n'), 2)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Formatting))
    if not result.wasSuccessful():
        sys.exit(1)
    print('gfp_spread_buckets formatting checks passed')
