#!/usr/bin/env python3
"""Check stable spread distances and scale values, including values above one."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
POOL = Path(__file__).resolve().parents[2] / 'gfp_pairwise_distances/mpr_only/in/mpr.gfp'


class Formatting(unittest.TestCase):
    def run_tool(self, *args, pool=POOL):
        result = subprocess.run([EXE, *map(str, args), str(pool)],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result.stdout

    def values(self, text, tag='DIST'):
        values = re.findall(r'^' + tag + r'<([^>]+)>$', text, re.MULTILINE)
        for value in values:
            self.assertRegex(value, r'^\d+(?:\.\d{1,4})?$')
        return values

    def test_normal_selection(self):
        text = self.run_tool()
        self.assertEqual(self.values(text), ['1', '0.4038', '0.3214', '0.1923', '0.113'])
        ids = [record.split('PCN<')[1].split('>')[0] for record in text.strip().split('|') if record.strip()]
        self.assertEqual(ids, ['aspirin', 'decane', 'caffeine', 'benzene', 'paracetamol'])
        self.assertEqual(self.values(text, 'SCALE'), [])

    def test_previously_selected(self):
        with tempfile.TemporaryDirectory() as tmp:
            selected = Path(tmp) / 'selected.gfp'
            selected.write_text(POOL.read_text().split('|')[0] + '|\n')
            text = self.run_tool('-A', selected)
            values = self.values(text)
            self.assertEqual(len(values), 5)
            self.assertEqual(values[0], '0.4038')
            self.assertIn('0', values)

    def test_scales_and_scaled_distances_above_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            pool = Path(tmp) / 'scaled.gfp'
            pool.write_text(POOL.read_text().replace('\n|', '\nWEIGHT<3.123456>\n|'))
            text = self.run_tool('-p', 'WEIGHT', pool=pool)
            self.assertEqual(self.values(text, 'SCALE'), ['3.1235'] * 4)
            values = self.values(text)
            self.assertEqual(len(values), 5)
            self.assertGreater(float(values[1]), 1)
            # Rounding the displayed scale must not affect the distance calculation.
            self.assertAlmostEqual(float(values[1]), 0.403846145 * 3.123456, delta=0.000051)
            selected = Path(tmp) / 'selected.gfp'
            selected.write_text(pool.read_text().split('|')[0] + '|\n')
            text = self.run_tool('-p', 'WEIGHT', '-A', selected, pool=pool)
            values = self.values(text)
            self.assertEqual(len(values), 5)
            self.assertGreater(float(values[0]), 1)

    def test_precomputed_distances(self):
        with tempfile.TemporaryDirectory() as tmp:
            pool = Path(tmp) / 'precomputed.gfp'
            pool.write_text(POOL.read_text().replace('\n|', '\nPREV<0.1234567>\nNNSMI<C>\nNNID<prior>\n|'))
            text = self.run_tool('-N', 'PREV', pool=pool)
            self.assertEqual(self.values(text)[0], '0.1235')


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Formatting))
    if not result.wasSuccessful():
        sys.exit(1)
    print('gfp_spread formatting checks passed')
