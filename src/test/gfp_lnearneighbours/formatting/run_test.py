#!/usr/bin/env python3
"""Check stable TDT formatting and keep protobuf distances unrounded."""
from pathlib import Path
import re
import subprocess
import sys
import unittest

EXE = str(Path(sys.argv[1]).resolve())
POOL = Path(__file__).resolve().parents[2] / 'gfp_pairwise_distances/mpr_only/in/mpr.gfp'


class Formatting(unittest.TestCase):
    def run_tool(self, *args):
        result = subprocess.run([EXE, '-p', str(POOL), *map(str, args), str(POOL)],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result.stdout

    def distances(self, text, tag='DIST'):
        values = re.findall(r'^' + re.escape(tag) + r'<([^>]+)>$', text, re.MULTILINE)
        for value in values:
            self.assertRegex(value, r'^\d+(?:\.\d{1,3})?$')
            self.assertGreaterEqual(float(value), 0)
            self.assertLessEqual(float(value), 1)
        return values

    def test_bounded_lists_and_small_distances(self):
        text = self.run_tool('-h', '-n', 2)
        self.assertEqual(len(self.distances(text)), 10)
        self.assertEqual(text.count('|\n'), 5)
        self.assertIn('DIST<0.113>\n', text)
        self.assertEqual(len(re.findall(r'^PCN<[^>]+ 2>$', text, re.MULTILINE)), 5)

    def test_self_matches_and_unfilled_slots(self):
        text = self.run_tool('-n', 1)
        self.assertEqual(self.distances(text), ['0'] * 5)
        text = self.run_tool('-h', '-n', 10)
        self.assertEqual(len(self.distances(text)), 20)
        self.assertEqual(len(re.findall(r'^PCN<[^>]+ 4>$', text, re.MULTILINE)), 5)
        text = self.run_tool('-h', '-n', 2, '-T', '0.0001')
        self.assertEqual(self.distances(text), [])
        self.assertEqual(len(re.findall(r'^PCN<[^>]+ 0>$', text, re.MULTILINE)), 5)

    def test_custom_distance_tags_and_smiles_suppression(self):
        for args in (('-D', 'SCORE'), ('-B', 'DIST=SCORE')):
            text = self.run_tool('-h', '-n', 2, *args)
            self.assertEqual(len(self.distances(text, 'SCORE')), 10)
            self.assertNotIn('DIST<', text)
        text = self.run_tool('-h', '-n', 2, '-B', 'nosmiles', '-B', 'nbrtag=COUNT')
        self.assertNotIn('$SMI<', text)
        self.assertEqual(text.count('COUNT<2>\n'), 5)
        self.assertEqual(len(self.distances(text)), 10)

    def test_deferred_sorting(self):
        normal = self.run_tool('-h', '-T', 1)
        sorted_at_end = self.run_tool('-h', '-T', 1, '-B', 'bignn')
        self.assertEqual(sorted(self.distances(normal)), sorted(self.distances(sorted_at_end)))
        self.assertEqual(len(self.distances(sorted_at_end)), 20)

    def test_protobuf_keeps_numeric_precision(self):
        text = self.run_tool('-h', '-n', 2, '-k')
        self.assertEqual(len(text.splitlines()), 5)
        values = re.findall(r'\bdist: ([0-9.eE+-]+)', text)
        self.assertEqual(len(values), 10)
        self.assertTrue(any(abs(float(value) - round(float(value), 3)) > 1e-5 for value in values))
        self.assertNotIn('DIST<', text)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Formatting)
    result = unittest.TextTestRunner().run(suite)
    if not result.wasSuccessful():
        sys.exit(1)
    print('gfp_lnearneighbours formatting checks passed')
