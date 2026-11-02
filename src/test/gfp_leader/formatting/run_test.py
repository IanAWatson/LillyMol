#!/usr/bin/env python3
"""Check stable cluster distances, discard reports and unrounded protobufs."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
POOL = Path(__file__).resolve().parents[2] / 'gfp_pairwise_distances/mpr_only/in/mpr.gfp'


class Formatting(unittest.TestCase):
    def run_tool(self, *args):
        result = subprocess.run([EXE, '-t', '0.99', *map(str, args), str(POOL)],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result.stdout

    def distances(self, text):
        values = re.findall(r'^DIST<([^>]+)>$', text, re.MULTILINE)
        for value in values:
            self.assertRegex(value, r'^\d+(?:\.\d{1,4})?$')
            self.assertGreaterEqual(float(value), 0)
            self.assertLessEqual(float(value), 1)
        return values

    def test_cluster_members(self):
        text = self.run_tool()
        self.assertEqual(self.distances(text), ['0.113', '0.1923', '0.3214', '0.4038'])
        self.assertEqual(text.count('PCN<'), 5)
        self.assertIn('CLUSTER<0>\nCSIZE<5>\n', text)
        self.assertEqual(text.count('|\n'), 1)

    def test_sorting_and_singletons(self):
        self.assertEqual(self.run_tool('-r'), self.run_tool())
        text = self.run_tool('-t', '0.00001')
        self.assertEqual(self.distances(text), [])
        self.assertEqual(text.count('CSIZE<1>\n'), 5)

    def test_previously_selected_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            selected = Path(tmp) / 'selected.gfp'
            selected.write_text(POOL.read_text().split('|')[0] + '|\n')
            report = Path(tmp) / 'discarded.tdt'
            text = self.run_tool('-A', selected, '-a', '0.2', '-L', report)
            distances = self.distances(report.read_text())
            self.assertEqual(sorted(distances), ['0', '0.113', '0.1923'])
            self.assertNotIn('PCN<aspirin>', text)
            self.assertNotIn('PCN<paracetamol>', text)
            self.assertNotIn('PCN<benzene>', text)
            self.assertEqual(text.count('PCN<'), 2)

    def test_protobuf_keeps_numeric_precision(self):
        text = self.run_tool('-k', '.')
        values = re.findall(r'\bdist: ([0-9.eE+-]+)', text)
        self.assertEqual(len(values), 4)
        self.assertTrue(any(abs(float(v) - round(float(v), 4)) > 1e-6 for v in values))
        self.assertNotIn('DIST<', text)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Formatting))
    if not result.wasSuccessful():
        sys.exit(1)
    print('gfp_leader formatting checks passed')
