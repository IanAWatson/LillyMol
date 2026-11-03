#!/usr/bin/env python3
"""Verify stable full, triangular and rectangular distance matrices."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
POOL = Path(__file__).resolve().parents[2] / 'gfp_pairwise_distances/mpr_only/in/mpr.gfp'
NAMES = ['aspirin', 'paracetamol', 'benzene', 'caffeine', 'decane']


class Formatting(unittest.TestCase):
    def run_tool(self, *args):
        result = subprocess.run([EXE, *map(str, args), str(POOL)],
                                capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result.stdout

    def matrix(self, text, names=NAMES):
        rows = [line.split() for line in text.splitlines()]
        self.assertEqual([row[0] for row in rows], names)
        for row in rows:
            for value in row[1:]:
                self.assertRegex(value, r'^\d+(?:\.\d{1,4})?$')
                self.assertGreaterEqual(float(value), 0)
                self.assertLessEqual(float(value), 1)
        return [row[1:] for row in rows]

    def test_full_matrix(self):
        matrix = self.matrix(self.run_tool())
        self.assertEqual(matrix[0], ['0', '0.113', '0.1923', '0.3214', '0.4038'])
        for i in range(5):
            self.assertEqual(len(matrix[i]), 5)
            self.assertEqual(matrix[i][i], '0')
            for j in range(5):
                self.assertEqual(matrix[i][j], matrix[j][i])
        self.assertEqual(self.run_tool('-e'), self.run_tool())

    def test_similarities(self):
        distances = self.matrix(self.run_tool())
        similarities = self.matrix(self.run_tool('-m'))
        for i in range(5):
            self.assertEqual(similarities[i][i], '1')
            for j in range(5):
                self.assertAlmostEqual(float(distances[i][j]) + float(similarities[i][j]), 1, delta=0.000101)

    def test_triangular_and_tab_header(self):
        full = self.matrix(self.run_tool())
        triangle = self.matrix(self.run_tool('-l'))
        self.assertEqual(triangle, [row[:i] for i, row in enumerate(full)])
        text = self.run_tool('-h', '-g', 'tab')
        header, *rows = text.splitlines()
        self.assertEqual(header.split('\t'), ['Name', *NAMES])
        self.assertEqual(self.matrix('\n'.join(rows)), full)
        self.assertTrue(all(len(row.split('\t')) == 6 for row in rows))

    def test_comparator_set(self):
        with tempfile.TemporaryDirectory() as tmp:
            comparator = Path(tmp) / 'comparator.gfp'
            comparator.write_text(POOL.read_text().split('|')[0] + '|\n')
            header, data = self.run_tool('-C', comparator).split('\n', 1)
            self.assertEqual(header.split(), ['ID', *NAMES])
            rows = self.matrix(data, names=['aspirin'])
            self.assertEqual(rows, self.matrix(self.run_tool())[:1])

    def test_legacy_fast_output(self):
        full = self.matrix(self.run_tool())
        fast = self.run_tool('-f')
        squeezed = self.run_tool('-f', '-f')
        for text in (fast, squeezed):
            rows = [line.split() for line in text.splitlines()]
            self.assertEqual([row[0] for row in rows], NAMES)
            for i, row in enumerate(rows):
                self.assertEqual(len(row), 6)
                for j, value in enumerate(row[1:]):
                    self.assertAlmostEqual(float(value), float(full[i][j]), delta=0.0051)
        self.assertIn(' .11', fast)


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Formatting))
    if not result.wasSuccessful():
        sys.exit(1)
    print('gfp_distance_matrix formatting checks passed')
