#!/usr/bin/env python3
"""Exercise file, command-line and random-pair distance formatting."""
import itertools
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
FIXTURES = Path(__file__).resolve().parents[1] / 'mpr_only'
POOL = FIXTURES / 'in/mpr.gfp'
PAIRS = FIXTURES / 'in/pairs'


class Formatting(unittest.TestCase):
    def run_tool(self, *args, pool=POOL, success=True):
        result = subprocess.run([EXE, '-p', str(pool), *map(str, args)],
                                capture_output=True, text=True, timeout=15)
        self.assertGreaterEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result.stdout

    def test_default_and_explicit_four_places(self):
        expected = (FIXTURES / 'out/stdout').read_text()
        self.assertEqual(self.run_tool(PAIRS), expected)
        self.assertEqual(self.run_tool('-o', 4, PAIRS), expected)

    def test_precision_and_integer_zero(self):
        reference = [float(line.split()[2]) for line in self.run_tool('-o', 7, PAIRS).splitlines()]
        for digits in (1, 2, 4, 7, 12):
            lines = self.run_tool('-o', digits, PAIRS).splitlines()
            self.assertEqual(len(lines), len(reference))
            for line, expected in zip(lines, reference):
                value = line.split()[2]
                self.assertRegex(value, rf'^\d+(?:\.\d{{1,{digits}}})?$')
                self.assertLessEqual(abs(float(value) - expected), 0.5 * 10 ** -digits + 1e-7)
            self.assertEqual(lines[-1].split()[2], '0')
        self.run_tool('-o', 0, PAIRS, success=False)

    def test_command_line_pairs_and_extra_tokens(self):
        args = []
        for line in PAIRS.read_text().splitlines():
            args.extend(('-K', ','.join(line.split())))
        self.assertEqual(self.run_tool(*args), self.run_tool(PAIRS))
        self.assertEqual(self.run_tool('-a', '-K', 'aspirin,paracetamol,extra'),
                         'aspirin paracetamol 0.113 extra\n')

    def test_random_pairs_use_the_same_formatter_and_valid_indices(self):
        identifiers = re.findall(r'^PCN<([^>]+)>', POOL.read_text(), re.MULTILINE)
        args = ['-o', 3]
        for lhs, rhs in itertools.permutations(identifiers, 2):
            args.extend(('-K', lhs + ',' + rhs))
        expected = {(lhs, rhs): distance for lhs, rhs, distance in
                    (line.split() for line in self.run_tool(*args).splitlines())}
        lines = self.run_tool('-r', 200, '-o', 3).splitlines()
        self.assertEqual(len(lines), 200)
        for line in lines:
            lhs, rhs, distance = line.split()
            self.assertEqual(distance, expected[(lhs, rhs)])

    def test_random_pairs_reject_a_pool_with_one_valid_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            pool = Path(tmp) / 'single.gfp'
            pool.write_text(POOL.read_text().split('|')[0] + '|\n')
            self.run_tool('-s', 2, '-r', 1, pool=pool, success=False)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Formatting)
    result = unittest.TextTestRunner().run(suite)
    if not result.wasSuccessful():
        sys.exit(1)
    print('gfp_pairwise_distances formatting checks passed')
