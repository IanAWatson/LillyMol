#!/usr/bin/env python3
"""Check portable formatting across the activity profiler's output modes."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())


class Formatting(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        records = [('a', '100', '10:1 20:1'), ('b', '010', '10:1 30:1'),
                   ('c', '001', '20:1 30:1'), ('d', '110', '40:1')]
        self.pool = self.root / 'pool.gfp'
        self.pool.write_text(''.join(f'PCN<{name}>\nFP<{fp}>\nNC<{nc}>\n|\n'
                                     for name, fp, nc in records))
        self.activities = self.root / 'activity.txt'
        self.activities.write_text('ID activity\na -2.5\nb 0.00012\nc 1.25\nd 12.5\n')
        self.classes = self.root / 'classes.txt'
        self.classes.write_text('ID class\na low\nb high\nc low\nd high\n')

    def run_tool(self, *args, classification=False, fingerprints=('FP,ascii', 'NC,nc,ascii'),
                 success=True):
        activity = self.classes if classification else self.activities
        command = [EXE, '-E', str(activity), '-E', 'skiphdr']
        for fp in fingerprints:
            command.extend(('-F', fp))
        if classification:
            command.append('-C')
        result = subprocess.run([*command, *map(str, args), str(self.pool)],
                                capture_output=True, text=True, timeout=15)
        self.assertGreaterEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result.stdout

    def assert_fixed_numbers(self, text, digits):
        count = 0
        for token in text.split():
            if re.fullmatch(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', token):
                count += 1
                self.assertRegex(token, rf'^-?\d+(?:\.\d{{1,{digits}}})?$')
                self.assertNotEqual(token, '-0')
                if '.' in token:
                    self.assertFalse(token.endswith('0'), token)
        self.assertGreater(count, 0)

    def test_fixed_and_sparse_activity_output(self):
        result = self.run_tool()
        self.assertIn('bit 0 set 2 ave 5 not set 2 ave 0.6251 diff 4.3749', result)
        self.assertIn('bit 10 set 2 btw -2.5 and 0.0001 ave -1.2499 diff 4.0625', result)
        self.assertIn('bit 40 set 1 btw 12.5 and 12.5 ave 12.5', result)
        self.assert_fixed_numbers(result, 4)
        self.assertEqual(result, self.run_tool('-M', 'prec=4'))

    def test_signed_differences_sparse_complement_and_table(self):
        table = self.root / 'table.txt'
        result = self.run_tool('-j', '-g', '-T', table)
        self.assertIn('Nset 2 1.25 to 12.5 ave 6.875 diff 8.1249', result)
        self.assertIn('diff -12.9166', result)
        self.assert_fixed_numbers(result, 4)
        text = table.read_text()
        self.assertIn('40 hit 1 ave 12.5 diff 9.6875', text)
        self.assertIn('10 hit 2 ave -1.2499 diff -4.0625', text)
        self.assert_fixed_numbers(text, 4)

    def test_precision_for_all_continuous_statistics(self):
        for digits in (1, 2, 7, 12):
            table = self.root / 'table.txt'
            result = self.run_tool('-j', '-g', '-M', f'prec={digits}', '-T', table)
            self.assert_fixed_numbers(result, digits)
            self.assert_fixed_numbers(table.read_text(), digits)
        self.assertIn('0.00012', self.run_tool('-M', 'prec=7'))
        self.run_tool('-M', 'prec=0', success=False)

    def test_classification_log_ratios(self):
        result = self.run_tool(classification=True)
        self.assertIn(' NB -0.585', result)
        self.assertIn(' NB 0.415', result)
        self.assert_fixed_numbers(result, 4)
        result = self.run_tool('-M', 'prec=7', classification=True)
        self.assertIn('-0.5849625', result)
        self.assert_fixed_numbers(result, 7)
        self.run_tool('-E', self.root / 'missing_classes.txt', classification=True, success=False)

    def test_fixed_only_and_sparse_only_modes(self):
        for fingerprints in (('FP,ascii',), ('NC,nc,ascii',)):
            for classification in (False, True):
                result = self.run_tool(fingerprints=fingerprints, classification=classification)
                self.assertTrue(result.strip())
                self.assert_fixed_numbers(result, 4)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Formatting)
    result = unittest.TextTestRunner().run(suite)
    if not result.wasSuccessful():
        sys.exit(1)
    print('gfp_profile_activity_by_bits formatting checks passed')
