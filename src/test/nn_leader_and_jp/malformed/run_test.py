#!/usr/bin/env python3
"""Reject malformed TDT neighbour counts before writing beyond the arrays."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())


class Malformed(unittest.TestCase):
    def run_tool(self, text, mode, success=False):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'input.nn'
            path.write_text(text)
            result = subprocess.run([EXE, *mode, str(path)], capture_output=True,
                                    text=True, timeout=15)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            # Positive exit status, rather than aborting or segfaulting.
            self.assertGreater(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, '')
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def test_excess_distances(self):
        for mode in (['-h', '-t', '0.2'], ['-T', '-h', '-t', '0.2'], ['-J', '.,1']):
            for neighbours in (0, 1, 3):
                with self.subTest(mode=mode, neighbours=neighbours):
                    # All neighbours refer to the parent, so every ID resolves.
                    # The extra distance must be rejected before the array write,
                    # including the zero-length allocation case.
                    text = '$SMI<C>\nPCN<m0>\n'
                    text += 'PCN<m0>\nDIST<0>\n' * neighbours
                    text += 'DIST<0.1>\n|\n'
                    result = self.run_tool(text, mode)
                    self.assertIn('too many distances', result.stderr)
                    self.assertIn(f'expected {neighbours}', result.stderr)
                    self.assertIn("'m0'", result.stderr)

    def test_missing_distances(self):
        for mode in (['-h', '-t', '0.2'], ['-T', '-h', '-t', '0.2'], ['-J', '.,1']):
            for actual in (0, 1):
                with self.subTest(mode=mode, actual=actual):
                    text = ('$SMI<C>\nPCN<m0>\nPCN<m0>\nPCN<m0>\n'
                            + 'DIST<0>\n' * actual + '|\n')
                    result = self.run_tool(text, mode)
                    self.assertIn('expected 2 neighbours', result.stderr)
                    self.assertIn(f'read 2 identifiers and {actual} distances', result.stderr)

    def test_empty_neighbour_lists_remain_valid(self):
        text = '$SMI<C>\nPCN<a>\n|\n$SMI<CC>\nPCN<b>\n|\n'
        for mode in (['-h', '-t', '0.2'], ['-T', '-h', '-t', '0.2'], ['-J', '.,1']):
            result = self.run_tool(text, mode, success=True)
            self.assertEqual(re.findall(r'^PCN<([^>]+)>', result.stdout, re.MULTILINE),
                             ['a', 'b'])
            self.assertEqual(re.findall(r'^CSIZE<(\d+)>', result.stdout, re.MULTILINE),
                             ['1', '1'])


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Malformed))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp malformed input checks passed')
