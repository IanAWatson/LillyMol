#!/usr/bin/env python3
"""Run the Butina edge-case suite against independent static initial counts."""
from pathlib import Path
import subprocess
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from butina.run_test import Butina, EXE, recount

class StaticButina(Butina):
    def run_tool(self, path, args, serialized):
        return super().run_tool(path, [*args, '-O', 'static'], serialized)

    def check(self, *args, **kwargs):
        return super().check(*args, **kwargs, static_order=True)

    def test_static_and_dynamic_really_differ(self):
        rows = [[(1, .125), (2, .125)], [(0, .125), (2, .125)],
                [(0, .125), (1, .125)], [(0, .125), (1, .125)], [(5, .125)], []]
        static = recount(rows, .25, [4] * 6, static_order=True)
        self.assertEqual(static, [[0, 1, 2], [3], [4, 5]])
        self.assertNotEqual(static, recount(rows, .25, [4] * 6))
        for serialized in (False, True):
            self.check(rows, .25, serialized)

    def test_parallel_thresholds_with_shared_reverse_index(self):
        rows, caps = self.random_rows(92, 64)
        for serialized in (False, True):
            path, ids = self.write(rows, caps, serialized)
            stem = self.directory / 'static'
            self.run_tool(path, ['-h', '-M', 'col=2', '-t', '1,0,.25,0',
                                 '-B', str(stem)], serialized)
            for i, threshold in enumerate((1, 0, .25, 0)):
                groups = recount(rows, threshold, caps, column=True, static_order=True)
                self.assert_clusters(Path(f'{stem}{i}.ldr').read_text(),
                                     [[ids[j] for j in group] for group in groups])

    def test_explicit_dynamic_matches_default(self):
        rows, caps = self.random_rows(87, 32)
        for serialized in (False, True):
            path, _ = self.write(rows, caps, serialized)
            args = ['-h', '-t', '.25']
            default = Butina.run_tool(self, path, args, serialized)
            explicit = Butina.run_tool(self, path, [*args, '-O', 'dynamic'], serialized)
            self.assertEqual(default.stdout, explicit.stdout)

    def test_bad_order_option(self):
        rows, caps = self.random_rows(9, 16)
        path, _ = self.write(rows, caps, False)
        for args in (['-O', 'static'], ['-T', '-O', 'bad']):
            result = subprocess.run([EXE, '-t', '.25', *args, str(path)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('requires -T', result.stderr)

if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(StaticButina))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp static Butina checks passed')
