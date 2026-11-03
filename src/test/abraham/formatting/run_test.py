#!/usr/bin/env python3
"""Check stable signed descriptors, pipeline output, and unchanged fingerprints."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

EXE = str(Path(sys.argv[1]).resolve())
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
MODELS = ROOT / 'data/queries/abraham'
MOLECULES = HERE / 'in/molecules.smi'


class Formatting(unittest.TestCase):
    def run_tool(self, *args, source=MOLECULES):
        result = subprocess.run([EXE, '-l', '-F', str(MODELS / 'Abraham'),
                                 '-P', str(MODELS / 'Alpha2H'), '-g', 'all',
                                 *map(str, args), str(source)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result.stdout

    def table(self, text, smiles=False):
        lines = [line.split() for line in text.splitlines()]
        self.assertEqual(len(lines), 7)
        start = 2 if smiles else 1
        for row in lines[1:]:
            self.assertEqual(len(row), len(lines[0]))
            for value in row[start:]:
                self.assertRegex(value, r'^-?\d+(?:\.\d{1,4})?$')
                self.assertNotEqual(value, '-0')
        return lines

    def test_signed_and_large_descriptors(self):
        rows = self.table(self.run_tool())
        self.assertEqual(rows[0], ['Name', 'abr_r2', 'abr_pi2h', 'abr_beta2h',
                                   'abr_beta2o', 'abr_logl', 'abr_alpha2h', 'abr_Vx', 'abr_Logp'])
        values = {row[0]: row[1:] for row in rows[1:]}
        self.assertEqual(values['methane'], ['0.248', '0.277', '0.071', '0.064',
                                            '0.13', '0.003', '24.95', '0.6657'])
        self.assertEqual(values['acetic_acid'][-1], '-1.3742')
        self.assertEqual(values['benzene'][-2:], ['71.64', '2.1104'])

    def test_pipeline_modes(self):
        ordinary = self.table(self.run_tool())
        pipeline = self.run_tool('-Y', 'wpipe')
        rows = self.table(pipeline, smiles=True)
        for expected, row in zip(ordinary[1:], rows[1:]):
            self.assertEqual(row[1:], expected)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'pipeline.txt'
            source.write_text(pipeline)
            for extra in ([], ['-Y', 'wpipe']):
                text = self.run_tool('-Y', 'rpipe', *extra, source=source)
                for original, actual in zip(rows[1:], [r.split() for r in text.splitlines()[1:]]):
                    # Existing columns are retained verbatim and new values match.
                    expected = original if extra else original[1:]
                    self.assertEqual(actual, expected + original[2:])

    def test_fingerprints_unchanged(self):
        # Captured before the descriptor formatting change.
        expected = (HERE / 'in/fingerprints.tdt').read_text()
        self.assertEqual(self.run_tool('-J', 'NCABR'), expected)

    def test_per_molecule_flushing(self):
        self.assertEqual(self.run_tool('-Y', 'flush'), self.run_tool())


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Formatting))
    if not result.wasSuccessful():
        sys.exit(1)
    print('abraham formatting checks passed')
