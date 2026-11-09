"""Run with python3 test/dicer_fragment_distribution_test.py."""
import importlib.util
from pathlib import Path
import tarfile
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / 'data/dicer_fragments/build_distribution.py'
spec = importlib.util.spec_from_file_location('fragment_distribution', SCRIPT)
distribution = importlib.util.module_from_spec(spec)
spec.loader.exec_module(distribution)


class DistributionTest(unittest.TestCase):
    def test_selection_sorting_and_reproducibility(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / 'DICER_1_1.textproto').write_bytes(
                b'iso: ATT smi: "[1NH3]" nat: 1 n: 2\n'
                b'iso: ATT smi: "[1CH4]" nat: 1 n: 9\n')
            (directory / 'DICER_1_3.textproto').write_bytes(b'nat: 3 n: 1\n')
            (directory / 'DICER_2_1.textproto').write_bytes(
                b'nat: 5 n: 100\nnat: 2 n: 1\nnat: 2 n: 4\n')
            (directory / 'DICER_2_2.textproto').write_bytes(b'nat: 3 n: 3\n')
            archive = directory / 'fragments.tar.xz'
            summary = distribution.build(directory, archive, 2)
            self.assertEqual(sum(count for _, count, _ in summary), 4)
            first = archive.read_bytes()
            with tarfile.open(archive) as tar:
                self.assertEqual(tar.getnames(), ['FRAG_1_1.textproto', 'FRAG_2_1.textproto'])
                self.assertEqual(tar.extractfile('FRAG_2_1.textproto').read(),
                                 b'nat: 2 n: 4\nnat: 2 n: 1\n')
                self.assertIn(b'n: 9', tar.extractfile('FRAG_1_1.textproto').readline())
            distribution.build(directory, archive, 2)
            self.assertEqual(first, archive.read_bytes())

    def test_bad_input_preserves_existing_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            archive = directory / 'fragments.tar.xz'
            archive.write_bytes(b'original')
            for bad in (b'nat: 2 n: 4\n', b'nat: 1\n'):
                (directory / 'DICER_1_1.textproto').write_bytes(bad)
                with self.assertRaises(ValueError):
                    distribution.build(directory, archive)
                self.assertEqual(archive.read_bytes(), b'original')

    def test_empty_selection_preserves_existing_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            archive = directory / 'fragments.tar.xz'
            archive.write_bytes(b'original')
            (directory / 'DICER_1_2.textproto').write_bytes(b'nat: 2 n: 1\n')
            with self.assertRaises(ValueError):
                distribution.build(directory, archive, 1)
            self.assertEqual(archive.read_bytes(), b'original')


if __name__ == '__main__':
    unittest.main()
