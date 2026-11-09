"""Run with python3 test/uncompress_and_install_test.py."""
import lzma
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

INSTALLER = Path(__file__).resolve().parents[1] / 'uncompress_and_install.sh'


@unittest.skipUnless(shutil.which('xz'), 'requires xz')
class SidechainInstallTest(unittest.TestCase):
    def run_installer(self, root):
        return subprocess.run(['bash', str(INSTALLER)],
                              env=dict(os.environ, LILLYMOL_HOME=str(root)),
                              capture_output=True, text=True)

    def test_reinstallation_updates_output_and_retains_archive(self):
        with tempfile.TemporaryDirectory(prefix='sidechain installation ') as tmp:
            root = Path(tmp)
            data = root / 'data'
            data.mkdir()
            archive = data / 'chembl_sidechains.textproto.xz'
            output = data / 'chembl_sidechains.textproto'
            for contents in (b'first collection\n', b'updated collection\n'):
                compressed = lzma.compress(contents)
                archive.write_bytes(compressed)
                result = self.run_installer(root)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(output.read_bytes(), contents)
                self.assertEqual(archive.read_bytes(), compressed)
                self.assertEqual(sorted(p.name for p in data.iterdir()),
                                 [output.name, archive.name])

    def test_corruption_preserves_previous_output_and_removes_temporary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root / 'data'
            data.mkdir()
            archive = data / 'chembl_sidechains.textproto.xz'
            # A truncated stream may emit bytes before reporting failure.
            archive.write_bytes(lzma.compress(b'new collection\n')[:-8])
            output = data / 'chembl_sidechains.textproto'
            output.write_bytes(b'existing collection\n')
            result = self.run_installer(root)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(output.read_bytes(), b'existing collection\n')
            self.assertTrue(archive.exists())
            self.assertEqual(len(list(data.iterdir())), 2)

    def test_corruption_does_not_create_final_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root / 'data'
            data.mkdir()
            (data / 'chembl_sidechains.textproto.xz').write_bytes(b'invalid xz')
            result = self.run_installer(root)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((data / 'chembl_sidechains.textproto').exists())
            self.assertEqual(len(list(data.iterdir())), 1)


if __name__ == '__main__':
    unittest.main()
