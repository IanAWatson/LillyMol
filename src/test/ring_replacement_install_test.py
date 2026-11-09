"""Run with python3 test/ring_replacement_install_test.py."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

INSTALLER = Path(__file__).resolve().parents[2] / 'data/ring_replacement/install.py'
spec = importlib.util.spec_from_file_location('ring_install', INSTALLER)
ring_install = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ring_install)


class RingInstallTest(unittest.TestCase):
    def collection(self, directory):
        manifest = {}
        for name, data in [('rings_5Al.smi', b'aliphatic\n'),
                           ('rings_5Ar.smi', b'aromatic\n'),
                           ('rings_56.smi', b'unspecified\n')]:
            (directory / (name + '.gz')).write_bytes(gzip.compress(data, mtime=0))
            manifest[name] = dict(bytes=len(data), records=1,
                                  sha256=hashlib.sha256(data).hexdigest())
        (directory / 'manifest.json').write_text(json.dumps(manifest))

    def test_platform_names_reinstallation_and_cleanup(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self.collection(directory)
            (directory / 'unrelated.smi').write_text('keep')
            for system, aliphatic, aromatic in [('Linux', '5A', '5a'),
                                                ('Darwin', '5Al', '5Ar'),
                                                ('Linux', '5A', '5a')]:
                self.assertEqual(ring_install.install(directory, system), 3)
                self.assertEqual((directory / f'rings_{aliphatic}.smi').read_bytes(), b'aliphatic\n')
                self.assertEqual((directory / f'rings_{aromatic}.smi').read_bytes(), b'aromatic\n')
                self.assertEqual(len(list(directory.glob('rings_*.smi'))), 3)
                self.assertEqual(len(list(directory.glob('*.gz'))), 3)
            # Existing output must be refreshed, rather than treated as proof
            # that all files are already installed.
            (directory / 'rings_5A.smi').write_text('stale')
            (directory / 'rings_56.smi').unlink()
            ring_install.install(directory, 'Linux')
            self.assertEqual((directory / 'rings_5A.smi').read_bytes(), b'aliphatic\n')
            self.assertTrue((directory / 'rings_56.smi').exists())
            self.assertEqual((directory / 'unrelated.smi').read_text(), 'keep')

    def test_corrupt_source_leaves_existing_outputs_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            self.collection(directory)
            ring_install.install(directory, 'Linux')
            before = {p.name: p.read_bytes() for p in directory.glob('*.smi')}
            (directory / 'rings_5Ar.smi.gz').write_bytes(gzip.compress(b'wrong\n'))
            with self.assertRaises(ValueError):
                ring_install.install(directory, 'Linux')
            self.assertEqual(before, {p.name: p.read_bytes() for p in directory.glob('*.smi')})

    def test_manifest_cannot_escape_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / 'manifest.json').write_text(json.dumps({'../outside.smi': {}}))
            with self.assertRaises(ValueError):
                ring_install.install(directory, 'Linux')


if __name__ == '__main__':
    unittest.main()
