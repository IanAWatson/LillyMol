"""Compatibility with metadata written by the pre-migration Python modules."""
import contextlib
import importlib
import io
import json
from pathlib import Path
import pickle
import unittest
from unittest import mock

from google.protobuf import text_format
from google.protobuf.runtime_version import VersionError
from lillymol_protos import check_runtime


class MetadataTest(unittest.TestCase):
    def test_legacy_metadata(self):
        fixtures = json.loads((Path(__file__).parent / 'testdata/legacy_metadata.json').read_text())
        for fixture in fixtures:
            with self.subTest(module=fixture['module']):
                module = importlib.import_module('lillymol_protos.' + fixture['module'])
                cls = getattr(module, fixture['message'])
                binary = bytes.fromhex(fixture['binary_hex'])
                message = cls.FromString(binary)
                expected = text_format.Parse(fixture['textproto'], cls())
                self.assertEqual(message, expected)
                self.assertEqual(cls.FromString(message.SerializeToString()), expected)
                self.assertEqual(pickle.loads(pickle.dumps(message)), message)
                self.assertEqual(cls.__module__, 'lillymol_protos.' + fixture['module'])

    def test_runtime_error_message(self):
        output = io.StringIO()
        with mock.patch.object(check_runtime.importlib, 'import_module',
                               side_effect=Exception('incompatible runtime')), contextlib.redirect_stderr(output):
            self.assertEqual(check_runtime.main(), 1)
        self.assertIn('incompatible runtime', output.getvalue())
        self.assertIn('copy_script_protos.sh', output.getvalue())

    def test_runtime_check(self):
        self.assertEqual(check_runtime.main(), 0)

    def test_runtime_update_instructions(self):
        errors = [
            VersionError('runtime is older than generated code'),
            ModuleNotFoundError("No module named 'google.protobuf'", name='google.protobuf'),
            ImportError("cannot import name 'runtime_version' from 'google.protobuf'"),
        ]
        for error in errors:
            with self.subTest(error=error):
                output = io.StringIO()
                with mock.patch.object(check_runtime.importlib, 'import_module', side_effect=error), \
                        mock.patch.object(check_runtime.sys, 'executable', '/tmp/env with spaces/bin/python'), \
                        contextlib.redirect_stderr(output):
                    self.assertEqual(check_runtime.main(), 1)
                self.assertIn("'/tmp/env with spaces/bin/python' -m pip install --upgrade protobuf",
                              output.getvalue())
                self.assertIn("uv pip install --python '/tmp/env with spaces/bin/python' --upgrade protobuf",
                              output.getvalue())
                self.assertNotIn('copy_script_protos.sh', output.getvalue())


if __name__ == '__main__':
    unittest.main()
