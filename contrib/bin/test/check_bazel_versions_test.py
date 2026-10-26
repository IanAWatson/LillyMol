"""Standard-library-only tests; no Bazel installation or network required."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
from urllib.error import HTTPError

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check_bazel_versions as checker


class VersionTest(unittest.TestCase):
    def test_bazel_order(self):
        ordered = ['1.0-alpha', '1.0-alpha.1', '1.0-alpha.beta',
                   '1.0-beta', '1.0-beta.2', '1.0-beta.11', '1.0-rc.1',
                   '1.0', '1.0.0', '1.0.0.bcr.2', '1.0.0.bcr.10', '1.0.1', '1.10']
        self.assertEqual(sorted(reversed(ordered), key=checker.version_key), ordered)
        self.assertEqual(checker.version_key('1.0+build.2'), checker.version_key('1.0'))
        # Numeric spelling is Bazel's final identifier tie-breaker.
        self.assertLess(checker.version_key('1.01'), checker.version_key('1.1'))
        self.assertLess(checker.version_key('1.9'), checker.version_key('1.a'))
        self.assertLess(checker.version_key('9999'), checker.version_key(''))

    def test_dates_and_prereleases(self):
        ordered = ['2024-07-02', '2024-07-02.bcr.1', '2025-08-12',
                   '2025-11-05', '2025-11-05.bcr.1']
        self.assertEqual(sorted(reversed(ordered), key=checker.version_key), ordered)
        for version in ordered + ['20260817.0', '5.0.1.bcr.2', '1.2.3+build']:
            self.assertFalse(checker.is_prerelease(version), version)
        for version in ['2.0-rc.1', '1.0-alpha', '2025-11-05.rc.1', '2025-99-99']:
            self.assertTrue(checker.is_prerelease(version), version)

    def test_invalid_versions(self):
        for version in ['v1/2', '1..2', '1.0-rc..1', '1.0+', '18446744073709551616']:
            with self.subTest(version=version), self.assertRaises(checker.CheckError):
                checker.version_key(version)


class CheckerTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.module = self.directory / 'MODULE.bazel'
        self.registry = self.directory / 'registry'
        (self.registry / 'modules').mkdir(parents=True)

    def metadata(self, name, versions, **kwargs):
        path = self.registry / 'modules' / name / 'metadata.json'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(dict(versions=versions, **kwargs)))
        return path

    def run_checker(self, *args):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            code = checker.main([str(self.module), *args])
        return code, output.getvalue()

    def test_read_module_without_execution(self):
        self.module.write_text('''# bazel_dep(name="comment", version="9.9")
module(name="test", version="1")
bazel_dep(
    name = "real_module",  # a multiline call with aliases
    repo_name = "alias", version = "1.2", dev_dependency = True,
)
new_local_repository(name="not_a_module", build_file_content="""
bazel_dep(name="inside_string", version="9.9")
""")
bazel_dep(name="local", version="1")
local_path_override(module_name="local", path="../local")
''')
        deps, overrides = checker.read_module(self.module)
        self.assertEqual(deps, {'real_module': '1.2', 'local': '1'})
        self.assertEqual(overrides, {'local': 'local_path_override'})

    def test_dynamic_and_included_declarations_rejected(self):
        for text in ['bazel_dep(name="x", version=VERSION)',
                     'include("other.MODULE.bazel")',
                     'bazel_dep(name="x", **options)',
                     'def configure():\n    bazel_dep(name="x", version="1")']:
            self.module.write_text(text)
            with self.subTest(text=text), self.assertRaises(checker.CheckError):
                checker.read_module(self.module)

    def test_local_report_and_no_changes(self):
        self.module.write_text('''bazel_dep(name="lib", version="1.9")
bazel_dep(name="re2", version="2024-07-02")
bazel_dep(name="local", version="1")
local_path_override(module_name="local", path="../local")
''')
        self.metadata('lib', ['1.9', '1.10', '2.0-rc.1', '1.11'],
                      yanked_versions={'1.11': 'bad release'})
        self.metadata('re2', ['2024-07-02', '2025-11-05.bcr.1'])
        original = self.module.read_bytes()
        code, output = self.run_checker('--registry', str(self.registry))
        self.assertEqual(code, 0, output)
        self.assertIn('1.10', output)
        self.assertNotIn('2.0-rc.1', output)
        self.assertIn('2025-11-05.bcr.1', output)
        self.assertIn('local_path_override', output)
        self.assertIn('2 updates available; 1 skipped; 0 errors', output)
        self.assertEqual(self.module.read_bytes(), original)
        code, output = self.run_checker('--registry', str(self.registry), '--include-prereleases')
        self.assertEqual(code, 0, output)
        self.assertIn('2.0-rc.1', output)

    def test_yanked_missing_and_newer_versions(self):
        metadata = {'versions': ['1.0', '1.1', '1.2'],
                    'yanked_versions': {'1.0': 'broken', '1.2': 'also broken'}}
        result = checker.check_dependency('lib', '1.0', metadata)
        self.assertEqual(result.status, 'YANKED')
        self.assertEqual(result.latest, '1.1')
        self.assertIn('broken', result.notes)
        result = checker.check_dependency('lib', '1.5', metadata)
        self.assertEqual(result.status, 'DECLARED_NEWER')
        self.assertIn('not listed', result.notes)
        result = checker.check_dependency('lib', '1.1+build', metadata)
        self.assertEqual(result.status, 'CURRENT')
        self.assertNotIn('not listed', result.notes)
        result = checker.check_dependency('lib', '1', {'versions': ['2-rc.1']})
        self.assertEqual(result.status, 'NO_ELIGIBLE_VERSION')

    def test_missing_and_malformed_metadata(self):
        self.module.write_text('bazel_dep(name="missing", version="1")\n'
                               'bazel_dep(name="bad", version="1")\n'
                               'bazel_dep(name="good", version="1")\n')
        self.metadata('bad', ['1']).write_text('{invalid JSON')
        self.metadata('good', ['1'])
        code, output = self.run_checker('--registry', str(self.registry))
        self.assertEqual(code, 1)
        self.assertIn('2 errors', output)
        self.assertIn('CURRENT', output)
        with self.assertRaises(checker.CheckError):
            checker.Registry(self.directory / 'absent')

    def test_https_metadata(self):
        self.module.write_text('bazel_dep(name="lib", version="1.9")')
        with mock.patch.object(checker, 'urlopen') as fetch:
            fetch.return_value.__enter__.return_value.read.return_value = b'{"versions":["1.9","1.10"]}'
            code, output = self.run_checker()
            self.assertEqual(code, 0, output)
            self.assertIn('UPDATE_AVAILABLE', output)
            self.assertEqual(fetch.call_args.args[0].full_url,
                             'https://bcr.bazel.build/modules/lib/metadata.json')
            self.assertEqual(fetch.call_args.kwargs['timeout'], 15)
        with mock.patch.object(checker, 'urlopen', side_effect=HTTPError('url', 404, 'missing', {}, None)):
            code, output = self.run_checker()
            self.assertEqual(code, 1)
            self.assertIn('HTTP error 404', output)


if __name__ == '__main__':
    unittest.main()
