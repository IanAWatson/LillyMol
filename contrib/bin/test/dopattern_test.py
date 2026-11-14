#!/usr/bin/env python3
"""Command execution, discovery and temporary-file handling for dopattern."""
from pathlib import Path
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'dopattern.rb'
RUBY = shutil.which('ruby')

@unittest.skipUnless(RUBY, 'ruby is required')
class CommandFiles(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.work = self.root / 'working directory'
        self.work.mkdir()

    def run_tool(self, *args, env=None):
        return subprocess.run([RUBY, str(SCRIPT), *args], cwd=self.work,
                              capture_output=True, text=True, env=env, timeout=15)

    def readonly(self):
        if os.geteuid() == 0:
            self.skipTest('root can write despite directory mode bits')
        self.work.chmod(0o555)
        self.addCleanup(self.work.chmod, 0o755)

    def command_path(self, stderr):
        match = re.search(r" < (.+)'\n", stderr)
        self.assertIsNotNone(match, stderr)
        return Path(shlex.split(match[1])[0])

    def test_parallel_and_xargs_are_mutually_exclusive(self):
        for options in [('-parallel', '2', '-xargs', '2'),
                        ('-xargs', '2', '-parallel', '2')]:
            with self.subTest(options=options):
                result = self.run_tool('-o', '1', *options, 'echo %')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('-parallel and -xargs are mutually exclusive', result.stderr)
                self.assertEqual(result.stdout, '')
                self.assertEqual(list(self.work.iterdir()), [])

    def test_writable_xargs_and_cleanup(self):
        result = self.run_tool('-v', '-o', '3', '-xargs', '2', 'echo %')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(sorted(result.stdout.splitlines()), ['1', '2', '3'])
        self.assertNotIn('Warning:', result.stderr)
        path = self.command_path(result.stderr)
        self.assertEqual(path.parent, self.work)
        self.assertFalse(path.exists())

    def test_readonly_xargs_falls_back_and_keeps_cwd(self):
        self.readonly()
        result = self.run_tool('-v', '-o', '2', '-xargs', '2', 'pwd')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), [str(self.work)] * 2)
        self.assertIn('Warning:', result.stderr)
        self.assertIn('trying /tmp', result.stderr)
        path = self.command_path(result.stderr)
        self.assertEqual(path.parent, Path('/tmp'))
        self.assertFalse(path.exists())

    def test_readonly_dry_parallel_retains_commands(self):
        self.readonly()
        result = self.run_tool('-dry', '-o', '2', '-parallel', '2', 'echo %')
        self.assertEqual(result.returncode, 0, result.stderr)
        path = self.command_path(result.stderr)
        self.addCleanup(path.unlink, missing_ok=True)
        self.assertEqual(path.parent, Path('/tmp'))
        self.assertEqual(path.read_text(), 'echo 1\necho 2\n')
        self.assertEqual(path.stat().st_mode & 0o777, 0o600)

    def test_writable_dry_runs_have_distinct_files(self):
        paths = []
        for _ in range(2):
            result = self.run_tool('-dry', '-o', '1', '-xargs', '1', 'echo %')
            self.assertEqual(result.returncode, 0, result.stderr)
            paths.append(self.command_path(result.stderr))
        self.assertNotEqual(paths[0], paths[1])
        self.assertTrue(all(p.read_bytes() == b'echo 1\0' for p in paths))

    def test_readonly_array_script_uses_same_fallback(self):
        self.readonly()
        result = self.run_tool('-dry', '-o', '2', '-cluster.seq', 'echo %')
        self.assertEqual(result.returncode, 0, result.stderr)
        match = re.search(r"commands in '(.+)'", result.stderr)
        self.assertIsNotNone(match, result.stderr)
        path = Path(match[1])
        self.addCleanup(path.unlink, missing_ok=True)
        self.assertEqual(path.parent, Path('/tmp'))
        self.assertIn('echo ${SGE_TASK_ID}', path.read_text())

    def test_asynchronous_submitter_can_keep_reading_file(self):
        # A local stub records the submitted filename; no scheduler is contacted.
        args_file = self.root / 'submitted'
        submitter = self.root / 'submit'
        submitter.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > ' + shlex.quote(str(args_file)) + '\n')
        submitter.chmod(0o755)
        result = self.run_tool('-o', '2', '-cluster', '-submit', str(submitter), 'echo %')
        self.assertEqual(result.returncode, 0, result.stderr)
        path = Path(args_file.read_text().strip())
        self.assertTrue(path.is_file())
        self.assertEqual(path.read_text(), 'echo 1\necho 2\n')

    def test_both_locations_fail_cleanly(self):
        mock = self.root / 'deny_tempfiles.rb'
        mock.write_text("require 'tempfile'\nmodule DenyTemporaryFiles\n"
                        "  def create(*args, **kwargs, &block)\n"
                        "    raise Errno::EACCES, 'test permission failure'\n"
                        "  end\nend\nTempfile.singleton_class.prepend(DenyTemporaryFiles)\n")
        env = dict(os.environ, RUBYOPT='-r' + str(mock))
        result = self.run_tool('-o', '1', '-xargs', '1', 'echo %', env=env)
        self.assertEqual(result.returncode, 3, result.stderr)
        self.assertIn('Warning:', result.stderr)
        self.assertIn('Cannot create temporary command file in /tmp', result.stderr)
        self.assertNotIn('.rb:', result.stderr)

    def test_xargs_preserves_shell_quotes_backslashes_and_newlines(self):
        for command in ("printf '%s\\n' 'a b' 'c\\d' @",
                        "printf '%s\\n' 'first\nsecond' @",
                        "printf '%s\\n' \"double quoted\" @"):
            serial = self.run_tool('-o', '1', '-f', '@', command)
            parallel = self.run_tool('-o', '1', '-f', '@', '-xargs', '2', command)
            self.assertEqual(serial.returncode, 0, serial.stderr)
            self.assertEqual(parallel.returncode, 0, parallel.stderr)
            self.assertEqual(parallel.stdout, serial.stdout)
        self.assertEqual(list(self.work.glob('DOPATTERNTMP*')), [])

    def test_xargs_reports_job_failure_and_cleans_up(self):
        result = self.run_tool('-o', '1', '-xargs', '2', 'false')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(list(self.work.glob('DOPATTERNTMP*')), [])

    def test_keep_going_still_reports_failure(self):
        for option in ('-k', '-ks'):
            result = self.run_tool('-o', '3', '-f', '@', option,
                                   'echo @; test @ -ne 2')
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(result.stdout.splitlines(), ['1', '2', '3'])
            if option == '-ks':
                self.assertEqual(result.stderr, '')

    def test_parallel_status_and_cleanup_with_local_stub(self):
        tools = self.root / 'tools'
        tools.mkdir()
        parallel = tools / 'parallel'
        parallel.write_text('#!/bin/sh\nexit 7\n')
        parallel.chmod(0o755)
        env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ['PATH'])
        result = self.run_tool('-o', '1', '-parallel', '2', 'echo %', env=env)
        self.assertEqual(result.returncode, 7, result.stderr)
        self.assertEqual(list(self.work.glob('DOPATTERNTMP*')), [])

    def test_sync_cleanup_and_status_with_local_stub(self):
        submitter = self.root / 'submit helper'
        for status in (0, 9):
            submitter.write_text(f'#!/bin/sh\nexit {status}\n')
            submitter.chmod(0o755)
            result = self.run_tool('-o', '1', '-sync', '-submit', str(submitter), 'echo %')
            self.assertEqual(result.returncode, status, result.stderr)
            self.assertNotIn('unlink', result.stderr)
            self.assertEqual(list(self.work.glob('DOPATTERNTMP*')), [])

    def test_rxdir_size_checks_use_the_scanned_directory(self):
        inputs = self.root / 'inputs'
        inputs.mkdir()
        (inputs / 'full.smi').write_text('C full\n')
        (inputs / 'empty.smi').touch()
        # Opposite contents in cwd catch both false inclusion and false exclusion.
        (self.work / 'full.smi').touch()
        (self.work / 'empty.smi').write_text('C other\n')
        for selector in (['-rx', r'(.+)\.smi$'], ['-suffix', 'smi']):
            result = self.run_tool('-dry', *selector, '-rxdir', str(inputs), '-s', 'echo %')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, 'echo full\n')
        (inputs / 'directory').mkdir()
        result = self.run_tool('-dry', '-subdir', '-rxdir', str(inputs), 'echo %')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, 'echo directory\n')

    def test_array_sequence_rejects_gaps_duplicates_and_wrong_order(self):
        for sequence in ('2,4', '1,1', '2,1', '0,1', '01,02', 'a,b'):
            result = self.run_tool('-do', sequence, '-cluster.seq', 'echo %')
            self.assertEqual(result.returncode, 3, result.stderr)
            self.assertIn('requires the ordered sequence', result.stderr)
            self.assertEqual(result.stdout, '')
            self.assertEqual(list(self.work.iterdir()), [])
        result = self.run_tool('-a', '2', '-o', '4', '-cluster.seq', 'echo %')
        self.assertEqual(result.returncode, 3, result.stderr)

    def test_array_script_uses_the_requested_placeholder_and_width(self):
        result = self.run_tool('-dry', '-do', '1,2,3', '-f', '@', '-w', '2',
                               '-cluster.seq', 'echo @')
        self.assertEqual(result.returncode, 0, result.stderr)
        path = Path(re.search(r"commands in '(.+)'", result.stderr)[1])
        for task in (1, 2, 3):
            env = dict(os.environ, SGE_TASK_ID=str(task))
            run = subprocess.run(['bash', str(path)], capture_output=True, text=True,
                                 env=env, timeout=15)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertEqual(run.stdout, f'{task:02d}\n')

    def test_array_submission_propagates_status_without_contacting_scheduler(self):
        mock = self.root / 'fake_array_submission.rb'
        mock.write_text("module FakeArraySubmission\n"
                        "  def system(*args)\n"
                        "    if args.first.to_s.start_with?('/etc/cluster-setup.sh && qsub ')\n"
                        "      super('bash', '-c', 'exit 11')\n"
                        "    else\n      super\n    end\n  end\nend\n"
                        "Object.prepend(FakeArraySubmission)\n")
        env = dict(os.environ, RUBYOPT='-r' + str(mock))
        result = self.run_tool('-o', '2', '-cluster.seq', 'echo %', env=env)
        self.assertEqual(result.returncode, 11, result.stderr)

if __name__ == '__main__':
    unittest.main()
