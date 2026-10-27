"""Run with a model-building Python environment and staged lillymol_protos.

Example: PYTHONPATH="$LILLYMOL_HOME/lib" python model_protos_integration_test.py
Requires the same pandas, scikit-learn, XGBoost and absl packages as the scripts.
"""
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
BIN = ROOT / 'contrib/bin'


class ModelScriptsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='model-protos-')
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.descriptors = self.directory / 'descriptor data.txt'
        self.activity = self.directory / 'activity data.txt'
        self.descriptors.write_text('ID f0 f1\n' + ''.join(
            f'mol{i} {i} {i % 3}\n' for i in range(12)))
        self.activity.write_text('ID activity\n' + ''.join(
            f'mol{i} {i / 2}\n' for i in range(12)))
        self.env = dict(os.environ, LILLYMOL_HOME=str(ROOT), PYTHON=sys.executable)
        self.env['PYTHONPATH'] = str(ROOT / 'lib') + os.pathsep + self.env.get('PYTHONPATH', '')

    def run_command(self, args):
        result = subprocess.run([str(arg) for arg in args], cwd=self.directory,
                                env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result

    def test_xgboost_build_and_python_score(self):
        model = self.directory / 'xgb model'
        self.run_command([BIN / 'xgbd_make.sh', '--activity', self.activity,
                          '--mdir', model, '--n_estimators', '3', '--max_depth', '2',
                          '--nthreads', '1', '--min_points', '1', self.descriptors])
        result = self.run_command([sys.executable, BIN / 'xgbd/xgbd_evaluate.py',
                                   '--mdir', model, self.descriptors])
        lines = result.stdout.strip().splitlines()
        header = lines.index('Id XGBD_activity')
        self.assertEqual(len(lines[header:]), 13)
        self.assertTrue((model / 'model_metadata.dat').is_file())

    @unittest.skipUnless((ROOT / 'bin' / platform.system() / 'xgboost_model_evaluate').is_file(),
                         'C++ XGBoost evaluator is not installed')
    def test_xgboost_build_and_ruby_cpp_score(self):
        model = self.directory / 'cpp_model'
        descriptors = self.directory / 'cpp_input.txt'
        descriptors.write_text(self.descriptors.read_text())
        self.run_command([BIN / 'xgbd_make.sh', '--activity', self.activity,
                          '--mdir', model, '--n_estimators', '3', '--max_depth', '2',
                          '--nthreads', '1', '--min_points', '1', descriptors])
        self.env['PATH'] = str(ROOT / 'bin' / platform.system()) + os.pathsep + self.env['PATH']
        result = self.run_command([BIN / 'xgbd_evaluate.sh', '-mdir', model, descriptors])
        self.assertEqual(sum(line.startswith('mol') for line in result.stdout.splitlines()), 12,
                         result.stdout + result.stderr)

    def test_random_forest_build_and_ruby_python_score(self):
        model = self.directory / 'rf model'
        self.run_command([BIN / 'rf_make.sh', '--activity', self.activity,
                          '--mdir', model, '--n_estimators', '3', '--max_depth', '2',
                          '--n_jobs', '1', self.descriptors])
        result = self.run_command([BIN / 'rf_evaluate.sh', '-mdir', model, self.descriptors])
        self.assertEqual(len(result.stdout.strip().splitlines()), 13)
        self.assertTrue(result.stdout.startswith('Id RF_activity\n'))

    def test_xgboost_ruby_does_not_require_python(self):
        self.env['PYTHON'] = '/no/such/python'
        result = self.run_command([BIN / 'xgbd_evaluate.sh'])
        self.assertIn('Scores an xgboost descriptor model', result.stderr)

    def test_ruby_reads_legacy_metadata(self):
        fixture = ROOT / 'src/lillymol_protos/testdata/legacy_metadata.json'
        code = '''require 'json'
require File.join(ENV.fetch('LILLYMOL_HOME'), 'lib/lillymol_protos/xgboost_model_pb')
require File.join(ENV.fetch('LILLYMOL_HOME'), 'lib/lillymol_protos/random_forest_model_pb')
classes = {'xgboost_model_pb2' => XgboostModel::XGBoostModel,
           'random_forest_model_pb2' => RandomForestModel::RandomForestModel}
JSON.parse(File.read(ARGV.fetch(0))).each do |fixture|
  klass = classes[fixture.fetch('module')]
  next unless klass
  message = klass.decode([fixture.fetch('binary_hex')].pack('H*'))
  abort 'Wrong response' unless message.response == 'potency'
  abort 'Wrong columns' unless message.name_to_col.to_h == {'f0' => 0, 'f1' => 1}
  abort 'Wrong parameters' unless message.parameters.max_depth == 3
  abort 'Round trip failed' unless klass.decode(klass.encode(message)) == message
end
'''
        self.run_command(['ruby', '-e', code, fixture])


if __name__ == '__main__':
    unittest.main()
