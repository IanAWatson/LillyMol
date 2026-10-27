"""Exercise generated bindings using Bazel's matching protobuf runtime."""
import unittest

from test import lillymol_tests_pb2


class TestSchema(unittest.TestCase):
    def test_round_trip(self):
        tests = lillymol_tests_pb2.TestCases(executable='fileconv', difftool='diff')
        case = tests.test.add(name='example', args=['-S - input.smi'])
        case.shell_script.script_name = 'run.rb'
        case.shell_script.args.append('--check')
        copy = lillymol_tests_pb2.TestCases.FromString(tests.SerializeToString())
        self.assertEqual(copy, tests)
        self.assertEqual(copy.test[0].shell_script.script_name, 'run.rb')

    def test_optional_presence(self):
        case = lillymol_tests_pb2.TestCase()
        self.assertFalse(case.HasField('non_zero_rc_expected'))
        case.non_zero_rc_expected = False
        self.assertTrue(case.HasField('non_zero_rc_expected'))


if __name__ == '__main__':
    unittest.main()
