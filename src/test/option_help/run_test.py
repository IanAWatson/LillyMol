#!/usr/bin/env python3
"""Check that tool-specific help works before inputs or other settings are read.

Run directly with --bindir src/bazel-bin/Molecule_Tools, or via
bazel test -c opt //Molecule_Tools:option_help_test.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


def command_line():
    parser = argparse.ArgumentParser(description=__doc__)
    default = Path(os.environ.get("TEST_SRCDIR", ".")) / os.environ.get(
        "TEST_WORKSPACE", "_main"
    ) / "Molecule_Tools"
    parser.add_argument("--bindir", type=Path, default=default)
    return parser.parse_args()


class OptionHelpTest(unittest.TestCase):
    def test_help_without_required_options_or_inputs(self):
        options = json.loads(Path(__file__).with_name("options.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            for tool, flags in options.items():
                executable = self.bindir / tool
                self.assertTrue(executable.is_file(), f"Missing executable {executable}")
                for flag in flags:
                    # Check both orderings of repeated options. Neither a bad
                    # qualifier nor a missing input should prevent useful help.
                    for values in (["help"], ["invalid-help-test", "help"],
                                   ["help", "invalid-help-test"]):
                        with self.subTest(tool=tool, flag=flag, values=values):
                            args = [str(executable)]
                            for value in values:
                                args.extend([f"-{flag}", value])
                            if len(values) > 1:
                                args.append(str(Path(tmp) / "missing.smi"))
                            result = subprocess.run(
                                args, cwd=tmp, stdin=subprocess.DEVNULL,
                                capture_output=True, text=True, timeout=15,
                            )
                            self.assertEqual(result.returncode, 0, result.stderr)
                            self.assertGreater(len(result.stderr.strip()), 40,
                                               "Missing detailed option help")
                            help_flag = "m" if (tool, flag) == ("tsubstructure", "n") else flag
                            self.assertIn(f"-{help_flag}", result.stderr)
                            self.assertEqual(result.stdout, "")
                            self.assertEqual(list(Path(tmp).iterdir()), [],
                                             "Help must not create output files")


if __name__ == "__main__":
    OptionHelpTest.bindir = command_line().bindir.resolve()
    unittest.main(argv=[__file__])
