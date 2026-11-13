These tests check tool-specific `-<option> help` requests without supplying
required queries, configuration, or input files. Each option is also checked
with an invalid qualifier before and after `help`, together with a missing
input file. All requests must exit successfully, write detailed help to stderr,
and create no output files.

Run from `src`:

```shell
bazel $(./output_user_root.sh) test -c opt //Molecule_Tools:option_help_test
```

Or use an existing directory of executables:

```shell
python3 test/option_help/run_test.py --bindir bazel-bin/Molecule_Tools
```

`options.json` lists the executable names and option letters. When adding a tool,
also add its executable to the test target's `data` in `Molecule_Tools/BUILD`.
The test target lives in that package so it can use the tools' private build
targets.
