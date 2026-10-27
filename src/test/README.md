# Executable regression tests

The regression tests live in `src/test`. The repository-root `test` path is a
permanent compatibility symlink to `src/test`, so existing commands and
`$LILLYMOL_HOME/test/...` references continue to work.

From `src`, build the tool under test with optimised settings, then run its
regression cases:

```shell
bazel build -c opt //Molecule_Tools:fileconv //Molecule_Tools:same_structures
ruby test/run_all_test.rb -bindir "$PWD/bazel-bin/Molecule_Tools" fileconv
```

The runner discovers cases beside itself and infers the repository root when
`LILLYMOL_HOME` is unset. Both `test/run_all_test.rb` (from `src`) and
`test/run_all_test.rb` through the repository-root compatibility link work.

## Test configuration bindings

`lillymol_tests.proto` describes the test configuration. Bazel generates its
Python and Ruby bindings using the protobuf compiler selected by `MODULE.bazel`:

```shell
bazel build //test:generated_protos
./copy_test_protos.sh
```

Run these commands from `src`. The copy step places `lillymol_tests_pb2.py` and
`lillymol_tests_pb.rb` beside the runner. It compares file contents rather than
timestamps, so a stale file copied from another machine is replaced. Generated
bindings are ignored by Git. `build_linux.sh` performs both steps automatically,
including in Python-only builds; a system `protoc` is not used for this schema.

Individual generation targets are `//test:lillymol_tests_py_proto` and
`//test:lillymol_tests_ruby_proto`. Check the Python bindings against Bazel's
protobuf runtime with:

```shell
bazel test //test:lillymol_tests_proto_test --test_output=errors
```

Programs run outside Bazel still need a compatible installed protobuf runtime.
The executable-test runner needs Ruby's `google-protobuf` gem. For Python,
use a protobuf runtime compatible with the generated code; copying bindings
does not install or update that runtime.
