# Standalone-script protocol buffers

This package builds and installs model metadata bindings for scripts that do
not need LillyMol's compiled Python extension. From `src`:

```shell
bazel build //lillymol_protos:generated_protos
./copy_script_protos.sh
```

The generated Python modules and Ruby bindings go in `../lib/lillymol_protos`.
An optional argument to `copy_script_protos.sh` selects a different `lib`
directory. The build script stages the default installation on every platform,
including builds without the compiled Python bindings or C++ XGBoost evaluator.

The schemas remain authoritative in `xgboost/` and `Utilities/General/`.
Separate script proto targets apply `import_prefix = "lillymol_protos"` to the
same source files. This preserves existing C++ includes and gives the generated
Python message classes importable module names, including for pickling. The
original proto packages, message names and field numbers are unchanged.

`package_init.py` is installed as `__init__.py`. Keeping the source directory
without a package initializer lets the installed package in `lib` take priority
even when Python is started from `src`.

Use `PYTHONPATH=/path/to/lib python -m lillymol_protos.check_runtime` to check a
particular Python environment. That environment needs a protobuf runtime that
supports the generated code. When protobuf is missing or incompatible, the
check prints a `python -m pip install --upgrade protobuf` command using the
exact interpreter that failed, and the equivalent
`uv pip install --python /path/to/python --upgrade protobuf` command.
No Python packages are installed by staging.

Run the metadata compatibility tests with:

```shell
bazel test //lillymol_protos:metadata_test --test_output=errors
```

The fixture contains binary and text metadata captured from the former contrib
bindings. Tests also verify Python pickling and the generated module names.
With a model-building environment, run the existing regression suite and the
build/evaluation integration checks from the repository root:

```shell
export PYTHONPATH="$PWD/lib:$PWD/contrib/bin${PYTHONPATH:+:$PYTHONPATH}"
python -m unittest xgbd.xgbd_make_test
python contrib/bin/test/model_protos_integration_test.py
```

The integration checks build small models and exercise Python scoring, the
Ruby-to-Python random-forest evaluator, legacy Ruby metadata loading and paths
containing spaces. Ruby checks require the `google-protobuf` gem.
