# Python
If you wish to build the python bindings, you will need a recent version of
python. Development was done with python3.11 but should work with other versions.

## Required Python Packages
You will need to install
```
pip install pybind11 absl-py protobuf numpy
apt install python-dev
```
Make sure that python-dev and libblas-dev are installed at system level.

```
sudo apt install python-dev libblas-dev
```

## Optional Python Packages
If you wish to use the xgboost QSAR model building tools in LillyMol,
```
pip install xgboost scikit-learn matplotlib and pandas.
```

If you wish to use the molecular property profile tool in LillyMol, also
```
pip install scipy numpy
```

Things seem to work seamlessly in virtualenv and `uv`.

Note that with the default build (below) Python bindings are not built,
but 'make all' will build python related targets.

## Wheel packaging
After building the pybind modules, LillyMol can stage those Bazel-built shared
objects into a Python wheel. The wheel infrastructure lives in
`${LILLYMOL_HOME}/python`.

```
cd ${LILLYMOL_HOME}/src
bazel build //pybind:all
./copy_shared_libraries.sh ../lib

cd ${LILLYMOL_HOME}/python
./scripts/stage_wheel_files.sh

python -m pip install build wheel setuptools
python -m build --wheel
 OR
uv install build wheel setuptools
uv build --wheel

```

The generated `python/prebuilt`, `python/build`, and `python/dist` directories
are build artifacts. Commit only the wheel-building infrastructure, not the
staged binaries or generated wheel.

