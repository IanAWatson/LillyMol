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

