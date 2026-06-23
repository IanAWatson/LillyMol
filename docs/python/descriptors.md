# Molecular Descriptors

The `IWDescr` Python class provides the molecular descriptors computed by
[`iwdescr.sh`](https://github.com/IanAWatson/LillyMol/blob/master/docs/Molecule_Tools/iwdescr.md).
It is intended for applications that need descriptor values directly in Python
rather than a descriptor file.

`IWDescr` always computes the complete descriptor set. The Python interface does
not provide options for enabling or disabling individual descriptor families.

## Requirements

The environment variable `LILLYMOL_HOME` must point to the LillyMol installation
or source tree. `IWDescr` uses it to locate the standard charge and
donor/acceptor queries:

```
export LILLYMOL_HOME=/path/to/LillyMol
```

NumPy must also be installed:

```
pip install numpy
```

For DataFrame output, pandas is also required:

```
pip install pandas
```

Construction fails if `LILLYMOL_HOME` is not defined or the standard query files
cannot be loaded.

## Basic Usage

```python
from lillymol import Molecule
from lillymol_tools import IWDescr

mol = Molecule()
if not mol.build_from_smiles("CCO ethanol"):
    raise ValueError("Invalid SMILES")

iwdescr = IWDescr()
names = iwdescr.feature_names()
values = iwdescr.process(mol)
```

`feature_names()` returns a `list[str]`. `process()` returns a one-dimensional
NumPy array with `dtype=float32`.

The names and values use the same ordering:

```python
name_to_column = {name: column for column, name in enumerate(names)}

amw = values[name_to_column["amw"]]
mxdst = values[name_to_column["mxdst"]]

print(f"AMW {amw:.2f}, longest path {int(mxdst)}")
```

The feature ordering is stable for the lifetime of an `IWDescr` object. Fetch
the names once and reuse the mapping for each molecule:

```python
iwdescr = IWDescr()
name_to_column = {
    name: column for column, name in enumerate(iwdescr.feature_names())
}

for mol in molecules:
    values = iwdescr.process(mol)
    print(values[name_to_column["amw"]])
```

All descriptors are computed by every call to `process()`, even when the
application uses only a few values.

## Batch Processing

For collections of molecules, `process_list()` is more convenient than calling
`process()` in a loop. It accepts a list of `Molecule` objects and returns
results for all of them in one call.

### NumPy Array (default)

```python
iwdescr = IWDescr()
molecules = [...]  # list of Molecule objects

X = iwdescr.process_list(molecules)
```

`X` is a two-dimensional NumPy array of `dtype=float32` with shape
`(n_molecules, n_descriptors)`. Rows correspond to molecules in the order they
were supplied; columns correspond to descriptors in the same order as
`feature_names()`. This layout is directly compatible with scikit-learn
estimators and other tools that expect a feature matrix.

### Pandas DataFrame

Pass `as_dataframe=True` to receive a `pandas.DataFrame` instead:

```python
df = iwdescr.process_list(molecules, as_dataframe=True)
```

The DataFrame uses descriptor names as column labels and molecule names
(from `Molecule.name()`) as the row index, so individual rows and columns can
be accessed by name:

```python
print(df["amw"])             # one descriptor across all molecules
print(df.loc["ethanol"])     # all descriptors for one molecule
```

pandas is only imported when `as_dataframe=True` is requested. If pandas is not
installed and this option is used, a standard `ImportError` is raised with a
clear message.

### Choosing between the two forms

Use `process_list()` (NumPy) when feeding descriptors directly into a
machine-learning pipeline, where a plain array is what the next step expects.

Use `process_list(..., as_dataframe=True)` when exploring results interactively,
writing to CSV, or working with downstream code that benefits from named columns
and a named index:

```python
# Quick exploration
df = iwdescr.process_list(molecules, as_dataframe=True)
df.describe()
df[["amw", "nrings"]].hist()
df.to_csv("descriptors.csv")
```

## Missing Values

Some descriptors are not defined for every molecule. Undefined values are
returned as NumPy `NaN` values:

```python
import numpy

values = iwdescr.process(mol)
defined = numpy.isfinite(values)
```

The same applies to entries in the array or DataFrame returned by
`process_list()`. Rows for molecules whose calculation fails will contain NaN
for the affected descriptors. In a DataFrame these show up naturally in pandas
summary statistics and can be handled with the usual `df.dropna()` /
`df.fillna()` methods.

## Descriptor Definitions

The descriptor names and their definitions are documented in the
[`iwdescr` descriptor table](https://github.com/IanAWatson/LillyMol/blob/master/docs/Molecule_Tools/iwdescr.md#descriptors).

## Molecule and Threading Behaviour

Descriptor calculation may modify the supplied `Molecule`. In particular,
legacy chirality calculations may remove chiral-centre information. Make a copy
before calling `process()` or `process_list()` if the original molecules must
remain unchanged.

An `IWDescr` instance must not be used concurrently from multiple threads.
Separate instances can be used independently.

## Performance

The Python interface is typically about 10–15% slower than descriptor
calculation through the C++ executable. This is relatively small overhead for a
Python binding because the descriptor calculation remains in C++ and each call
returns the complete result as one NumPy array.

For large collections, prefer `process_list()` over a Python-level loop calling
`process()`. Both execute the same C++ code per molecule, but `process_list()`
allocates the output array once and avoids the per-call Python/C++ boundary
overhead of the loop.

Create one `IWDescr` object and reuse it rather than constructing a new object
for every molecule or batch.

## Complete Example

See
[`contrib/examples/descriptor_computation.py`](https://github.com/IanAWatson/LillyMol/blob/master/contrib/examples/descriptor_computation.py)
for a command-line example that reads molecules and retrieves descriptors by
name.

A batch example using `process_list`:

```python
from lillymol import Molecule
from lillymol_tools import IWDescr

smiles_list = [
    ("CCO",           "ethanol"),
    ("c1ccccc1",      "benzene"),
    ("CC(=O)O",       "acetic_acid"),
    ("c1ccccc1O",     "phenol"),
]

molecules = []
for smi, name in smiles_list:
    mol = Molecule()
    mol.build_from_smiles(f"{smi} {name}")
    molecules.append(mol)

or

molecules = slurp("/path/to/file.smi")

iwdescr = IWDescr()

# As a NumPy array for ML pipelines
X = iwdescr.process_list(molecules)
print(X.shape)   # (4, n_descriptors)

# As a DataFrame for exploration
df = iwdescr.process_list(molecules, as_dataframe=True)
print(df[["amw", "nrings", "aromring"]])
```
