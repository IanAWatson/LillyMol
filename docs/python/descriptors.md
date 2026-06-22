# Molecular Descriptors

The `IWDescr` Python class provides the molecular descriptors computed by
[`iwdescr.sh`](/docs/Molecule_Tools/iwdescr.md). It is intended for applications
that need descriptor values directly in Python rather than a descriptor file.

`IWDescr` always computes the complete descriptor set. The Python interface does
not provide options for enabling or disabling individual descriptor families.

## Requirements

The environment variable `LILLYMOL_HOME` must point to the LillyMol installation
or source tree. `IWDescr` uses it to locate the standard charge and
donor/acceptor queries:

```shell
export LILLYMOL_HOME=/path/to/LillyMol
```

NumPy must also be installed:

```shell
pip install numpy
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

## Missing Values

Some descriptors are not defined for every molecule. Undefined values are
returned as NumPy `NaN` values:

```python
import numpy

values = iwdescr.process(mol)
defined = numpy.isfinite(values)
```

## Descriptor Definitions

The descriptor names and their definitions are documented in the
[`iwdescr` descriptor table](/docs/Molecule_Tools/iwdescr.md#descriptors).

## Molecule and Threading Behaviour

Descriptor calculation may modify the supplied `Molecule`. In particular,
legacy chirality calculations may remove chiral-centre information. Make a copy
before calling `process()` if the original molecule must remain unchanged.

An `IWDescr` instance must not be used concurrently from multiple threads.
Separate instances can be used independently.

## Performance

The Python interface is typically about 10–15% slower than descriptor
calculation through the C++ executable. This is relatively small overhead for a
Python binding because the descriptor calculation remains in C++ and each call
returns the complete result as one NumPy array.

For large collections, create one `IWDescr` object and reuse it rather than
constructing an object for every molecule.

## Complete Example

See
[`contrib/examples/descriptor_computation.py`](/contrib/examples/descriptor_computation.py)
for a command-line example that reads molecules and retrieves descriptors by
name.
