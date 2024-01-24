# Fingerprints in LillyMol

Eventually LillyMol python bindings will include support for the
GFP fingerprint framework - a multi-fingerprint composite fingerprint
used extensively at Lilly.

In the meantime, certain fingerprints can be constructed as counted byte
numpy arrays. These can be used for similarity calculations.

In LillyMol Molecules and Fingerprints are completely separate entities.
Fingerprints can be generated from molecules, but once generated, they
retain no connection to the Molecule that formed the fingerprint.

For efficiency, this initial implementation returns fingerprints as
numpy vectors, containing unsigned byte values. These can be
processed efficiently. The function `tanimoto` will compute a pair-wise
similarity between two such objects.

Speed is reasonable. On CPU from 2017
```
model name      : Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz
```
reading 100k molecules takes about 1.7 seconds.

Building 100k EC type fingerprints takes 6.3 seconds

Building 100k Linear fingerprints takes 11.2 seconds

Building 100k Atom Pair fingerprints takes 8.8 seconds

All numbers include the 1.7 seconds to read the smiles.

See the script [fingerprint_nn.py](fingerprint_nn.py) for
an exmaple that finds all nearest neighbours - which
will be infeasible for large input files.

Similarities are computed via the `tanimoto` function that
computes a Tanimoto similarity between two unsigned byte
arrays.
