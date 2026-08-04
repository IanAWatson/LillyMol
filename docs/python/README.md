# LillyMol Python

Python bindings for LillyMol, built with [pybind11](https://github.com/pybind/pybind11).

These started as a prototyping convenience and have become fast enough to use in
anger. LillyMol is a C++ cheminformatics toolkit that has always been built around
speed, and most of that survives the trip into python, because almost every call
does real work in C++ rather than shuffling python objects. Where we have measured
against RDKit on the same task, LillyMol python has run 3 to 7 times faster - see
[Performance](#performance).

The bindings cover the Molecule, substructure searching, reactions, fingerprints
and similarity, chemical standardisation, and a growing number of the LillyMol
command line tools. Names and signatures are still occasionally adjusted, so
check the reference if something does not behave as you expect.

## Contents

| Document | What is in it |
| -------- | ------------- |
| [Build.md](Build.md) | Python packages you need, building, and wheel packaging |
| [LillyMolPython.md](LillyMolPython.md) | The API reference - Molecule, Atom, Bond, Set_of_Atoms and Ring methods, substructure searching, reactions, GFP |
| [Tools.md](Tools.md) | LillyMol command line tools exposed to python - QED, Lipinski, MedchemWizard, unique molecules, Position3D |
| [descriptors.md](descriptors.md) | Molecular descriptors, singly and in batches, numpy and pandas forms |
| [fingerprints.md](fingerprints.md) | Fingerprints and similarity |
| [tsubstructure.md](tsubstructure.md) | Substructure searching as the `tsubstructure` tool does it |
| [structure_database.md](structure_database.md) | Looking molecules up in a structure database |
| [synthetic_precedent.md](synthetic_precedent.md) | Synthetic precedent from a database of known molecules |
| [4-pyridone.py](4-pyridone.py), [reaction_example.py](reaction_example.py) | Worked examples |

## Installing and verifying

See [Build.md](Build.md) for the full story. The short version, for working
inside a LillyMol checkout:

```
cd ${LILLYMOL_HOME}/src
bazel build -c opt //pybind:all
./copy_shared_libraries.sh ../lib
```

then run python through the wrapper, which sets `PYTHONPATH` and
`LD_LIBRARY_PATH` for you:

```
${LILLYMOL_HOME}/run_python.sh my_script.py
```

Build `-c opt`. A default `fastbuild` is enormously slower and is not what any of
the timings here describe.

Check it works:

```python
from lillymol import *

m = MolFromSmiles('c1ccccc1CC(=O)N')
print(m.natoms(), m.aromatic_ring_count(), m.unique_smiles())
# 10 1 O=C(N)Cc1ccccc1
```

## First steps

`MolFromSmiles` will look familiar to RDKit users, and returns `None` on failure.
Diagnostics go to stderr; nothing but the `None` arrives in python.

```python
from lillymol import *

m = MolFromSmiles('c1ccccc1')
bad = MolFromSmiles('c1ccc')          # None, with a complaint on stderr
mols = MolFromSmiles(['C', 'CC', 'C1CC1'])   # a list form also exists
```

Reading files is usually what you want. `ReaderContext` infers the file type from
the name:

```python
import lillymol_io

with lillymol_io.ReaderContext('/path/to/file.smi') as reader:
    for mol in reader:
        print(mol.name(), mol.natoms())
```

There is no `None` molecule from a reader. A connection table error stops the read
unless you have raised the allowed error count. Note also that in LillyMol a
molecule with a bad valence is still a valid molecule - call `valence_ok()` if you
care.

Molecules are always mutable. There is no distinction between a molecule you may
edit and one you may not, which removes a whole category of ceremony, at the cost
of making it your job to know when a function you called has changed its argument.

## The lazy Molecule

**This is the one thing to understand before writing anything substantial.**

A LillyMol `Molecule` computes nothing it was not asked for. Ring membership,
fragment membership, aromaticity and canonical order are all derived properties,
computed on demand and thrown away when you edit the molecule. That is why
building a molecule from smiles and asking only for `natoms()` is nearly free.

The consequence is that `Atom`, `Bond` and `Ring` objects hold no pointer back to
their `Molecule`, and so cannot know anything molecule-derived. A freshly parsed
benzene will tell you its bonds are not aromatic and not in a ring:

```python
benzene = MolFromSmiles('c1ccccc1')
bond = next(iter(benzene[0]))
bond.is_aromatic(), bond.nrings()      # (False, 0)   <- not yet computed

benzene.compute_aromaticity_if_needed()
bond.is_aromatic(), bond.nrings()      # (True, 1)
```

Note what that first line does **not** do: it does not raise. You get `False` and
`0`, which are perfectly plausible answers, and wrong. Nothing will tell you.

So if you are going to interrogate `Atom`, `Bond` or `Ring` objects, force the
perception you need first. Which call to use:

| You need | Call | Also gives you |
| -------- | ---- | -------------- |
| `Bond.nrings()`, `Bond.IsInRing()` | `mol.ring_membership()` | nothing else |
| `Bond.is_aromatic()`, `Ring.is_aromatic()` | `mol.compute_aromaticity_if_needed()` | ring membership as well |

Both are cheap once done - a second call is a flag test. Two traps worth naming:

- **`mol.nrings()` is not a force call.** It returns the ring count and does not
  push ring membership down onto the bonds. Neither does `natoms()`,
  `number_fragments()`, or anything else that merely counts.
- Any molecule-level aromaticity question, `mol.is_aromatic(0)` for instance,
  forces aromaticity as a side effect, so code that happens to ask one first
  appears to work. Do not rely on it; say what you mean.

None of this applies if you ask the `Molecule` rather than its parts.
`mol.is_aromatic(atom)`, `mol.nrings(atom)`, `mol.ring_bond_count(atom)` and
friends always give the right answer, because the molecule perceives whatever it
needs. **Asking the Molecule is both safer and faster.** Which brings us to:

## Working efficiently

`Atom` and `Bond` objects are convenient and, in python, expensive - each one you
touch has to be wrapped in a python object. The same question asked of the
`Molecule` costs one call and no wrapper.

```python
m.atomic_number(3)      # prefer this
m[3].atomic_number()    # same answer, more work
```

This is not a micro-optimisation. In a real script that walks every atom's
neighbours, replacing bond iteration with `connections()` - one call returning a
list of atom numbers - took the neighbour-walking cost from **166 us to 9 us per
molecule**, and the whole script from 12.7 s to 8.9 s.

```python
# One call per atom, no wrapper objects. bond_type_between_atoms returns the
# BondType directly, again avoiding a python Bond.
for a in range(m.natoms()):
    for nbr in m.connections(a):
        if m.bond_type_between_atoms(a, nbr) == BondType.DOUBLE_BOND:
            ...
```

Use `mol.bonds()` when you genuinely want every bond once - iterating atoms and
then their bonds visits each bond twice. Use `bond_between_atoms(a1, a2)` when you
need a real `Bond` to ask about ring membership or direction, and
`bond_type_between_atoms(a1, a2)` when the type is all you want. Both return
`None` when the atoms are not bonded.

## Translating RDKit code, and writing new code with an LLM

Translating an existing RDKit program to LillyMol python is now a realistic thing
to hand to an LLM, and we have done it successfully on non-trivial programs. Two
things to keep in mind.

**LillyMol is deliberately narrower than RDKit.** RDKit is vast and tries to cover
everything; LillyMol concentrates on core cheminformatics and is unapologetic
about it. There is no depiction or 2D coordinate generation, no conformer
generation, no force fields, and nothing like RDKit's full descriptor catalogue.
If the program you are translating leans on any of that, the translation will stop
being a translation. Check what the program actually uses before starting.

**Do not expect identical numbers.** Both toolkits are correct; they answer
slightly different questions. There is no universally accepted definition of
aromaticity, and the same is true of several other everyday quantities. Real
divergences we have measured:

| Quantity | Why the two differ |
| -------- | ------------------ |
| Aromaticity | No agreed definition. Expect the great majority of rings to agree and a small tail not to |
| H-bond donors and acceptors | LillyMol's defaults are faithful to Lipinski's paper - donors are a count of N-H and O-H *hydrogens*, so an NH2 is two. RDKit counts matching atoms, and refines which heteroatoms qualify. `rdkit_num_h_donors()` and `rdkit_num_h_acceptors()` are provided for comparison |
| Rotatable bonds | Differing treatments of conjugated linkages, sulfonamides and terminal groups |
| Molecular weight | LillyMol has no table of isotopic masses, so `amw()` erases isotopic labels and warns. `[37C]OC` weighs the same as `COC` |
| Canonical smiles | Both are canonical, and they are different canonical forms. Never compare signature *strings* across toolkits - compare the groupings they induce |

If you are having an LLM write a tool from scratch, you have a choice. Writing it
against LillyMol alone is the simplest path. Writing it against both and comparing
is more work but tells you something worth knowing: where the two agree you can
believe the answer, where they disagree you have found a definition worth
understanding rather than an average worth taking, and you get a performance
comparison for free.

When comparing, compare the right thing. Canonical smiles differ between
toolkits by construction, so equality of output strings is the wrong test. Compare
the partition - do both tools put the same molecules in the same buckets, and find
the same number of distinct answers.

## Performance

Measured on the same input, same machine, `-c opt`:

| Task | RDKit | LillyMol python |
| ---- | ----- | --------------- |
| Multi-parameter optimisation score, 20k molecules | 9.13 s | 1.37 s |
| Position-resolved chemotype extraction, 25k molecules | 27.3 s | 8.4 s |

The second is graph traversal and canonicalisation heavy, the first is descriptor
heavy, which is roughly why the ratios differ. Neither number is a benchmark
suite, and your workload may sit anywhere.

The cost you can control is the number of times you cross the C++/python boundary.
A trivial bound call - `natoms()`, `atomic_number()` - costs on the order of 150 ns,
almost all of it the crossing rather than the work. Anything that hands back a
python object costs more. Favour calls that answer a whole question at once, and
prefer the `Molecule` over its parts.

## Threads, and substructure searching

Substructure searching releases the GIL for the duration of the C++ work, so
python threads can search concurrently. What you get depends entirely on which
entry point you use. Measured on 32 cores, one query, 4000 drug sized molecules
per thread. All three below are `TSubstructure`, so they are directly comparable.

| Approach | 1 thread | 4 threads | 16 threads |
| -------- | -------- | --------- | ---------- |
| `substructure_search(mol)`, one call per molecule | **55k mol/s** | 193k (3.5x) | 183k (3.3x) |
| `substructure_search(list_of_molecules)` | 39k mol/s | 85k (2.2x) | 96k (2.5x) |
| `substructure_search(list_of_smiles)` | 47k mol/s | 182k (3.9x) | **565k (12.1x)** |

Read that table before choosing.

**One call per molecule stops improving at about three cores.** Each call releases
the GIL, but the loop around it is python and reacquires the GIL every iteration,
so you end up bounded by how fast a single python loop can run. Adding cores past
three or four gains nothing. If that is fast enough it is also the simplest thing
to write, and the fastest single-threaded option, since it neither copies
molecules nor reparses smiles.

**Passing a list of `Molecule` objects is the one to avoid.** It is the slowest of
the three and scales worst - beaten by simply calling once per molecule, at every
thread count. pybind11 converts the python list into a `std::vector<Molecule>`
before the search starts, copy constructing every molecule, and that happens with
the GIL held so it cannot overlap between threads. The same conversion is why
there is no batch `label_matched_atoms`: it would label copies and discard them.

**Passing a list of smiles is what scales.** All the work, parsing included,
happens inside one C++ call with the GIL released, so threads genuinely overlap.
It does more work per molecule than the alternatives - a smiles parse is 14.3 us
against 5.4 us to copy a Molecule - and still wins by a wide margin once there are
cores to spend. It peaked at 668k mol/s on 24 threads in a separate run; past
that, contention on this 32 core machine took it back down.

Whichever you use, the safety contract is the same, and the GIL was previously
hiding it:

- **Give each thread its own query.** `Substructure_Query` and `TSubstructure` are
  thread compatible, not thread safe; a search accumulates match state in the
  query. Building a query from smarts costs about 7 us, so per-thread queries are
  free in practice.
- **Give each thread its own molecules.** A search *writes* to the molecule it
  searches, forcing ring and aromaticity perception on it. Two threads searching
  one shared `Molecule` is a data race, which rules out the natural-looking
  "one molecule, many queries in parallel" pattern.

## Reference

The API reference is [LillyMolPython.md](LillyMolPython.md). Every method named in
its tables is checked against the built module by
[`src/pybind/lillymol_doc_test.py`](/src/pybind/lillymol_doc_test.py), so a
renamed binding fails a test rather than quietly misleading a reader. Run it with

```
cd ${LILLYMOL_HOME}/src
PYTHONPATH=bazel-bin/pybind python3 pybind/lillymol_doc_test.py
```

That test does not check prose or examples. If you change behaviour, the
surrounding text still needs a human.
