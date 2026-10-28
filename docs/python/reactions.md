# Reactions in Python

LillyMol reactions support in-place transformations, generation of separate
products at individual matching sites, and combinatorial enumeration with stored
sidechain reagents. The active nanobind API is available directly from `lillymol`:

```python
from lillymol import MolFromSmiles, Reaction, ReactionIterator
from lillymol import Set_of_Atoms, SidechainMatchConditions, FileType
```

The older `lillymol_query` and `lillymol_reaction` modules are compatibility shims;
importing the query module first is no longer required.

## Running the examples

The [example directory](../../src/nanobind/examples/reactions/) contains four
standalone scripts and their input files. Each script checks its expected
products and can be run from any working directory. From an installed checkout,
with a virtual environment matching the staged LillyMol Python ABI:

```shell
"$LILLYMOL_HOME/run_python.sh" "$LILLYMOL_HOME/src/nanobind/examples/reactions/in_place.py"
```

Replace `in_place.py` with `per_embedding.py`, `explicit_reagents.py`, or
`enumerate_reagents.py` to run the other workflows. The scripts below are shown
in full, including the assertions exercised by their Bazel tests.

## Constructing a reaction

`Reaction.read(fname)` and `Reaction.read_proto(fname)` read a reaction
**textproto** file. Both use the reaction reader, rather than the inherited
substructure-query reader. Each returns `True` on
success and `False` on failure. `construct_from_textproto(text)` accepts the
textproto as a string; `construct_from_smirks(text)` accepts SMIRKS. Both return
booleans. Check construction before using the reaction.

File-based construction supplies the reaction filename for relative references.
String-based textproto construction uses the current directory as its base.
For SMIRKS, atom maps establish correspondence between left and right sides;
leaving-group handling also depends on the reaction's SMIRKS configuration.
The simple SMIRKS example below changes a mapped atom's element without removing
atoms. The amide example uses explicit textproto removal and joining operations.

`name()` returns the reaction name/comment, and `set_name(text)` changes it.
A named reaction is represented as, for example,
`Reaction(sidechains=1, name='Amide coupling')`. An unnamed reaction is represented
as `Reaction(sidechains=1)`.

## In-place transformations

Use `in_place_transformations(molecule)` when all matching sites should be
transformed in the same molecule. It returns a Python boolean, not a count of
transformations. The example replaces both chlorine atoms, then checks that a
second call finds no remaining matches.

This example reads [chlorine_to_fluorine.textproto](../../src/nanobind/examples/reactions/chlorine_to_fluorine.textproto).
The `atom: 0` in `change_element` means the first atom of the scaffold query.
For reactions with sidechains, this method uses the first stored reagent from
each sidechain; it does not enumerate a reagent library. Use the iterator
workflow below for enumeration.

[Run this example](../../src/nanobind/examples/reactions/in_place.py).

```python
"""Replace both chlorine atoms in one molecule; check the resulting structure."""
from pathlib import Path

from lillymol import MolFromSmiles, Reaction


def main():
    reaction = Reaction()
    reaction_file = Path(__file__).with_name("chlorine_to_fluorine.textproto")
    assert reaction.read(str(reaction_file))
    molecule = MolFromSmiles("ClCCCCl dichloropropane")
    assert reaction.in_place_transformations(molecule)
    assert molecule.unique_smiles() == MolFromSmiles("FCCCF").unique_smiles()
    assert not reaction.in_place_transformations(molecule)  # No chlorine remains.
    print(reaction)
    print(molecule.smiles(), molecule.name())


if __name__ == "__main__":
    main()
```

## Multiple reactive sites

A molecule with several reactive sites can generate several products, and the
choice of how to handle those sites should be explicit. Both `trxn` and the
Python bindings provide controls for this. In Python, distinguish the scaffold's
matches from the sidechain reagents' matches:

- For scaffold sites, loop over `substructure_search_matches()` and generate a
  separate product for each embedding, as shown below. In contrast,
  `in_place_transformations()` applies the transformations to the same molecule.
- For sidechain sites, set
  `conditions.set_make_new_reagent_for_each_hit(1)` **before loading reagents**
  to store a separate alternative for each reagent embedding. The reaction
  iterator then enumerates those alternatives. Without suitable match conditions,
  a reagent with multiple matches can cause loading to fail.

If the intended workflow requires a single reactive site per molecule, filtering
the input first can simplify enumeration. For example:

```shell
tsubstructure -s '1[OH]-C=O' -m single_acid -o smi acid.smi
```

The leading `1` requires exactly one match to the acid query. Matching molecules
are written to `single_acid.smi`; molecules with no matching group or more than
one are excluded. This pattern identifies a protonated carboxylic acid, so the
input's protonation state matters. Use a query consistent with the corresponding
reaction's scaffold or sidechain query: one acid group does not imply that the
molecule has only one reactive group of every kind.

A query embedding is not always a chemically distinct reactive site. Symmetry
and alternative matches within the same functional group can yield equivalent
products. Match controls can reduce redundant embeddings; product deduplication
using `unique_smiles()` can remove duplicate structures. Limiting a search to
one match with `set_max_matches_to_find(1)` selects one match; it does **not**
establish that the molecule has only one reactive site.

## One product per scaffold embedding

Use explicit embeddings when each matching site should generate a separate
product from the original scaffold. This example constructs a SMIRKS reaction
and changes one of two chlorine sites at a time.

`substructure_search_matches(scaffold)` returns a list of lists of atom indices,
or `None` when there is no match. Pass each match directly to
`perform_reaction`; an existing `Set_of_Atoms` embedding is also accepted. The
product is a new `Molecule`, or `None` on
failure. Do not modify the scaffold between finding and using its embeddings:
those indices refer to that scaffold's atom numbering.

Different embeddings can produce identical molecules. This method does not
implicitly deduplicate products; use `unique_smiles()` as a key if desired.

[Run this example](../../src/nanobind/examples/reactions/per_embedding.py).

```python
"""Generate separate products from the two chlorine sites of one scaffold."""
from lillymol import MolFromSmiles, Reaction


def main():
    reaction = Reaction()
    assert reaction.construct_from_smirks("[Cl:1]>>[F:1]")
    reaction.set_name("One chlorine to fluorine")
    scaffold = MolFromSmiles("ClCCC(Cl)C dichloropentane")
    original = scaffold.unique_smiles()
    matches = reaction.substructure_search_matches(scaffold)
    assert matches is not None
    assert len(matches) == 2
    products = []
    for match in matches:
        product = reaction.perform_reaction(scaffold, match)
        assert product is not None
        product.set_name(scaffold.name() + " " + reaction.name())
        products.append(product.unique_smiles())
        print(product.smiles(), product.name())
    expected = {MolFromSmiles(s).unique_smiles()
                for s in ("FCCC(Cl)C", "ClCCC(F)C")}
    assert set(products) == expected
    assert scaffold.unique_smiles() == original
    assert reaction.substructure_search_matches(MolFromSmiles("CCO")) is None


if __name__ == "__main__":
    main()
```

## Explicitly supplied reagents

For a single combination of reagents, `perform_reaction([scaffold, reagent, ...])`
is convenient. The scaffold must be first, followed by one reagent for each
sidechain in sidechain order. The result is a `Molecule`, or `None` on failure.
This overload uses the first scaffold embedding; use explicit embeddings for
control over the reaction site.

The example constructs a reaction from the contents of
[amide.textproto](../../src/nanobind/examples/reactions/amide.textproto), then joins
acetic acid and methylamine. The scaffold query numbers carbonyl carbon as atom
0, oxygen as atom 1, and hydroxyl oxygen as atom 2. It removes atom 2 and joins
scaffold atom 0 to sidechain query atom 0 (the amine nitrogen).

This convenience overload clears stored sidechain reagents before loading its
arguments and clears them after successful product generation. Use a separate
reaction object for this workflow and for preloaded-library enumeration. The
input reagent structures are copied into the reaction and remain unchanged.

[Run this example](../../src/nanobind/examples/reactions/explicit_reagents.py).

```python
"""Form one amide from an explicitly supplied scaffold and sidechain reagent."""
from pathlib import Path

from lillymol import MolFromSmiles, Reaction


def main():
    reaction = Reaction()
    text = Path(__file__).with_name("amide.textproto").read_text()
    assert reaction.construct_from_textproto(text)
    assert reaction.number_sidechains() == 1
    scaffold = MolFromSmiles("CC(=O)O acetic_acid")
    amine = MolFromSmiles("CN methylamine")
    original = [m.unique_smiles() for m in (scaffold, amine)]
    product = reaction.perform_reaction([scaffold, amine])
    assert product is not None
    assert product.unique_smiles() == MolFromSmiles("CC(=O)NC").unique_smiles()
    assert [m.unique_smiles() for m in (scaffold, amine)] == original
    assert reaction.number_sidechains_with_reagents() == 0
    print(product.smiles())


if __name__ == "__main__":
    main()
```

## Enumerating stored reagent libraries

Load sidechain reagents once, then use `ReactionIterator` to select each stored
combination. Sidechains and reagent indices are zero-based. A sidechain is a
reaction component; it can hold many alternative reagent molecules.

The example loads [amines.smi](../../src/nanobind/examples/reactions/amines.smi)
into sidechain 0, processes two acids, and skips a third molecule that does not
match. It checks that the four expected products are generated.

Compute matches inside the scaffold loop. Reset the iterator for each matching
scaffold and call `increment()` after processing each reagent combination.
`reagent_names(iterator)` returns names in sidechain order. With several
sidechains, the iterator traverses combinations of their stored reagents. Keep
the reaction and its reagent libraries alive and unchanged during iteration;
recreate the iterator after adding or removing reagents.

[Run this example](../../src/nanobind/examples/reactions/enumerate_reagents.py).

```python
"""Enumerate two acids against two stored amines; skip a nonmatching scaffold."""
from pathlib import Path

from lillymol import (MolFromSmiles, Reaction, ReactionIterator,
                     SidechainMatchConditions)


def main():
    directory = Path(__file__).parent
    reaction = Reaction()
    assert reaction.read(str(directory / "amide.textproto"))
    conditions = SidechainMatchConditions()
    count = reaction.add_sidechain_reagents(
        0, str(directory / "amines.smi"), conditions)
    assert count == 2
    assert reaction.sidechain_name(0, 0) == "methylamine"
    assert reaction.number_sidechains_with_reagents() == 1
    iterator = ReactionIterator(reaction)
    products = []
    for smiles in ("CC(=O)O acetic_acid", "CCC(=O)O propionic_acid", "CCC no_acid"):
        scaffold = MolFromSmiles(smiles)
        matches = reaction.substructure_search_matches(scaffold)
        if matches is None:
            continue
        iterator.reset()
        while iterator.active():
            names = reaction.reagent_names(iterator)
            for match in matches:
                product = reaction.perform_reaction(scaffold, match, iterator)
                assert product is not None
                products.append(product.unique_smiles())
                print(product.smiles(), scaffold.name(), *names)
            iterator.increment()
    expected = {MolFromSmiles(s).unique_smiles() for s in
                ("CC(=O)NC", "CC(=O)NCC", "CCC(=O)NC", "CCC(=O)NCC")}
    assert len(products) == 4
    assert set(products) == expected
    assert reaction.remove_all_reagents() == 1  # Number of sidechains, not reagents.
    assert reaction.number_sidechains_with_reagents() == 0


if __name__ == "__main__":
    main()
```

## Loading and matching sidechain reagents

```python
conditions = SidechainMatchConditions()
count = reaction.add_sidechain_reagents(0, "amines.smi", conditions)
count = reaction.add_sidechain_reagents(
    0, "amines.dat", conditions, file_type=FileType.SMI
)
```

`file_type` is keyword-only and defaults to `None`, which triggers inference
from the filename. An explicit type overrides the suffix. The return value is
the number of stored reagents added by this call, excluding previous loads;
`0` indicates that no reagents were added or loading failed. This is a reagent
count, not the number of sidechain components. Multiple embeddings of a reagent
can yield several stored alternatives when requested by the match conditions.

`add_sidechain_reagent(sidechain, molecule, conditions)` adds a single molecule
and returns a boolean. The reaction stores its own copy. Loading can fail if a
reagent does not match the sidechain query or has unacceptable multiple matches.
A failed load is not transactional: earlier reagents from the file can remain.

`SidechainMatchConditions` controls reagent matching; it does not configure
scaffold matching. A reaction inherits `SubstructureQuery`, so scaffold controls
such as `set_find_unique_embeddings_only(1)` and `set_max_matches_to_find(n)` apply
to the reaction itself. Symmetry and matching settings affect how many products
are generated. Neither matching nor reaction application performs chemical
standardisation automatically.

### SidechainMatchConditions Methods

| Method | Effect |
| --- | --- |
| `set_make_new_reagent_for_each_hit(value)` | Store separate alternatives for multiple reagent embeddings. |
| `set_max_matches_to_find(n)` | Limit matches searched for each reagent. |
| `set_strip_reagents_to_largest_fragment(value)` | Reduce reagent copies to their largest fragment. |
| `set_ignore_not_reacting(value)` | Skip nonmatching reagents instead of failing the load. |
| `set_find_unique_embeddings_only(value)` | Request unique reagent embeddings. |
| `set_one_embedding_per_start_atom(value)` | Limit reagent embeddings by their starting atom. |
| `set_ignore_symmetry_related_matches(value)` | Ignore symmetry-related reagent matches. |

## Reaction Methods

This table summarizes the entry points used in the examples. Always check
`None` or the documented success result before accessing a product.

| Method | Result and behavior |
| --- | --- |
| `read(fname)` | `bool`; read a textproto reaction file. |
| `read_proto(fname)` | `bool`; alias for the reaction-specific `read(fname)`. |
| `construct_from_textproto(text)` | `bool`; construct from a textproto string. |
| `construct_from_smirks(text)` | `bool`; construct from SMIRKS. |
| `name()` | Reaction name/comment. |
| `set_name(text)` | Assign the reaction name/comment. |
| `number_sidechains()` | Number of sidechain components, including empty ones. |
| `number_sidechains_with_reagents()` | Number of sidechain components holding reagents. |
| `sidechain_name(sidechain_number, reagent_number)` | Reagent name, or `None` for an unavailable reagent. |
| `add_sidechain_reagents(sidechain, fname, match_conditions, *, file_type=None)` | Number of reagents added by this call; `0` on failure or no additions. |
| `add_sidechain_reagent(sidechain, mol, match_conditions)` | `bool`; add a molecule to one sidechain. |
| `reagent_names(iterator)` | Names of the currently selected reagents, in sidechain order. |
| `remove_all_reagents()` | Delete stored reagents; return the number of sidechain components. |
| `substructure_search_matches(mol)` | List of scaffold embeddings, or `None` if unmatched. |
| `in_place_transformations(mol)` | `bool`; modify the input molecule at its matching sites. |
| `perform_reaction(scaffold, embedding)` | New `Molecule` or `None`; accepts a list of atom indices or `Set_of_Atoms`. |
| `perform_reaction(scaffold, embedding, iterator)` | New `Molecule` or `None`; accepts a list or `Set_of_Atoms` and uses the iterator's reagent combination. |
| `perform_reaction(reagents)` | New `Molecule` or `None`; scaffold first, then sidechain reagents; uses first scaffold embedding. |

Other overloads accept an output `Molecule`, a single sidechain molecule, or a
list of sidechain molecules. These differ in matching requirements, return types,
and stored-reagent handling. In particular, the `(scaffold, embedding, sidechains)`
and `(scaffold, sidechains)` forms clear stored reagents after product generation;
the latter requires exactly one scaffold embedding. `perform_reaction_to_list`
returns a list, can return a partial list on failure, and uses only the first
scaffold embedding despite its name. Prefer the explicit-embedding iterator
workflow when every scaffold site and every stored reagent combination matters.

In-place operations modify the input and are not transactional on failure.
A false result alone does not distinguish no match from a transformation error;
search first if that distinction matters. Returning a product does not imply
that it has been standardised, deduplicated, or checked for your downstream
chemical constraints.

## ReactionIterator Methods

| Method | Result and behavior |
| --- | --- |
| `initialise(reaction)` | Initialize the iterator from the reaction's stored reagents. The constructor also accepts a reaction. |
| `active()` | Whether a reagent combination remains. |
| `reagent(sidechain)` | Zero-based selected reagent index for a sidechain. |
| `increment()` | Advance to the next combination; use this rather than C++ `iter++`. |
| `reset()` | Return to the first combination; use only with populated sidechains. |

## Building and testing

Build and stage the nanobind modules as described in
[LillyMol Python](LillyMolPython.md#using-lillymol-python). From `src`, run:

```shell
bazel test -c opt --config=native \
  //nanobind:reaction_in_place_example_test \
  //nanobind:reaction_per_embedding_example_test \
  //nanobind:reaction_explicit_reagents_example_test \
  //nanobind:reaction_enumerate_reagents_example_test
```

These tests execute the same scripts linked above, using the Bazel-built binding.
