# medchemwizard

## Overview

`medchemwizard` generates analogues by applying a configured set of medicinal
chemistry transformations to each input molecule. The standard LillyMol
configuration contains more than 200 reactions under
`${LILLYMOL_HOME}/data/MedchemWizard`. Some transformations are additive, such
as adding a methyl group, while others perform molecular rearrangements or
isostere-style replacements.

The tool writes every product that passes filtering. It does not rank products,
estimate synthetic accessibility, or choose a preferred transformation. The
output is intended for downstream filtering, scoring, database lookup, or
inspection.

The original transformation set was derived from medicinal chemistry transform
work from Abbott Laboratories, but the original literature reference is not
currently known. The implementation should therefore be treated as a practical
LillyMol transformation library rather than a formally cited reproduction.

## Quick Start

Apply the standard reaction set and append the transformation name to each
product record:

```shell
medchemwizard \
  -R "${LILLYMOL_HOME}/data/MedchemWizard/REACTIONS" \
  -W space input.smi > products.smi
```
Specifying the default reaction file is built into the convenience wrapper medchemwizard.sh
found in ${LILLYMOL_HOME}/contrib/bin

```shell
medchemwizard.sh -W space input.smi > products.smi
```

Input is any LillyMol-supported molecule format. SMILES input is commonly used:

```text
CCO ethanol
```

Output is SMILES followed by the product name. With `-W space`, the reaction
name is appended to the parent name as another space-separated token:

```text
CCCO ethanol add_propyl
```

The exact product set depends on the reaction files in the selected
configuration.

## Reaction Configuration

`-R <fname>` specifies a file containing the reactions to apply. This option is
required by the command-line tool. The standard file is:

```text
${LILLYMOL_HOME}/data/MedchemWizard/REACTIONS
```

Each non-comment record in that file names a reaction file in the same
directory. Records beginning with `PROTO:` name textproto/json reaction files;
other records name the older LillyMol reaction format.

A custom configuration can be built by making a different reaction-list file
that contains only the reaction files desired for a particular workflow. This is
often preferable to applying the full standard set when only a narrow class of
changes is wanted.

## Product Names

Use `-W <sep>` to append the reaction name to each product name:

```shell
medchemwizard.sh -W space input.smi > products.smi
```

The separator is passed through LillyMol's standard character-name conversion,
so `space` is a convenient way of producing a three-token output suitable for
text processing:

```text
product_smiles parent_id reaction_name
```

Without `-W`, products retain the parent molecule name.

## Recursion

By default, each transformation is applied only to the input molecule. Use
`-D <depth>` to recursively apply transformations to generated products:

```shell
medchemwizard.sh \
  -D 1 -W space input.smi > recursive_products.smi
```

Recursive generation can increase output size very rapidly. Size filters,
uniqueness filters, and maximum-hit limits are usually advisable when using
`-D`.

## Controlling Product Count

Use atom-count limits to discard products outside the useful size range:

```shell
medchemwizard.sh -m 8 -M 45 -W space input.smi > products.smi
```

Use `-x <n>` to limit how many embeddings of any one reaction are processed per
molecule:

```shell
medchemwizard.sh -x 5 -W space input.smi > products.smi
```

By default, the tool reports when reaction matches are truncated. Add `-x quiet`
to suppress those warnings.

Use `-U each` to suppress duplicate products generated from the same parent
molecule, or `-U all` to suppress duplicate products across the entire run:

```shell
medchemwizard.sh -U each -W space input.smi > products.smi
```

Uniqueness is based on LillyMol unique SMILES after any requested product
postprocessing.

## Preprocessing and Product Filtering

Input molecules can be preprocessed with common LillyMol options:

```shell
medchemwizard.sh -l -c -g all -W space input.smi > products.smi
```

Important options are:

| Option | Meaning |
| --- | --- |
| `-l` | Reduce input molecules to the largest fragment before transformation. |
| `-c` | Remove all chirality from input molecules. |
| `-g ...` | Apply standard LillyMol chemical standardisation to input molecules. |
| `-y` | Apply the same preprocessing to products before uniqueness and atom-count checks. |
| `-V .` | Discard products with invalid valence. |
| `-V <fname>` | Discard products with invalid valence and write them to `<fname>`, adding `.smi` if needed. |

`-A`, `-E`, and `-i` are the standard LillyMol aromaticity, element, and input
type options. If no aromaticity option is specified, Daylight aromaticity is
used.

## Protecting Parts of the Molecule

Use `-s <smarts>` or `-q <query>` to mark atoms that must not be changed. Any
reaction embedding that touches a protected atom is skipped. This is useful when
a known pharmacophore, binding motif, or other important substructure should be
kept fixed while transformations are explored elsewhere.

```shell
medchemwizard.sh \
  -s '[CX3](=O)[OH]' -z i -W space input.smi > products.smi
```

If protected-atom queries are specified, the default behaviour is to fail when
none of them match an input molecule. Use `-z i` or `-z ignore` when it is
acceptable for a molecule to contain no protected atoms.

The protection currently applies to atoms in the reaction embedding. It is not a
full bond-change or radius-based protection model. If neighbouring atoms must
also be protected, include them in the SMARTS or query definition.

## Option Summary

| Option | Meaning |
| --- | --- |
| `-R <fname>` | Reaction-list file. Required for the CLI. |
| `-D <depth>` | Recursively transform generated products to `depth`; default is no recursion. |
| `-x <n>` | Use no more than `n` matches to any reaction in a molecule. |
| `-x quiet` | Suppress warnings when matches are truncated. |
| `-U each` | Suppress duplicate products within each parent molecule. |
| `-U all` | Suppress duplicate products across the whole run. |
| `-m <natoms>` | Discard products with fewer than `natoms` atoms. |
| `-M <natoms>` | Discard products with more than `natoms` atoms. |
| `-W <sep>` | Append reaction names to product names using separator `sep`. |
| `-q <query>` | Protect atoms matching a LillyMol query file. |
| `-s <smarts>` | Protect atoms matching SMARTS. |
| `-z i`, `-z ignore` | Ignore molecules that do not match any protected-atom query. |
| `-V .` | Discard products with bad valences. |
| `-V <fname>` | Write discarded bad-valence products to a SMILES file. |
| `-y` | Postprocess generated products as input molecules are processed. |
| `-l` | Reduce input molecules to largest fragment. |
| `-c` | Remove chirality from input molecules. |
| `-g ...` | Chemical standardisation. |
| `-i <type>` | Input type. |
| `-A ...` | Aromaticity specification. |
| `-E ...` | Element specification. |
| `-v` | Verbose output, including per-reaction hit counts. |

## Python Interface

`MedchemWizard` is also available from the `lillymol_tools` Python module. The
Python object stores the reaction set and options, and `process(mol)` returns a
list of generated `Molecule` objects. The input molecule is copied internally
and is not modified.

```python
from lillymol import MolFromSmiles
from lillymol_tools import MedchemWizard

mol = MolFromSmiles("CCO ethanol")
if mol is None:
    raise ValueError("invalid SMILES")

wizard = MedchemWizard()
wizard.initialise_from_environment()
wizard.set_append_names(True)
wizard.set_name_separator(" ")
wizard.set_max_atoms(25)

products = wizard.process(mol)
for product in products:
    print(product.smiles(), product.name())

print(wizard.stats())
```

`initialise_from_environment()` reads:

```text
${LILLYMOL_HOME}/data/MedchemWizard/REACTIONS
```

For custom reaction sets, use `read_reactions(fname)`.

```python
wizard = MedchemWizard()
wizard.read_reactions("/path/to/REACTIONS")
```

Available Python methods are:

| Method | Meaning |
| --- | --- |
| `initialise_from_environment()` | Load the standard reaction set from `LILLYMOL_HOME`. |
| `read_reactions(fname)` | Load a reaction-list file. |
| `number_reactions()` | Return the number of loaded reactions. |
| `process(mol)` | Return generated products as a list of `Molecule` objects. |
| `add_do_not_change_smarts(smarts)` | Protect atoms matching SMARTS. |
| `add_do_not_change_query(fname)` | Protect atoms matching a LillyMol query file. |
| `set_ignore_do_not_change_queries_not_matching(value)` | Allow molecules with no protected-query match. |
| `set_max_depth(value)` | Set recursive transformation depth. |
| `set_max_hits(value)` | Set maximum matches processed for each reaction. |
| `set_max_atoms(value)` | Set maximum product atom count. |
| `set_min_atoms(value)` | Set minimum product atom count. |
| `set_unique_within_molecule(value)` | Suppress duplicate products per parent. |
| `set_unique_across_all_molecules(value)` | Suppress duplicate products across the object lifetime. |
| `set_append_names(value)` | Append reaction names to product names. |
| `set_name_separator(sep)` | Set the product-name separator used with appended names. |
| `set_discard_bad_valences(value)` | Discard products with invalid valence. |
| `stats()` | Return processing counters as a Python dictionary. |

Protected-atom behaviour mirrors the command-line tool:

```python
wizard.add_do_not_change_smarts("[CX3](=O)[OH]")
wizard.set_ignore_do_not_change_queries_not_matching(True)
products = wizard.process(mol)
```

`stats()` includes counters such as `molecules_read`, `molecules_produced`,
`duplicate_molecules_suppressed`, `bad_valences_discarded`,
`molecules_not_matching_do_not_change_queries`, and
`embeddings_rejected_for_changing_protected_atoms`.
