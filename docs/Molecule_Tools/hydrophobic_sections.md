# hydrophobic_sections

`hydrophobic_sections` describes connected regions of hydrophobic atoms. It can
also describe hydrophilic regions and write molecules with isotopic labels on
selected atoms.

```shell
hydrophobic_sections molecules.smi > descriptors.txt
hydrophobic_sections -G '[#6]' -G '[F,Cl,Br,I]' molecules.smi
hydrophobic_sections -G '[#6]' -G minsize=4 -G iso=7 -S labelled molecules.smi
```

Without `-G`, the existing simple hydrophobicity rules are used. Supplying `-G`
queries replaces those rules; add `-G simple` to combine queries with the rules.
The simple rules exclude carbon atoms with nonaromatic multiple bonds to
heteroatoms, including aldehyde and imine carbons, even when those carbons bear
hydrogen. Nonaromatic carbons attached only to carbon or hydrogen are included,
including substituted alkene and alkyne carbons without attached hydrogens.
Aromatic carbons bearing halogens remain connected to the hydrophobic ring section.
An aromatic carbon attached to an external carbon remains eligible even when
that external carbon belongs to a polar group such as an imine, carbonyl, or
nitrile. The external carbon is classified separately and remains excluded when
multiply bonded to a heteroatom.
Explicit queries can still select atoms excluded by the simple rules.

Repeated queries add atoms to the selection. All atoms in each matching embedding
are selected, so use a single-atom query when only one atom should be selected.
Query files can be supplied with `-G q:filename` or `-G Q:filename`.

Use `-L` for hydrophilic atom queries, or `-x` to treat atoms outside the retained
hydrophobic regions as hydrophilic. `-L minsize=N` sets the minimum hydrophilic
region size; `-G minsize=N` sets the minimum hydrophobic region size. Both default
to one. Smaller regions are excluded from atom counts, section descriptors,
isotope labels, and inter-section calculations. With `-x`, the complement is
formed after the hydrophobic minimum-size filter.

Output includes atom counts and fractions, section sizes, aromatic atom counts,
and other section descriptors. `HPOmxarm` and `HPOavarm` report the maximum and
mean aromatic atom counts per hydrophobic section; the corresponding `HPI`
columns describe hydrophilic sections. Earlier versions incorrectly counted
nonaromatic atoms in these columns. Earlier versions also applied minimum-size
settings only to inter-section calculations.

The reusable C++ `hydrophobic::Sections` class in
`Molecule_Tools/hydrophobic_sections_lib.h` accepts a molecule, an atom-selection
mask, and an optional minimum section size. It returns section membership,
section sizes, total retained atoms, and the largest section size. Atom-selection
rules and preprocessing remain the caller's responsibility. The command-line
tool retains its existing selection rules and descriptor column order.

Neutral, two-connected aromatic sulfur (for example, in thiophene) is hydrophobic
under the simple rules and is excluded from the simple hydrophilic selection.
This rule requires atomic number 16, two connections, aromaticity, and formal
charge zero. It does not include charged aromatic sulfur or oxidised sulfur.

## Restricting analysis to input isotope labels

Use `-I iso` to consider only atoms with a nonzero input isotope, or `-I noiso`
to consider only atoms with no input isotope. These options cannot be combined.
After capturing the selection, all input isotopes are cleared, before queries
or simple rules are applied. Output isotopes come only from `-G iso=` and
`-L iso=`. Without `-I`, existing isotope handling is unchanged.

```shell
hydrophobic_sections -I iso -G simple -G minsize=2 -G iso=1 -S regions input.smi
```

Ignored atoms cannot belong to either type of section, including hydrophilic
sections selected by `-x`. They break section connectivity, and minimum sizes
apply to the resulting restricted sections. Molecular context for atom
classification is retained: ignoring an oxygen does not turn its carbonyl
carbon into a hydrocarbon carbon.

Atom fractions and the size-to-molecule ratios use the number of eligible
atoms as their denominator. Ignored atoms do not contribute to the neither
category. If no atoms are eligible, all descriptors are zero. Distances and
terminal-section detection retain their existing whole-molecule definitions;
ignored atoms are not deleted from the molecular graph.

## Numeric output

Values within `1e-5` of an integer are written as integers. Other descriptor
values use `Fraction_as_String` with three decimal places and trailing zeros
removed, including values greater than one. This replaces the former
three-significant-digit formatting and avoids platform-dependent scientific
notation and formatting differences.
