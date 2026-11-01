# dicer_mcs

`dicer_mcs` finds shared fragments in dicer textproto output and labels their
corresponding atoms in the parent molecules. It searches the fragments supplied
by dicer; it does not perform an unrestricted maximum common substructure search.

```shell
numbered_smiles -I -n 0 -S numbered molecules.smi
dicer -B proto -I inc=100 -k 2 numbered.smi > fragments.textproto
dicer_mcs -m 2 fragments.textproto > common.smi
```

Perform structure-changing preprocessing before numbering. Parent isotopes must
contain each value from zero through `natoms - 1` exactly once; the atom with
label zero has no explicit isotope in SMILES. Dicer adds 100 at broken-bond
attachment points. The reader removes those increments to recover the parent
mapping and stores cut locations separately as fragment atom numbers. Repeated
increments at the same atom are supported. All isotopes are removed before
computing fragment unique SMILES and canonical atom ranks: attachment locations
do not affect either fragment identity or the correspondence labels written on
parents.

The fixed attachment increment makes original atom labels of 100 or greater
ambiguous. Parents with more than 100 atoms are therefore rejected explicitly.
Fragments must contain fewer than 100 atoms. For parents below that size, the
fragment-size check is skipped; the exceptional 100-atom-parent path checks
parsed fragment sizes, rather than trusting optional proto atom counts.

## Selecting a common fragment

`-m N` discards fragments smaller than `N`. Each fragment's support is the number
of distinct input parent records containing it, regardless of how many times it
occurs within a parent. The winner has the greatest support, then the greatest
atom count, then the lexicographically smallest unique SMILES. A shared fragment
need not occur in every input parent. Only parents containing the winning
fragment are written; stderr includes its SMARTS and support.

Use `-s N` to require support in at least `N` parents, or `-s NN%` for a
whole-number percentage from 1 through 100. Percentage thresholds round up to a
parent count, using all input parents as the denominator, including parents
without surviving fragments. For example, `-s 50%` with three parents requires
two. Counts larger than the input set are allowed but cannot yield a result.

```shell
dicer_mcs -m 2 -s 80% fragments.textproto > common.smi
dicer_mcs -m 2 -s 5 fragments.textproto > common.smi
```

With `-s`, the winner is the **largest** fragment meeting the support threshold;
equal-sized fragments are ranked by support and then unique SMILES. Without
`-s`, the original support-first policy remains. `-s 100%` requires a common
fragment in every parent, while `-s 1` selects the largest available fragment.
No qualifying fragment produces a diagnostic and nonzero exit status.
The threshold affects global selection only; pairwise `-Y` comparisons retain
their existing behavior.

No surviving fragments produces a diagnostic and nonzero exit status, rather
than an empty successful result. A parent without fragments is otherwise valid
and does not prevent other parents from contributing.

## Attachment-independent common cores

A bond cut in one parent can correspond to a hydrogen in another. Fragments with
the same isotope-free structure therefore share support even when different
atoms were cut, or no atoms were cut in one of the fragments. This intentionally
broadens matching compared with earlier versions, which included attachment
labels in fragment identity. Element, aromaticity, bond order, and the existing
stereochemical representation are retained.

Generated SMARTS require at least the number of connections present within the
fragment at every atom. Extra external substituents are permitted anywhere;
attachment annotations from one exemplar do not constrain other parents. A
single-atom fragment has no degree restriction. For example, an ethane fragment
produces `[CD>0][CD>0]`, which can match ethane or a carbon-carbon single bond in
a larger molecule.

Attachment metadata is retained for future reporting and expansion work.
Pairwise expansion continues to examine all matched atoms, without requiring a
cut site in both parents. For symmetric fragments, isotope-free canonical ranks
provide a correspondence, but alternative symmetry-related mappings are not
searched to optimise expansion.

## Pairwise output

`-Y pairs.smi` additionally writes pairwise matches. Parents with identical names
are skipped. Each pair reports the largest matching supplied fragment, with
SMILES, parent name, `LHS` or `RHS`, matched atom count, and matched fraction of
that parent. All occurrences of each isotope-free fragment are retained.

`-P UST:Z` (or another LillyMol atom typing specification) enables greedy atom
expansion for pairwise matches. The tool first identifies the largest shared
fragment size, then tries every occurrence pairing for all shared fragment
identities of that size. Each trial starts with cleared labels. The winner has
the greatest expanded atom count; ties use fragment unique SMILES, then the
lexicographically smallest list of `(lhs atom number, rhs atom number)` pairs.
Occurrence lists are ordered by their parent atom numbers in canonical fragment
rank order, making the chosen seed reproducible when fragment records are
reordered. Without `-P`, all trials have the seed size and the same deterministic
tie-breaking applies.

Only the winning seed pair is written; if its expansion adds atoms, another
pair is written with the expanded labels, counts, and fractions. Smaller shared
seeds are not tried, even if one might yield a larger expansion. No disconnected
unions of seeds are formed. Each expansion remains greedy with no backtracking
or alternative symmetry mappings, and does not change global fragment selection.
Global output uses the first deterministically ordered occurrence of its selected
fragment. A pair without a shared fragment emits no pairwise records.

The number of trials grows as the product of the occurrence counts in the two
parents for each largest shared fragment identity. Highly repetitive molecules
can therefore increase run time. There is currently no cap on trial counts.

Invalid minimum sizes, invalid isotope mappings, unsupported sizes, and output
open failures return nonzero status.

The executable regression suite uses `tsubstructure` alongside `dicer_mcs` to
verify that generated SMARTS match contributing parents and reject incompatible
atom and bond types. Build/install both tools when running these tests.
