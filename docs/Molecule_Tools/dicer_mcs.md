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
mapping, then uses isotope 1 to distinguish attachment points in fragment
identity. Repeated increments at the same atom are supported.

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

## Pairwise output

`-Y pairs.smi` additionally writes pairwise matches. Parents with identical names
are skipped. Each pair reports the largest matching supplied fragment, with
SMILES, parent name, `LHS` or `RHS`, matched atom count, and matched fraction of
that parent. Equal-sized candidates use unique SMILES to break ties.

`-P UST:Z` (or another LillyMol atom typing specification) enables greedy atom
expansion for pairwise matches. The seed pair is written first; if expansion
adds atoms, another pair is written with the expanded labels, counts, and
fractions. This heuristic has no backtracking and does not change the global
fragment selection. A pair without a shared fragment emits no pairwise records.

Invalid minimum sizes, invalid isotope mappings, unsupported sizes, and output
open failures return nonzero status.
