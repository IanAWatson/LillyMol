# dicer_mcs

`dicer_mcs` finds shared fragments in dicer textproto output and labels their
corresponding atoms in the parent molecules. It searches the fragments supplied
by dicer plus each eligible whole parent as a zero-cut candidate; it does not
perform an unrestricted maximum common substructure search.

The maximum common substructure (MCS) problem is known to be hard and existing solutions
are notoriously slow. `dicer_mcs` is not an MCS tool. It is a scalable way of
exploring common fragments in a set of molecules, using `dicer` generated fragments.
By default, `dicer` does not break ring bonds, so candidates requiring cuts
through a ring are generally unavailable. Capping can nevertheless expose
individual ring atoms as fragment endpoints. Other bonds excluded by dicer's
cutting rules can also prevent useful common fragments from being generated.
The combination of `dicer` and `dicer_mcs` is not an MCS solution.

There are also judgement issues. Across a set of molecules the actual MCS
might be a Carbon atom, but this is not useful. Across a set of primary amines,
the MCS might be `N-C`. But there might be larger fragments, that occur less often
that are interesting/informative. The `-R` option to `dicer_mcs` presents a table
of fragments, atom counts, and contributing-parent coverage. Coverage counts
parents that emit a candidate, rather than all parents matching its SMARTS.
This table can be examined to determine the most interesting/useful fragments
in the dataset.

```shell
numbered_smiles -I -n 0 -S numbered molecules.smi
dicer -B bscb -B nbamide -B brcb -B proto -k 2 -I capboth -I inc=100 numbered.smi > fragments.textproto
dicer_mcs -m 2 -R fragments.tsv fragments.textproto > common.smi
```

The algorithm does not require size-sorted input. However, sorting the starting
molecules by increasing atom count **before running dicer** can make the report
more useful: each fragment's first contributing parent will then be one of its
smallest contributing examples. Use [msort](msort.md), for example:

```shell
msort -k natoms molecules.smi > sorted.smi
numbered_smiles -I -n 0 -S numbered sorted.smi
dicer -B bscb -B nbamide -B brcb -B proto -k 2 -I capboth -I inc=100 numbered.smi > fragments.textproto
dicer_mcs -m 2 -R fragments.tsv fragments.textproto > common.smi
```

Apply any structure-changing preprocessing before sorting as well if the sort
should reflect the molecules actually analysed. Sorting is for choosing readable
examples; it does not broaden candidate generation or change support counts.

Perform structure-changing preprocessing before numbering. Parent isotopes must
contain each value from zero through `natoms - 1` exactly once; the atom with
label zero has no explicit isotope in SMILES. Dicer adds 100 at broken-bond
attachment points. The reader removes those increments to recover the parent
mapping and stores cut locations separately as fragment atom numbers. Repeated
increments at the same atom are supported. All isotopes are removed before
computing fragment unique SMILES and canonical atom ranks: attachment locations
do not affect either fragment identity or the correspondence labels written on
parents.

Element caps from dicer's `-I cap` and `-I capboth` modes carry the isotope
marker 10000 in addition to their existing payload. For example, 10117 means a
cap whose payload is 117. The reader first records the cap flag and subtracts
10000, then decodes attachment increments and parent numbering as usual. The
flag is retained separately from attachment metadata, canonical ranks, and
fragment identity. It does not split support or appear in correspondence labels
or highlighted examples. Marked caps must be terminal, with payloads below
10000; cap markers are not valid on the original parent record. Legacy fragments
without markers remain readable, but their cap provenance cannot be recovered
reliably. Only explicitly marked caps receive atomic-number SMARTS as described below;
ordinary fragment atoms retain their atomic symbols.

The fixed attachment increment makes original atom labels of 100 or greater
ambiguous. Parents with more than 100 atoms are therefore rejected explicitly.
Fragments must contain fewer than 100 atoms. For parents below that size, the
fragment-size check is skipped; the exceptional 100-atom-parent path checks
parsed fragment sizes, rather than trusting optional proto atom counts.

## Whole parents as zero-cut fragments

Each parent with fewer than 100 atoms is automatically included as an implicit
fragment formed by making no cuts, subject to `-m`. This also applies when dicer
emitted no fragments for that parent. Original atom labels provide the mapping;
the implicit fragment has no attachment points. If dicer already emitted the
whole structure, it is not added again. Support is still counted once per parent.

This allows the whole of a smaller molecule to match a fragment of a larger one.
For example, ethylamine (`NCC`) can now contribute its complete structure and
match the `NCC` fragment from tert-butylamine (`NC(C)(C)C`), even if dicer only
emitted `N` and `CC` from ethylamine.

A 100-atom parent can still contribute smaller explicit dicer fragments, but is
not added whole: the existing fragment limit is fewer than 100 atoms. With `-v`,
this omission is reported. Parents larger than 100 atoms remain unsupported by
the fixed isotope encoding.

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
fragment in every parent, while `-s 1` selects the largest available fragment
(including eligible whole parents).
No qualifying fragment produces a diagnostic and nonzero exit status.
The threshold affects global selection only; pairwise `-Y` comparisons retain
their existing behavior.

No surviving fragments produces a diagnostic and nonzero exit status, rather
than an empty successful result. A parent without emitted dicer fragments is
otherwise valid and contributes its implicit whole-parent candidate when that
candidate meets the size requirements.

## Parents without the selected candidate fragment

Use `-X filename` to write non-contributing parents as SMILES with their original
names and without atom-number isotopes. A `.smi` suffix is appended unless the
filename already ends in `.smi`; for example, `-X missing` writes `missing.smi`.
Collision checks use this final filename.

```shell
dicer_mcs -s 80% -X missing.smi fragments.textproto > common.smi
```

Membership is determined by whether a parent contributes the selected candidate,
through an emitted dicer fragment or its implicit whole structure, not by a
substructure search. A parent can contain a smaller core but be written to `-X`
because dicer's bond-breaking rules, cut limit, or fragment filters prevented
that fragment from being generated.
Parents without any surviving fragments are included in `-X` as well.

Records retain input order. If every parent contributes, the file is empty.
If no fragment survives or meets minimum support, all parents are written to
`-X` and the command still returns nonzero, with its existing no-result
diagnostic. Invalid input prevents analysis and is not treated as a non-match.

`-X` requires a separate file (`-` is rejected because stdout carries matching
parents). Input files cannot be used as `-X` or `-Y` outputs, and the two output
paths must differ; aliases through symlinks and hard links are checked. Output
open/write failures return nonzero status.

## Fragment report

Use `-R filename` for a regular, headered report of every candidate surviving
`-m`, including whole-parent candidates and candidates below minimum support.
The filename is used as specified, without adding a suffix. The default
separator is a tab; `-R sep=tab` makes this explicit. Other named separators are
parsed by LillyMol's `char_name_to_char`, for example `comma`, `space`, `vbar`,
and `semic` (semicolon), or supply a literal single character.

```shell
dicer_mcs -s 80% -R fragments.tsv fragments.textproto > common.smi
dicer_mcs -s 80% -R fragments.csv -R sep=comma fragments.textproto > common.smi
```

| Column | Meaning |
|---|---|
| `SMILES` | Isotope-free canonical candidate structure. |
| `NATOMS` | Candidate atom count. |
| `NPARENTS` | Number of distinct contributing input parents. |
| `FPARENTS` | `NPARENTS` divided by all input parents, formatted to four decimal places with trailing zeros removed. |
| `FIRST_PARENT` | Full name of the first contributing parent in input order. |
| `SMARTS` | Relaxed common-core query, allowing external substituents. |
| `NWHOLE` | Parents whose whole structure contributes this candidate. |
| `NCUT` | Parents with at least one occurrence marked with a cut attachment. |
| `SELECTED` | One for the selected global core, zero otherwise. |
| `ATTACHMENT_SMILES` | One exemplar with cut attachment atoms labelled isotope 1. |
| `EXAMPLE_SMILES` | Complete `FIRST_PARENT` structure, with isotope 1 on all atoms of the exemplar fragment and zero elsewhere. |

`NWHOLE` and `NCUT` are independent per-parent counts, not occurrence counts.
Repeated cut fragments within a parent increase neither `NPARENTS` nor `NCUT`.
These columns allow whole-parent-only candidates to be distinguished from cores
also produced by fragmentation, without removing useful candidates from the
report. The attachment SMILES uses the deterministically first occurrence in
`FIRST_PARENT`. Its labels describe that exemplar only; attachment patterns may
differ in other occurrences and are not part of candidate identity.

`EXAMPLE_SMILES` highlights one deterministic occurrence in `FIRST_PARENT`,
using the stored atom correspondence without another substructure search. With
capped fragments, the original parent atoms represented by caps are highlighted;
the example retains their original substituents and bonds. Whole-parent
candidates highlight every atom. The example is made from a copy with the input
numbering isotopes cleared, so report generation does not alter other output.
Coverage is contributing-parent coverage, not unrestricted substructure-match
coverage. Sorting the starting molecules by atom count as described above is
optional and helps select smaller examples.

Rows are sorted by decreasing `NPARENTS`, decreasing `NATOMS`, then increasing
canonical SMILES. Fields containing the separator, quotes, or line breaks are
quoted with double quotes; embedded quotes are doubled. Thus the same quoting
convention works for TSV, CSV, and the other separators. A quote or a line-break
character cannot itself be the separator.

If no candidate meets minimum support, the report is still written with all
`SELECTED` values zero; the command retains its nonzero no-result status. If no
candidates survive `-m`, the report contains only its header. `-R` cannot share
an input file or the final `-X`/`-Y` output path, including file aliases, and `-`
is rejected because stdout is reserved for labelled parents.

The default stderr output now reports the number of unique candidates and the
selected core instead of dumping every fragment. `-v` retains the compact
best-per-atom-count summary. Structural matches in stdout and pairwise output
are unaffected by requesting a report.

## Label examples for a chosen fragment

After examining the `-R` table, copy a fragment from its `SMILES` column and
select its contributing parents with `-F`:

```shell
dicer_mcs -F 'NCC' fragments.textproto > NCC_examples.smi
```

The supplied SMILES must be connected and contain fewer than 100 atoms. It is
canonicalized without isotopes, so equivalent spellings such as `CCN` select
the same identity. This is an exact candidate-identity
lookup, not a substructure search of the original parents. It considers all
recorded occurrences and the implicit whole-parent candidate. A parent that
contains NCC but did not contribute it as a candidate is not selected.

The default writes each qualifying parent once, in input order. All atoms in
all qualifying occurrences receive isotope 1; other isotopes are cleared. The
original parent name is retained and `nmatch=N` is appended. Occurrences are
counted as **distinct sets of original parent atoms**, so repeated fragment
records, capped/uncapped forms and symmetry-related mappings of the same atom
set count once. Overlapping occurrences are allowed. The highlighted union can
contain more atoms than the requested fragment, and `nmatch` preserves the fact
that the result is ambiguous even when labels do not distinguish the overlaps.

To write one parent per distinct occurrence instead, add `-F separate`:

```shell
dicer_mcs -F 'NCC' -F separate fragments.textproto > NCC_occurrences.smi
```

Each record labels just one occurrence and appends both `nmatch=N` and
`match=i/N`. Occurrences are ordered by their sorted parent atom numbers,
independently of the order of fragment records.

### Require particular attachment positions

A positive isotope in the supplied SMILES means that atom must be an attachment
point in a contributing occurrence. These requirements do not change the
isotope-free structural identity being selected. For example:

```shell
# Require the middle carbon to be an attachment point.
dicer_mcs -F 'N[1CH2]C' fragments.textproto > middle_attachment.smi
# Require the terminal carbon to be an attachment point.
dicer_mcs -F 'NC[1CH3]' fragments.textproto > terminal_attachment.smi
```

Unlabelled atoms remain unconstrained and can also be attachment points. Any
positive isotope value acts as the same requirement; it is not an original atom
number, required isotope value, or required number of cuts. Use ordinary valid
SMILES with the appropriate bracketed hydrogen counts.

Requirements are checked against each recorded fragment's decoded attachment
metadata, before deduplicating its parent atom set. Equivalent atoms can exchange
positions when testing the requirements; the chosen representative from the
report does not limit which occurrences qualify. Numbering isotopes alone do
not identify attachment points: the usual `-I inc=100` output supplies that
information. Cap markers are decoded separately. When a cap qualifies, its
original parent atom is highlighted, retaining the parent's actual bonds and
substituents. Attachment requirements can reduce coverage relative to the
report row, which combines all attachment patterns.

### Streaming and command options

This mode processes one parent at a time, without retaining the dataset or
ranking fragments. Multiple input files are processed in command-line order;
`-` reads textproto records from standard input. It writes SMILES to stdout and
reports contributing-parent and distinct-occurrence counts on stderr. No
contributing parent gives nonzero exit status. Invalid records and input/output
errors also give nonzero status; already written records remain partial output.

Specify one fragment SMILES, optionally `-F separate`. Explicit selection cannot
be combined with `-m`, `-s`, `-P`, `-R`, `-X` or `-Y`; these belong to the global
analysis mode. The normal `-A`, `-E` and `-v` options remain available. Output
labels are uniformly isotope 1, including when different positive isotopes were
used to specify attachment requirements.

## Attachment-independent common cores

A bond cut in one parent can correspond to a hydrogen in another. Fragments with
the same isotope-free structure therefore share support even when different
atoms were cut, or no atoms were cut in one of the fragments. This intentionally
broadens matching compared with earlier versions, which included attachment
labels in fragment identity. Candidate identity retains element, aromaticity,
bond order, and the existing stereochemical representation.

Generated SMARTS require at least the number of connections present within the
fragment at every atom. Extra external substituents are permitted anywhere;
attachment annotations from one exemplar do not constrain other parents. A
single-atom fragment has no degree restriction. `D>0` is omitted for singly
connected atoms because the query bond already requires a connection. For example,
an uncapped ethane fragment produces `[C][C]`.

Only explicitly marked **cap atoms** use atomic-number SMARTS, such as `[#6]`
and `[#7]`, allowing either aromatic or aliphatic atoms of that element. Ordinary
fragment atoms retain their symbols: a genuine terminal primary-amine nitrogen
is `[N]`, not `[#7]`. This distinction uses the stored 10000 isotope marker;
it cannot be recovered reliably from older unmarked cap files.

For example, an NCC candidate with a genuine nitrogen, genuine internal carbon,
and marked carbon cap produces `[N][CD>1][#6]`. Its terminal cap can match the
aromatic attachment carbon of benzylamine, while the internal carbon must still
be aliphatic. The report SMARTS describes the chosen occurrence in `FIRST_PARENT`.
When occurrences map to exactly the same parent atom subset, the one with fewer
caps is preferred, retaining more genuine atom specifications. A nitrogen that
actually is a cap still uses `[#7]`.

This change affects generated SMARTS in both stderr and the report, not fragment
identity, support counts, selection or pairwise comparison. Dummy and
non-periodic-table elements retain symbolic matching.

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
fragment size with an induced occurrence correspondence, then tries every
occurrence pairing for all shared fragment identities of that size. Each trial starts with cleared labels. The winner has
the greatest expanded atom count; ties use fragment unique SMILES, then the
lexicographically smallest list of `(lhs atom number, rhs atom number)` pairs.
Occurrence lists are ordered by their parent atom numbers in canonical fragment
rank order, making the chosen seed reproducible when fragment records are
reordered. Without `-P`, all trials have the seed size and the same deterministic
tie-breaking applies.

Only the winning seed pair is written; if its expansion adds atoms, another
pair is written with the expanded labels, counts, and fractions. Smaller shared
seeds are not tried once a valid seed size is found, even if one might yield a
larger expansion. If all occurrence pairings of a larger seed violate induced
connectivity, the search proceeds to smaller shared seeds. No disconnected
unions of seeds are formed. Each expansion remains greedy with no backtracking
or alternative symmetry mappings, and does not change global fragment selection.
Global output uses the first deterministically ordered occurrence of its selected
fragment. A pair without a shared fragment emits no pairwise records.

The number of trials grows as the product of the occurrence counts in the two
parents for each largest shared fragment identity. Highly repetitive molecules
can therefore increase run time. There is currently no cap on trial counts.

### Induced connectivity

Every pairwise seed and expansion must preserve the graph induced by the
matched parent atoms. For every two corresponding atom pairs, a bond must be
present in both parents or absent in both. Existing bonds must have the same
type; aromatic bonds are compared as aromatic regardless of their particular
Kekule single/double assignment. Thus an expansion cannot silently introduce or
lose a ring closure, or change the order of a bond to an already matched atom.

Seed validation matters when dicer has broken a ring bond while retaining both
endpoint atoms. Such a fragment can have the same fragment SMILES in two
parents whose induced connectivity differs. Incompatible occurrence pairings
are skipped, including without `-P`. If none are valid at a given seed size,
smaller shared fragments are considered; no valid seed gives no pairwise output.

During expansion, incompatible neighbours are skipped and the next available
neighbour is tried. The search is still greedy: accepted additions are not
backtracked, and no alternative symmetry mappings are generated. These checks
concern bond connectivity and type, not a new stereochemical matching policy.
Global fragment support and selection remain unchanged; the induced constraint
applies to pairwise correspondences, not to the set-wide fragment lookup.

Invalid minimum sizes, invalid isotope mappings, unsupported sizes, and output
open failures return nonzero status.

The executable regression suite uses `tsubstructure` alongside `dicer_mcs` to
verify that generated SMARTS match contributing parents and reject incompatible
atom and bond types. Build/install both tools when running these tests.

## Developer notes

See [candidate coverage and support limitations](dicer_mcs_developer_notes.md)
for an investigation of small amines and a proposed fragment-extension approach.
