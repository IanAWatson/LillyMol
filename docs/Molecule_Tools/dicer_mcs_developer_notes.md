# Developer note: candidate coverage and support in dicer_mcs

## Observation

An investigation of 100 small primary amines from Enamine exposed a substantial
difference between a fragment being present as a substructure of a parent and
being available to `dicer_mcs` as an explicitly generated fragment. The input
was sorted by molecule size, which made the first discrepancies easy to inspect.
These results describe this particular dataset, rather than a benchmark or a
general estimate of coverage.

The original fragmentation command was:

```shell
dicer -B bCCD3 -B nbamide -B brcb -B proto -k 2 -z i -c -g all \
  -I inc=100 -M 16 -M maxnr=9 /tmp/p2_num.smi > /tmp/frags.textproto
dicer_mcs -X /tmp/Xfile -m 2 -v /tmp/frags.textproto
```

A substructure search for `[NH2]-[CR0]-[CR0]` matched 94 of the 100 parents.
The query requires two single bonds and two acyclic carbons, but does not require
those carbons to be saturated or free of other substituents.

The following counts came from separate runs on the same numbered input:

| Fragmentation settings | Parents emitting NC | Parents emitting NCC |
| --- | ---: | ---: |
| Original command | 66 | 61 |
| Add `-B bscb`, retain `-k 2` | 90 | 73 |
| Add `-B bscb`, increase to `-k 4` | 95 | 75 |

These are distinct-parent counts for exact isotope-free fragment identities,
including the implicit whole-parent candidates added by `dicer_mcs`. They are
not counts obtained by searching every original parent for the candidate.
The comparison with the 94 query matches concerns NCC; that query does not
provide a corresponding reference count for NC.

## Why propylamine initially did not provide NCC

The bond-selection qualifiers have different meanings:

* `-B brcb` allows ring–chain bonds to break.
* `-B bCCD3` allows C–C single bonds with at least one endpoint having more than
  two connections.
* `-B bscb` allows C–C single bonds generally, subject to other exclusions in
  the bond-selection rules.

Propylamine, `NCCC` (ID:EN300-20505), has neither a ring nor a highly connected
carbon. The original command emitted ammonia and propane, but no NCC fragment.
Adding `-B bscb` makes the required C–C cut available. Increasing `-k` alone
cannot make a bond eligible for cutting.

## Why 19 discrepancies remain with bscb and four cuts

In input order, the first parent hits lacking an emitted NCC fragment included:

| Record (one-based) | Parent | Required operation to isolate NCC |
| --- | --- | --- |
| 3 | `NCC=C` | Cut a double bond |
| 4 | `NCC#C` | Cut a triple bond |
| 14 | `NCC(=C)C` | Cut a double bond, as well as remove another substituent |
| 29 | `NCCF` | Cut a C–F bond |

All 19 remaining discrepancies were explained by these policies: 15 parents
would require cutting a multiple bond and four would require cutting C–F
bonds. Thus the 94 parent query hits and 75 generated NCC occurrences are
consistent with the fragmentation rules. Ring-attached amines also explain
some missing NC/NCC fragments, but are excluded by the acyclic-carbon query
and are not among these 19 discrepancies.

The hard-coded rules consider single, non-ring bonds and exclude bonds to
halogens. Protecting C–F bonds was intended to avoid splitting CF3 groups,
which is generally desirable. The consequence is broader: a fluorinated
carbon cannot contribute as an unadorned carbon in a smaller candidate when
isolating that candidate requires removing fluorine. Similarly, preserving
multiple bonds keeps an unsaturated carbon tied to its multiple-bond partner.
These are sensible fragmentation policies, but restrict the common
substructures the algorithm can discover.

## Implications for the MCS approach

There are two independent limitations:

1. **Candidate availability.** A common substructure must first be generated
   somewhere as a candidate. Ordinary bond-cut fragments cannot represent every
   useful common atom subset under the chosen cutting rules. Limits on the
   number of simultaneous cuts impose further restrictions.
2. **Support measurement.** The current support calculation counts parents
   contributing an exact candidate fragment, rather than all parents containing
   that candidate as a substructure. Even a candidate generated elsewhere may
   be undercounted in parents from which it cannot be isolated.

Removing attachment isotopes from fragment identity and including whole parents
as zero-cut candidates address other sources of missed correspondence, but do
not remove either limitation above. More cuts can help only where the necessary
bonds are eligible. No unrestricted MCS guarantee follows from exhaustive
combinations of a restricted set of cuts.

Without `-s`, selection maximizes generated-fragment support before atom count,
so NC wins in these runs. With `-s`, the tool maximizes size among candidates
meeting the support threshold. For example, `-s 60%` selected NCC with support
61 in the original run. Changing the selection objective does not change how
support is measured.

One possible extension is to generate candidates with dicer, then measure their
support by searching all original parents. That could improve support estimates
without more extensive fragmentation, but adds computation and still cannot
recover candidates that were never generated. Query semantics would also need
to be explicit: the reported SMARTS can contain attachment-derived minimum
connectivity constraints, and those affect which parents match.

## Idea to revisit: extend fragments across a cut

A proposed direction is to add an atom to a generated fragment corresponding to
the atom at the other end of a broken bond, using the same atom type as that
opposite endpoint. This would allow a retained fragment to represent some
context on the discarded side without retaining that entire component.

For example, cutting `N-CH2 | CH=CH2` already gives an N–C fragment. Adding the
carbon immediately across the cut could make an N–C–C candidate available
without breaking the C=C bond in the original parent. Similarly, cutting
`N-CH2 | CH2-F` could expose the adjacent carbon without cutting C–F. This is
an extension of the atom subset represented by a fragment, not a proposal to
make multiple bonds or protected C–F bonds generally breakable.

This is a design idea, not implemented behavior. Before implementing it,
resolve the following:

* Define whether the added endpoint is a real mapped parent atom in a query
  graph, or a synthetic cap in a chemically valid molecule. An ordinary carbon
  cap can acquire hydrogens or lose unsaturation, so these representations do
  not automatically have equivalent atom types or matching semantics.
* Specify which atom properties are retained. Atomic number alone, an existing
  LillyMol atom type, aromaticity, formal charge, and original unsaturation can
  lead to different equivalence classes. If original unsaturation is required,
  the extended unsaturated endpoint may intentionally fail to match a saturated
  endpoint; that would not reproduce the permissive NCC query above.
* Preserve the original crossing bond type and the added atom's parent mapping.
  Mark genuinely synthetic atoms separately so they cannot be mistaken for
  numbered parent atoms or inflate MCS size and support.
* Define what happens with multiple boundary bonds, repeated additions of the
  same opposite atom, and edges between newly included and already retained
  atoms. Induced matching requires accounting for all parent edges within the
  resulting atom subset, not just the edges used to construct it.
* Decide whether expansion is one hop only, whether all cuts are extended at
  once or in combinations, and whether an added endpoint can itself expand.
  These choices directly affect candidate count and runtime.
* Retain a canonical identity and deduplicate by parent, irrespective of which
  cut or original fragment generated the extended candidate.
* Decide whether existing and extended candidates share one support definition,
  or whether candidate generation is separated from parent-wide verification.
  Extension alone does not establish full substructure support.

A useful first experiment would retain the current algorithm as the default
and make one-hop endpoint extension optional. Test propylamine, allylamine,
propargylamine, fluoroethylamine, CF3-containing examples, branched amines, and
ring-attached amines. Check candidate identity, original-atom correspondence,
parent support, and runtime separately. The 100-amine investigation provides a
small motivating example; it should not dictate broader matching semantics.

## Initial element-cap experiment

Dicer now has an optional `-I cap` mode implementing element caps on recursive
chain cuts. It emits capped versions only, retains caps during recursion, uses
ordinary C/N for originally aromatic endpoints, and can combine with
`-I inc=100` to preserve original parent numbering. No synthetic element markers
are used. Existing size checks see the capped components without separate
cap accounting.

On the same 100 amines, with `-B bscb -k 4` and otherwise the original settings,
NC support increased from 95 to 100, but NCC support decreased from 75 to 62.
There were 137 unique candidate identities instead of 160. The anticipated
NCC candidates are generated for the fluorinated and unsaturated examples,
but compulsory capping also converts other bare NCC fragments into larger or
branched candidates. This demonstrates that the change trades one candidate
family for another rather than monotonically extending the original set.
The default support-first selection still chooses NC.

## Recording capped and uncapped forms together

The optional `-I capboth` mode records both forms after one uniqueness check
on the retained original atom subset. Caps are removed from a temporary copy
for the uncapped record; the capped component continues through recursion.
The extra work is copying, cap removal, filtering and canonicalization rather
than a second recursive traversal. Both forms retain the normal proto layout
and numbered-parent correspondence. Intermediate selections of some caps are
not enumerated.

A comparison using `-B bscb` and the same remaining settings gave:

| Mode | Cut limit | Candidate identities | NC support | NCC support | NCCC support |
| --- | ---: | ---: | ---: | ---: | ---: |
| Uncapped | 2 | 160 | 90 | 73 | 45 |
| Capped only | 2 | 137 | 100 | 62 | 19 |
| Both forms | 2 | 160 | 100 | 95 | 51 |
| Uncapped | 4 | 160 | 95 | 75 | 45 |
| Capped only | 4 | 137 | 100 | 62 | 19 |
| Both forms | 4 | 160 | 100 | 95 | 51 |

The unchanged number of candidate identities does not imply unchanged coverage:
existing identities now occur in more parents. All 94 parents matching
`[NH2]-[CR0]-[CR0]` contribute NCC in the dual-form run. The 95th contributor is
cyclopropylmethylamine, ID:EN300-27211 (`NCC1CC1`), where the ring carbon becomes
the cap. It does not match the acyclic-carbon query, but contributes a valid
NCC atom subset under the cap equivalence.

With both forms, the default support-first objective still chooses NC at 100
parents. Requiring `-s 95%` instead selects NCC at 95 parents. These observations
support the usefulness of dual-form recording for this dataset without changing
the broader limitations of a fragmentation-based candidate search.

## Cap provenance without a proto change

Element caps are now distinguished by a reserved isotope marker of 10000 added
to their existing isotope payload. This replaces the initial experiment's
unmarked caps. The reader strips the marker before decoding the attachment
increment, retaining a separate set of cap atom indices. Fragment identity and
canonical correspondence are still based on the isotope-free structure.
For example, 10117 decodes to cap provenance, attachment increment 100 and
parent atom 17. Both writers and readers share the marker constant in
`dicer_lib.h`. Capping modes reject overlapping input/payload values; the reader
rejects malformed marked caps and markers on parent records. These checks do
not introduce a new proto field.

Unmarked experimental files still work for structural analysis but cannot
reliably distinguish caps from ordinary terminal atoms. Generated SMARTS now use atomic numbers only for explicitly marked caps;
genuine atoms retain their element symbols and aromaticity constraints. When
otherwise equivalent occurrences map to the same parent subset, prefer fewer
caps for the exemplar. This replaces the initial all-terminal-atom relaxation.
