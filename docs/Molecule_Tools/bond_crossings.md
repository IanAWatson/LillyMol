# bond_crossings

Scores crossed bonds in **existing 2D molecular coordinates**, providing a
measure of depiction congestion suitable for ranking molecules or descriptor
analysis. It does not generate or optimise a depiction, or compute a graph's
minimum crossing number. No optional coordinate-generation library is needed.

```shell
bond_crossings molecules.sdf > scores.txt
bond_crossings -l -z i molecules.sdf > scores.txt
bond_crossings -d molecules.sdf > detailed_scores.txt
```

If there is no need to retain the 2D generated coordinates
```shell
make_2d_coordinates -o smi3d -p best -e -n -w -c -v -S - input.smi | \
bond_crossings -z i -l -v -i smi - > scores.txt
```
but note that generation of coordinates is much slower than detecting
crossed bonds.

The output is a headered, space-separated descriptor table with one row per
molecule. Identifiers use the first token of the molecule name; unnamed records
receive `molecule_N`. Floating point output uses stable four-decimal formatting
with trailing zeros removed. Calculations use their original precision.

For example, sort by decreasing cubic score (column 12):

```shell
tail -n +2 scores.txt | sort -k12,12nr > ranked_rows.txt
```

## What is counted

Bonds are straight centre-line segments in XY, including double, triple, aromatic,
and wedge/hash bonds as one segment each. A crossing is a proper intersection of
the interiors of two segments. Two bonds sharing an atom are never a crossing.
Each crossing bond pair is counted once, even when several pairs cross at the
same location. The number of distinct participating bonds is reported separately.

At least one member of a pair must be a ring bond. **Chain–chain crossings are
never counted**, including in the optional pair diagnostics. Thus valid acyclic
inputs always have zero crossing counts and score, although their distance
descriptors are still calculated. Actual crossings in a monocyclic drawing are
counted; the program does not assume the supplied drawing is optimal.

Fused and bridged rings share a system. Rings joined only at a spiro atom remain
separate. Ring bonds are assigned from the rings, avoiding the ambiguity of a
single atom label at a shared spiro atom. Per-system atom counts include that
shared atom in each system; `bc_nratoms` counts it only once.

A ring–chain crossing is assigned to the system containing the crossed ring
bond, regardless of where the chain is attached. Crossings between ring bonds
of different systems are reported separately.

## Scoring and normalisation

For system i, let c_i be its within-system ring–ring crossings plus its
ring–chain crossings. Let I be the number of between-system ring–ring crossings.
The score is

```
score = sum(c_i ** 3) + I
```

One, two, and three crossings in one system contribute 1, 8, and 27. One crossing
in each of two systems contributes 2. Between-system crossings contribute a
linear penalty of one each, without assigning them to either local system.
Every properly crossed eligible pair therefore makes the total score positive.
Degeneracies and near misses do not change the primary score.

| Column | Definition |
| --- | --- |
| `bc_natoms` | Atoms in the molecule after optional `-l`; explicit hydrogens are retained |
| `bc_nrs` | Ring systems, keeping spiro-only connections separate |
| `bc_nratoms` | Distinct ring atoms |
| `bc_crossings` | All eligible proper crossing pairs |
| `bc_ring_ring` | Ring–ring crossing pairs within the same system |
| `bc_ring_chain` | Ring–chain crossing pairs |
| `bc_inter_system` | Ring–ring crossing pairs between different systems |
| `bc_crossed_bonds` | Distinct bonds participating in proper crossings |
| `bc_nrs_crossed` | Systems with c_i > 0; excludes between-system crossings |
| `bc_max_crossings` | Largest c_i, or zero when there are no systems |
| `bc_score` | Cubic score plus the between-system penalty |
| `bc_score_per_atom` | Score / atom count |
| `bc_score_per_rs` | Score / ring-system count |
| `bc_score_per_ring_atom` | Score / distinct ring-atom count |
| `bc_worst_density` | Maximum c_i ** 3 / atoms in system i |
| `bc_min_2d` | Closest pair of distinct atoms, including bonded pairs |
| `bc_min_2d_per_atom` | Closest distance / atom count |
| `bc_max_2d` | Largest pairwise XY distance |
| `bc_max_2d_per_atom` | Largest XY distance / atom count |
| `bc_topological_diameter` | Maximum shortest-path distance in bonds |
| `bc_max_2d_per_topology` | Largest XY distance / topological diameter |

A ratio with a zero denominator is zero. Pairwise distances are zero for a
single atom. Distance measures retain the input coordinate units and change
with uniform scaling; crossing counts and scores do not. A large simple
substituent can dilute whole-molecule normalisations, while `bc_worst_density`
retains the local ring-system measurement.

## Input policy

Input must be nonempty, connected, and have finite coordinates with
`highest_coordinate_dimensionality() == 2`. This is LillyMol's coordinate test,
not the molfile's 2D/3D header flag. It requires Z coordinates to be zero and may
reject axis-aligned or coordinate-free records such as isolated methane.

* `-l` calls `reduce_to_largest_fragment_carefully()` before validation and
  scoring, preserving the selected fragment's coordinates.
* Disconnected inputs otherwise fail with a message suggesting `-l`.
* `-z i` skips records failing only the dimensionality check. The skipped count
  is reported on stderr, even without `-v`. Other invalid inputs still fail.
* `-i type` selects the input type explicitly. Otherwise suffixes determine it.
  A single `-` reads SDF from stdin by default.
* `-A ...` supplies the usual aromaticity controls; `-v` reports molecule counts.

An invalid record causes nonzero exit status. Rows already written before that
record remain partial output and should not be taken as a complete result.

## Optional diagnostics

`-d` appends the following columns:

| Column | Definition |
| --- | --- |
| `bc_degenerate_bonds` | Bonds whose XY length is at most the geometric tolerance |
| `bc_overlaps` | Eligible non-adjacent bond pairs overlapping collinearly |
| `bc_near_misses` | Eligible non-crossing, non-touching pairs closer than the near-miss threshold |
| `bc_atoms_on_ring_bonds` | Non-endpoint atom/ring-bond pairs where the atom lies on the bond interior |
| `bc_coincident_atoms` | Distinct atom pairs separated by at most the tolerance |
| `bc_min_nonbonded_2d` | Closest pair of nonbonded atoms, or zero if none exists |

The reference length is the upper median of positive bond lengths, or 1 if none
exists. Geometric tolerance is `1e-7` times that reference length. Near misses
use 0.1 times the reference length by default; `-t fraction` changes this and
requires `-d`. A zero fraction disables near-miss detection. Degenerate segments
are excluded from pair intersection tests. Endpoint touches are not proper
crossings or near misses; atoms on ring-bond interiors have their own diagnostic.

## Library and validation

The reusable implementation is `Depict/bond_crossings.h`, target
`//Depict:bond_crossings`. `Analyse` returns a status and resets its `Result` on
every call. `Result::systems()` exposes the individual system sizes, ring counts,
and crossing counts for callers needing more than the CLI's summary. Ring
perception and distance caches may be populated, but coordinates and molecular
connectivity are not changed.

```shell
bazel build -c opt //Molecule_Tools:bond_crossings
bazel test -c opt //Depict:bond_crossings_test
ruby test/run_all_test.rb -bindir "$PWD/bazel-bin/Molecule_Tools" bond_crossings
```

Tests use explicit synthetic coordinates, including crossed rings, chain-only
crossings, substituent crossings, fused/bridged cages, spiro systems, separate
systems, rigid transforms, and optional diagnostics. No depiction engine is
needed to build or test the analysis.
