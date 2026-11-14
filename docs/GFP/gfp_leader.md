# gfp_leader

A clustering tool also known as sphere exclusion.

Text `DIST<...>` values use platform-independent formatting rounded to four
decimal places, with trailing zeros omitted. This applies to cluster output and
the `-L` report of molecules excluded by previously selected molecules (`-A`).
Clustering thresholds and sorting use the original distances. Textproto output
(`-k .`) retains the original numeric precision.

## Use Case
The most common use case for `gfp_leader` is when you have a set of molecules
which can be ordered by some desirability function, and when you need to select
a limited number of desirable, yet diverse candidates.

The desirability measure might be things like

. A docking score
. One or model model scores
. Medchem demerits
. Heavy atom count
. Other...

There is almost always some relative desirability score that can be assigned to
a set of molecules, even if just heavy atom count in order to favor smaller molecules.

## Algorithm
A radius must be specified. The guarantee is that no two selected molecules will
be closer than this distance. What distance to use will be a function of what
fingerprints are being used. For the standard `gfp` set, typical distances might
be in the range of 0.1 to 0.2, although there is no 'correct' value. It depends on
how many molecules you have, how many you must select, and how confident you are in
the scoring function.

If the scoring function is very good, you might be comfortable with a larger radius,
confident that the model is giving accurate scores to non leader molecules in
each cluster. If on the other hand you had a very weak desirability, heavy atom
count, you might want a shorter radius, since the 'model' is not very predictive
of the other molecules in each cluster. There are many factors to consider, and
choice of threshold cannot be strictly prescribed.

Some have criticized leader for the need to have an arbitrary distance. More
expensive clustering algorithms may not properly emphasise the most desirabile
molecules, and/or have their own tuning parameters.

It is also possible that any given choice of the threshold will yield too many,
or too few clusters. You may need to run multiple times with different thresholds
in order to get a desirable outcome.

Leader scales well. Worst case would be if every molecule ended up in its own
cluster, a very small radius and no duplicates. That would take `n * (n - 1) / 2`
calculations. The other extreme would be if every molecule ended up in the same
cluster, which would take `n - 1` calculations. So it is always better than `n^2`.

Leader starts with the order set of candidates. The first item is marked as selected
and then all subsequent items that are within `threshold` of that first item are
placed in the first cluster. It then returns to the beginning of the list and finds
the first unselected item. Subsequent items that are within `threshold` of that 
leader are placed in its cluster. This continues until all items have been selected,
or until the number of items needed have been selected.

## HowTo
starting with a smiles file, a likely sequence of commands might be
```
gfp_make file.smi > file.gfp
gfp_leader -t 0.15 -v file.gfp > file.ldr
```

If you have a large file, you may want to run `gfp_leader_tbb` in order to 
run multi-threaded. If you are using the standard fingerprints `gfp_make -STD` then
you can use `gfp_leader_standard` which is also multi-threaded.

To explore different distance thresholds for either `gfp_leader` or
`gfp_leader_standard`, consider [nn_leader_and_jp](nn_leader_and_jp.md).
Compute nearest neighbours once, then reuse that file to cluster at several
thresholds without repeating the fingerprint comparisons:

```shell
nn_leader_and_jp -t 0.15,0.2,0.25 -B clusters neighbours.nn
```

This runs one thread per threshold and writes `clusters0.ldr`, `clusters1.ldr`
and `clusters2.ldr`, in threshold order. Add `-r` if the neighbour file contains
TFDataRecord serialized protos. Generate distance-limited neighbour lists that
cover at least the largest threshold you want to explore; limiting the number
of neighbours can omit eligible members and change the clustering. Preserve
the desired input order, since it determines leader priority.

## DISE-style input ordering

If you are looking for Directed Sphere Exclusion (DISE), its reference-derived
priority can be supplied by ordering the input before ordinary leader selection.
The exclusion rule itself is unchanged. To calculate distances to a set of
reference molecules, see the [gfp_lnearneighbours sorting recipe](gfp_lnearneighbours.md#sorting-by-proximity-to-reference-molecules-dise-style-ordering),
which uses the nearest `DIST` value and `tdt_sort` to rank candidate records.
Use the resulting ranking to order molecules before fingerprint generation.

Closest-to-any-reference ordering is a useful alternative to sorting by a
medicinal desirability score. It differs from the original DISE lexicographic
ordering, where subsequent reference similarities only break ties in the first.
This applies to input prioritisation for both `gfp_leader` and
`gfp_leader_standard`.

## Options

The following options are supported by the base, serial, version.
```
Performs leader clustering on a set of fingerprints
Usage <options> <input_file>
 -C <number>      maximum number of clusters to find
 -t <dis>         specify distance threshold
 -t col=<nn>      threshold is column <nn> of the name field
 -t tag=<TAG>     threshold for each molecule in dataitem <TAG>
 -H <TAG>         threshold for each molecule in dataitem <TAG>
 -m <number>      maximum cluster size
 -M <tag>         max cluster size for each molecule in <TAG>
 -M col=nn        max cluster size is column <nn> of the name field
 -S <TAG>         score tag
 -S col=nn        score is column <nn> in the name field
 -I <TAG>         specify identifier tag
 -r               sort clusters by distance from leader
 -E <dataitem>    specify pool object dataitems to be echo'd (default $SMI and PCN)
 -E ALL           echo all dataitems from the pool file
 -A <file>        file(s) of previously selected molecules - discard all within threshold
 -a <dist>        use <dist> as the threshold when comparing against the -A file
 -L <file>        write fingerprints discarded by -A file(s)
 -D ...           miscellaneous options, enter '-D help' for info
 -s <number>      specify max pool size (not needed)
 -F ...           gfp options, enter '-F help' for details
 -V ...           Tversky specification, enter '-V help' for details
 -v               verbose output
```

### -C ncluster
The number of clusters to form. By default, clustering will continue until
all fingerprints have been placed into a cluster.

### -t distance
The threshold used for cluster formation. Note that a constant threshold is
used for all selections. This may, or may not be what you want.

### -t tag=TAG
Specify a per-molecule threshold. For the most desirable molecules you may want
to use a small threshold, so only the very close neighbours get clustered together
with the most desirable molecules. Then use a larger threshold with less desirable
molecules.
```
$SMI<C>
PCN<methane>
THRESHOLD<0.11>
|
```
This is the same as `-H THRESHOLD`. It is more common to place the threshold
in a column of the name, see below.

### -t col=col
Instead of a TDT tag, per-molecule thresholds can be specified as part of the name
field
```
$SMI<C>
PCN<methane 0.11>
|
```
In the example above, `-t col=2` means that the threshold for each molecule if
found in column 2 of the name.

### -m max
Maximum cluster size. Once a cluster has \<max\> members, stop adding
items to that cluster. 

### -M TAG
Specify a per molecule maximum cluster size via a TDT tag.
```
$SMI<C>
PCN<methane>
MAX_CLUSTER_SIZE<50>
...
```
In this case, 'methane', if it is selected as a leader, can only form a
cluster with as many as 50 items. More commonly specifying a maximum
cluster size is done via a column in the name.
```
$SMI<C>
PCN<methane 50>
...
```
and then specifying `-M col=2`.

### -S tag -S col=nn
Select the next leader by highest score, from a TDT tag or a one-based column
in the identifier. Without `-S`, input order determines leader priority.
The former `gfp_leader -R` score-distance adjustment has been replaced by
small-cluster reassignment, described below.

### -A fname
Specify a file of fingerprints containing molecules that are to be treated as
having been previously selected. For example if running two different selection
methods, a docking score and a QSAR model, you can either

* Combine the scores into a composite measure and sort by that.
* Run leader with one method, then again with the other method with the first
set marked as already selected.

A workflow that needs to select 1000 moleules, based on two different scoring
methods might look something like.
```
sort_by_docking_score file.smi > file.sorted1.smi
gfp_make.sh file.sorted1.smi > file.sorted1.gfp
gfp_leader -C 500 -t 0.2 -v file.sorted1.gfp > file.sorted1.ldr

# Fetch the selected molecules, without neighbours.
nplotnn -n 0 file.sorted1.ldr > file.sorted1.sel

# The selections from the first method are now available. Fetch those
# fingerprints from the previously generated file, or regenerate.

fetch_tdt_quick -c 2 file.sorted1.sel file.sorted1.gfp > file.sorted1.sel.gfp
sort_by_qsar_score file.smi > file.sorted2.smi
gfp_make file.sorted2.smi > file.sorted2.gfp
gfp_leader -A file.sorted1.sel.gfp -C 500 -t 0.2 -v file.sorted2.ldr
nplotnn -n 0 file.sorted2.ldr file.sorted2.ldr.smi
```
The two files, 'file.sorted1.sel' and 'file.sorted2.sel' contain the
selections.

### -a threshold
By default, unselected items are discarded if they are within the
global threshold (-t) of any previously selected item (-A). That
can be changed via the `-a` option.

### -r
By default, when a cluster is written, it will be in the order of the
molecules in the input file - which was ordered by desirability. In some
cases it may be desirable to have cluster members sorted by distance
from the leader.

### -s \<size\>
The number of fingerprints in the input. Seldom needed. By default the tool counts
the fingerprints, and then allocates internal arrays. Once upon a time this
may have been slow.

## Reassigning small clusters

`gfp_leader_standard -t 0.2 -R 1 fingerprints.gfp` first forms the usual
leader clusters, then dissolves singletons into larger clusters. `-R 2` dissolves
clusters containing one or two members, counting the leader. The same `-R <size>`
option is available in `gfp_leader`. For precomputed neighbour lists,
[nn_leader_and_jp](nn_leader_and_jp.md#reassigning-small-clusters) also supports
`-R`, subject to the distances available in those lists.

Bare `-R <n>` assigns each dissolved member to its nearest surviving leader.
Use `-R <n>:any` to search all original members of the larger clusters and
join the cluster containing the closest member instead. For example:

```shell
gfp_leader_standard -t 0.2 -R 1 fingerprints.gfp > nearest_leader.tdt
gfp_leader_standard -t 0.2 -R 1:any fingerprints.gfp > nearest_member.tdt
```

The same syntax applies to `gfp_leader`. The cheaper nearest-leader search is
the default; `:any` explicitly selects the previous nearest-member behaviour.
Destinations and eligible targets are fixed from the initial clustering;
newly assigned members do not attract subsequent assignments. Equal distances
favour the earliest eligible target in input order. Members
of a dissolved cluster may join different clusters. If no larger cluster exists,
the tool warns and retains the original clustering.

Reassignment can exceed the initial distance threshold and any initial cluster
size constraint. It preserves destination leaders, and output `DIST` values
remain distances to those leaders, rather than to the closest member. Empty
clusters disappear and output cluster numbers are consecutive. Both ordinary
and `-k .` textproto output support reassignment. Molecules excluded by `-A`,
or not clustered when `-n` (`gfp_leader_standard`) or `-C` (`gfp_leader`) stops
processing, remain excluded.

The final search compares every dissolved member with surviving leaders, or
with all original members of larger clusters when `:any` is specified. `gfp_leader_standard` uses the OpenMP thread count selected
with `-h`; the `gfp_leader` final search is serial. All fingerprints remain
available; neither implementation squeezes its pool.
When no members move, the existing output path is used. After reassignment,
member lists are built in a single pass to avoid repeated scans of the pool.

Without `-R`, `gfp_leader` retains its write-as-you-go output. With `-R`, it
retains membership and defers output until reassignment is complete. Existing
`-E` echo selection, `-r` distance sorting and `-k` output options still apply.
The final search uses the configured distance metric, including Tversky, with
the eligible destination target as the first argument. The initial `-X` early-abandon
cutoff does not apply. An item without any comparable eligible target stays
in its original cluster, with a warning. If a partially dissolved source cluster
loses its leader, the first remaining input member becomes its leader.

`gfp_leader -R` cannot be combined with `-D sbc` or `-D sbcr`: those options
discard members during within-cluster sub-clustering. The initial maximum
cluster sizes (`-m`/`-M`) do not constrain the reassignment pass.

## Paired scans in gfp_leader_standard

`gfp_leader_standard -p 2 -h 8 -t 0.2 fingerprints.gfp` can process two
successive leaders in a single candidate scan. If the next two unassigned
molecules are farther apart than the threshold, both will become leaders; each
remaining candidate is compared with both, giving the first leader priority.
If they belong to the same sphere, the tool falls back to the ordinary scan.

`-p 1` uses the default one-leader path. Both modes preserve leader order,
cluster membership, distances, `-n` limits (including odd limits), and `-R`
reassignment. Paired scans can reduce memory traffic and parallel-region
synchronisation, but also calculate some second distances that separate scans
would avoid. The benefit depends on threshold, dataset, and thread count;
benchmark representative inputs before adopting it as a default.

Before active-index compaction, a local benchmark on `rand50k.gfp`, using an
optimised native build, gave the
following median wall times over three runs (including reading and writing):

| Threshold | Threads | `-p 1` | `-p 2` |
| --- | --- | --- | --- |
| 0.20 | 4 | 12.48 s | 10.30 s |
| 0.20 | 8 | 10.47 s | 9.07 s |
| 0.35 | 4 | 5.46 s | 4.49 s |
| 0.35 | 8 | 4.50 s | 3.74 s |

All output hashes agreed between modes and thread counts. These measurements
are specific to that dataset and machine; on a 10k subset at threshold 0.5,
paired scans did not consistently improve elapsed time.

## Active candidate indices in gfp_leader_standard

Both ordinary and paired leader scans use an input-ordered array of candidate
indices. Molecules excluded by `-A` are left out when the array is built.
Selected slots are marked inactive during the parallel scan. When approximately
one quarter of the remaining range is inactive, the indices are compacted
between scans, preserving input order. This reduces scans of already assigned
items and helps distribute useful work among OpenMP workers.

Fingerprint objects, original molecule indices, cluster assignments and output
metadata stay in their original storage. Compaction therefore also works with
both `-R <n>` and `-R <n>:any`; the final pass can still access every original
fingerprint. Leader selection, `-n` limits, distances and output ordering are
unchanged. No extra command-line option is required.

A local before/after comparison on `rand50k.gfp` used median wall times over
three alternating runs of each executable, including input and output. At eight
threads the results were:

| Threshold | Leaders per scan | Previous scans | Active indices |
| --- | --- | --- | --- |
| 0.20 | 1 | 10.06 s | 9.38 s |
| 0.20 | 2 | 8.60 s | 8.55 s |
| 0.35 | 1 | 4.35 s | 3.99 s |
| 0.35 | 2 | 3.73 s | 3.57 s |

Output hashes matched in every comparison. Gains were generally clearer on the
larger file than on the 10k file, where results were mixed; the paired path had
less remaining scanning overhead to remove. These are measurements for one
machine and dataset, rather than expected speedups for all workloads.
