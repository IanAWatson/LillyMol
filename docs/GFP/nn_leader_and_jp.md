# nn_leader_and_jp

Clusters an existing nearest-neighbour file using sphere exclusion (leader),
Taylor–Butina, or Jarvis–Patrick clustering.

```shell
nn_leader_and_jp -t 0.2 neighbours.nn > clusters.tdt
nn_leader_and_jp -T -h -t 0.2 neighbours.nn > butina.tdt
nn_leader_and_jp -J 10,3 neighbours.nn > jp.tdt
```

Use `-r` for TFDataRecord serialized nearest-neighbour protos instead of TDT
input. `-a` writes SMILES with identifiers instead of TDT records.

The `-t` cutoff accepts the inclusive range `[0, 1]`, whether specified as a
single value, a CSV list, or repeated `-t` options. Zero admits only zero-distance
neighbours; one admits every neighbour present in the input list. Neighbours at
the threshold are eligible; neighbours beyond it remain available for later clusters. Each threshold in a multi-threshold
run starts with all molecules available.

Leader and Taylor–Butina member distances use platform-independent formatting
with four decimal places and trailing zeros omitted, in both output modes.
Clustering comparisons retain the original precision. Jarvis–Patrick members
have no distance to a designated leader, so their output omits distances.

Multiple thresholds, for example `-t 0.2,0.3 -B clusters`, run concurrently,
one thread per threshold, and write `clusters0.ldr`, `clusters1.ldr`, and so on.
File numbers follow the threshold order on the command line, regardless of
completion order. Each file is identical to a separate single-threshold run,
with the same distance formatting. This applies to both leader and Taylor–Butina
clustering; Jarvis–Patrick parameter runs remain sequential.

The loaded neighbour data is shared read-only. Each threshold has its own
selection array and output stream, so memory for selection state grows with the
number of thresholds. A single threshold uses the original serial path to
standard output. Verbose reports are protected from interleaving, but their
completion order may vary. A worker failure gives a nonzero exit status after
all started workers have finished.

Leader and Taylor–Butina require a global distance threshold via `-t`.
Per-molecule threshold tags (`-H`) are no longer supported.

Use `-m <number>` for a global maximum cluster size, or `-M col=<col>` to
read the maximum from a column of each molecule's name. Columns are numbered
from one, including the identifier token. Tag-based `-M` specifications are
rejected; only `col=<col>` is supported. Limits include the leader. A global
`-m` limit must be at least two; a per-molecule limit of one is allowed and
produces a singleton. A per-molecule limit overrides the global `-m` limit.
Use `-h` to exclude self-neighbours from both leader and Taylor–Butina
clustering, including the neighbour counts used to select Butina centres.
This does not exclude zero-distance neighbours with different identifiers.
The option does not change Jarvis–Patrick neighbour-list semantics.

Input order defines molecule priority: put the highest-scoring or most desirable
molecules first before generating the neighbour file. With multiple input files,
priority follows the command-line file order and then record order within each
file. Explicit scoring options (`-S`) are no longer supported.

Leader clustering chooses the first remaining molecule. Random leader selection
(`-e`) is no longer supported. Taylor–Butina chooses the molecule with the largest number
of unselected neighbours within the threshold; ties go to the earliest molecule
in the input. Input priority therefore resolves ties without overriding the
neighbour-count criterion.

Taylor–Butina maintains these counts after each cluster using reverse neighbour
links and an indexed priority queue. It therefore retains dynamic reordering,
including for asymmetric neighbour lists, without recounting the entire pool.
The reverse index is shared between threshold workers; counts and priorities
are private to each worker. This adds memory linear in the eligible edges and
pool size, in exchange for reducing the work when many clusters are formed.

A single input file can be read without `-s`. Multiple input files require
`-s <number>`, which sets one overall limit, rather than a limit per file.
Reading stops at that number even in the middle of a file; later files are not
opened. A limit larger than the available data simply processes everything,
without a warning.

When the limit is reached, neighbours whose identifiers are absent from the
retained pool are omitted. Clustering therefore uses only the retained subset;
it is not equivalent to clustering the complete collection and filtering its
output afterwards. This also reduces Jarvis–Patrick shared-neighbour counts,
and `-J` list-size limits apply to the remaining neighbour lists. Unknown
neighbours remain errors when the pool limit has not been reached.

Neighbour distances must be finite and in `[0, 1]`; NaN and infinities are rejected
when the data is loaded.

Malformed input is an error, not an end-of-file marker. Truncated TFRecord
frames, checksum failures, invalid protobufs, missing identifiers, and incomplete
TDT records cause a nonzero exit status instead of clustering a partial dataset.
TDT records must contain SMILES and identifier fields and end with `|`; identifiers
must be nonempty. As with normal records, data beyond an explicit `-s` limit is
not examined.

Distance thresholds are checked against the maximum stored neighbour distance in each
input file. A threshold above that maximum produces a warning, but clustering
continues using the stored neighbours. Files without neighbour distances also
produce a warning. This is only a sanity check: the observed maximum is not the
generation cutoff, and count-limited lists cannot establish complete distance
coverage. The check uses records read under `-s`, before self or excluded-subset
neighbours are removed. Jarvis–Patrick does not use this distance check.

## Reassigning small clusters

```shell
nn_leader_and_jp -h -t 0.2 -R 1 neighbours.nn > reassigned.tdt
nn_leader_and_jp -r -T -h -t 0.2,0.3 -R 2 -B reassigned neighbours.rio
```

`-R <size>` enables a final pass for leader or Taylor–Butina clustering. Size
includes the leader: `-R 1` reassigns singletons, while `-R 2` also dissolves
pairs. Each member independently joins the cluster containing its nearest
**stored** neighbour among original members of clusters larger than this limit.
Only the initial larger clusters are destinations; newly reassigned members do
not attract further items. A pair may split between two destinations. Equal
distances favour the earliest neighbour in molecule input order, regardless of
its position among equal-distance entries in a neighbour list.

The final pass examines stored distances beyond the initial `-t`. Resulting
clusters can exceed both the original threshold and `-m`/`-M` size constraints.
Use neighbour lists with enough distance coverage for this search; this tool
cannot calculate missing distances from fingerprints. Members with no stored
neighbour in a destination cluster retain their original membership, with a
warning. If no larger clusters exist, the original clustering is retained.
Molecules left unclustered by `-C` are neither reassigned nor destinations.

Destination leaders remain unchanged. A partially dissolved source cluster
whose leader moved uses its first remaining input member as its new leader.
`DIST` continues to mean distance to the final leader. Either stored direction
can supply that distance; if neither exists, the member is written without a
`DIST` field (or without a distance column in `-a` output). A warning reports
the number of missing leader distances. The attracting member's distance is
never substituted for an unavailable leader distance.

With `-R`, output is deferred until the final pass, and empty clusters are
removed with consecutive output cluster numbers. Without `-R`, write-as-you-go
output is unchanged. Membership and reassignment state are private to each
threshold worker. `-R` is incompatible with Jarvis–Patrick (`-J`). Verbose
cluster-size statistics describe the final clustering.
