# nn_leader_and_jp review follow-up

Address correctness issues before introducing the proposed optimisations. This
checklist records the review and small-input reproductions from November 2026.
The existing regression cases passed despite the bugs below.

## Addressed

- [x] Remove `-H`, which dereferenced an uninitialised per-item threshold array.
  Remove that array and the threshold-index arguments throughout the clustering
  code. Only global thresholds specified with `-t` remain supported.
- [x] Restrict `-M` to `col=<col>`. Tag-based limits were accepted but never read.
  Reject those specifications rather than silently ignoring them.

- [x] Remove `-S` and all explicit score storage and parsing. Input order now
  defines priority, with earlier molecules preferred. This eliminates duplicate
  score-file completeness bugs and unsupported score tags with proto input.
- [x] Remove the disabled `NEED_TO_SORT_EVERYTHING` implementation, which
  referenced the retired scores and contained an uninitialised index and an
  invalid comparator.

- [x] Correct capped cluster sizes, headers and returned counts. The limit
  includes the leader; counting and writing use the same neighbour prefix.
  Name-column limits of one yield singletons, and `-h` applies consistently.
  Preserve the uncapped self-neighbour output convention when `-h` is absent.

- [x] Fix TDT parent identifier extraction to start immediately after the tag,
  preserving the first character. Test `-M col=1` with multi-digit and
  single-digit names, and single-character identifiers in leader, Butina and
  JP modes, using default and custom identifier tags.

- [x] Remove `-e`, the random-leader selection path and RNG initialisation.
  This also eliminates the bias in the previous random-start/next-available
  selection scheme. All leader selections now follow input priority.

- [x] Check TDT array bounds before writing neighbour IDs or distances, and
  require both final counts to equal the allocated neighbour count. Move
  debug array reads after validation. Regression tests cover excess distances
  (including zero-neighbour records), missing distances, and valid empty
  neighbour lists in all three modes.

- [x] Require `-s` for multiple input files and enforce one overall pool limit
  in both passes. Do not open later files once the limit is reached. Oversized
  limits process all available data silently. When the explicit limit is
  reached, omit unresolved neighbours and compact indices/distances together
  in place; otherwise unresolved IDs remain errors. Subset clustering, including
  JP shared-neighbour counts, uses only retained neighbours. Tests cover both
  input formats and all clustering modes.

- [x] Separate record, EOF and error outcomes in first-pass readers. TFRecord
  frame/read failures and protobuf parse failures now reject the whole run,
  including after a valid prefix or in a later file. Second-pass readers must
  fill every expected item rather than accepting a positive partial count.
- [x] Reject empty/malformed TDT identifiers and missing metadata or record
  terminators without assertions. Reject empty numeric distance fields before
  substring extraction. Empty SMILES remain allowed. Preserve `-s` behaviour:
  corrupt data beyond the explicit cutoff is not read or validated.

- [x] Accept the inclusive threshold range `[0, 1]` in single, CSV and
  repeated `-t` arguments. Reject empty values and invalid or out-of-range
  numbers. Tests check zero-distance matches, distance-one neighbours,
  isolated molecules and identical results across all three option forms.

- [x] Remove self neighbours once during ingestion when `-h` applies (leader
  and Taylor–Butina). Combine this with omitted-neighbour compaction, keeping
  indices/distances aligned without another allocation. Remove self checks
  from cluster formation. Mixed self/no-self lists no longer bias Butina's
  input-order tie breaking; JP's existing interpretation is preserved.
- [x] Reject non-finite distances during TDT/proto loading, before pruning or
  clustering. NaN previously passed proto range comparisons and could join
  molecules even at threshold zero. No finite-value checks are added to the
  clustering kernels.

## Correctness status

All correctness issues recorded in this review have been addressed.

## Behaviour decisions

- Input order defines priority. Ordinary leader clustering takes the first
  remaining molecule. Taylor–Butina maximises available neighbours and breaks ties by
  earliest input position. Preserve these rules in future parallel reductions.
- JP currently requires sufficient shared neighbours, without requiring the
  two items to occur in each other's neighbour lists. Preserve this variant
  during optimisation unless a deliberate behaviour change is agreed.
- Neighbour lists are assumed sorted by distance. Do not add repeated sorting
  or validation in the clustering kernels.
- Neighbour IDs are assumed distinct within a list. Malformed duplicate IDs
  can inflate counts and size headers; normal GFP-generated lists should not
  contain them. Any future validation should happen during ingestion, not in
  the clustering kernels.

## Optimisation work, after correctness

1. **Implemented:** JP union–find merging and one-pass member grouping,
   preserving component/member order by original input position (details below).
2. Investigate an inverted neighbour index to avoid all-pairs JP testing.
3. Maintain Taylor–Butina available-neighbour counts with reverse adjacency
   and an indexed priority queue, instead of recounting the pool each round.
   Account for directed/truncated neighbour lists.
4. Consider contiguous neighbour indices/distances with row offsets, using
   counts from the first pass to allocate and filling them during the second.
5. Separate immutable graph data, options and per-run state. Parallelise
   independent parameter runs, or pair testing and deterministic reductions,
   without racing on selection arrays, union state or output streams.

Benchmark realistic neighbour densities and cluster sizes. Parallelising the
current repeated scans alone will retain their unfavourable scaling.


## Parallel threshold implementation

Leader and Taylor–Butina threshold runs now use one `std::thread` per threshold.
Single-threshold runs remain serial. Each worker owns its selected array, cluster
statistics and output stream; the shared pool is const during clustering.
The dispatcher joins all started workers, propagates errors/exceptions and keeps
filenames tied to the original threshold indices. A dispatcher-owned mutex
protects verbose messages and complete reports, without locking neighbour loops.
No worker-count option or inner-loop parallelism is introduced.

Regression tests compare output files byte-for-byte with independent serial
runs for both formats/algorithms, caps, subsets, repeated thresholds, verbose
output and a failed output alongside successful workers. Remaining optimisation
ideas above are still proposals; JP runs are not parallelised by this change.


## JP union–find implementation

The pair-test criterion and all-pairs traversal are unchanged. Qualifying pairs
now join components using union by size and path halving, instead of scanning
and relabelling the entire pool. Already-connected pairs require no relabelling.
Once comparisons finish, a reverse input-order pass builds linked member lists
using array indices. Output visits each component at its earliest input member,
then follows its members in input order. The internal root need not be that
first member. This removes repeated output scans without sorting or copying
neighbour data; working storage remains linear in pool size.

A breadth-first graph reference tests components independently of union–find,
including late bridges, nonminimal roots, interleaved components, self/directed
lists, neighbour-list limits, multiple parameter runs, subsets, both input and
output formats. These tests pass against both the old and new executables.
The complete existing regression suite also passes.

Local synthetic timings (seconds, median of three runs, including reading and
output) compared the previous executable with this change. Output bytes matched
in every comparison. Merge-heavy inputs give every molecule the same single
neighbour, causing every pair to qualify; singleton inputs have empty lists.
These illustrate the removed costs rather than predict real-data speedups.

| Input | Molecules | Previous | Union–find | Speedup |
|---|---:|---:|---:|---:|
| Merge-heavy | 500 | 0.0708 | 0.0150 | 4.7x |
| Merge-heavy | 1,000 | 0.4777 | 0.0390 | 12.3x |
| Merge-heavy | 2,000 | 2.6732 | 0.1122 | 23.8x |
| Singletons | 2,000 | 0.1037 | 0.0998 | 1.0x |
| Singletons | 4,000 | 0.1282 | 0.1121 | 1.1x |
| Singletons | 8,000 | 0.2179 | 0.1550 | 1.4x |

JP remains serial. All-pairs neighbour-list testing is still quadratic and is
the next candidate for algorithmic or parallel work.
