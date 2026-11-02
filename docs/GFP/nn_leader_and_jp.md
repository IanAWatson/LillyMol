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

The `-t` cutoff is inclusive: neighbours at the threshold are eligible; neighbours
beyond it remain available for later clusters. Each threshold in a multi-threshold
run starts with all molecules available.

Leader and Taylor–Butina member distances use platform-independent formatting
with four decimal places and trailing zeros omitted, in both output modes.
Clustering comparisons retain the original precision. Jarvis–Patrick members
have no distance to a designated leader, so their output omits distances.

Multiple thresholds, for example `-t 0.2,0.3 -B clusters`, write
`clusters0.ldr`, `clusters1.ldr`, and so on, with the same distance formatting.
