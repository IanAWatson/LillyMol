`in/neighbours.tfrecord` is a TFDataRecord file of `nnbr::NearNeighbours` protos,
generated from the adjacent human-readable neighbour file:

```
nn2proto -T in/neighbours.tfrecord in/neighbours.nn
```

It covers multiple neighbours, empty lists, self-neighbours, distinct zero-distance
neighbours, and a missing query ID. Regenerate with the optimized LillyMol tool.


`in/neighbours.indexed.tfrecord` contains equivalent indexed lists. Its source,
`in/indexed.nn`, also includes empty records for all neighbour IDs, since indexed
references require a complete record-index-to-name mapping:

```shell
nn2proto -T /tmp/neighbours.named.tfrecord in/indexed.nn
nn_id_to_ndx -S in/neighbours.indexed.tfrecord /tmp/neighbours.named.tfrecord
```

The regression compares named/indexed continuous and classification tables,
descriptive output, summaries and textprotos. It checks both `-N` option orders,
invalid combinations and protection against overwriting the indexed input.
