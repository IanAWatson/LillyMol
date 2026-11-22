# DBF
## Distance Between Features
This tool supports pharmacaphore type distance profiling in either
2D (bond distance) or 3D (geometric distance), or both.

It is driven by queries that identify the features of interest. Then,
for each atom of interest, statistics about the bonds between, or distance
between those atoms are reported.

If there are two queries, then the distances 1-1, 1-2 and 2-2 distances will be reported.

For example if we start with benzene, and specify two features, both 'c', an aromatic
carbon
```
dbf -d -s c -s c -p 2d benzene.sdf
```
we get
```
Name dbf_n0 dbf_n1 dbf_minaa dbf_maxaa dbf_aveaa dbf_minab dbf_maxab dbf_aveab dbf_minbb dbf_maxbb dbf_avebb
benzene 6 6 1 3 1.8 1 3 1.8 1 3 1.8

```
The first two columns show that there were 6 occurrences of each feature in
the molecule.

The `-d` option is important because it enables accumulation of data across
bonded atoms - which may or may not be of interest. In this case we enable it
and find that aromatic carbon atoms in a benzene ring are between 1 and 3
bonds away from other aromatic carbon atoms. Note that the zero distance,
same matched atom, match is not reported.

Then follows the mean bond distance, 1.3 for the first matched atom. Since
all queries are the same the 'ab' and 'bb' columns show the same information.

If this is done in 3D, change `-p 2d` to `-p 3d` we get
```
Name dbf_n0 dbf_n1 dbf_3minaa dbf_3maxaa dbf_3aveaa dbf_3minab dbf_3maxab dbf_3aveab dbf_3minbb dbf_3maxbb dbf_3avebb
benzene 6 6 1.3830 2.7660 2.0640 1.3830 2.7660 2.0640 1.3830 2.7660 2.0640
```
where the topological distances are replaced by spatial distances.

## HOWTO
The following options are recognised
```
  -p 2d          compute topological distances
  -p 3d          compute spatial     distances
  -q <query>     identify features as queries
  -s <smarts>    identify features as smarts
  -e             implement the default positive, negative, acceptor, donor queries (no smarts needed)
  -k <nkeep>     report the <nkeep> shortest distances for each type
  -z first       take the first match when multiple matches to a query
  -z ignore      ignore any query that matches multiple times
  -z oknomatch   ignore molecules where the query doesn't match
  -b <number>    ignore atom pairs <number> bonds apart or closer
  -d             allow bonded features
  -T min,max,dx  profile each enabled distance type in bins of width dx (stderr)
  -C <file>      detect all spatial distances and write to <file>
  -r             compute spatial/topological distance ratios
  -W ...         various hard-coded queries, enter '-W help' for info
  -L <fname>     write labelled molecules to <fname> - last query match wins!
  -i <type>      input type
  -A <qualifier> Aromaticity, enter "-A help" for options
  -t ...         element transformation options, enter '-t help'
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -N <...>       Charge assigner specification, enter "-N help" for usage
  -H <...>       donor/acceptor options, enter '-H help' for details
  -u <string>    missing value specification
  -n             suppress normal output
  -l             reduce to largest fragment
  -E ...         standard element options, enter '-E help' for info
  -v             verbose output
```

The `-p` option controls what calculation is done. If both `-p 2d` and `-p 3d` are
specified, then another set of features that are the ratios of the 3d distances to 
number of bonds between the features.

The `-q` and `-s` options are used to specify the features of interest. Note that
the script `dbf.sh` has hard coded queries already included - charges and hydrogen
bond donors and acceptors.

## Distance histograms

Use `-T min,max,dx` to collect dataset-wide distance histograms, where `dx` is
spacing between bucket centers, rather than the number of buckets. Values are
assigned to the nearest center using the existing LillyMol histogram class. For example:

```shell
dbf -s '[#6]' -p ratio -T 0,20,1 molecules.sdf > distances.txt 2> histograms.txt
```

For each pair of features, stderr contains a labeled histogram for each enabled
distance type: topological (bond counts), spatial (coordinate units), and
spatial/topological ratio when requested. All use the same supplied range and
bucket width. Each accepted pair contributes one sample, so molecules with more
matching pairs contribute more samples. Counts include the final molecule and
are accumulated across all input files. Distances outside the histogram range
are omitted from the histogram; descriptor calculations are unchanged.
