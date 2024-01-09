# train_test_split_optimise

This tool uses a previous computed nearest neighbour dataset to create one or more train/test
splits that optimise the distances between the training and test splits.

## Algorithm
You specify the fraction of the dataset that will be in the training set via the `-f` option.
Each split starts with a random split of the data. Items are randomly shifted as needed in
order to get the precise number of items in each split.

Then follow some number of optimisation steps, specified with the `-o` option. At each step
two items, from opposite sides of the split, are selected for swapping their set membership.
The change in inter-set distances is computed, and if more favourable, the switch is accepted,
otherwise it is rejected and reversed.

## Details
The tool becomes viable froma run-time perspective by keeping track of the neighbours of each
item, so when two molecules are swapped from train to test, only the distances associdated with
those two molecules are recomputed.

Generally assume that we only care about relatively short distances, and all distances beyond
some threshold are just considered "substantially different" and are all considered the same.
This is not really necessary, but given the run-time advantages seems reasonable.

## HowTo
The process involves three steps.

1. Run gfp_nearneighbours_single_file to find all neighbours within a given distance
2. Convert the resulting nearest neighbour file to TFDataRecord format - eventually
this will not be necessary.
3. Run train_test_split_optimise on that data.

Specifically this might look like

```
gfp_nearneighbours_single_file -T 0.45 -v /tmp/all.gfp > /tmp/all.nn
```
where we select a distance of 0.45. For typical gfp fingerprints, this is a long
distance, and perhaps a shorter value would be better. On the 3100 molecule
set I have been using for testing, this takes 2 seconds.

To convert the text nearest neighbour file to TFDataRecord serialized protos
```
nn2proto -v -T all.tfdata all.nn
```
This is fast.

Then to run train_test_split_optimise that might look like
```
train_test_split_optimise -S split -f 0.5 -n 10 -o 20000 -v all.tfdata
```
We are asking for 10 splits, `-n 10`, each with a 50/50 split between
train and test. In order to generate each split, a random split is generated
and 20k optimisation steps are taken. With this particular dataset, this
calculation takes less than 15 seconds.

The results are written as a number of split* files, as specified by the `-S split`
option combination. These come in pairs, a file of identifiers, and a corresponding
smiles file. For example
```
splitR0
splitE0
```
contain the identifiers for the tRain and tEst splits of the first split. The
corresponding .smi files contain the smiles. 

## Diagnostics
As part of the optimization, a series of `split_stats*.txt` files are created. These
might look like
```
Dist Count Cumulative
0.01 6 6
0.02 9 15
0.03 12 27
0.04 13 40
0.05 16 56
0.06 11 67
0.07 15 82
0.08 24 106
0.09 23 129
```
which shows that there are 6 instances of molecules with a distance within 0.01
that are separated across the train/split divide.
