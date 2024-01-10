# train_test_split_optimise

This tool uses a previous computed nearest neighbour dataset to create one or more train/test
splits that optimise the distances between the training and test splits. When building models
these kinds of splits will often be the most difficult to predict.

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

1. Run `gfp_nearneighbours_single_file` to find all neighbours within a given distance
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
This is fast. Eventually `gfp_nearneighbours_single_file` will be able to
generate this file directly.

Then to run train_test_split_optimise that might look like
```
train_test_split_optimise -S split -f 0.5 -n 10 -o 200000 -r 10000 -v all.tfdata
```
We are asking for 10 splits, `-n 10`, each with a 50/50 split between
train and test. In order to generate each split, a random split is generated
and 200k optimisation steps are taken. With this particular dataset, this
calculation takes less than about 35 seconds on a 2017 consumer class cpu.

The results are written as a number of split* files, as specified by the `-S split`
option combination. These come in pairs, a file of identifiers, and a corresponding
smiles file. For example
```
splitR0
splitE0
```
contain the identifiers for the tRain and tEst splits of the first split. The
corresponding .smi files contain the smiles. 

Also written are a series of `split_stats*` files. These hold the distribution of
distances between train and test. At the start of optimization, a random split
is written to `split_stats_rand.txt` which might look like
```
Dist Count Cumulative
0 2 2
0.01 54 56
0.02 124 180
0.03 210 390
0.04 233 623
0.05 309 932
0.06 434 1366
0.07 437 1803
0.08 580 2383
0.09 607 2990
0.1 646 3636
0.11 736 4372
0.12 775 5147
0.13 932 6079
0.14 1047 7126
0.15 1194 8320
0.16 1268 9588
0.17 1399 10987
0.18 1500 12487
0.19 1703 14190
0.2 1875 16065
0.21 2046 18111
0.22 2246 20357
0.23 2465 22822
0.24 2755 25577
0.25 3016 28593
0.26 3328 31921
0.27 3656 35577
0.28 4025 39602
0.29 4391 43993
0.3 5094 49087
0.31 5420 54507
0.32 6048 60555
0.33 6732 67287
0.34 7682 74969
0.35 8629 83598
0.36 9833 93431
0.37 11418 104849
0.38 13076 117925
0.39 15013 132938
0.4 17502 150440
0.41 19809 170249
0.42 22381 192630
0.43 25434 218064
0.44 28736 246800
0.45 33066 279866
0.46 1718 281584
0.47 2184886 2466470
```
here we see that there are a great many instances of molecules spanning the train/test divide that
are at relatively close distances. After optimisation, one of the other split_stats files might
look like
```
Dist Count Cumulative
0.04 2 2
0.05 3 5
0.06 1 6
0.07 1 7
0.08 2 9
0.09 4 13
0.1 1 14
0.11 1 15
0.12 1 16
0.13 1 17
0.14 3 20
0.15 8 28
0.16 5 33
0.17 10 43
0.18 7 50
0.19 9 59
0.2 16 75
0.21 19 94
0.22 31 125
0.23 27 152
0.24 35 187
0.25 36 223
0.26 42 265
0.27 64 329
0.28 57 386
0.29 134 520
0.3 166 686
0.31 208 894
0.32 277 1171
0.33 360 1531
0.34 499 2030
0.35 679 2709
0.36 816 3525
0.37 1182 4707
0.38 1608 6315
0.39 2178 8493
0.4 3149 11642
0.41 4400 16042
0.42 6048 22090
0.43 7952 30042
0.44 10918 40960
0.45 14633 55593
0.46 848 56441
0.47 2410029 2466470
```
where we see that the number of pairs at short distances has decreased substantially. In this
particular case, the average of pairs *not* at the maximum distance is 0.3722 before optimisation,
and 0.4689 after. Originally there were 2.18M cross set pairs at the maximum distance, and after
optimisation, that rises to 2.41M.

The `-r 1000` option says to report progress every 1000 optimisation steps. A typical output might
look like
```
0 193000 score 115644689 cmp 113170036 accepted 2208 0.0114404 last successful 190190
0 194000 score 115639380 cmp 113170036 accepted 2208 0.0113814 last successful 190190
0 195000 score 115641857 cmp 113170036 accepted 2208 0.0113231 last successful 190190
0 196000 score 115636336 cmp 113170036 accepted 2208 0.0112653 last successful 190190
0 197000 score 115636609 cmp 113170036 accepted 2208 0.0112081 last successful 190190
0 198000 score 115635180 cmp 113170036 accepted 2208 0.0111515 last successful 190190
0 199000 score 115641424 cmp 113170036 accepted 2208 0.0110955 last successful 190190
Writing split 0 score 115646849 computed 115646849 diff 0
starting_score 113170036 score 115646849 improvement 2476813 across split 46.8876 at max 2406506
Split 0 accepted 2208 of 200000 steps
```
We see that no changes were accepted in the last 10 optimisations, so 200k optimisation
steps were not necessary for this split. Other splits did record successful switches close
to 200k steps, and one never made a change after 176k switches. Unpredictable. Run time
is attractive, so let it run.

## Other Options
Initially this was expected to be a potentially very long running tool, so a `-t` option
was added to give the opportunity to optimise each split for a fixed amount of time.

While the initial near neighbour determination might be done with a particular maximum
distance, there is conceptually no reason that same maximum distance needs to be used
for the optimisation. So as the dataset is read, nearest neighbours beyond the `-T`
option are discarded. While this seems reasonable, for reasons I do not yet understand
it breaks the internal optimisation shortcuts. Do not use.

