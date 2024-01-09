# gfp_make

gfp_make is the primary means by which structures are converted into fingerprint form.

A typical usage might be
```
gfp_make file1.smi > file1.gfp
gfp_make file2.smi > file2.gfp
gfp_* file1.gfp file2.gfp > result
```
where the result is some kind of comparison between 'file1' and 'file2'.

## Details
gfp_make is a fundamentally simple too, that converts a set of fingerprint specifications
given on the command line to a (possibly pipelined) invocation of one or more
tools that add fingerprints to an input stream.

For example if you invode with the default arguments, and add the -v (verbose) option,
it will show what commands are executed. That might look like
```
temperature -J MPR -E autocreate -g all -l file.smi |
        maccskeys -E autocreate -n -J FPMK -J LEVEL2=FPMK2 -f - |
        iwfp -E autocreate -J FPIW -f -
```
where the default fingerprint, -MPR -IW -MK -MK2, have been turned into three
program invocations. The tool [temperature](/docs/Molecule_Tools/temperature.md)
initiates processing, standardising the smiles, stripping to the largest
fragment, and then writing the fingerprint, with molecular properties, to the output.
The next stage of the pipeline is an invocation of [maccskeys](/docs/Molecule_Tools/maccskeys.md)
which is instructed to generate two fingerprints, the normal, and the level 2 fingerprint
which accounts for number of times set. The resulting stream is then passed to
[iwfp](/docs/Molecule_Tools/iwfp) which adds the linear fingerprint to the
stream. gfp_make is responsible for constructing these command pipelines, which
can be of arbitrary complexity.

### Implementation
The current gfp_make is a perl script, that was first initiated in the 1990's and
has been in continuous use ever since. It has grown considerably and is now too
complex. Work is underway on a ruby alternative that promises to be more
flexible and maintainable. The design has however proven to be remarkably
adaptable and useful, with gfp_make now also serving as the core of svmfp models.

## Fingerprints
Numerous tools within LillyMol have been
implemented so they are compatible with the expectations of gfp_make.

Primarily this involves two functionalities.

1. Read a molecule and generate a fingerprint.
2. Read an already formed fingerprint stream froms stdin and insert one or more extra fingerprints.

We see that in the above example, temperature is satisfying the first requirement,
while maccskeys and iwfp the second. The order in which various tools are added
to the pipeline is hard coded within the logic of gfp_make, and is not meaningful.

`gfp_make` supports a great many fingerprint types, most of which are generated
by executables available within LillyMol. Some of those are desribed here.

## Atom Typing
The concept of the shape/kind of fingerprint, and the atom typing used to describe
the atoms are orthogonal concepts. Therefore, for many fingerprints, the shape
of the fingerprint can be combined with a separate atom typing specification. Generally
these will be of the form
```
-EC3:AY
```
In this case, an Extended Connectivity (Morgan) type fingerprint (radius 3) is generated. The
atom typing is `UST:AY` (see [atom typing](/docs/Molecule_Lib/atom_typing.md). In
this particular case, the type comprises two atomic properties, aromaticity
and compressed atomic number.

A linear fingerprint might be specified as
```
-IW8:AHY1
```
which is a linear path fingerprint, max length 8, and with an atom typing
that is
* aromaticity
* hcount
* compressed atomic number
which then undergoes 1 round of Morgan type shell expansion. In effect this
fingerprint descrives a tube rather than a path. Similarly the Atom Pair
fingerprint can have atom types that have undergone round(s) of shell expansion
thereby describing dumbell shaped motifs
```
-AP10:CPY1
```
is an atom pair fingerprint, up to 10 bonds apart, and with the atom property
being a combination of
. number of connections
. pi electrons
. compressed atomic number
which then undergoes one round of shell expansion.

See [atom typing](/docs/Molecule_Lib/atom_typing.md)

