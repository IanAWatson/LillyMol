# common_names

## Purpose
There are two utilities that deal with duplicate molecules, `common_names`
and `unique_molecules`.

If all you care about is eliminating duplicates as quickly and efficiently
as possible, use `unique_molecules`.

The purpose of `common_names` is to allow identification of which 
molecules are duplicates.

For example, given the file 'input.smi' containing
```
C methane
C carbon
CC ethane
```
the command
```
common_names -S common -v input.smi
```
would result in the file 'common.smi' containing
```
C methane:carbon
CC ethane
```

## HOWTO
The help message for `common_names` is
```
Groups identical molecules and writes groups with concatenated names.
common_names -c -l -g all -S output -v input.smi
 -a             compare graph forms - add 2nd -a option to include H count, recommended.
 -c             exclude chirality information.
 -x             exclude directional bonds.
 -l             strip to largest fragment.
 -I             remove isotopes before storing.
 -T ...         standard element transformation options, enter '-T help'.
 -g ...         chemical standardisation options

 -D <separator> separator for when storing duplicate entries, default ':'.
 -f             single pass operation, smiles output only.
 -y             write first name and count of smiles only.
 -s <size>      maximum number of molecules to process. Will count if not specified.
 -r <number>    report progress every <number> molecules processed.
 -S <name>      specify name for output.
 -X ...         miscellaneous and obscure options.
 -i <type>      specify input file type.
 -o <type>      specify output file type.
 -E ...         standard element options.
 -A ...         standard aromaticity options.
 -v             verbose output.
```
## Background
Molecules can be compared for equality in different ways. Several
command line options govern what components of the input molecules
go into the comparison.

## Large Datasets
A simple invocation of common_names on 100k molecules
```
common_names -S common rand100k.smi
```
takes about 5 seconds on relatively old hardware. If you have larger datasets
consider using ${LILLYMOL_HOME}/contrib/bin/unique_molecules_parallel.sh
which can invoke common names via
```
unique_molecules_parallel.sh -thr 16 -exe common_names -S common -c -l ... big.smi
```
On a modern server class machine with 16 cores, this has processed about
24M molecules in less than two minutes.

## Options
The first group of options control changes to the molecules controlling
unique smiles generation.

It is important to note that duplicate groups are represented by the
first molecule encountered in the input, see the discussion below with
the -a option.

### -c
Remove chirality before comparing. Be cognisant of the accuracy
of chirality information you may have.

### -x
Remove cis-trans bonds before comparing. This is highly recommended
since this information is frequently, incomplete and/or not accurate.
In addition, LillyMol cannot reliably generate unique smiles for
certain cis-trans bonding arrangements - it has never been a 
priority given the general unreliability of the input data.

### -l
Strip counterions before comparing.

### -I
Remove any isotopic information before comparing.

### -T
Frequently it is desirable to consider the heavy halogens to be
equivalent. Use the -T option to enable this.
```
-T 'I=Cl' -T 'Br=Cl'
```
enables such a transformation.

There is an important caveat with this, and other operations that change
the molecule. In all cases, what shows up in the output file (-S) is the
smiles of the *first* molecule that is part of the identified group.
For example, given
```
BrC BrC
ClC ClC
FC FC
IC IC
```
running
```
common_names -T I=Cl -T Br=Cl
```
yields
```
BrC BrC:ClC:IC
FC FC
```
We see that the smiles for the group is the BrC, which is the first
representative of the group encountered in the file. It was changed
to 'ClC' in order to generate the unique smiles used for identifying
the group, but in the output, its original structure is preserved.
This behaviour is intentional and preserves the original
representation of the first molecule encountered in a duplicate group.
An alternative design would be to write the transformed structure used
for comparison, but common_names currently preserves the original
input representation.
In the future, perhaps this behaviour should be controlled by an option, which if 
activated, would write the transformed molecules instead.

### -a
Compare molecules in their graph form. This involves major
changes to the molecule: removing chirality, isotopes and
charges, all bonds become single bonds. The unique smiles of
this 'transformed' molecule are then compared.
In this way, benzene and cyclohexane
are identical, with 'smiles' of 'C1CCCCC1'.
Add a second `-a` option and the hydrogen count of
the starting molecule is included with the unique smiles, so 
benzene and cyclohexane would no longer be the same.
```
C1CCCCC1H6 benzene
C1CCCCC1H12 cyclohexane
```
### -D \<separator\>
`common_names` generates an output file that has molecules found to
be duplicates grouped together. So if the input contained
```
C methane
C CH4
```
the output file would contain
```
C methane:CH4
```
where the `:` character is set by the -D option.

## What is a duplicate?
In summary, this is how the various options that control
unique smiles generation influence comparisons.

| Option Set | Compared Features               |
| ---------- | ------------------------------- |
| none       | exact unique smiles             |
| -c         | chirality ignored               |
| -c -x      | chirality and cis/trans ignored |
| -l         | only consider largest fragment  |
| -I         | isotopes ignored                |
| -a         | transformed to graph            |
| -a -a      | graph + hydrogen count          |
| -T ...     | elements transformed            |

### -f
Single-pass mode. By default common_names makes two passes over
the input files, applying structure transformations and building
a hash of unique_smiles and a second pass that re-reads the
input and does output. With the -f option in effect, the first
pass accumulates more information, requiring more memory, but
now including the smiles of the molecule.
This avoids the second pass. Only smiles output is supported,
and some optons that affect molecule output are not applicable.
Use this option if dealing with files that are not large and you
want smiles output. One advantage of the -f option is that it
supports stdin as input
```
generate molecules ... | common_names -g all -l -c -S - - | ...
```
### -s \<size\>
If not specified, common names will first count the number of
molecules in the input file(s).  It does this because it maintains
certain fixed-size internal data structures, and needs to know, in
advance, how many molecules it needs to process.  This should be
modernised.

The value specified can be larger than the actual number of molecules in the
input file(s), but cannot be smaller.

It is not used when using the `-f` option, since that uses
dynamic arrays internally.

### -r \<number\>
Report progress every \<number\> molecules processed. This can be
helpful when processing large datasets. Again, for truly
big datasets, use unique_molecules_parallel.sh as noted above.

### -y
Only works in the absence of the `-f` option. Rather than concatenate
all names found, just write the first one and a count. So, our
methane example would be written as
```
C methane 2
```
This option can significantly reduce memory usage when many
duplicates are present because only the first identifier and a
count are retained.

### -R \<rxn\>
Apply one or more reactions to molecules before doing the comparison.
Probably better to run trxn before invoking common_names.
Only old style reaction files are supported.

### -S \<stem\>
Specify a file name stem for the output. This is mandatory. By
default smiles are written, so '-S common' will result in the
output file being 'common.smi'. Using '-S -' for stdout does not
work here.

### -X
Several options are available via the -X option.

#### -X num

In the `-S` output file, include a count of the number of instances of each
structure.

#### -X proto=\<fname\>
Write the results as textproto form to \<fname\>. This might look like
```
smiles: "[N+]1(=CN(C)C=C1)CCCC" key: "CCCC[n+]1c[n](C)cc1" id: "CHEMBL3182180" id: "CHEMBL3184676" id: "CHEMBL3559949" 
```
Where `id` is a repeated field, holding the id's of all the molecules matching
the key.

The other options are common across different tools.

## Typical Workflows

### Find exact duplicates
```bash
common_names -S common input.smi
```

### Exact duplicates, with chemical standardisation and largest fragment.
```bash
common_names -g all -l -S common input.smi
```

### Ignore stereochemistry
```bash
common_names -c -x -S common input.smi
```

### Compare graph forms
```bash
common_names -a -S common input.smi
```

### Compare graph forms and retain hydrogen count
```bash
common_names -a -a -S common input.smi
```

### Generate a textproto duplicate report
```bash
common_names -X proto=dups.textproto -S common input.smi
```
