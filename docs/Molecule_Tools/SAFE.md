# SAFE
Following the very insightful work of [Noutahi et. al](https://arxiv.org/pdf/2310.10773.pdf)
LillyMol includes a utility for converting from smiles to SAFE form. This implementation uses bond breaking
rules based on dicer fragment formation rules.

## HOWTO
```
mol2SAFE file.smi > file.sfae.smi
```
will generate the SAFE representation for each input molecule.

## Performance
Performance seems reasonable. Processing 2.18M molecules (max 50 heavy
atoms) from Chembl, takes less than 2 minutes, or 19k per second on hardware
```
model name      : Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz
```
from 2017. This rises to over 3 minutes if fragment statistics are accumumalated
with the `-S` option below.

## Usage

The following options are recognised
```
Transform smiles to SAFE representations.
Input must be a smiles file.
 -p             write the parent molecule as well
 -c             remove chirality
 -l             strip to largest fragment
 -y             re-use ring numbers in generated ring numbers
 -I <iso>       place a fixed isotope on all attachment points
 -P <atype>     atom typing specification. If specified the isotope will be the atom type.
 -S <fname>     write fragment statistics in dicer_data::DicerFragment textproto form
 -M ...         constraints on max fragment size
 -M <n>         fragments can contain no more than <n> atoms
 -M maxnr=<n>   the maximum number of non-ring atoms that can be in a fragment
 -g ...         standardisation options
 -z             ignore connection table errors on input
 -T ...         Element transformations - suggest '-T Br=Cl -T I=Cl'
 -t             test the unique smiles of the parent against that formed from the SAFE representation.
                This only works if chirality is removed and isotopes NOT added.
 -v             verbose output
```
## Options
### -p
Write the parent molecule in addition to the SAFE form. Output might look like
```
C1=CN=C(C=N1)NS(=O)(=O)C1=CC=CC=C1 CHEMBL596295
N%10S%11(=O)=O.C%111=CC=CC=C1.C1=CN=C%10C=N1
```
This can be helpful during debugging and examiniation.

### -c
Remove chirality. During fragmentation much of the chirality in the molecule will be
destroyed. Unless you are sure of the chiral information associated with input molecules
it might be best to just discard all chiral information.

### -l
Discard likely salts. Generally desirable.

### -y
By default, the tool will not re-use ring numbers, but it is likely desirable to 
minimise the number of ring opening and closing numbers used. Use the `-y` option to enable
re-use of ring numbers. This represents a small run-time penalty, while possibly making
the resulting strings more useful.

### -I
Specify an isotopic label to place on the attachment points. This will be desirable if
you wish to store SAFE strings as smiles - the SAFE form can easily be reconstructed
from an isotopic smiles.

### -P
specify [atom typing](/docs/Molecule_Lib/atom_typing.md). If specified, the atom type
number will be applied as an isotopeisotope. For example
```
mol2SAFE -I UST:AY -p file.smi
```
might result in
```
O=C1N(NC(=O)C=CC(=O)O)C(=NN1)CC CHEMBL1988240
O%10.C%11=[3001C]%12.[3001C]%10%12=O.C%13C.N%14[3001C]%11=O.O=C1[6044N]%14[3038C]%13=NN1
``` 
The atom type is the atom type of the atom from which the atom with the isotope was
detached, This can then facilitate subsequent precise reconstruction based on SMART
fragments.

### -S <fname>
Write a summary of the fragments encountered. This is a dicerdata::DicerFragment
textproto form. Typical entries might look like
```
iso: ATT smi: "C12=C[1CH]=[1CH]C=C1CCC(=O)C2" par: "CHEMBL192267" n: 1 
iso: ATT smi: "[1NH]1C(=O)C=NC(=O)[1CH2]1" par: "CHEMBL286918" n: 1 
iso: ATT smi: "O1[1CH2]CCC2=CC=CC=C12" par: "CHEMBL157323" n: 127 
iso: ATT smi: "C12(COC1)[1NH]CCOC2" par: "CHEMBL4751373" n: 1 
iso: ATT smi: "[1CH]1=C[1CH]=CC2=C1OC1=C2CNCC1" par: "CHEMBL3696983" n: 143 
iso: ATT smi: "[1CH2]1CC2N(CCC2)CC1" par: "CHEMBL4951195" n: 10 
iso: ATT smi: "[1CH2]1CC2COC(=S)N2C1" par: "CHEMBL3947263" n: 1 
iso: ATT smi: "[1CH]1=C2C3(OCCN2N=C1)COC3" par: "CHEMBL3954592" n: 1 
iso: ATT smi: "[1NH2+]1CCC2(CC1)C1=C(C=CC=C1)CCN2" par: "CHEMBL1183307" n: 4 
iso: ATT smi: "C1=C2C3C(C=C1)[1CH]=[1CH]N3[1CH2]CO2" par: "CHEMBL4439452" n: 1 
```
Where which was generated from
```
mol2SAFE -I 1 -S /tmp/Sfile -y -M maxnr=10 -M 16 -c -v chembl.smi
```

### -M
This only applies with the -S option is specified.
Unless constrained, mol2SAFE will generates and accumulates arbitrary sized fragments,
which may not be desirable. We can limit the size of the fragments written via the `-S`
option with the `-M` option. Specifying a number imposes a limit on the overall
size of fragments accumulated. But what is often found in cases like this is that
having a larger number in order to accommodate certain ring systems is desirable,
this may then enable processing of long chain fragments. So the construct `-M maxnr=nn`
is also available, which allows a constraint on the number of non ring atoms. I have
found the following combination useful.
```
-M maxnr=10 -M 16
```
In this case the largest fragment contains 16 atoms, but within that constraint,
fragments containing more than 10 non-ring atoms will be discarded.

This also has the very desirable effect of limiting the number of fragments
generated. When doing frgament replacement, it is usually desirable to use
well precedented fragments as replacements.

### -g
Apply chemical standardisation to the incoming molecules. Generally this will be
desirable since the default fragmentation scheme built into mol2SAFE works best
with LillyMol standardised molecules.

For example, it is assumed
that all Nitro groups are represented in their `O=N=O` form, rather than O=[N+]-[O-]`. If presented
with the latter form, the `N-O` bond will be broken. Avoid this by adding the `-g all` option,
which will transform the molecules into LillyMol standardised form
[standardise](/docs/Molecule_Lib/chemical_standardisation.md).

Note too that one of the transformations
performed is to remove explicit Hydrogen atoms, which will also likely be desirable.

### -z
Ignore connection table errors in the input. LillyMol may be unable to read
certain aromatic smiles, due to differing rules about aromaticity. By default,
processing ceases once a connection table error is encountered. If the `-z`
option is used, the unreadable smiles is ignored.

### -T 
Element transformation options. In most QSAR datasets it will not be
advantageous to differentiate among the heavy halogens, so considering
I, Br and Cl to be equivalent will often be beneficial
```
-T I=Cl -T Br=Cl
```
does this

### -t
This is a testing option. In order to ensure correctness of the result,
the resulting SAFE form can be interpreted as a molecule, and that can
be compared with the starting molecule. If the SAFE transformation has been 
implemented correctly, it should yield the same molecule. For the
Chembl dataset being used for testing, this is the case.

Note that in order for this to work, chirality must be removed and
no isotopes applied.

```
mol2SAFE -t -v file.smi
```
should do it.
### -v
Verbose output. If this is in effect, statistics on overall
data generation will be reported, otherwise the tool operates silently.
