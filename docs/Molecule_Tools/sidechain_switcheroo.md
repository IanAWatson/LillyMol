# sidechain_switheroo
This is a simple *de-novo* tool that can quickly generate large numbers
of molecules.

The tool works by first removing sidechains from molecules, and then re-attaching
those sidechains across all molecules in the input.

For this reason, it works by reading a set of input molecules into memory. Unlike
most LillyMol tools, it cannot work on a semi-infinite infinite stream of molecules
read from stdin.

A query is used to identify the bonds separating the 'scaffolds' from the
'sidechains'. By default this is `[aD3]-!@{a<10}[R0]`. This means that only
aromatic substituents are identified - by default. Also, substituents
are defined as having fewer than 10 atoms. Adjust to taste.

Any smarts can be specified via the -s option, with the first matched atom
assumed to be an atom in the scaffold, the second matched atom being an
atom in the sidechain, with the bond between them being broken.

During testing the following smarts were found to work well.
```
sidechain_switcheroo -s '[aD3]-!@{a<10}[R0]' -s '[CG0D3R0]-!@{a<10}[R0]' -h -z i -V -x 3 -v file.smi > new.smi
```
This breaks bonds between an aromatic atom, as well as substituents attached to 
three connected, fully saturated, non ring carbon atoms. Again, adjust to taste.

Note however that given just 10 random Chembl molecules as input, the command above
generated 5900 molecules in about 0.4 seconds. This tool needs to be used carefully
or else numbers may explode.  

The `-h` option adds Hydrogen as a substituent, which means that some of the
molecules generated will have lost a substituent at a site where previously there
was a substituent.

The `-x` option is important. This specifies the maximum number of scaffold sites
that are simultaneously altered. This is one of the most significant sources of the 
combinatorics problems with this tool. By default, without the -x option, single
sidechains will be moved from one site to another - inter or intra molecular
transfers. The `-x` option specifies how many simultaneous transfers can be done.
The random 10 Chembl molecules used for testing, generate 14 unique sidechains, but
the number of ways one can select 3 from 14 is large (2184).

Output is straightforward, sidechains will be removed from one location and
moved to other sites - either in the same molecule or in other molecules.
