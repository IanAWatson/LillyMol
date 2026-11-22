# numbered_smiles

## Purpose
`numbered_smiles` is a tool intended for use by expert. It consumes a structure
file and writes output that will have either isotopic labels, or atom map numbers
corresponding to whatever kind of numbering is specified on the command line. This
can be very useful for investigating complex structural issues.

In its most simple form, given ethane, `CC`
```
numbered_smiles ethane.smi
```
will generate a file `ethane_num.smi` that contains `C[1CH3] ethane`, where atom
0 has been assigned isotope 0, no isotope, and atom 1 has isotope 1.

Without `-S`, each input gets its own output beside the input file. Only the
filename's final extension is removed before `_num` and the output-format suffix
are added. For example, `run.v2/sample.v3.smi` produces
`run.v2/sample.v3_num.smi`, or `run.v2/sample.v3_num.sdf` with `-o sdf`.

An explicit `-S` combines inputs under that output stem; `-S -` writes to stdout.
Before processing, planned output names are checked against every input and
against each other. Conflicting names produce an error before outputs are opened.

`numbered_smiles` is used for the `-N` option to vf.

## Options
The help message is
```
  -n atom        number by atom number (the default)
  -n symm        number by symmetry class
  -n canon       number by canonical order
  -n ring        number by ring system
  -n <number>    number atoms starting at a non-negative <number>
  -e <step>      atom-number increment; multiplier for other numbering modes (default 1)
  -P <atype>     number by atom typing (no -n option needed)
  -I             overwrite any existing isotopic information
  -r             stably sort atoms by their assigned labels
  -H             make all implicit Hydrogens explicit
  -l             reduce to the largest fragment before processing
  -b ...         write bond symmetry information to stderr
  -m             overwrite atom maps with numbering; preserve isotopes
  -S <stem>      specify output file name stem
  -K ...         standard smiles control options
  -i <type>      specify input type
  -A <qualifier> Aromaticity, enter "-A help" for options
  -v             verbose output
```


## Numbering and sorting

Atom-number labels are `start + atom_index * step`. The default start is zero,
so the first atom normally has no visible isotope. Starting values must be
non-negative and the step must be positive. For example:

```shell
numbered_smiles -n 1 -e 3 -S - input.smi
```

For `CCO`, this writes `[1CH3][4CH2][7OH]`. In other numbering modes, `-e`
multiplies the computed value: `-n z -e 3` labels carbon 18 and oxygen 24.
Values exceeding the isotope or atom-map range are rejected rather than wrapped.

Existing isotopes are preserved unless `-I` is specified. In atom-map mode
(`-m`), all atom maps are overwritten, isotopes are preserved, and `-I` has no
effect. For example, `-m -n 1 -e 3` assigns maps 1, 4, 7, regardless of input
isotopes or existing maps.

`-r` sorts the atom table by the labels actually assigned or preserved, in
ascending order. Equal labels retain their input order, including zero ring-system
labels on non-ring atoms. A SMILES traversal need not follow atom-table order;
SDF output shows that order directly.

The previously accepted `-q` and `-s` options had no implementation and have
been removed. Use `tsubstructure` for query-based selective labeling.

Use `-l` to retain the largest fragment before numbering or other molecule
processing. Fragment selection uses `Molecule::reduce_to_largest_fragment_carefully()`.
Without `-l`, all fragments are retained.

```shell
numbered_smiles -l -n 1 -S numbered salts.smi
```

Missing or unreadable input files, fatal molecule-reading errors, invalid options,
and output failures produce a nonzero exit status. Buffered output is checked
before switching output files and before successful completion. Input-error
tolerance remains available through the usual `-i` directives, such as `-i ICTE`.
