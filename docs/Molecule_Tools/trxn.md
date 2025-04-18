# trxn

`trxn` makes molecules based on a reaction specification
and one or more sets of input molecules. It can operate on a
single molecule - performing some kind of intra-molecular
change, or to perform a reaction involving any number of
different sidechains.

While `trxn` does support `smirks` as a reaction specification,
the `-K` option, it works best when dircted by a reaction file.
A reaction file precisely describes the changes to be performed
on each atom matched by a query. For example in an acid + amine
reaction, the smarts for the acid might be `[OH]-C=O`. The
reaction file will
need to specify that the bond between matched atoms 0 (`[OH]`) and
1 (`C`) is to be broken, and that matched atom 0, the `[OH]` is 
to be removed. While this may seem tedious at first, the
precise degree of control available becomes an advantage.

Whenever using `trxn` this correspondence between matched atom
number, starting at zero, is crucial. All directives for changes
are specified based on one or more matched atom numbers.

For example, an simple acid-amine reaction might look like

```
name: "Acid + amine"
scaffold {
  id: 0
  smarts: "[OH]-[CD3]=O"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "[NH2]-[CX4]"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
The `name` field is optional.

The order of the atoms in the smarts defines the atom numbering
for the directives that change the molecule. Again, atom numbers
start with 0.

In this case we chose to call the acid the
scaffold, and the amine the sidechain. That is entirely arbitrary, but it
does mean that in order to perform this reaction, the command line must
look something like
```
trxn -P acid_amine.rxn acid.smi amine.smi
```
where the order of the input files on the command line corresponds to
the ordering within the reaction file.

The same effect can be achieved by reversing the role of scaffold and
sidechain,
```
name: "Amine + Acid"
scaffold {
  id: 0
  smarts: "[NH2]-[CX4]"
}
sidechain {
  id: 1
  smarts: "[OH]-[CD3]=O"
  remove_atom: 0
  join {
    a1: 0
    a2: 1
    btype: SS_SINGLE_BOND
  }
}
```
and then running
```
trxn -P amine_acid.rxn amine.smi acid.smi
```
where again, the order of reagents in the reaction file must correspond
with what is given on the command line.

If there is a single sidechain, the same reagent is being reacted with
each scaffold, that can be specified in the reaction file via a
`reagent` directive in the sidechain. For example
a Buchwald reaction that just uses di-ethyl amine 'CCNCC' might be

```
name: "Buchwald"
scaffold {
  id: 0
  smarts: "Br-c"
  remove_atom: 0
}
sidechain {
  id: 1
  reagent: "CCNCC diethylamine"
  smarts: "[ND2](-[CX4])-[CX4]"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
and in this case `trxn` would be invoked with just the file of aryl Bromides
```
trxn -P buchwald_fixed.rxn aryl_bromine.smi
```

Note that with both of these reactions there is a single atom lost, which
is discarded via the `remove_atom` directive. If that was not specified,
those atoms would remain with their original bonding - likely leading
to a valence error. If you instread wish to preserve the lost atoms as
a separate fragment, you can break the bond between the lost atom and
the rest of the molecule.

The fixed reagent Buchwald reaction would be
```
name: "Buchwald"
scaffold {
  id: 0
  smarts: "Br-c"
  break_bond {
    a1: 0
    a2: 1
  }
}
sidechain {
  id: 1
  reagent: "CCNCC diethylamine"
  smarts: "[ND2](-[CX4])-[CX4]"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
which generates products that include disconnected Bromine atoms.

Once a bond has been broken, if that defines a fragment, rather
than atom, to be lost, the remove_fragment directive will remove all atoms
in the fragment defined by the matched atom.

Note that multiple atoms (or fragments) can be removed, since in 
the proto definition file, [reaction.proto](/src/Molecule_Lib/reaction.proto)
the `remove_atom` and `remove_fragment` directives are repeated fields. Removing
multiple atoms can be either
```
scaffold {
  remove_atom: 0
  remove_atom: 1
}
```
or
```
scaffold {
  remove_atom: [0, 1]
}
```
where the atoms to be removed are entered either as separate directives or
as a vector of matched atoms.

## Combinatorial Generation
`trxn` was originally developed to process combinatorial libraries
and as such, can be fast if doing that. For historical reasons,
it divides reagents into one called `scaffold` and all others
are referred to as `sidechains`. This is an artificial
distinction and is not necessary.

That said, there are some important differences between the
scaffold and the sidechains. For example, on startup, all sidechains
are read at once. During this time, their substructure query is run
in order to determine the atoms involved in the reaction. That
information is stored, so the sidechain can be cheaply re-used
against any number of scaffolds. The expensive substucture matching
is done once. On the other hand, the scaffold file is read
one molecule at a time.

In large libraries, this information retained with the sidechains
may cause memory problems. Split the library into chunks, which
is probably a good thing for cpu considerations anyway.

For each `scaffold` or `sidechain` specification in the reaction
specification, there must be a corresponding file on the command
line containing those structures. The first file will be the 
scaffold molecules, and subsequent files belong to each sidechain.
Any number of sidechains can be specified.

## Historical
The original reaction specification scheme for `trxn` was based on an
old file format from Cerius-2. Those reaction files are still recognised
but no new functionality is being added to them. They are not documented
here.

## Simple Reaction
The only mandatory component of a reaction is a scaffold directive. For example
the alchemists's dream reaction would be the unimolecular reaction
```
scaffold {
  id: 0
  smarts: "[Pb]"
  change_element {
    atom: 0
    element: "Au"
  }
}
```
which transforms Lead to Gold.

The important principle is that the smarts establishes an atom numbering.
In this case, there is only one atom in the smarts, and therefore matched
atom 0 is defined. In this case, matched atom 0 will be a Lead atom.

The `change_element` directive specifies that matched atom 0, the Lead atom,
should be changed to Gold.

The Bootlegger's reaction for making alcohol might (naievely) be implemented
by adding an Oxygen atom to a Ethane molecule
```
name: "bootlegger"
scaffold {
  id: 0
  smarts: "CC"
}
sidechain {
  id: 1
  reagent: "O"
  smarts: "[O]"
  join {
    a1: 0
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
Unfortunately if this is run with default settings
```
trxn -P bootlegger.rxn ethane.smi
```
the product will be `OCCO`. This is because the smarts for ethane `CC` matches
the ethane molecule two ways because of symmetry. Add the `-u` option to the
`trxn` invocation to suppress these duplicate atom matches. Or the `-k` option
which suppresses symmetry related matches such as these.

Alternatively we can add a directive to the scaffold message to only find unique
embeddings
```
name: "bootlegger"
scaffold {
  id: 0
  smarts: "CC"
  match_conditions {
    find_unique_embeddings: true
  }
}
sidechain {
  id: 1
  reagent: "O"
  smarts: "[O]"
  join {
    a1: 0
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
and then no command line options are needed. Both means accomplish the
same effect. Longer term, we hope to make most command line options
also settable from within the reacton file.

In order to hide evidence of a crime, the Bootlegger may wish to
denature their ethanol by converting it to methanol!
```
name: "ethanol to methanol"
scaffold {
  id: 0
  smarts: "[OH]-CC"
  break_bond {
    a1: 1
    a2: 2
  }
  remove_atom: 2
}
```
The bond between matched atoms 1 and 2, the two carbon atoms, is broken.
Then matched atom 2 is removed - otherwise the product would be a disconnected
molecule.

Formal charges are straightforward.
```
name: "charge acid"
scaffold {
  id: 0
  smarts: "[OH]-C=O"
  formal_charge {
    atom: 0
    formal_charge: -1
  }
}
```
The smarts defines atom numbering and the `formal_charge` directive specifies
which atom gets assigned the formal charge.

The `change_formal_charge` directive can be used to add or subtract a value
from an existing formal charge. Similarly the `change_isotope` message can
be used to change an existing isotopic value
```
name: "increment isotope"
scaffold {
  id: 0
  smarts: "a-[9]"
  change_isotope {
    atom: 1
    delta: 1
  }
}
```
will change any atom, singly bonded to an aromatic atom, that has isotope 9
to isotope 10 (delta is 1). Delta can also be negative and it is a fatal error
to form a negative isotopic value.


## Fixed Reagent
A sidechain may involve a single reagent that is added to each
molecule. You can either place that single molecule in a file, or
specify that reagent in the sidechain. In proto notation, that sidechain
might look like
```
sidechain {
  id: 2
  reagent: "[U]"
  smarts: "[U]"
  join {
    a1: 0
    a2: 0
    btype: SS_DOUBLE_BOND
  }
}
```
which adds a doubly bonded Uranium atom (!!) to matched atom 0
in the scaffold (a1). `trxn` should be able to discern that this
sidechain does not need a file on the command line.

`trxn` is necessarily complex, because reaction enumeration
can be very complex.

NOTE: It should be considered essential to review structures
generated by `trxn`. And don't just look at the top of the file,
view samples from throughout the file. Few of us can write a
reaction file that works first time. This is hard.

# Smiles Quirk
While trxn can be very efficient, there is a feature of the smiles
language that can lead to very efficient reaction enumeration. The
smiles 'C1.C1' is a valid smiles for ethane, and `C1CC2.C1CC2`
is the same as `C1CCCCC1`, etc. So if you were to start with
isotopically labelled reagents, you could use text processing
to convert those isotopes to ring opening/closing numbers, join
the reagents as text, and have `fileconv` read the resulting
smiles. This will always be faster than `trxn` although it does
require some potentially error prone preprocessing.

# HOWTO
The following options are recognised
```
  -r <file>      specify single reaction file
  -D <file>      specify ISIS reaction file
  -D rmfrag      remove small fragments from the results of ISIS reaction file reactions
  -D rmunmp      remove unmapped elements shown on LHS but absent on RHS
  -P <file>      reaction as proto file
  -K <smirks>    reaction specified by smirks string. F:<fname> means smirks is in file
  -z i           ignore molecules not reacting
  -z w           write molecules not reacting
  -Z             ignore sidechains not reacting
  -C <string>    append <string> to the name of all changed molecules
  -C ifmult      only append the -C string in the case of multiple scaffold matches
  -m <number>    the maximum number of scaffold reaction sites to process
  -m do=number   process site number <number> in the scaffold
  -m each        enumerate each scaffold hit separately
  -m RMX         ignore scaffolds that generate multiple substructure hits
  -X <symbol>    extract/remove all atoms of type <symbol>. No bonds changed
  -I             change isotopes to natural form in product molecules
  -M all         generate all regio-isomers from multiple sidechain matches
  -M do=number   process site number <number> in the sidechains
  -M mskip=text  append <text> to names where just one possible sidechain attachment chosen
  -M write=file  write non-reacting sidechains to <file>
  -M RMX         ignore any sidechains with multiple substructure matches
  -V <fname>     molecules with invalid valences ignored and written to <fname>
  -l             strip reagents to largest fragment
  -L             strip products to largest fragment
  -f             function as a TDT filter
  -n <xxx>       number assigner options, enter '-n help' for details
  -W <string>    token put between names of products (default " + ")
  -u             one embedding per start atom
  -d             suppress duplicate molecules - only checks current molecule
  -k             don't perceive symmetry equivalents in the scaffold
  -J ...         various special purpose options, enter '-J help' for details
  -E <symbol>    create an element with symbol <symbol>
  -E autocreate  automatically create new elements when encountered
  -o <type>      specify output file type(s)
  -S <string>    create output files with name stem <string>
  -A <qualifier> Aromaticity, enter "-A help" for options
  -g <qualifier> chemical standardisations, enter "-g help" for usage
  -i <type>      specify input file type
  -v             verbose output
```

# Reactions
The reaction specification can be any of

1. Legacy reaction file format `-r`
2. Text format Proto reaction file format `-P`
3. ISIS reaction file format `-D`
4. smirks `-K`.

The legacy reaction file format is well tested, but no new features are
added to it. The proto file format is newer, and therefore less tested,
but new features will only be supported here. ISIS reaction files have
been used extensively with `trxn` although no new developement has
occurred for some time. smirks support is provided, although not
extensively tested, but seems to work well.

People have asked why the need for more verbose reaction representations
than smirks. Primarily there are concepts that would be difficult to
express with smirks, increment an isotope for example, or use any of
the more complex substructure concepts enabled by query files. But the main
reason is that like smarts, smirks can be thought of as a write-only
language, sometimes being very difficult to deciper once written. In
smirks, the changes needed have to be implied by looking at the
state of the atoms on the lhs and rhs of the reaction, whereas 
the approach here requires very specific instructions for any change
made. Each approach has advantages and disadvantages. The clear disadvantage
of reaction files is their verbosity. But that disadvantage becomes an
advantage in clarity and control.

The other advantage of the approach taken is that the query that 
defines a reaction can be available to other tools. In most systems
the smirks is the only instantiation of the query defining a reaction,
and ends up being copied and pasted if used elsewhere. With `trxn`
queries can be put in query files, and used by any tool.

Note that `trxn` only supports one reaction at a time. Another tool,
`molecular_transformations` is designed to process multiple reactions.
`make_these_molecules` enumerates specific subsets from larger
enumerations, again with multiple reactions.

# Multiple Matches
One of the hardest things with any enumeration is the problem of
multiple matches. This can arise from one of two sources.

1. The reagent has multiple functional groups.
2. Your query to define the reagent is not specific enough.

While the second one can be fixed, the first one is harder. If there
are multiple functional groups, where reactivity at each site is 
equally likely, how should this be handled?

# Best Practice
Generally things will be easiest if you can pre-process all reagents
prior to starting `trxn`. Use `tsubstructure` with numeric qualifiers
to identify reagents that have the desired functionality
```
tsubstructure -s '1[OH]-C=O&&0[OH]-S=O' -g all ...
```
will look for molecules that contain a single carboxylic acid, and zero
occurrences of a sulf* acid.

If you can do all duplicate functionality processing *before* getting
to `trxn`, interacting with `trxn` will be much easier.

But there are of course times when multiple functional groups are
present and must be handled.

# Multiple Functional Groups
This is probably the hardest part of reaction enumeration, there are many things
that can go wrong.

`trxn` differentiates multiple functional groups in the scaffold and
in the sidechains. Multiple functional groups in the scaffold are
handled with the `-m` option, and Multiple functional groups in sidechain(s)
are handled with the `-M` option. The same `-M` directive applies to all
sidechains.

Note that there are options to matching that can turn what might otherwise
be multiple matches into unambiguous matches.

### -u
When embeddings are found, only find one embedding per start atom. Normally
a search for a nitro `C-N(=O)=O` would find two embeddings. Suppress this
with the `-u` option.

### -k
Do not find symmetry equivalent matches. If there is only one sidechain
this should be the default, but if there are multiple sidechains, it is
unclear what the best approach should be.

There are also options under `-J` that allow further refinement of how
to deal with multiple matches.

## Scaffold
The following options govern multiple functional groups in the scaffold.
```
  -m <number>    the maximum number of scaffold reaction sites to process
  -m do=number   process site number <number> in the scaffold
  -m each        enumerate each scaffold hit separately
  -m RMX         ignore scaffolds that generate multiple substructure hits
```

Note that the default behaviour is to react all scaffold sites at once,
generating a molecule that will have sidechains attached to each possible
site. This may be what you want, and is a good default.

### -m \<number\>
The maximum number of scaffold reaction sites to process. If there are more
than \<number\> embeddings of the scaffold query, the first \<number\> will
be used. Note that this will be an arbitrary choice, dependent on the
ordering of the atoms in the connection table. Probably not a great idea.

### -m do=number
If there are multiple matches in the scaffold, process match number
\<number\>. Since which embedding is found first is an arbitrary function
of what is in the connection table, this is probably not a great idea
either.

### -m each
Enumerate each scaffold hit separately. This is also a commonly
used approach. For each scaffold site, make a copy of the scaffold
and place the sidechain at just that one site.

### -m RMX'
Ignore scaffolds that generate multiple substructure hits. Again a
commonly used approach. The best practices section previously
recommended using `tsubstructure` to eliminate these problems
before coming to `trxn`.

## Sidechains
There are different options for specifying how multiple matched atoms in
the sidechain might influence the enumeration. The following options control
sidechain matching
```
  -M all         generate all regio-isomers from multiple sidechain matches
  -M do=number   process site number <number> in the sidechains
  -M mskip=text  append <text> to names where just one possible sidechain attachment chosen
  -M write=file  write non-reacting sidechains to <file>
  -M RMX         ignore any sidechains with multiple substructure matches
```
### -M all
Generate all regio-isomers from multiple sidechain matches. Since there
is only one attachment point in the scaffold, what this does internally
is to create copy of that sidechain molecule, that uses the other
mathed atoms as the reaction site.

### -M do=number
Process site number \<number\> in the sidechains. Again, this will be an
arbitrary choice. This might be useful for exploratory work, but not
for a production enumeration.

### -M mskip=text
Append <text> to names where just one possible sidechain attachment chosen.
This is used in conjunction with the `do=` directive. Again, for 
exploration only. This should move to the `-J` miscellaneous option.

### -M write=file
Write non-reacting sidechains to \<file\>. This is better done with
`tsubstructure` as a preprocessing step. If this is retained, it should
move to the `-J` miscellaneous option. Do not depend on this.

### -M RMX
Ignore any sidechains with multiple substructure matches. Again,
using `tsubstructure` for preprocessing might be simpler, although
this is certainly convenient.

In practice only three of these possibilities are used.

1. Leave out the `-M` option, sidechains are all single function.
2. `-M all` is used. Enumerating all possibilities is what is wanted.
3. `-M RMX` ignore sidechains with multiple matches.

But again, best practice would be to prepare reagents with
`tsubstructure` first, so that there is no need for the `-M` option,
or possibly `#2` above.

# Unreactive Molecules
The opposite of molecules not matching the reaction queries multiple
times is the case of molecules that do not match the queries at all.
While it might be desirable to eliminate these up front, they can
be conveniently handled with `trxn`.

## -z i
Ignore scaffold that do not match the scaffold query.
## -z w
Write the non matching scaffolds to the ouput. This is a convenient
way of selectively applying a transformation to some of the molecules
in a set.
## -Z
Ignore sidechains that do not match their query.

# Matched Atoms
Fundamental to the operation of `trxn` is the concept of a matched atoms.
Substructure queries (often smarts) define the reaction, and each substructure
query defines a numbering. For example if the query is for an acid,
`[OH]-C=O` then the `OH` is matched atom 0, the Carbon is matched atom 1
and the second Oxygen is matched atom 2. All directives in reaction files
depend on these atom numberings.

In addition, each scaffold and sidechain has its own substructure query,
which defines a numbering within that entity. So across all reacting
species, each atom can be defined by the id of the component, and the
matched atom number.

When performing operations across different components, there are conventions
that apply. The most common inter-component directive is `join`, that
makes a bond between a matched atom in the sidechain and a matched atom
in the scaffold. By convention, the first atom in a `join` directive is
assumed to be a matched atom number in the scaffold, and the second is in
the sidechain in which the directive is found.

So,
```
join {
  a1: 2
  a2: 1
  btype: SS_SINGLE_BOND
}
```
is a directive to create a single bond between matched atom 2 in the
scaffold and matched atom 1 in the current sidechain. This is equivalent to
```
join {
  c1 {
    component: 0
    atom: 2
  }
  a2: 1
  btype: SS_SINGLE_BOND
}
```
Note that rather than `a1`, which means the first atom (assumed to be in the scaffold),
we have `c1` which
means the first component in the bond. So if the join directive uses `a1` it is
assumed to be a matched atom number in the scaffold, which is component 0. Otherwise
it can be `c1` in which case it fully describes the component number and matched
atom number.
That is the InterParticleBond proto, which has a
oneof.  Clearly `a2` can be changed as well, so if the current
scaffold is `id` 2, is equivalent to
```
join {
  c1 {
    component: 0
    atom: 2
  }
  c2 {
    component: 2
    atom: 1
  }
  btype: SS_SINGLE_BOND
}
```
This is the most general specification of matched atoms within `trxn`.
So it becomes possible to create bonds between sidechains.

```
join {
  c1 {
    component: 1
    atom: 2
  }
  a1: 3
  btype: SS_DOUBLE_BOND
}
```
creates a double bond between matched atom 2 in component 1, and matched
atom 3 in the current component.

Generally directives are owned by individual scaffolds and sidechains. And
usually this is clearest. Remove matched atom 2 is unambiguous and easy, it
applies to the current component. But as might be inferred from the above,
there is nothing stopping component 1 from specifying an operation involving
components 0 and 2. 

But a more interesting design would be to define a series of reagents,
and then all directives that could possibly involve multiple components
would be owned by the reaction itself, and all
matched atoms would need to be qualified with a component. Directives
that were local to a component, `break_bond` would not be available
globally. But the current form has advantages too.

# Output Filters and Processing
Several options are available for altering and filtering the
molecules generated.

## -V \<fname\>
When molecules with bad valence are generated, write to \<fname\> and
not to the default output stream. Unfortunately, it is very common
to accidentally generate molecules with invalid valence and
examining and debugging is a common undertaking. Having them
written to a separate file aids rapid debugging.

If \<fname\> is 'NONE' then products with invalid valences are discarded,
but not written - no discard file is opened.

## -d
For the current scaffold, suppress any duplicates formed. Note that
each scaffold is reacted with all sidechains, so it is quite
possible that duplicate products will appear. For example if there
are halogen duplicates in the reagents, I and Br, I and Cl, Cl and Br.
These may react and all generate the same product.

Note that there is no provision for removing duplicate molecules
across scaffolds. For that, pipe the output of trxn to
`unique_molecules`.

## -I
Remove any isotopic labels in the products. Isotopes are frequently
used to mark reacting atoms. The reaction can specify that isotopes
are set back to zero, but it is also convenient to just have
`trxn` itself remove them when done.

## -J maxat=nn
Discard any product molecule having more than \<nn\> atoms. Combinatorial
libraries can generate large molecules. It may be desirable to suppress
them.

## -J rmph
Remove explicit hydrogen atoms from product molecules.

## -J rmhsqb
Remove unnecessary [] in product molecules. This is an area of
frustration and ambiguity. When molecules enter with isotopic
labels, those atoms necessarily have square brackets, which
among other things indicates that the Hydrogen count is known.
How should these atoms be processed during a reaction? That is
unclear. Most of the time we will want to simply discard
that hydrogens known information. Use this option.

## -J rmxhbv
Remove explicit hydrogens causing bad valences. If the product molecule
contains explicit Hydrogen atoms attached to atoms that have
a bad valence, try removing them one at a time until the
valence of the atom is OK. Note that it might remove all
the explicit hydrogen connections and still have a bad 
valence.

## -J coords
Like all LillyMol tools, trxn is 3-D capable. With this option in effect,
3-D smiles will be written if smiles output is requested.

## -J minpfs=\<n\>  -J maxpfs=\<n\>
Reactions can generate fragments that are too small (or too large). Sometimes
it is possible to construct the reaction queries so that undesirable
fragment sizes are not generated, but sometimes this is challenging. Instead,
these directives allow filtering of product molecules if they contain
fragments that violate these directives.

. minpfs=\<n\>  If the product contains a fragment with fewer than \<n\> atoms,
then discard the entire product molecule.
. maxpfs=\<n\>  If the product contains a fragment with more than \<n\> atoms,
then discard the entire product molecule.

## -L
If there are multiple fragments generated by the reaction, only retain
the largest. In the reaction specification, fragments can be removed,
but it can also be convenient to just have `trxn` remove small
fragments when done. Just make sure that your desired product molecule is
larger than anything else that might be produced!



# Output
In a rather unconventional way, `trxn` has a default output file, `trxn.smi`.
This is inconsistent with most LillyMol utilities, but is retained for
historical reasons. This can be changed with the `-S` directive. Note that
the output format can be specified with the `-o` option. If you are 
processing 3D molecules, you will likely need to add `-o sdf` in order
to have 3D data generated. As is usually the case, multiple `-o` options
can be specified, so `-o sdf -o smi` will generate two output files.
The `-J coords` option will result in smiles containing coordinates for more
compact 3d output.

## Names
By default, the name of each new molecule is the concatenation of the
names of the scaffold and sidechain(s) that generate the new molecule. The
separator between those components is ` + ` by default. That can be
changed with the `-W` option.

As an alternative, a sequential number can be assigned to each generated
molecule with the `-n` option. This is seldom used.

In the case where non reacting molecules are being written, the `-C` option
can be used to append a string to the names of those molecules that
have been changed.

The combination `-C ifmult` means only append the -C string in the
case of multiple scaffold matches. This way names that will be duplicates
can be quickly identified.

There is no 'right' way to name the products of an enumeration, but
trxn allows several possible approaches. The concatenation of names
default is usually adequate.

# Miscelaneous Options.
The `-J` option hides several special optional forms.
```
 -J wrscaf      write scaffolds  to the output stream
 -J wrsdch      write sidechains to the output stream
 -J blki        bonds lose Kekule identity
 -J appinfo     append text info from reagents to products
 -J onlyreact=RX only react molecules where the name matches RX
 -J maxat=nn    discard any product molecule having more than <nn> atoms
 -J rmncm       ignore multiple substructure matches involving non changing atoms
 -J rmovm       ignore multiple substructure matches involving     changing atoms
 -J exph        make implicit Hydrogen atoms explicit (changes reaction)
 -J exphR       make implicit Hydrogen atoms explicit on reagents (reaction not changed)
 -J rmph        remove hydrogen atoms from product molecules
 -J rcksm       when multiple scaffold hits present, re-check matches for activity
 -J numok       keep non-unique embeddings - default is unique embeddings only
 -J isonum      isotopically label atoms with their initial atom number
 -J msm=<s>     text designating multiple scaffold matches, use NONE to skip
 -J marvin      the input reaction file (-D) has come from Marvin
 -J keepatmn    retain any atom map numbers in output molecules
 -J larf        in smirks, if an atom is lost, remove the fragment
 -J rmhsqb      remove unnecessary [] in product molecules
 -J rmxhbv      remove explicit hydrogens causing bad valences
 -J minpfs=<n>  discard products with a fragment with < minpfs atoms
 -J maxpfs=<n>  discard products with a fragment with > maxpfs atoms
 -J mfpseparate write multi fragment products as separate molecules
```

## -J wrscaf
Write scaffolds to the output stream. This can be helpful when
a single stream of molecules is being procesed and you want to
have the starting form retained for comparison.

## -J wrsdch
Write sidechains to the output stream. Useful for debugging.

## -J blki
Bonds lose Kekule identity. By default LillyMol aromatic
forms can retain their single/double bond status. Sometimes
that is useful but mostly not. Turn that off with this option.

## -J appinfo
Append text info from reagents to products. LillyMol molecules
can have an array of text records - frequently read from either
an SDF file or a TDT. If this is specified, that extra information
is transferred during the reaction process.

## -J onlyreact=RX
Only react molecules where the name matches RX, a regular rexpression.
This can sometimes be convenient if the set of desired
reagents can be specified with a regular expression. For smiles however
using grep on the input file also works.

## -J maxat=nn
Discard any product molecule having more than \<nn\> atoms.

## -J rmncm
Ignore multiple substructure matches involving non changing atoms.
Multiple matches are the bane of enumeration. But if multiple
embeddings only involve non reacting atoms, the problem may
be less severe. This is an obscure option, but sometimes useful.

## -J rmovm
Ignore multiple substructure matches involving changing atoms. Obscure.

## -J exph
Make implicit Hydrogen atoms explicit (changes reaction).

## -J exphR
Make implicit Hydrogen atoms explicit on reagents (reaction not changed).
For this to work, you will need to make sure the reaction specification
anticipated explicit Hydrogen atoms.

## -J rmph
Remove explicit hydrogen atoms from product molecules.

## -J rcksm
When multiple scaffold hits present, re-check matches for activity. There
can be cases of multiple reacting sites where performing the first
change invalidates the query match that was discerned for another match.
This says to re-perform the substructure search after each site is
reacted. Note that things could go the other way too, a change might
induce multiple matches, but that is not checked.

## -J numok
Keep non-unique embeddings - default is unique embeddings only.

## -J isonum
Isotopically label atoms with their initial atom number. Useful for
debugging.

## -J msm=\<s\>
Text designating multiple scaffold matches, use NONE to skip.

## -J marvin
The input reaction file (-D) has come from Marvin. This enables
processing of some Marvin specific features.

## -J keepatmn
Retain any atom map numbers in output molecules.

## -J larf
In smirks, if an atom is lost, remove the fragment.

## -J rmhsqb
Remove unnecessary [] in product molecules.

## -J rmxhbv
Remove explicit hydrogens causing bad valences.

## -J minpfs, maxpfs
These are used to discard products where any fragment violates these atom
count constraints. This is most commonly used when using `trxn` to identify
substituents and fragments. The same effect could be accomplished
via post-processing with `fileconv` but this will be more straightforward.

## -J mfpseparate
If the product molecule contains multiple fragments, write them as
separate molecules. Note that each fragment must pass any constraints
specified by `minpfs` and/or `maxpfs`. The same effect could be 
accomplished with `mkfrag` to convert these multi-fragment molecules
to fragments, but that is much less efficient.

For example, running a simple bond breaking reaction across 2.2M molecules
took 6 minutes. Then using mkfrag to convert these to separate molecules
took an extra 6 minutes. By using the `-J mfpseparate` option on `trxn`
the total run time was 8 minutes.

## -J nomshmsg
Do NOT write 'hits in scaffold' messages for multiple scaffold query hits.
In large enumerations where this is expected, the messages are uninformative
and slow down processing.

## -J noschmsg
Do NOT write warning messages about no sidechain substructure matches.
Note however that this information can be very important, so only suppress
this is you are confident about the reaction being done.

## -J rpt=<n>
Report progress every <n> products written. By default trxn works silently
and can quickly generate large numbers of molecules. This can help keep
track of progress - helpful too if you know how many products to expect.

## 3D
While the primary purpose of `trxn` is manipulating connection
tables, it is fully 3d aware if the input data contains
3d structural information.

The queries that define the reaction define atom numbers that
are used as directives for breaking and forming bonds, etc..
Those same atom numbers are available for setting bond lengths,
bond angles and torsions. By convention, when a value is set
the first atoms in the specification are held still and the
subsequent atoms move. For example, a bond length is defined
by two atoms and a distance. The first atom will be held
stationary, and the second atom, and all atoms attached to it, will move along
the line defined by the bond. A bond angle
requires three atoms, and in that case, the first two atoms
are held stationary, and the third atom, and all atoms attached to it, move.
Similar for torsions.

For example

```
bond_length {
  a1: 3
  a2: 2
  distance: 1.45
}
```
would see matched atom 3 held in place, and matched atom 2, and all atoms
attached to it,
moved along the existing vector so that the separation becomes 1.45.
Note however that there has been no specification on bond angles or torsions, so this
is very unlikely to be a realistic geometry. Note that the
bond_length directive is a proto oneof, and so rather than
atom numbers components can be specified.
```
bond_length {
  c1 {
    component: 1
    atom: 3
  }
  a2: 2
  distance: 1.45
}
```
shifts the fragment defined by matched atom 2 in our component
so that it is 1.45 Angstroms away from matched atom 3 in component 1.

Bond angles and dihedral angles behave similarly.

### Coordinate Transfer
A common operation is to replace a given sidechain with a different one,
while preserving aspects of the existing positioning. In order to be solvable
a minimum of 3 matched atoms must be specified. If fewer than 3 matched
atoms are used, you will manually need to set the angles. Here is an example
of replacing a fragment, and manually specifying the geometry. 

The example describes an aromatic ring that has a substituent that is 
to be replaced with a different one. Within the scaffold, three matched
atoms are defined, and the exocyclic bond is broken, and the old fragment
removed.

The sidechain is defined by an isotopic label, and joined to the scaffold
in a typical 2d way. What makes this a 3D operation is the `coordinate_transfer`
directive, which says that matched atom 0 in the sidechain should go in the
same location as matched 2 in component 0 (the scaffold). That is the atom
that was removed.
Note that the `bond_length` directive is processed after the `coordinate_transfer`
directive.

Note that this example has only one `coordinate_transfer` directive, so the
orientation of the fragment is undefined - it preserves whatever
orientation it had. This may be desirable or not.
```
scaffold {
  id: 0
  smarts: "a:a!@A"
  break_bond {
    a1: 1
    a2: 2
  }
  remove_fragment: 2
}

sidechain {
  id: 1
  smarts: "[1]~[!#1]"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }

  coordinate_transfer {
    atoms {
      c1 {
        component: 0
        atom: 2
      }
      a2: 0
    }
  }

  bond_length {
    c1 {
      component: 0
      atom: 1
    }
    c2 {
      component: 1
      atom: 0
    }
    distance: 1.42
  }

  bond_angle {
    c1 {
      component: 0
      atom: 1
    }
    c2 {
      component: 1
      atom: 0
    }
    c3 {
      component: 1
      atom: 1
    }
    angle: 108
  }

  dihedral_angle {
    c1 {
      component: 0
      atom: 0
    }
    c2 {
      component: 0
      atom: 1
    }
    a3: 0
    a4: 1

    angle: 90
  }
}
```
This is very complex, but mostly seems to work. Note however that in some
test cases, some dihedral angles came out as -90 instead of 90. The cause
is unclear.

### 3D positioning
If you are starting with two fragments wit arbitrary spatial orientations and positions, and
wish to join them in a spatially reasonable way, one of the fragments will need to
be translated and rotated so the atoms in each fragment are "pointed towards" each
other in a plausible way. This is accomplished with the `align_3d` directive in
the InterParticleBond message, the `join` directive. For example to perform
an acid + amine reaction in 3d one could use
```
name: "acid_amine"
scaffold {
  id: 0
  smarts: "[OH]-C(=O)*"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "NC"
  join {
    a1: 1
    a2: 0
    align_3d: 1.32
  }
  bond_angle {
    c1 {
      component: 0
      atom: 1
    }
    a2: 0
    a3: 1
    angle: 108.0
  }
  dihedral_angle {
    c1 {
      component: 0
      atom: 3
    }
    c2 {
      component: 0
      atom: 1
    }
    a3: 0
    a4: 1
    angle: 145
  }
}
```
This not only positions the amine "pointing towards" the acid, 1,32 Angstroms away, 
but then once that has been done, sets various bond angles and dihedrals - in this case just one of
each, but those are repeated messages. 

This functionality is accomplished by the `Position3D` function, which operates on
fragments within the same molecule.

Within 2D `trxn` a common reaction is acid + amine. That can be accomplished
via reactions, as demonstrated above, or the underlying function can be called directly.
The scenario is to have an acid that is fixed in space, and we have a
list of 3D amines that are randomly positioned and oriented, and which each
need to be positioned near the acid for bond formation.

If doing this in python, a typical sequence of code might look like
```
# Read 3d acid molecule from somewhere. Save the number of atoms.
acid = read_from_somewhere()
atoms_in_acid = acid.natoms()

# Identify the atom in ‘acid’ that will form the bond -
# for simplicity we  assume the [OH] has been stripped off.
# This function presumably does some kind of substructure search,
# or maybe just look for an isotopically labelled atom, or...
acid_atom = identify_acid_join(acid)

# Read a list of amines from somewhere and iterate through them.
for amine in amines:
  # Which atom in this amine is the join atom - the N atom.
  amine_atom = identify_amine_atom(amine)

  # make a copy of the acid so we can add to it and change it
  # and leave `acid` ready for the next fragment.
  mcopy = copy.copy(acid)

  # Add the amine – no bonds created, two fragments in mcopy,
  # no change in `amine` atom positions.
  mcopy += amine

  # Orient the atoms in the fragment
  Position3D(mcopy, acid_atom, distance, atoms_in_acid + amine_atom)

  # Now that it is positioned, add a bond.
  mcopy.add_bond(acid_atom, atoms_in_acid + amine_atom, SINGLE_BOND)

  # Further processing on the newly formed amide `mcopy`.
```
Note that we did not really need to make a copy of `acid`. We could
have added each `amine` to `acid` and then when work began on the
next amine, `acid.resize(atoms_in_acid)` would remove all the atoms
that had been added by the previous amine, leaving the acid in its
original form. Depends on what subsequent processing might have
done with the newly formed molecule.

In the repo there is an example python script that was used for
testing `Position3D`, [test_position_3d](/contrib/bin/test_position_3d.py).
That tool looks for all breakable bonds in a
set of molecules, breaks each of those bonds, randomly positions
one of the fragments, and then re-assembles the starting molecule.

Given 10k random molecules from Chembl, this generates 76k 
broken and rejoined molecules and takes about 10 seconds.

That script also demonstrates how you can them perform a systematic
dihedral scan around the newly formed bond, to whatever degree of
granularity is needed. Of course if you have specified a dihedral
angle in the reaction, this is not necessary.

# Examples

There are a variety of reactions in the `data/MolecularVariants` directory
of this repo. Those are transformations describing possibly interesting
bio-isostere pairs.

## React an acid with an amine.
Note that the amine query is not terribly
precise. At least one hydrogen attached, fully saturated and no attached
heteroatoms.
```
scaffold {
  id: 0
  smarts: "[OH]-C=O"
  remove_atom: 0
}
sidechain {
  id: 1
  smarts: "[NHG0T0]"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
```
Note that this reaction works best with implicit Hydrogens. It should
still work if Hydrogens are explicit, but would leave an unattached
Hydrogen atom if run that way. It also works for secondary amines
because the `H` requirement in the smarts is interpreted
as 'at least 1 Hydrogen attached'.

## Making/Changing Bonds
There are two directives which can be used to add/alter a bond.

The `make_bond` directive was shown above. It works regardless
of whether there is an existing bond between those atoms or not.
If there is no bond there, one will be added. If there is an existing
bond present, the bond type may be changed. On the other hand the
`change_bond` directive only works if the atoms are already bonded.
If the atoms are not already bonded, it is a fatal error. This is
largely a usability issue, because while `make_bond` seems more
like "create a bond", `change_bond` is clearly associated with
changing an existing bond. Syntax of the two messages is the same,
a1, a2 and btype.

## Breaking Molecules
Reactions are just changes to a connection table. They can both form, change or
break bonds, add or remove atoms.

This reaction breaks off a CF3, converts an alcohol to an an aledhyde,
changes a phenol to aniline, and an ethyl substituent to ethylene. Yes this
is a silly reaction.
```
scaffold {
  id: 0
  smarts: "FC(F)(F)-C.C(=O)-[O].[OH]-a.[CD1]-[CD2H2]"
  remove_fragment: 0
  break_bond {
    a1: 1
    a2: 4
  }
  remove_atom: 7
  change_element {
    atom: 8
    element: "N"
  }
  make_bond {
    a1: 10
    a2: 11
    btype: SS_DOUBLE_BOND
  }
}
```
But once the pattern is understood, a great deal of flexibility
becomes available. See `reaction.proto` for other directives
that are available within a component. There is considerable
flexibility dealing with isotopes, but see also the isotopic
processing in both `fileconv` and `tsubstructure`.

## Repeated Fields
Where multiple values are allowed, standard proto text format repeated
field representations can be used. For example to remove a CF3 one could do
```
scaffold {
  id: 0
  smarts: "FC(F)(F)-C"
  remove_atom: 0
  remove_atom: 1
  remove_atom: 2
  remove_atom: 3
}
```
which is the same as
```
scaffold {
  id: 0
  smarts: "FC(F)(F)-C"
  remove_atom: [0, 1, 2, 3]
}
```
which is of course the same as the bond breakage and fragment removal
in the earlier example.


## Unfix Implicit Hydrogens
As mentioned, molecules that enter trxn with square brackets can be problematic.
Those can be dealt with via `trxn` command line options, or with directives.
In this case we are removing an isotope from a Carbon atom that has isotope 3.
```
scaffold {
  id: 0
  smarts: "[3C]"
  isotope {
    atom: 0
    isotope: 0
  }
}
```
If that runs on `[3C]C`, which is deficient in Hydrogens, it produces '[C]C` which
I think is correct. The input molecule did not have a proper hydrogen count,
changing the isotope has nothing to do with that.

If we add an `unfix_implicit_hydrogens` directive,
```
scaffold {
  id: 0
  smarts: "[3C]"
  isotope {
    atom: 0
    isotope: 0
  }
  unfix_implicit_hydrogens: 0
}
```
Now the output contains no square brackets. The unfix_implicit_hydrogens directive
was applied to matched atom 0.


## Inactive
Dealing with differential reactivity can be very complex.
Generally, best practice would be to concentrate such efforts
on the reagent preparation step, where a great deal of detailed
effort can be devoted to identifying the reacting atoms, possibly
ending up with a set of carefully isotopically labelled reagents.

While I believe that is a fundamentally better way to operate,
reaction files do have the ability to specify motifs that cause
the query match to not happen. Note that the same could be
achieved with the use of a query file in the reaction, as we will see.

We react an aniline with an acid, but not if there is an ortho
Fluorine - I have no idea if that is a thing, but it is a nice
simple example. First with no constraints on the aniline
```
scaffold {
  id: 0
  smarts: "[NH2]-c(:c):c"
}
sidechain {
  id: 1
  reagent: "OC=O"
  smarts: "OC=O"
  remove_atom: 0
  join {
    a1: 0
    a2: 1
    btype: SS_SINGLE_BOND
  }
}
```
But if an inactivating directive is added
```
scaffold {
  id: 0
  smarts: "[NH2]-c(:c):c"
  inactive  {
    query {
      smarts: "F-c:c-[NH2]"
    }
  }
}
```
the ortho-fluoro variant now reports no matches.

We could achieve the same thing by using a query file or with
an inline query. The query for the requirement looks like
```
query {
  smarts: "[NH2]-c(:c)(:c)"
  environment_no_match {
    attachment {
      attachment_point: [2, 3]
      btype: SS_SINGLE_BOND
    }
    smarts: "F"
  }
}
```
THis could be placed in a file, "aniline.qry" and incorporated as
```
  query_file: "aniline.qry"
```
Or it could also be incorporated directly into the reaction file as
```
scaffold {
  id: 0
  query {
    query {
      smarts: "[NH2]-c(:c)(:c)"
      environment_no_match {
        attachment {
          attachment_point: [2, 3]
          btype: SS_SINGLE_BOND
        }
        smarts: "F"
      }
    }
  }
}
```
The nested `query` directives are needed. While this nicely
demonstrates the flexibility of protos, I think in this case the
`query_file` directive would be better. That way the query is also
available to things like `tsubstructure`.

In summary, the `inactive` directive might be useful, but putting
such complexity into a query file, where it can be reused, should
generally be preferred.

## Keep Fragment
`fileconv` has extensive fragment handling means, but so does trxn.
You can specify matched atoms that specify fragments that are to 
be retained. For example if you wanted to keep all fragments with
more than 1 isotope 3 atom, that might be
```
scaffold {
  id: 0
  smarts: "([3C].[3C])"
  keep_fragment: 0
}
```

## Toggle Kekule Form
One of the significant limitations of LillyMol is the dependence
on Kekule forms. This can create problems when doing reactions.

It may be necessary to ensure that a reagent is in
a particular Kekule form. This can be accomplished by something like
```
toggle_kekule_form {
  a1: 0
  a2: 1
  btype: SS_SINGLE_BOND
}
```
If the current Kekule form is for there to be a single (and aromatic) bond between
matched atoms 0 and 1, this will be a no-op. Otherwise it will attempt to
rearrange the bonds in order to achieve the desired bonding. This usually
works fine, but can fail...

Whenever you encounter this problem, it is a low-point in the use of LillyMol.
There are however other advantages we gain by sticking to a requirement of
valid Kekule forms.


## No-Op Reaction
There have been circumstances where a reaction was needed, which did
nothing to its inputs, just reported a successful outcome. The
`noop_reaction` attribute specifies such a condition. It does not
work with `trxn`.

## Stereochemistry
Stereochemistry can be challenging in LillyMol reactions.

Currently chirality is not properly supported in smirks, so
```
trxn -K '[Cl:1][C:2][C:3].[I:4]>>[Cl:1][C@H:2]([C:3])[I:4]' ...
```
does not work. TODO:ianwatson implement this.

In 2023 rudimentary support for Cahn Ingold Prelog stereochemistry
specifications was added. As of writing, only the first shell of
atoms is examined. TODO:ianwatson implement the standard shell
expansion.

So, today, the easiest way to specify stereochemistry in a reaction
is via CIP style stereo designators. That might look something like
```
scaffold {
  id: 0
  smarts: "ClC(F)"
}
sidechain {
  id: 1
  reagent: "I"
  smarts: "I"
  join {
    a1: 1
    a2: 0
    btype: SS_SINGLE_BOND
  }
}
cip_stereo {
  atom {
    component: 0
    atom: 1
  }
  rs: CIP_S
}
```
which generates `I[C@H](Cl)F` as a product. Note that the `cip_stereo`
specification is an attribute of the reaction itself. This will usually
be the case, where a newly created centre must be given a specific
chirality.

If on the other hand, if the atoms involved can be localised to one of
the reagents, the CIP stereo specification can be a property of that location.

For example, we might have a set existing molecules, and need to enforce
a specific chirality, that might look like
```
scaffold {
  smarts: "C[CD3H](N)a"
  cip_stereo {
    atom: 1
    rs: CIP_R
  }
}
```
where the CIP stereo designator is fully contained within the scaffold.
Same for any 'sidechain' specificication. Again, until shell expansion is
implemented, this will remain challenging.

Stereochemistry likely spans multiple components, and therefore
is a property of the reaction, and not of any component. A simple, but
quite verbose, example is
```
scaffold {
  id: 0
  smarts: "C[CD3]([Cl])N"
}

reaction_stereo_center {
  center {
    atom {
      component: 0
      atom: 1
    }
  }
  top_front {
    atom {
      component: 0
      atom: 0
    }
  }
  top_back {
    atom {
      component: 0
      atom: 2
    }
  }
  left_down {
    atom {
      component: 0
      atom: 3
    }
  }
  right_down {
    implicit_hydrogen: true
  }
}
```
To specify a chiral center five atoms must be specified,
the center and either 3 or 4 attached atoms. These have names
indicative of their spatial positions. Each one can be either a
matched atom, in which case it must specify that component and
atom, or it can be an implicit_hydrogen. This is a proto oneof.

When give the smiles `CC([Cl])N` as input, it generates `C[C@H]([Cl])N`.

The `reaction_stereo_center` proto also has an `optional` directive
which means that the stereo center is created if possible, but ignored
otherwise.

Components have both `invert_chirality` and `remove_chirality` directives
which operate on their own matched atoms.

Newer versions of `trxn` support a `cip_stereo` directive which 
can use Cahn Ingold Prelog (CIP) notations to set a chiral centre.
```
cip_stereo {
  atom: 2
  rs: CIP_R
}
```
But note that support for Cahn Ingold Prelog type chirality is
currently only rudimentary in LillyMol. But this offers the
advantage of brevity. The current limitation is that shell
expansion is not implemented, so only the adjacent atoms
are considered when deciding on the chirality. 

## Stereo Preserving Atom Replacement
Unfortunately `trxn` performs its operations sequentially on a Molecule
object. The Molecule is smart, and when it detects that an atom that is
part of a chiral center has been removed, or unattached, it will
destroy that chiral center. Generally this is a good thing, but in
the case of a reaction, might not be what is desired. Frequently
the objective is to preserve the chirality of an atom, but to
make a change to one of its attachments.

To deal with that particular case the `replace_atom` directive is used.
Whereas join specifies two atoms between which a bond is to be
made, `replace_atom` is different. If the first atom in the replacement
is part of a chiral center, that chiral center is altered to include
the new atom. Consider the chiral molecule generated above, `C[C@H]([Cl])N`.
When reacted with
```
scaffold {
  id: 0
  smarts: "C[CD3]([Cl])N"
  remove_fragment: 2
}
sidechain {
  id: 1
  reagent: "Br"
  smarts: "*"
  replace_atom {
    c1 {
      component: 0
      atom: 2
    }
    a2: 0
  }
}
```
generates `C[C@H](Br)N`. Note there is an unfortunate quirk, `remove_atom: 2` does
not work in this case. Atom removals are done earlier in processing, whereas
fragment removals are done later. This would be difficult to alter.

## Inverting isotopes (uncommon).

A matched atom can have its isotopic information inverted.
```
name: "invert isotope"
scaffold {
  id: 0
  smarts: "a"
  invert_isotope {
    atom: 0
    isotope: 2
  }
}
```
Any matched atom with an existing isotope of 0, will be set to 2. Any
atom with an existing isotope will be set to 0.
