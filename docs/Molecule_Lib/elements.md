# Elements

LillyMol contains a periodic table data structure, and each
element is an Element object. Each Atom has a pointer to
the Element that defines that Atom. Every atom must have an
Element.

LillyMol is unusual in that the universe of elements is not
restricted to the existing periodic table. Any two letter
combination can be an element, being added to the periodic
table as needed. Once that functionality was in place, we
realised that there was nothing special about two character
elements, so the concept was extended to any string.

HOWTO:

Most LillyMol tools have a `-E` option which enables customisation
of the elements. The usage message is
```
  -E autocreate  automatically create new elements when encountered
  -E PTABLE=file use 'file' for element data
  -E O:Xx        element Xx is considered Organic
  -E anylength   elements can be of arbitrary length
  -E nsqb=xx     element <xx> needs square brackets
  -E <symbol>    automatically create \<symbol\>
```

### -E autocreate
Automatically create new elements when encountered. This applies to
two letter combinations only.
All combinations of the form `[A-Z][a-z]` are recognised as elements.
This also enables such things as `R1` as a valid element. But non
alphanumeric characters are not permitted.

### -e PTABLE=file
Use 'file' for element data. This is retained for historical reasons
only. No current use cases are anticipated.

### -E O:Xx
Element Xx is considered Organic. A common use might be `-E O:Se`.

### -E anylength
Elements can be of arbitrary length. So things like `Ala` and `Phe`
become valid elements. The restriction to alphanumeric characters
is also lifted. Using this implies `-E autocreate` so you do not
need to specify both options.

In bracketed SMILES, an arbitrary element name must start with a letter.
Leading digits specify an isotope: `[1Foo]` is isotope 1 of the pseudo-element
`Foo`, not an element called `1Foo`. The complete symbol is recognized, so
`[Cz]` remains the pseudo-element `Cz` rather than carbon followed by `z`.

The characters `@`, `+`, `-`, and `:` delimit SMILES atom attributes and cannot
be part of an element name in this representation. A trailing `H` or `H` followed
by a hydrogen-count digit is also an attribute. For example:

| SMILES atom | Element | Isotope | H count | Charge | Atom map |
| --- | --- | --- | --- | --- | --- |
| `[Foo]` | Foo | 0 | 0 | 0 | 0 |
| `[FooH]` | Foo | 0 | 1 | 0 | 0 |
| `[1FooH2+:17]` | Foo | 1 | 2 | +1 | 17 |
| `[Cz:20]` | Cz | 0 | 0 | 0 | 20 |

Consequently, a literal pseudo-element name ending in an apparent hydrogen-count
suffix cannot be represented unambiguously and is interpreted as the shorter
name plus that count. This does not depend on which names have already been
created. A standalone `[H]` remains the hydrogen element. Periodic-table atoms
retain their ordinary hydrogen counts, charge, chirality, isotope, and atom-map
attributes when arbitrary names are enabled. These annotations do not assign
valences to pseudo-elements.

If you need to perform substructure searches against these elements,
you will need to use the #{symbol} smarts extension. See the file
`Molecule_Tools/introduction.cc` for examples of how to do substructure
searches over peptides.

### -E nsqb=xx
Element <xx> needs square brackets. Use case is unclear, almost all
elements need square brackets. Mostly this gets used with Hydrogen.
We have encountered explicit Hydrogens written as either `H` or `[H]`.

### -E <symbol>
Automatically create \<symbol\>. Generally it is easier to just use
`autocreate` and allow LillyMol to create whatever otherwise unrecognised
elements it encounters. But for more control, the set of acceptable
extra elements can be specified via the `-E` option.

## Element_Matcher
Several tools use the Element_Matcher class to specify atoms to
be matched, based on their element. For example, the '-X' option
to fileconv removes atoms with a specified element. The Elements_to_Remove
object uses an Element_Matcher in order to identify the atoms to
remove. In order to see the functionality try
```
fileconv -X help
```
Other tools containing this functionality should respond similarly.
