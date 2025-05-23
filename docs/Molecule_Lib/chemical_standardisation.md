# Chemical Standardisation

# Background.

When processing molecules it is convenient to ensure
that all molecules processed are in the same form. For example
all acids should be either all protonated, or all charged. Trying
to deal with a mixture will lead to complex, error prone code.

Similarly representational issues, like how to represent a Nitro
group, are also made much easier if all molecules adhere to a
single representational form.

Tautomers, where a Hydrogen may appear on different atoms, are
another problem where standardization can ensure that the Hydrogens
are localized on a specific location.

All three cases are handled by what is called chemical standardization
in LillyMol, and is implemented via the `-g` option in almost
all tools.

It is important to achknowedge that this is a controversial
area, and different people will come to different conclusions
about what kinds of representations to use. From a cheminformatics
perspective, it matters less what the representation is, than
the fact that all instances of a functinal group get treated the same
way. So whether one prefers charge separated nitro groups, or
neutral forms with a five valent Nitrogen atom, matters less than
just the idea of them all being the same - whatever that canonical
form is. Same with tautomeric representations. The exact tautomer
might be unknown, but if all instances are forced into the
same form, even if it is wrong, that should help find patterns
in a dataset.

In LillyMol we prefer neutral atoms wherever possible. One of
the goals of chemical standardization is to neutralize charged
atoms wherever possible. Therefore five valent Nitrogen atoms
in Nitro groups.

Note that both nitro representations are "wrong". The charge
separated form, `O=[N+]-[O-]` implies that the oxygen atoms
are different. Also any detection of formal charges will be
rendered more complex because of the need to identify, and ignore,
these actually neutral forms. A five valent Nitrogen is not
good either. The smiles language does not have a good way of
representing groups like this, although perhaps `N(:O):O` could
be used.

## HOWTO
Entering `-g help` to most LillyMol tools yields the usage message
```
 -g nitro       transform nitro groups to N(=O)=O
 -g n+o-        transform charge separated [N+]-[O-] (includes nitro)
 -g n+n-        transform charge separated [N+]-[N-] to N=N
 -g s+c-        transform [S+]-[C-] to S=C
 -g all+-       transform all [X+]-[Y-] to X=Y
 -g xh          remove hydrogens
 -g amine       change all amines
 -g o-          protonate all O- groups
 -g n-          protonate all N- groups
 -g nrmch       remove all hydrogens except those to chiral centres
 -g covm        break covalent bonds between Oxygen and Na,K
 -g isolc       assign formal charges to isolated Na, K, .. and Halogens
 -g guan        convert guanidines to -N-C(=N)-N form
 -g Rguan        convert Ring type guanidines to -N-C(=N)-N form
 -g azid        convert charge separated azids [N-]=[N+]=N to N#N=N
 -g msdur       convert misdrawn ureas, O-C(=N)-N to O=C(-N)-N
 -g msdsa       convert misdrawn sulfonamides, O=S(O)=N to O=S(=O)N
 -g fcor        for converting back from corina mangled structures
 -g ehlast      move all explicit Hydrogen atoms to last in the connection table
 -g fmrk        reverse transformations applied to .mrk files
 -g fwih        fix obviously wrong implicit hydrogen settings
 -g imidazole   convert imidazoles to have nH near cD3
 -g charged_imidazole   convert charged imidazoles to have n+ near cD3
 -g pyrazole    convert pyrazoles to have nH near electron withdrawing
 -g triazole    convert triazoles to have nH near electron withdrawing
 -g tetrazole   convert tetrazoles to have nH near attachment
 -g ltlt        convert lactim to lactam form (non ring)
 -g ltltr       convert lactim to lactam form (ring)
 -g isoxazole   convert Hydroxy isoxazoles to O= forms
 -g arguan      aromatic "gauanidines" - adjacent to =O, better name needed
 -g pyrazolone    convert pyrazolone to keto form
 -g aminothazole  convert -N=c1scc[nH]1 to -[NH]c1sccn1
 -g keto_enol     convert enol to keto forms (no adjacent heteroatoms)
 -g 4-pyridone    convert 4 pyridol to pyridone form
 -g surea         convert S-C(=N)-N to S=C(-N)-N
 -g 124-triazine  convert S-C(=N)-N to S=C(-N)-N
 -g enol-fused    convert [O,S;D1]-c(:n):[aD3x3] to O=C form
 -g isotope       convert all isotopic atoms to non isotopic forms
 -g to2ap         convert exocyclic N=C:[nH] to N-c:n (InChI generates these)
 -g all         ALL the above standardistions
 -g rvnitro     convert O=N=O nitro groups to charge separated
 -g rvnv5       convert all 5 valent N atoms to charge separated
 -g APP=<xxx>   append 'xxx' to the name of changed molecules
 -g APP=EACH    append the reason for each change
 -g EXT:/path/to/textproto   file of external standardisation(s) standardisation::Standardisation proto
                                 one proto per line, eg: smarts: "[ND1H1]=[CD3]-[OD1H]" smiles: "N-C=O"
```

## Discussion
Many of these are straightforward, but some have proven to be
extremely difficult. The lactam/lactim transformation can exist in
aromatic systems where multiple instances of such groups can
coexist. At this stage, the tool converts the simple cases, but
gives up in the face of more complex, possibly linked groups.

### Activation
The usual way of using chemical standardisation is via `-g all` which turns
on all default transformations. It is also possible to activate individual
transformations only. For example, to activate only the transformation of
Nitro groups to charge neutral forms
```
-g nitro
```
will do that. To activate all transformations, except the nitro transformation use
negation
```
-g all -g -nitro
```
after all transformations have been activated, specific transformations can be
disabled by prepending a '-' to the directives.

### External
While it has been desirable to handle many cases in c++, there are often ad-hoc
cases where on a project specific basis, certain groups need to be transformed, in
on way or another.

The Standardisation [proto](/src/MoleculeLib/standise.proto) allows specification of
a smarts to identify the atoms to be changed, and a smiles that specifies how the
bonds should be changed. For example, a transformation to charged acids might look
like
```
smarts: "[OH][C,S]=O" smiles: "[O-]-*=O" name: "acid"
```
in this case a formal charge is placed on the Oxygen atom.

Note that the smiles must be a valid smiles, but even if the smarts
specifies atoms that must be in a ring, the smiles does not necessarily
need that. For example
```
smarts: "[NR0]=c1:[cD2]:[cD2]:[nH]:c:c1" smiles: "N-C=CC=N" name: "para-amino"
```
specifies a full ring in the query (7 query atoms), but the smiles only matches the first 
5 matched atoms, with no ring information.

Multiple protos can be in a file and the transformations added to the standardisation via
```
-g EXT:/path/to/file.textproto
```
