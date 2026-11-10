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
 -g keto_enol     convert O/S enol to keto forms (no adjacent heteroatoms)
 -g 4-pyridone    convert 4-pyridol/thiol to pyridone/thione forms
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

### Lactim/lactam and sulfur analogues

`-g ltlt` converts eligible non-ring `N=C-OH` and `N=C-SH` groups to their
amide and thioamide forms. The terminal oxygen or sulfur must have one
connection and the central carbon must have three connections. Nitrogen can
be two-connected, or a neutral terminal `=NH`. Explicit `[OH]`, `[SH]`, and
`[NH]` hydrogen specifications are supported.

| Input | Standardised form |
| --- | --- |
| `CC(O)=NC` | `CC(=O)NC` |
| `CC(S)=NC` | `CC(=S)NC` |
| `CC(O)=N` | `CC(=O)N` |
| `CC(S)=N` | `CC(=S)N` |
| `NC(S)=N` | `NC(=S)N` |

Terminal `=NH` support broadens equivalence relative to older releases, which
left these forms unchanged. Charged terminal imines and specified C=N
stereochemistry remain excluded from this transformation. Two-connected
central carbons, such as `OC=N`, also remain excluded. Other transformations
activated by `-g all` can independently affect charges or related groups.

`-g ltltr` handles eligible ring analogues, including oxygen and sulfur forms.
Both transformations are included in `-g all`. The separate `-g keto_enol`
option handles carbon-based `C=C-OH`/`C=C-SH` groups and is not the option for
these nitrogen-adjacent groups.

### Other oxygen/sulfur tautomer rules

The following transformations now also handle the sulfur analogues. This
broadens standardised equivalence relative to earlier versions.

| Option | Example sulfur transformation |
| --- | --- |
| `-g keto_enol` | `CC(S)=CC` → `CC(=S)CC` |
| `-g 4-pyridone` | `Sc1ccncc1` → `S=C1C=CNC=C1` |
| `-g oxopyrimidine` | `S=c1cc[nH]c(C)n1` → `S=c1ccnc(C)[nH]1` |

General keto/enol processing retains its conservative restrictions: the central
carbon is non-ring, its other neighbors are carbons, and interacting groups
or nearby unsaturation can prevent conversion. Existing C=S groups participate
in the same interference checks as C=O groups, including mixed oxygen/sulfur
molecules. Ethers and thioethers are not OH/SH tautomers and are excluded.

The default keto preference can be changed for handlers that support it:

```shell
fileconv -g keto_enol -g ENOL -S - molecules.smi
```

For example, `CC(=S)CC` becomes `CC(S)=CC`. `-g KETO` restores the default
preference. The named `4-pyridone` and `oxopyrimidine` transformations retain
their own fixed canonical choices. `4-pyridone` remains opt-in and is not
enabled by `-g all`; use `-g all -g 4-pyridone` to combine it with the default
transformations.

Oxopyrimidine processing rejects competing C=O **or C=S** groups attached to
the hydrogen-bearing nitrogen; converting one group must not silently favor
an incompatible tautomer of another. The existing restrictions on complex
ring environments remain in place.

Hydrogen counts specified in brackets, such as `[OH]`, `[SH]`, and `[CH2]`, are
refreshed on atoms whose hydrogen count changes during these conversions. This
also fixes earlier oxygen cases that could retain an obsolete hydrogen count
and produce an invalid-valence product.

### Isolated amidines and pyrimidinetriones

`-g amidine` chooses the terminal-imine form of eligible isolated amidines:

```text
CC(N)=NC  ->  CC(=N)NC
C(N)=NC   ->  C(=N)NC
```

The initial rule requires a terminal nitrogen and an unsubstituted or
carbon-substituted amidine carbon. When the imine nitrogen has a substituent,
that substituent must be a saturated, non-aromatic carbon. Ring nitrogens,
shared/conjugated amidine systems, N-aryl and N-acyl groups, N-N/N-O
substituents, and specified C=N stereochemistry are excluded. With `-g all`,
normal amine neutralisation runs first, so an eligible protonated amidine
also reaches its neutral canonical form in one pass. `-g amidine` alone does
not neutralise charged groups.

`-g pyrimidine` now handles isolated pyrimidine-2,4,6-trione systems and their
oxygen/sulfur analogues. This option previously had no implemented conversion.
It considers the whole ring together, avoiding competing local proton shifts:

```text
N1=C(O)C=C(O)N=C1O  ->  N1C(=O)CC(=O)NC1=O
```

Partially keto starting forms, N-alkyl substitution, substitution at position
5, and mixed oxygen/sulfur systems are supported. The conversion conserves
hydrogen count. Fused/spiro systems, charged participating atoms, explicit
hydrogen attachments, and a specified chiral centre at position 5 are excluded.
The rule is skipped with the `ENOL` preference; it does not implement a general
enol-preferred pyrimidine canonicalisation.

Both options are included in `-g all` and broaden the equivalence of structures
relative to earlier versions. They can be disabled with
`-g all -g -amidine -g -pyrimidine`. These are deliberately restricted rules;
general heteroatom-rich and ring-conjugated tautomer systems remain outside
their scope.

For enol-preferred processing, independent eligible keto groups are now all
converted in one pass. An indexing error previously caused some groups to be
left for a subsequent pass. Interacting groups remain excluded.

Use `tstandardise` to check random atom orderings and repeated processing:

```shell
tstandardise -g all -u -n 100 -v -F failures.smi molecules.smi
```

`-u` removes marked centres that LillyMol does not consider actually chiral.
This avoids failures caused by such annotations on symmetric rings. Omit it
when those annotations are part of what you want to test.

Inspect the reported failure count and `failures.smi`; the tester can return
success even when it reports inconsistencies. A failure can also originate in
SMILES/aromaticity/stereochemistry handling rather than a standardisation rule.
Use `-f` to stop with a nonzero exit status on a failure. Invalid starting
valence is distinguished from invalid valence introduced by standardisation;
the latter is counted as a failure and written to the failure file.

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

The Standardisation [proto](../../src/Molecule_Lib/standardise.proto) allows specification of
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
6 matched atoms, with no ring information.

Multiple protos can be in a file and the transformations added to the standardisation via
```
-g EXT:/path/to/file.textproto
```

An external directive can be used on its own or together with built-in options.
Each proto must occupy one line in the external file.

Instead of inline SMARTS, use a query file:

```text
query_file: "amide.textproto" smiles: "N-C=O" name: "amide"
```

Relative query paths are resolved against the external standardisation file's
location. Environment variables such as `${LILLYMOL_HOME}` are expanded first.
Both textproto queries and legacy MSI queries are accepted; their contents
identify the format, independently of filename suffix.

External transformations require one query component with usable embeddings.
Template atom `i` corresponds to initial query atom number `i`; a template may
cover a prefix of the query, but every template atom must exist in the query.
Template bonds must already exist in the target: this facility changes bond
orders and formal charges, and does not add or remove atoms or bonds. Invalid
mappings are diagnosed and skipped before any part of that embedding is changed.
The query must constrain the target sufficiently for the imposed bonds and
charges to give chemically valid products. Test external rules with
`tstandardise` before applying them to a large collection.
