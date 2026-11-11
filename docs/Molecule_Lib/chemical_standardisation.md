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
 -g nitroso      convert neutral aliphatic C-N=O groups to C=N-OH
 -g phosphorus   convert neutral P-OH tautomers to P(=O)H
 -g beta-ketoester canonicalise acyclic beta-ketoester enols
 -g beta-diketone canonicalise beta-diketone enols to diketo form
 -g ascorbate    canonicalise isolated ascorbate rings to lactone/enediol form
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

### Aliphatic nitroso groups

`-g nitroso` converts eligible aliphatic C-nitroso groups to oximes. For
example, `CCN=O` becomes `CC=NO`. It is enabled by `-g all`, broadening
default tautomer equivalence, and can be disabled with `-g all -g -nitroso`.
The oxime is the fixed choice, independently of `KETO` or `ENOL`.

Nitrogen and oxygen remain neutral throughout the conversion. The rule does
not introduce charge-separated representations. Persistent charges elsewhere,
such as a pyridinium group, are retained by this transformation.

The carbon must be neutral, non-aromatic and fully saturated, with at least
one implicit hydrogen. Its other neighbours must be carbon or hydrogen. The
nitroso nitrogen must be neutral and two-connected, with one single bond to
carbon and one double bond to a neutral terminal oxygen. Aromatic nitroso
compounds, nitro groups, N-nitroso compounds, and additional heteroatom
substituents on the reacting carbon are excluded from this initial rule.

The hydrogen moves implicitly from carbon to oxygen. Isotope labels and
explicit hydrogen atoms are retained. If no implicit carbon hydrogen is
available, the conversion is skipped. The ordinary explicit-hydrogen removal
in `-g all` runs beforehand, while isotopic hydrogen atoms remain explicit.

A conversion touching a marked chiral centre is skipped. Unaffected
stereocentres are retained, and the newly formed C=N geometry is unspecified.
Existing oxime geometry is preserved. Processing is gated by the existing
nitrogen and oxygen counts. Reports include nitroso-to-oxime counts, and
`-g APP=EACH` appends `STD:nitroso` when this rule changes a molecule.

### Phosphorus OH tautomers

`-g phosphorus` converts neutral, three-coordinate P(III) groups bearing an OH
substituent to the P(=O)H form. For example, `P(O)(O)C(N)C` becomes
`P(=O)(O)C(N)C`. LillyMol supplies the implicit hydrogen on phosphorus after
the conversion, preserving molecular formula.

The rule is enabled by `-g all` and can be disabled with
`-g all -g -phosphorus`. The existing phosphorus count guards both the processing
check and the rule invocation. Molecules without phosphorus incur no additional
atom scan from this rule. The P=O choice is independent of `KETO`/`ENOL`.

This first implementation supports carbon and oxygen substituents on neutral
phosphorus. It skips marked phosphorus stereochemistry and charged neighbours.
A hydrogen attached explicitly to oxygen is not moved. Other implicit OH groups
on the same phosphorus remain eligible, including when an isotopic hydrogen is
attached elsewhere. `-g all` removes ordinary explicit hydrogens before applying
this rule. Oxygen isotope labels are retained; the lowest oxygen isotope is
chosen if multiple OH groups are eligible. Existing P=O groups are left as read.

Reports include phosphorus tautomer counts. `-g APP=EACH` appends
`STD:phosphorus` when this transformation changes a molecule.

### Acyclic beta-ketoester enols

`-g beta-ketoester` converts eligible `RO-C(=O)-CH=C(OH)-R` motifs to
`RO-C(=O)-CH2-C(=O)-R`. For example, `COC(=O)C=C(O)C` becomes
`COC(=O)CC(C)=O`. The ester bonds are retained. This broadens default
equivalence because the rule is included in `-g all`. Disable it with
`-g all -g -beta-ketoester`; the beta-diketone rule remains independently
selectable. The ester rule is skipped with the `ENOL` preference.

The ester carbonyl, central carbon and enol carbon must be acyclic. The central
carbon must have one hydrogen and no additional substituent in the starting
enol. The ester oxygen must be neutral and bonded to another carbon. Alkyl,
aryl and vinyl esters are supported. Carboxylic acids, anhydrides, amides,
thioesters, lactones and cyclic beta-ketoesters are excluded from this initial
rule. Enolisation at the ester end is also outside its scope.

Formula, isotope labels and unaffected stereocentres are preserved. Conversions
that touch a marked chiral centre or remove specified alkene geometry are
skipped. The rule precedes beta-diketone processing so a newly formed ketone
can enable another local enol conversion in the same standardisation pass.
Reports include beta-ketoester counts, and `-g APP=EACH` appends
`STD:beta-ketoester` when this rule changes a molecule.

### Beta-diketone enols

`-g beta-diketone` converts eligible `O=C–CH=C(OH)–` motifs to the
diketo form. The rule applies to acyclic motifs and isolated non-aromatic
carbocycles. For example, `CC(=O)C=C(O)C` becomes `CC(=O)CC(C)=O`.

In an isolated carbocycle, the central carbon may also have a carbon
substituent, including an aryl group. For example,
`OC1=C(c2ccccc2)C(=O)CCC1` becomes `O=C1C(c2ccccc2)C(=O)CCC1`.
Both eligible enol orientations converge to the same diketo form. This
extends the equivalence of structures under default standardisation.

The rule requires a ketone at one end and an unprotected enol at the other.
It skips substituted acyclic central carbons, specified alkene geometry,
and fused rings. A conversion that would modify an atom with a marked
chiral centre is skipped. Existing centres at unaffected atoms are retained.
Newly formed stereocentres remain unspecified, with no `Chiral_Centre`
object created. It is enabled by `-g all`, skipped with the `ENOL`
preference, and can be disabled with `-g all -g -beta-diketone`.

### Ascorbate rings

`-g ascorbate` chooses the conventional lactone/enediol form for a restricted
family of isolated ascorbic-acid-like rings. It is included in `-g all` and
therefore broadens default equivalence relative to earlier versions. Disable it
with `-g all -g -ascorbate`.

```text
OCC(O)C1OC(O)=C(O)C1=O  ->  OCC(O)C1OC(=O)C(O)=C1O
OCC(O)c1oc(O)c(O)c1O    ->  OCC(O)C1OC(=O)C(O)=C1O
```

The aromatic trihydroxyfuran tautomer is converted to the same lactone. Eligible
diketo forms are also handled. The central enediol oxygen may carry a phosphate
ester, which is retained:

```text
OCC(O)C1OC(O)=C(OP(=O)(O)O)C1=O
  -> OCC(O)C1OC(=O)C(OP(=O)(O)O)=C1O
```

The ring must contain one oxygen and four carbons, with three consecutive
oxygen-substituted carbons. The remaining ring carbon must have one carbon
substituent and one hydrogen in the target. The carbonyl and terminal enediol
oxygens must be unprotected; the central oxygen can be OH or O-P. Other ether
and ester substitutions are excluded. Ring atoms and participating exocyclic
oxygens must be neutral; phosphate charges on other oxygens are preserved by
this rule. Other options in `-g all` can neutralise those charges independently.

The transformation preserves hydrogen count and isotope labels. Reduced or
oxidised forms with different hydrogen counts, fused/spiro rings, and sulfur
analogues are excluded. Unchanged stereocentres are preserved; a conversion
that changes a bond at a marked chiral centre is skipped. An aromatic starting
form does not specify the configuration of the newly saturated ring carbon,
so no configuration is assigned. Explicit hydrogen atoms on oxygens whose
protons would move are excluded. Ordinary explicit hydrogens are removed
beforehand by `-g all`; isotopic hydrogens are not moved by this rule.

This is a fixed canonical choice, independent of `-g KETO` or `-g ENOL`.
It does not enable carbonyl migration in open-chain sugars or cyclic polyols.
Verbose reports include ascorbate counts, and `-g APP=EACH` appends
`STD:ascorbate` when this rule changes a molecule.

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
