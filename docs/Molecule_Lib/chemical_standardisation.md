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
 -g biguanide   canonicalise eligible adjacent guanidines together
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
 -g amidoxime    convert neutral acyclic N-hydroxyamidines to amino-oximes
 -g cyclic-amidine canonicalise isolated neutral cyclic amidines
 -g nitroso      convert neutral aliphatic C-N=O groups to C=N-OH
 -g phosphorus   convert neutral P-OH tautomers to P(=O)H
 -g beta-enaminone convert eligible acyclic beta-ketoimines to enaminones
 -g beta-hydrazone prefer acyclic beta-keto hydrazones (included in all)
 -g beta-enehydrazine prefer acyclic beta-enehydrazines (alternative preference)
 -g alpha-ketoamide canonicalise simple acyclic alpha-ketoamide enols
 -g alpha-ketoester canonicalise simple acyclic alpha-ketoester enols
 -g alpha-ketoacid canonicalise simple acyclic alpha-ketoacid enols
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

### Coordinated biguanide standardisation

`-g biguanide`, included in `-g all`, handles eligible pairs of guanidine
centres joined through a two-connected nitrogen. It chooses a neutral form
with one C=N bond to that bridging nitrogen and one terminal imine:

```text
CN(C)C(=N)NC(=N)N -> CN(C)C(N)=NC(=N)N
```

This broadens default tautomer equivalence. Disable it with
`-g all -g -biguanide` to retain the individual guanidine rules. The preference
is independent of `KETO`/`ENOL`. Reports count changed biguanide units, and
`-g APP=EACH` appends `STD:biguanide`. Chemical-standardisation verbose output (for example,
`fileconv -v -v`) distinguishes units handled from units
changed, including canonical units protected from individual processing.

The coordinated rule runs before the individual guanidine rules. Both carbon
centres are removed from their shared candidate list whenever a unit is
handled, including when it already has the preferred bonds. This prevents the
individual rules from undoing the coordinated choice on either half.

The initial scope requires two non-ring guanidine carbons, a non-ring
2-connected bridging nitrogen, and ordinary amino/imino valences. External
nitrogens may carry saturated carbon or aromatic carbon substituents, including
membership in an external saturated ring. Neutral and protonated N precursors
are accepted; positive core charges are neutralised together with bond changes
and hydrogen reassignment. Charges elsewhere are unaffected by this rule.
Explicit H/D neighbours, participating chiral annotations, specified geometry,
heteroatom substitutions on the external nitrogens, N-acyl/conjugated substituents, four-connected
nitrogens, and larger overlapping guanidine chains are excluded.

At least one centre must have an unsubstituted terminal nitrogen for the
terminal C=N bond. When both bridge orientations are possible, the bridge
double bond is placed on the side with more external carbon substituent bonds.
Equal, nonzero substitution counts are outside this initial rule: resolving
those cases requires a more general comparison of their substituents. Parent
biguanide is supported; isotope labels distinguish its otherwise equivalent
ends and terminal nitrogens where necessary. Existing double bonds between
chemically equivalent atoms are retained where possible to avoid spurious
changes. Unsupported units remain subject to the existing individual rules.

The choice of a bridge-imine form is supported by the experimental and
computational evidence discussed in
[Fraczkiewicz et al., 2024](https://doi.org/10.1002/minf.202400088).
The orientation tie-breaks are canonicalisation conventions, not predictions
of solution populations for every substituted biguanide. No E/Z geometry or
new chiral-centre annotations are assigned. Heavy-atom isotope labels and
remote stereochemistry are retained.

### Completing lactim conversion after fused-ring enol changes

`-g enol-fused` can expose a ring lactim that was excluded during initial
perception. When two OH/SH lactims share one nitrogen, both are initially
excluded to avoid competing conversions. Converting one with `enol-fused`
can make the other eligible for `-g ltltr`.

For example, CHEMBL1899721 previously needed two standardisation passes:

```text
OC1=NC(O)=CC2=NN(C(N)=C12)c1ccccc1
  -> O=C1NC(O)=CC2=NN(C(N)=C12)c1ccccc1
  -> O=C1NC(=O)C=C2NN(C(N)=C12)c1ccccc1
```

With both rules enabled (as in `-g all`), an actual `enol-fused` change now
triggers fresh lactim perception and one completion of the existing ring
lactim rule after the normal transformations. This preserves the initial
snapshot for other rules and avoids repeating the entire standardisation
pipeline. The additional perception is limited to molecules changed by
`enol-fused`; it does not copy the molecule. Disabling either rule retains
that rule's usual scope. Completion is reported as `STD:LTLTR` with
`-g APP=EACH` and counted by the ring lactim rule.

The fused-ring enol conversion also clears bracket-specified hydrogen counts
on atoms whose bonds change. This prevents a labelled `[18OH]` becoming an
invalid `[18OH]=C` instead of `[18O]=C`. Heavy-atom isotope labels are retained.

### Guanidine charge and valence handling

`-g guan`, included in `-g all`, relocates eligible guanidine C=N bonds.
When a positively charged iminium nitrogen loses its double bond, charge and
implicit-hydrogen information are now updated consistently in both the
terminal-nitrogen shortcut and general conversion path. For example:

```text
[N+]1(=C(N)N)CCCCC1 -> N=C(N)N1CCCCC1
```

Previously, the shortcut could leave a three-single-bond `[N+]` with no
hydrogen. This introduced invalid valence, and the later amine neutralisation
could not repair it because that rule requires a hydrogen on the charged
nitrogen. Resetting the hydrogen counts also handles bracket-specified NH
counts correctly. Explicit isotope-labelled hydrogen atoms are retained.

Charge counters are updated for the actual charged nitrogens, even when the
general resolver has permuted its nitrogen variables. If that resolver selects
a three-connected nitrogen to carry the neutral C=N bond, the unresolved group
is left unchanged rather than introducing invalid valence. This guard does
not provide a new canonicalisation strategy for those substituted groups.

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

### Amidoximes

`-g amidoxime` chooses the neutral amino–oxime representative of eligible
N-hydroxyamidines:

```text
CC(=N)NO       -> CC(N)=NO
ONC(=N)C(N)C   -> ON=C(N)C(N)C
```

A hydrogen moves from the hydroxylamine nitrogen to the terminal imine
nitrogen. The OH hydrogen stays on oxygen. This option is included in `-g all`
and broadens default tautomer equivalence. Disable it with
`-g all -g -amidoxime`. The preference is independent of `KETO`/`ENOL` because
it selects a nitrogen tautomer without changing the OH bond.

The initial scope requires a neutral, acyclic `R-C(=NH)-NH-OH` group, with a
carbon atom for R, a terminal imine NH and a terminal OH. Aromatic carbon
substituents are allowed, but the participating atoms cannot be aromatic or
in rings. Formamidoxime (R = H), guanidine analogues (R = N), N-substituted
nitrogens, O-alkyl/O-acyl groups, and charged participating atoms are excluded
by the isolated rule. With `-g all`, ordinary charge neutralisation and explicit
hydrogen removal run first, so eligible charged precursors can also reach the
neutral representative in one pass.

Explicit hydrogen atoms attached to either nitrogen or oxygen prevent the
isolated conversion, protecting isotopic hydrogen. Heavy-atom isotope labels
and remote stereochemistry are retained. Marked chiral centres on participating
atoms and specified geometry on the group prevent conversion. No C=N geometry
or chiral-centre annotations are assigned to the product. Multiple eligible
groups are processed in the same pass.

The preference agrees with IR/NMR evidence for amino–oxime structures in
solution reported by [Bell, Nambury and Bauer, 1964](https://doi.org/10.1021/jo01033a016).
Calculations by [Tavakol and Arshadi, 2009](https://doi.org/10.1007/s00894-008-0435-4)
also favour amidoxime over imino–hydroxylamine. This is a supported
canonicalisation preference, not a prediction of populations for arbitrary
substituents, pH, or solvent. Charge-separated aminonitrone forms are outside
this rule.

Reports count amidoxime conversions separately. `-g APP=EACH` appends
`STD:amidoxime` when this rule changes a molecule.

### Locally resolved cyclic guanidines

`-g Rguan`, included in `-g all`, also handles neutral guanidines in isolated,
otherwise saturated five- and six-membered rings containing two nitrogens:

```text
N1C(=Nc2ccccc2)NCC1 -> N1=C(Nc2ccccc2)NCC1
```

The exocyclic nitrogen must have one alkyl or aryl carbon substituent.
N-acyl and other unsaturated aliphatic external branches are excluded.
Specified geometry, marked stereocentres on the ring or reacting nitrogen,
and explicit hydrogen neighbours on reacting nitrogens prevent this additional
conversion. Existing handling of other ring guanidines is retained.

The representative has an endocyclic C=N bond. The ring nitrogen is chosen
by a local lexicographic comparison: lower nitrogen isotope first, then lower
connectivity of the adjacent ring carbon, then lower isotope on that carbon.
The comparison does not depend on the incoming double-bond position or the
nitrogen hydrogen counts. Both incoming endocyclic forms are examined, so
an existing C=N may move to the other ring nitrogen.

If these properties tie and both neighbouring ring carbons are matching CH2
groups, the two choices are equivalent; an existing endocyclic bond is retained.
More substituted ties remain unchanged. There is no shell expansion or full
canonical ranking. This is a canonicalisation convention rather than a
prediction of solution tautomer populations. Isotopes and remote stereochemistry
are retained, no charges are introduced, and changed nitrogen hydrogen
specifications are reset.

Disable the ring rule with `-g all -g -Rguan`. Conversions use the existing
ring-guanidine counter and `STD:Rguanidine` annotation with `-g APP=EACH`.

### Simple alpha-ketoamides

`-g alpha-ketoamide` selects the keto form of eligible neutral, acyclic
alpha-ketoamide enols:

```text
NC(=O)C(O)=C             -> NC(=O)C(=O)C
CN(C)C(=O)C(O)=Cc1ccccc1 -> CN(C)C(=O)C(=O)Cc1ccccc1
```

Included in `-g all`, this broadens default tautomer equivalence. Disable it
with `-g all -g -alpha-ketoamide`. The rule is skipped with `ENOL`.
Primary, secondary and tertiary amides are accepted, including amides whose
nitrogen belongs to a saturated ring. The amide nitrogen must be neutral,
with only single bonds and hydrogen/carbon substituents. Hydroxamates,
hydrazides, imides and related additional N-acyl/N-thioacyl groups are excluded.
The amide carbonyl is retained unchanged.

The reacting chain must be acyclic, with a terminal enol OH and a CH2 or CH
carbon at the other end of the double bond. That carbon can have a saturated
or aromatic carbon substituent. Additional aliphatic conjugation, adjacent
carbonyl branches, cyclic reacting chains and more substituted double bonds
are excluded. Specified geometry, marked chiral centres on reacting atoms,
and explicit hydrogen neighbours on the enol OH or reacting carbon prevent
conversion. Heavy-atom isotope labels and remote stereochemistry are retained;
no charges or new chiral-centre annotations are introduced. The chosen form
is a canonicalisation preference, not a prediction of solution populations.

Reports count this rule separately; `-g APP=EACH` appends
`STD:alpha-ketoamide` to changed molecules.

### Simple alpha-ketoesters

`-g alpha-ketoester` selects the keto form of eligible neutral, acyclic
alpha-ketoester enols:

```text
COC(=O)C(O)=C          -> COC(=O)C(=O)C
CCOC(=O)C(O)=Cc1ccccc1 -> CCOC(=O)C(=O)Cc1ccccc1
```

This rule is included in `-g all`, broadening default tautomer equivalence.
Disable it with `-g all -g -alpha-ketoester`. It is skipped with `ENOL`.
The ester must be an ordinary alkyl or aryl ester. The reacting chain must
be acyclic, with a terminal enol OH and a CH2 or CH carbon at the other end
of the double bond; that carbon can have a saturated or aromatic carbon
substituent. Anhydrides, enol esters, cyclic reacting chains, additional
aliphatic conjugation and adjacent carbonyl branches are excluded.

Specified geometry, marked chiral centres on reacting atoms, and explicit
hydrogen neighbours on the enol OH or reacting carbon prevent conversion.
Heavy-atom isotopes and remote stereochemistry are preserved. The ester
carbonyl is unchanged, no charges are introduced, and no new chiral-centre
annotations are created. This is a canonicalisation preference, not a
prediction of solution tautomer populations.

Conversions are reported separately; `-g APP=EACH` appends
`STD:alpha-ketoester` to changed molecules.

### Simple alpha-ketoacids

`-g alpha-ketoacid` chooses the keto representative of eligible neutral,
acyclic alpha-ketoacid enols:

```text
O=C(O)C(O)=Cc1ccccc1 -> O=C(O)C(=O)Cc1ccccc1
O=C(O)C(O)=C         -> O=C(O)C(=O)C
```

This option is included in `-g all` and broadens default tautomer equivalence.
Disable it with `-g all -g -alpha-ketoacid`. The rule is skipped with `ENOL`;
other enabled keto/enol rules retain their existing behaviour.

The general `keto_enol` rule excludes unsaturated neighbouring groups. Those
safeguards remain intact; the specialised rule recognises a carboxylic acid
adjacent to the enol explicitly. The initial scope requires a neutral COOH,
a terminal enol OH, and a non-ring C=C whose other carbon is CH2 or CH attached
to a saturated or aromatic carbon. This includes pyruvic-acid and
phenylpyruvic-acid forms. Substituted enol carbons, cyclic reacting chains,
additional aliphatic conjugation and most adjacent carbonyl branches are
excluded, as are replacement of the COOH by an ester or amide, sulfur
replacements, and hydrates.

One additional carbonyl branch is supported explicitly:

```text
O=C(O)C(O)=CC(=O)NR2 -> O=C(O)C(=O)CC(=O)NR2
```

Here `NR2` denotes a neutral amide nitrogen with hydrogen and/or carbon
substituents, including cyclic amines. Hydroxamates, hydrazides and imides
are excluded. This broadens the equivalence classes formed by `-g all`;
the amide carbonyl remains unchanged. Existing beta-diketone
handling remains responsible for its own more complex cases.

The isolated rule accepts only neutral participating atoms. With `-g all`,
carboxylate/enolate neutralisation can run first so eligible anionic precursors
reach the neutral keto representative in one pass. No charges are introduced.
The OH hydrogen transfers to the carbon, and bracket-specified hydrogen counts
on changed atoms are reset. Explicit hydrogen neighbours on the enol OH,
carboxyl OH, or reacting carbon prevent conversion, protecting isotope-labelled
hydrogen. Heavy-atom isotope labels and remote stereochemistry are retained.
Specified geometry and marked chiral centres on participating atoms prevent
conversion; no new stereocentre annotations are created.

The representative is a canonicalisation preference. For phenylpyruvic acid,
[solution NMR work](https://doi.org/10.1515/znc-1976-1-204) reports keto favoured
in water, while a [DMSO study](https://doi.org/10.2116/analsci.34.235) measured
93.5% enol at room temperature. The rule does not predict solvent-specific
populations and does not equate hydrated and non-hydrated structures.

Reports count conversions separately; `-g APP=EACH` appends
`STD:alpha-ketoacid` when the rule changes a molecule.

### Acyclic beta-enaminones

`-g beta-enaminone` chooses the conjugated keto–enamine representative of
eligible acyclic beta-ketoimines:

```text
CC(=O)CC(C)=NCC -> CC(=O)C=C(C)NCC
```

Both sides retain a ketone C=O. A hydrogen moves from CH2 to nitrogen, forming
the conjugated `O=C-C=C-NH` chain. Experimental solution NMR supports this
choice for closely related compounds, including an ethylene-linked
bis-enaminone with the same connectivity as the archive example
CHEMBL3211628/CHEMBL4752297. See [Eshghi et al., 2014](https://doi.org/10.1016/j.crci.2014.01.004)
and the acyclic N-aryl study by
[Darugar et al., 2021](https://doi.org/10.1016/j.molliq.2021.116035).
This is a canonicalisation preference, not a prediction of solution populations
for arbitrary substituents, solvents, or pH.

The rule is included in `-g all`, broadening default equivalence. Disable it
with `-g all -g -beta-enaminone`. It is skipped with `ENOL`. The initial scope
requires an acyclic ketone–CH2–C=N chain, a carbon substituent on the imine
carbon, and a carbon substituent on the neutral, two-connected imine nitrogen.
Non-aromatic carbon substituents at the ends of the reacting chain and at
nitrogen must be saturated; aromatic substituents are allowed. N-acyl,
N-hydrazone, cyclic reacting chains, acids/esters/amides in place of the ketone,
and substituted CH2 linkers are excluded. If both sides of the imine carbon
provide an eligible CH2–ketone group, conversion is skipped rather than making
an atom-order-dependent choice.

Specified geometry on the reacting bonds and marked chiral centres at
participating atoms prevent conversion. Explicit hydrogen neighbours on the
CH2 linker are excluded, preserving isotopic hydrogen. Heavy-atom isotope
labels and remote stereochemistry are retained. No double-bond geometry or
chiral-centre annotations are assigned to the product. Ordinary explicit
hydrogens can first be removed by `-g all`.

Reports count conversions separately; `-g APP=EACH` appends
`STD:beta-enaminone` for changed molecules.

### Acyclic beta-keto hydrazones

`-g beta-hydrazone` chooses the keto–hydrazone representative of eligible
hydrazine-substituted enamines:

```text
CC(=O)C=C(C)NNC(=O)c1ccncc1 -> CC(=O)CC(C)=NNC(=O)c1ccncc1
```

A hydrogen moves from the nitrogen next to the alkene onto the intervening
carbon, replacing C=C–NH with C–C=N. No charges or stereocentre annotations
are introduced. `-g beta-hydrazone` is included in `-g all`. The aim is consistent grouping
for modelling; this choice does not predict the dominant solution tautomer.
Disable it with `-g all -g -beta-hydrazone`.

The alternative `-g beta-enehydrazine` chooses the reverse representative:

```text
CC(=O)CC(C)=NNC(=O)c1ccncc1 -> CC(=O)C=C(C)NNC(=O)c1ccncc1
```

To use this preference with ordinary standardisation:

```shell
fileconv -g all -g -beta-hydrazone -g beta-enehydrazine -S standardised input.smi
```

The two options are mutually exclusive. Enabling both produces an error before
molecules are read; validation uses the final command-line configuration, so
`-g beta-enehydrazine -g all -g -beta-hydrazone` is also valid. C++ callers
configuring the object directly can use `ok()` to validate the final state;
`process()` refuses to modify molecules while both preferences are active.
Both directions are skipped with `ENOL`.

The same structural restrictions apply in both directions. A second eligible
CH2–ketone branch at the hydrazone carbon makes the reverse choice ambiguous;
such structures are left unchanged in both preferences, including the
corresponding enehydrazine forms. Defined geometry is preserved rather than
assigned or discarded during conversion. Newly generated double bonds have
unspecified geometry.

Solution preferences in related hydrazone systems vary. See
[Ahlbrecht and Henk, 1976](https://doi.org/10.1002/cber.19761090438)
for beta-ketoester/beta-ketonitrile derivatives and
[Sergienko et al., 2022](https://new.ras.ru/upload/iblock/0e3/jo4dmafqp2gg7krd57kz5l3xo1uer33z.pdf)
for related beta-dicarbonyl acylhydrazone solid structures. These do not
establish a universal solution preference for the ketone subset handled here.
The N-carbon enaminone preference described above cannot simply be extended
to N–N systems.

The initial scope requires an acyclic ketone, an unsubstituted CH between the
ketone and enamine carbon, a carbon substituent on the enamine carbon, and two
neutral NH nitrogens. The distal NH must have a carbon substituent; alkyl,
aryl, and acyl substituents are supported. Non-aromatic carbon substituents
at the two ends of the carbon chain must be saturated, excluding competing
conjugated branches. Acids, esters, amides in place of the ketone, thiones,
cyclic reacting chains, substituted linker carbons, and substituted or terminal
NH2 nitrogens are excluded. This is not a general imine/enamine rule.

Specified geometry on the reacting chain and marked chiral centres at
participating atoms prevent conversion. Heavy-atom isotope labels and remote
stereocentres are preserved. Explicit hydrogen attachments on either nitrogen
or the CH linker prevent conversion, so isotopic hydrogen is not transferred.
As usual, `-g all` can first remove ordinary explicit hydrogen atoms and
neutralise charged amines; the isolated rule only accepts neutral groups.

Reports count these conversions separately. `-g APP=EACH` appends
`STD:beta-hydrazone` or `STD:beta-enehydrazine`, according to the selected
preference, when the rule changes a molecule.

### Isolated cyclic amidines

`-g cyclic-amidine` chooses an endocyclic C=N bond and an exocyclic amino
nitrogen in isolated five- and six-membered cyclic amidines:

```text
N=C1NCCC1   ->  NC1=NCCC1
N=C1NCCCC1  ->  NC1=NCCCC1
```

This rule is included in `-g all` and broadens default tautomer equivalence.
Disable it with `-g all -g -cyclic-amidine`. Its fixed representative is also
used with the `ENOL` preference. It is separate from the non-ring `amidine`
rule, which prefers a terminal imine.

The ring must contain one nitrogen and otherwise saturated carbon atoms,
except for the amidine carbon. The starting ring nitrogen must be NH and the
exocyclic imine nitrogen terminal. Single-bonded substituents on the remaining
ring carbons are allowed. Fused/spiro rings, N-substituted groups, shared
amidines, conjugated or N-acyl systems, specified imine geometry, and marked
chiral centres at participating atoms are excluded. Remote stereocentres and
isotope labels on participating heavy atoms are retained; no new stereocentre
annotations or charges are created.

Explicit hydrogen atoms bonded to either participating nitrogen prevent this
conversion, protecting isotopic hydrogen. With `-g all`, ordinary explicit
hydrogens are removed and amine neutralisation runs before this rule, so
eligible protonated precursors can also reach the neutral representative in
one pass. `-g cyclic-amidine` alone does not neutralise charged groups.

Reports count cyclic-amidine conversions separately. `-g APP=EACH` appends
`STD:cyclic-amidine` when the rule changes a molecule.

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
