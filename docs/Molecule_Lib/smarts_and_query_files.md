# When SMARTS becomes hard to maintain

SMARTS is useful for compact local structural patterns. As requirements grow,
a single expression can become difficult to read, review and change. LillyMol
query files let us name and separate concepts such as rings, environments,
substituents and regions, while retaining SMARTS for the local patterns where
it is convenient.

These examples compare executable queries with their stated intent. They are
not a claim that a long SMARTS is impossible to shorten, or that character count
alone measures quality. The practical question is whether another person can
understand the requirements and change one without accidentally changing another.

## Run the examples

All query files, input molecules, expected results and validation scripts are
in [smarts_query_examples](smarts_query_examples/). For example:

```shell
tsubstructure \
  -q PROTO:${LILLYMOL_HOME}/docs/Molecule_Lib/smarts_query_examples/two_donors.textproto \
  -m matches molecules.smi
```

Run the complete validation from any directory:

```shell
python3 ${LILLYMOL_HOME}/docs/Molecule_Lib/smarts_query_examples/validate_examples.py \
  --tsubstructure ${LILLYMOL_HOME}/bin/Linux/tsubstructure
```

Use `bin/Darwin` on a Mac. The validator uses the LillyMol command-line tool;
it does not require LillyMol Python bindings. To also sanity-check the portable
SMARTS with RDKit, run it in a Python environment containing RDKit and add
`--rdkit`. RDKit is used only for this independent comparison:

```shell
python3 ${LILLYMOL_HOME}/docs/Molecule_Lib/smarts_query_examples/validate_examples.py \
  --tsubstructure ${LILLYMOL_HOME}/bin/Linux/tsubstructure --rdkit
```

The validator runs in a temporary directory, so shared SMARTS-list imports must
resolve relative to their query files rather than the current directory. It
also checks expected matching identifiers, the generated portable SMARTS, and
isotope counts for the region-selection example.

## 1. An isolated aromatic six-membered ring with at least two donors

First define what “donor” means. There is no universal electron-donating-group
primitive in SMARTS. Here we deliberately use an illustrative structural list:

- Neutral hydroxyl oxygen.
- Neutral divalent oxygen attached to an aliphatic tetrahedral carbon, representing
  alkoxy groups.
- Neutral, trivalent nitrogen with three total connections (including hydrogens),
  excluding attachment to the
  specified carbon/sulfur/phosphorus multiple-bond groups and to a cyano carbon.

The nitrogen valence restriction matters: LillyMol can represent nitro groups
as neutral `N(=O)=O`, which must not be mistaken for an amine. Alkyl groups and
other possible donors are deliberately absent. Neither this list nor the later
withdrawing-group list should be treated as a general electronic-effects model.

“Isolated” means that the ring shares no bonds or atoms with another ring.
It can still be connected to another ring through an ordinary bond, as in a
biphenyl. The ring can contain aromatic heteroatoms. For benzene rings alone, use `c`
throughout the portable ring pattern and add `heteroatom_count: 0` to the ring
`base` message in the proto query.

### Standard SMARTS

One explicit construction enumerates the ortho, meta and para relationships of
two qualifying positions. Every ring atom is aromatic, in a six-membered ring,
and a member of exactly one perceived ring:

```text
[$([a;r6;R1]1(-[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]):[a;r6;R1](-[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]):[a;r6;R1]:[a;r6;R1]:[a;r6;R1]:[a;r6;R1]1),$([a;r6;R1]1(-[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]):[a;r6;R1]:[a;r6;R1](-[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]):[a;r6;R1]:[a;r6;R1]:[a;r6;R1]1),$([a;r6;R1]1(-[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]):[a;r6;R1]:[a;r6;R1]:[a;r6;R1](-[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]):[a;r6;R1]:[a;r6;R1]1)] two_donors
```

The trailing `two_donors` is the query name in the `.smt` file, not part of the
SMARTS. The SMARTS itself has **724 characters**, and the donor definition occurs
six times. The two selected positions are distinct, but extra donors are allowed.

These are three relative arrangements, not three separate chemical definitions.
Changing which groups count as donors requires updating all six copies. Changing
the required number of donors requires changing the positional construction.
[generate_smarts.py](smarts_query_examples/generate_smarts.py) makes that expansion
reproducible; the expanded files do not require the generator at search time.

### LillyMol query file

```textproto
name: "Isolated aromatic six-membered ring with at least two donor positions"
query {
  ring_specifier {
    ring_size: 6
    aromatic: true
    fused: 0
    spiro_fusion_count: 0
    base {
      environment: ">1a-[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]"
    }
  }
}
```

Here `>1` counts matching ring positions. In a ring environment the first atom
of the environment belongs to the ring; subsequent atoms describe the attached
group. The donor definition occurs once. Changing the threshold to at least
three changes `>1` to `>2`.

The queries agree on all 19 examples in
[two_donors.smi](smarts_query_examples/two_donors.smi): nine match. Positive
examples include catechol, resorcinol, hydroquinone, mixed donors, a three-donor
ring, and a heteroaromatic ring. Negative examples include phenol, a fused
aromatic system, one donor on each of two separate rings, and a phenol whose
other substituent is an acylated amine, sulfonamide, ammonium, cyanamide, nitro
or acetoxy group.

This is a comparison of molecule-level presence. The portable SMARTS returns
its root atoms; the ring-only proto is a ring predicate. Their embedding counts
and returned atoms are not asserted to be identical.

## 2. Exactly three donor positions and no withdrawing groups on that ring

Now refine the question without redefining the chemistry. We require exactly
three donor-bearing ring positions, and forbid an illustrative list of directly
attached withdrawing groups: nitrile, nitro, carbonyl attached through its carbon,
and CF3. Both nitro charge representations are included for input portability.
Other substituents, such as methyl, are allowed.

A fully expanded standard SMARTS is provided in
[three_donors_no_ewg.smt](smarts_query_examples/three_donors_no_ewg.smt). It is
**3,094 characters**. Three positional classes cover the three selected donor
positions; every remaining position must fail the donor test, and every ring
position must fail the withdrawing-group test. It is executable standard SMARTS,
not a proposed new syntax or an abbreviation.

In the LillyMol query we separate the ring from the attachment requirements:

```textproto
name: "Exactly three donor positions and no withdrawing groups on one isolated ring"
query {
  ring_specifier {
    ring_size: 6
    aromatic: true
    fused: 0
    spiro_fusion_count: 0
    base { set_global_id: 1 }
  }
  # Return the six ring atoms; the first must belong to a qualifying ring.
  smarts: "[/IWgid1a]1:a:a:a:a:a:1"
  environment {
    attachment {
      attachment_point: [0, 1, 2, 3, 4, 5]
      btype: SS_SINGLE_BOND
    }
    smarts_file: "donors.textproto"
    hits_needed: 3
    max_env_matches_per_anchor: 1
  }
  environment_no_match {
    attachment {
      attachment_point: [0, 1, 2, 3, 4, 5]
      btype: SS_SINGLE_BOND
    }
    smarts_file: "withdrawers.textproto"
  }
}
```

The six atoms in the main SMARTS are one ring. `set_global_id` and `/IWgid1`
connect that ring to the isolation/aromaticity constraints. The attachment list
allows an environment at any of its six positions. `hits_needed: 3` makes the
count exact; `max_env_matches_per_anchor: 1` ensures we count ring positions,
rather than counting multiple matches at one position.

The chemical definitions live in two reusable files:

[donors.textproto](smarts_query_examples/donors.textproto):

```textproto
# Neutral hydroxyl, alkoxy, and selected amine-like nitrogen.
smarts: "[O;X2;+0;H1]"
smarts: "[O;X2;+0]-[C;X4]"
smarts: "[N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)]"
```

[withdrawers.textproto](smarts_query_examples/withdrawers.textproto):

```textproto
# Illustrative list, not a general classification of electronic effects.
smarts: "C#N"
smarts: "N(=O)=O"
smarts: "[N+](=O)[O-]"
smarts: "C(=O)[#6,O,N]"
smarts: "C(F)(F)F"
```

Five of 14 molecules match. The tests cover all three positional classes, mixed
donors, an allowed methyl substituent, too few or too many donors, each forbidden
group, a fused system and donors split across rings. Adding a new withdrawing
group edits one shared list, rather than each repeated negative recursive test.

## 3. Exactly two complete acyclic substituents, each containing 3–8 atoms

This requirement concerns the **entire substituent**, rather than the existence
of a short path inside it. An isolated aromatic six-ring must have exactly two
external connections. Each attached component must contain 3–8 atoms and no rings:

```textproto
name: "Isolated aromatic six-ring with exactly two acyclic 3-to-8-atom substituents"
query {
  ring_specifier {
    ring_size: 6
    aromatic: true
    fused: 0
    spiro_fusion_count: 0
    base {
      ncon: 2
      substituent {
        min_natoms: 3
        max_natoms: 8
        nrings: 0
        hits_needed: 2
      }
    }
  }
}
```

`ncon: 2` limits the ring's total external connections. `hits_needed: 2` requires
two qualifying substituents. Branches and heteroatoms count towards their full
sizes. A nine-atom chain fails even though it contains many three-atom paths;
a cyclopropyl-containing substituent fails even when its total size is allowed.

Four of ten molecules match, including branched and heteroatom-containing
substituents. The negatives check the lower and upper size boundaries, a third
substituent, a cyclic sidechain, a fused central ring and only one substituent.

There is no standard SMARTS primitive meaning “the complete attached connected
component has 3–8 atoms.” One could enumerate bounded atom/bond topologies, or
match candidate anchors and measure their substituents in code. This example
uses the query-file concept directly; no portable SMARTS equivalent is supplied.

## 4. An amine with a reusable exclusion list

A local nitrogen test is compact. A growing list of disqualifying attachments
makes the same test cumbersome when repeated at several sites. Here a neutral
trivalent nitrogen is excluded if attached to any group in
[amine_exclusions.textproto](smarts_query_examples/amine_exclusions.textproto):

```textproto
# Groups disqualifying a nitrogen when attached through the first atom.
smarts: "C=O"
smarts: "C=S"
smarts: "C=N"
smarts: "S(=O)=O"
smarts: "P=O"
smarts: "C#N"
smarts: "N=O"
```

The one-site query is:

```textproto
name: "amine with reusable exclusions"
query {
  smarts: "[N;X3;v3;+0]"
  environment_no_match {
    attachment { attachment_point: 0 btype: SS_SINGLE_BOND }
    smarts_file: "amine_exclusions.textproto"
  }
}
```

The two-site query reuses the same list at both nitrogens:

```textproto
name: "diamine with reusable exclusions"
query {
  smarts: "[N;X3;v3;+0]-[C;X4]-[C;X4]-[N;X3;v3;+0]"
  environment_no_match {
    attachment { attachment_point: 0 btype: SS_SINGLE_BOND }
    smarts_file: "amine_exclusions.textproto"
  }
  environment_no_match {
    attachment { attachment_point: 3 btype: SS_SINGLE_BOND }
    smarts_file: "amine_exclusions.textproto"
  }
}
```

Both queries refer to the same file; no recursive exclusions need to be copied
between the nitrogen sites or between queries. Relative imports resolve against
the query-file directory. `${LILLYMOL_HOME}` can also be used in imported paths.

Portable equivalents are included as [amine.smt](smarts_query_examples/amine.smt)
and [diamine.smt](smarts_query_examples/diamine.smt), respectively 86 and 187
characters. These are not particularly long SMARTS, which is useful context:
query files offer naming and reuse even before expressions become enormous.

Eight of 19 molecules have at least one qualifying nitrogen. Only two have the
two qualifying nitrogens linked through `N-C-C-N`. Molecules with one blocked
nitrogen still match the one-site query, but fail the two-site query. Other
negatives include amide, urea, thiourea, amidine, sulfonamide, phosphoramidate,
cyanamide, nitramine, nitrosamine and quaternary ammonium examples.

## 5. Identify and label the largest qualifying inter-ring region

This is a selection operation as well as a structural search. Find non-ring
scaffold regions connecting exactly two rings, require at least two region atoms,
then retain those with the greatest atom count:

```textproto
name: "Largest qualifying inter-ring region; retain ties"
query {
  inter_ring_region {
    min_natoms: 2
    ring_connections: 2
    selection: MAX_ATOMS
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
  # Put every selected atom into one embedding for isotope labelling.
  compress_embeddings: true
}
```

The region constraints are applied before selection. `MAX_ATOMS` keeps all tied
largest qualifying regions. `set_global_id` makes their atoms available to the
main SMARTS; `compress_embeddings` combines them for labelling:

```shell
tsubstructure \
  -q PROTO:${LILLYMOL_HOME}/docs/Molecule_Lib/smarts_query_examples/largest_region.textproto \
  -j 1 -j same -m labelled molecules.smi
```

The tests verify four labelled atoms when choosing a four-atom linker over a
two-atom linker, four labelled atoms across two tied two-atom linkers, and two
labelled scaffold atoms when a pendant methyl must be excluded. A one-atom linker
and a molecule without rings fail.

SMARTS can identify candidate paths and local patterns. Selecting the largest
connected region from all qualifying regions normally calls for postprocessing.
The query file expresses that operation directly. It also makes the treatment
of ties explicit; it does not arbitrarily select the first region encountered.

## 6. Rooted atom constraints

`matched_atom_must_be` is another way to move recursive atom tests out of the
main SMARTS. The main query describes the scaffold or relationship of interest;
a separate block describes the acceptable structural roles of a particular
matched atom. Its SMARTS are **rooted**: their first atom must be the specified
matched atom, not just any atom somewhere in the molecule.

Positive entries in one block are alternatives: at least one must match. Entries
starting with `!` are exclusions: any one matching rejects the embedding. If the
block contains only exclusions, the specified atom must avoid all of them.
The leading `!` here is part of LillyMol's message-level convention, rather than
an ordinary standalone SMARTS negation.

### A selected cyclic amine on an aminoethoxyaryl scaffold

A realistic library-selection question is: find an aminoethoxyaryl motif whose
chain nitrogen is part of a piperidine, morpholine or piperazine skeleton, while
excluding selected carbonyl and other attachments at that same nitrogen.
This is a structural selection rule, not a prediction of protonation or basicity.

The portable recursive-atom version is:

```text
[N;X3;v3;+0;$(N1CCCCC1),$(N1CCOCC1),$(N1CCNCC1);!$(N-C=[O,S,N]);!$(N-S=O);!$(N-P=O);!$(N-C#N)]-[CH2]-[CH2]-[O;X2]-[a] cyclic_aminoether
```

The corresponding query separates the motif from the nitrogen's permitted and
forbidden contexts:

```textproto
name: "Aminoethoxyaryl with a selected cyclic amine at the chain nitrogen"
query {
  # Query atom 0 is the chain nitrogen, not a nitrogen elsewhere in the molecule.
  smarts: "[N;X3;v3;+0]-[CH2]-[CH2]-[O;X2]-[a]"
  matched_atom_must_be {
    atom: 0
    # At least one positive pattern must match, rooted at that same nitrogen.
    smarts: "N1CCCCC1"  # Piperidine skeleton.
    smarts: "N1CCOCC1"  # Morpholine skeleton.
    smarts: "N1CCNCC1"  # Piperazine skeleton.
    # Any one of these rooted matches rejects the embedding.
    smarts: "!N-C=[O,S,N]"
    smarts: "!N-S=O"
    smarts: "!N-P=O"
    smarts: "!N-C#N"
  }
}
```

Query atom zero is the chain nitrogen. The three positive entries describe its
allowed ring skeletons, including substituted versions. The negative entries
then exclude selected undesirable roles even if the ring skeleton matched.
A lactam can contain the `N1CCCCC1` skeleton: the separate exclusion makes clear
why it is nevertheless rejected.

Five of 15 molecules pass. The positives include all three ring classes, a
substituted piperidine and a piperazine acylated on its **other** nitrogen.
Negatives include lactam, thiolactam and cyclic amidine anchors; five- and
seven-membered rings; an acyclic amine; a charged nitrogen; and a longer chain.
Another negative contains a morpholine elsewhere in the molecule but has the
wrong nitrogen at the aminoethoxyaryl anchor. Finding the right group somewhere
in the molecule is not sufficient.

The validator labels exactly the five main-query atoms in every positive.
The other atoms of the cyclic-amine context are not added to the embedding.

### Small acyl groups on a secondary amide

For an amide-focused library, suppose the allowed acyl groups are acetyl,
glycolyl and fluoroacetyl. The main query describes a secondary amide with an
aliphatic carbon next to its carbonyl. The rooted alternatives describe that
specific carbon:

```textproto
name: "Secondary amide with acetyl, glycolyl or fluoroacetyl acyl group"
query {
  # Atom order is N=0, carbonyl C=1, O=2, alpha carbon=3.
  smarts: "[N;X3;v3;H1;+0]-[C](=[O])-[C;X4]"
  matched_atom_must_be {
    atom: 3
    smarts: "[CH3]"
    smarts: "[CH2]-[O;X2;H1;+0]"
    smarts: "[CH2]-F"
  }
}
```

Compare the portable form:

```text
[N;X3;v3;H1;+0]-[C](=[O])-[C;X4;$([CH3]),$([CH2]-[O;X2;H1;+0]),$([CH2]-F)] small_acyl_amide
```

The branch oxygen is query atom 2, so the alpha carbon is atom 3. This is one
reason to retain short main SMARTS and comment their atom ordering: the index
refers to the query, not a target-molecule atom number or an atom-map label.

Five of 14 molecules pass. Negatives include propionyl, branched acyl,
difluoroacetyl, methoxyacetyl and chloroacetyl groups, plus primary and tertiary
amides and a urea. A propionamide carrying an alcohol on its N-side chain still
fails: the hydroxyl must be connected to the **acyl alpha carbon**, not merely
present elsewhere. All positives have exactly four main-query atoms labelled.

### One rejection list applied to both matched nitrogens

The earlier diamine example can also be written with a single matched-atom block:

```textproto
name: "Two amine sites with one shared rooted rejection list"
query {
  smarts: "[N;X3;v3;+0]-[C;X4]-[C;X4]-[N;X3;v3;+0]"
  matched_atom_must_be {
    atom: [0, 3]
    smarts: "!N-C=O"
    smarts: "!N-C=S"
    smarts: "!N-C=N"
    smarts: "!N-S(=O)=O"
    smarts: "!N-P=O"
    smarts: "!N-C#N"
    smarts: "!N-N=O"
  }
}
```

Every atom listed in `atom: [0, 3]` must satisfy the same rooted restrictions.
The atom list is an AND requirement, not a choice of which atom to inspect.
The two accepted molecules and all negatives agree with the earlier environment
version and with the portable [diamine SMARTS](smarts_query_examples/diamine.smt).
Isotope labelling returns just the four atoms of the `N-C-C-N` motif.

This avoids duplicating an atomic recursive expression at both nitrogen sites.
Unlike `environment` and `environment_no_match`, `matched_atom_must_be` does not
currently have a `smarts_file` field: the sharing here is within one block.
Use the environment form when a shared external list across queries is the
main requirement.

### Choosing between the constructs

Use `matched_atom_must_be` when the condition is naturally about the structural
role of an already matched atom. Its rooted context can revisit atoms already
in the main embedding. It constrains an embedding without enlarging it.

Use an environment when the condition is naturally about additional groups
attached to matched atoms, particularly when counting those groups or importing
shared lists. Environments ordinarily search unmatched context atoms; this is
a different matching rule, not just another spelling of the rooted test.

Neither construct eliminates the need to define the chemistry carefully. They
make it easier to explain which condition belongs to which matched atom, and
to test alternatives and exclusions independently of the scaffold.

## What was checked

The examples were validated with the current LillyMol `tsubstructure` and RDKit
2024.09.5. RDKit independently parsed and searched the six distinct portable SMARTS used by seven cases;
all matching identifier sets agreed with LillyMol and the expected results.
The two feature-specific examples use LillyMol alone. In total the validator
performs 23 checks, plus the isotope-count assertions.

| Example | Input molecules | Matching molecules | Portable SMARTS checked in RDKit |
| --- | ---: | ---: | --- |
| At least two donors | 19 | 9 | Yes |
| Exactly three donors, no withdrawing group | 14 | 5 | Yes |
| Two complete acyclic substituents | 10 | 4 | No equivalent supplied |
| One qualifying amine | 19 | 8 | Yes |
| Two qualifying amines | 19 | 2 | Yes |
| Largest inter-ring region | 5 | 3 | Selection uses LillyMol |
| Selected cyclic amine at an aminoether anchor | 15 | 5 | Yes |
| Small acyl groups on a secondary amide | 14 | 5 | Yes |
| Diamine with shared matched-atom restrictions | 19 | 2 | Yes |

The comparisons concern these fixtures, not a proof of equivalence for every
possible molecular graph. Aromaticity, ring perception and charge normalization
can differ between toolkits. `R1` refers to toolkit-perceived ring membership,
not a ring-basis-independent graph-theoretic definition. These fixtures do not
run chemical standardization; RDKit may normalize the neutral nitro input to its
charge-separated representation when reading it.

The examples intentionally separate chemical policy from query mechanics. The
most useful demonstration is often a maintenance exercise: change the donor
list, change “at least” to “exactly,” add a size limit, or reuse the exclusions
at another site, and inspect which parts of the query have to change.

See [Substructure searching](substructure.md) for SMARTS extensions and
[Textproto query files](substructure_proto.md) for the message-based query syntax.
