# Formal Charges
This directory contains queries that define formal charges applied to molecules.

These rules are from Fred Bruns, who worked at Lilly.

The rules actually applied are in the file [queries](queries). If you wish to use a subset
of the rules implemented, simply remove what you don't wish to use from that file.

These queries can inter-operate with the default Hydrogen Bond donor acceptor assignment
queries [hbonds](/data/queries/hbonds/AAREADME.md)


## Selected rules and limitations

Aromatic sulfonamides are deprotonated when the aromatic ring bears a direct
halogen or a methoxy group. These are separate queries so that both kinds of
activation can be found in the same molecule. The methoxy query requires
`[OD2-0]-[CD1H3]`; ethoxy, fluoromethyl and trifluoromethyl groups do not satisfy
that rule. The heteroaryl sulfonamide query separately recognises adjacent
aromatic oxygen, sulfur or unsubstituted pyridine-like nitrogen.

The amidine exclusions apply symmetrically to the two nitrogens: an aryl group
on the central carbon together with an aryl group on either nitrogen prevents
charging. Moving the imine double bond between the nitrogens therefore does
not change that exclusion. N-phenylacetamidines remain eligible for charging.

The phosphorus query permits at most two matches per invocation, across the
whole molecule. This is not a limit on the final accumulated charge: running
charge assignment repeatedly on a phosphate with three acidic oxygens can
remove another proton on a subsequent invocation. Molecules containing several
phosphorus centres are also subject to the same molecule-wide match limit.
