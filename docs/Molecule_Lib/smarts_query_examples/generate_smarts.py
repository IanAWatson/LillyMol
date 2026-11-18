#!/usr/bin/env python3
"""Generate portable SMARTS counterparts for the documented queries."""
from pathlib import Path

DONOR = '[$([O;X2;+0;H1]),$([O;X2;+0]-[C;X4]),$([N;X3;v3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N)])]'
WITHDRAWER = '[$(C#N),$(N(=O)=O),$([N+](=O)[O-]),$(C(=O)[#6,O,N]),$(C(F)(F)F)]'
EXCLUSIONS = ['C=O', 'C=S', 'C=N', 'S(=O)=O', 'P=O', 'C#N', 'N=O']


def ring(donors, exact=False):
    atoms = []
    for position in range(6):
        primitive = 'a;r6;R1'
        if exact:
            primitive += ';!$(a-' + WITHDRAWER + ')'
            if position not in donors:
                primitive += ';!$(a-' + DONOR + ')'
        atom = '[' + primitive + ']'
        if position == 0:
            atom += '1'
        if position in donors:
            atom += '(-' + DONOR + ')'
        atoms.append(atom)
    return ':'.join(atoms) + '1'


def patterns():
    # Dihedral symmetry leaves three positional classes for two substituents,
    # and three classes for three substituents on an otherwise generic hexagon.
    two = '[' + ','.join('$(' + ring({0, p}) + ')' for p in (1, 2, 3)) + ']'
    three = '[' + ','.join('$(' + ring(set(p), exact=True) + ')'
                          for p in ((0, 1, 2), (0, 1, 3), (0, 2, 4))) + ']'
    amine = '[N;X3;v3;+0;' + ';'.join('!$(N-' + item + ')' for item in EXCLUSIONS) + ']'
    diamine = amine + '-[C;X4]-[C;X4]-' + amine
    cyclic_amine = ('[N;X3;v3;+0;$(N1CCCCC1),$(N1CCOCC1),$(N1CCNCC1);'
                    '!$(N-C=[O,S,N]);!$(N-S=O);!$(N-P=O);!$(N-C#N)]'
                    '-[CH2]-[CH2]-[O;X2]-[a]')
    small_acyl = '[N;X3;v3;H1;+0]-[C](=[O])-[C;X4;$([CH3]),$([CH2]-[O;X2;H1;+0]),$([CH2]-F)]'
    return {'cyclic_aminoether': cyclic_amine, 'small_acyl_amide': small_acyl,
            'two_donors': two, 'three_donors_no_ewg': three,
            'amine': amine, 'diamine': diamine}


if __name__ == '__main__':
    root = Path(__file__).resolve().parent
    for name, smarts in patterns().items():
        (root / (name + '.smt')).write_text(smarts + ' ' + name + '\n')
