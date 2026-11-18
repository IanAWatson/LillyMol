#!/usr/bin/env python3
"""Generate deliberately complex recursive SMARTS equivalent to an aspirin query.

Run via mkaspirin.sh (or LillyMol's run_python.sh). Chemistry perception uses
LillyMol; generation itself needs only the Python standard library.
"""
import argparse
from dataclasses import dataclass
import random
import sys

ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O'
CONTEXTS = ('atom', 'bond', 'radius2', 'molecule')


@dataclass(frozen=True, slots=True)
class AtomInfo:
    atomic_number: int
    degree: int
    hydrogens: int
    aromatic: bool
    ring_bonds: int
    fused_size: int

    def predicate(self, rng, extensions):
        """Keep the baseline constraints and optionally repeat implied facts."""
        ring = self.ring_bonds > 0
        body = (f'0#{self.atomic_number};' + ('a' if self.aromatic else 'A') +
                f';H{self.hydrogens};D{self.degree};X{self.degree + self.hydrogens}' +
                f';R{int(ring)};x{self.ring_bonds};+0')
        if ring:
            body += ';r6'
        variants = ['', f'#{self.atomic_number}', f'D{self.degree}',
                    f'H{self.hydrogens}', f'X{self.degree + self.hydrogens}',
                    'a' if self.aromatic else 'A', 'C,c' if self.atomic_number == 6 else 'O']
        if extensions:
            variants += [f'D>{self.degree - 1}', f'D<{self.degree + 1}',
                         f'X>{self.degree + self.hydrogens - 1}',
                         f'H<{self.hydrogens + 1}']
            if self.hydrogens:
                variants.append(f'H>{self.hydrogens - 1}')
            # Always exercise a relational constraint when extensions are enabled.
            body += f';D>{self.degree - 1}'
            if ring:
                directive = rng.choice([f'/IWfss{self.fused_size}',
                                        f'/IWfss>{self.fused_size - 1}',
                                        f'/IWfss<{self.fused_size + 1}'])
                body = directive + ';' + body
                variants += ['r>5', 'r<7']
        variant = rng.choice(variants)
        return body + (';' + variant if variant else '')


@dataclass(frozen=True, slots=True)
class RingToken:
    number: int
    bond: str
    opening: bool


@dataclass(slots=True, eq=False)
class AtomNode:
    atom: int
    depth: int
    predicate: str
    alternatives: list | None = None


class AspirinModel:
    """LillyMol perceives the fixed graph; templates retain source atom indices."""
    def __init__(self):
        try:
            from lillymol import MolFromSmiles
        except ImportError as error:
            raise ValueError('LillyMol bindings are required; use mkaspirin.sh with '
                             'a Python matching lib/lillymol.so.soabi') from error
        self.molecule = MolFromSmiles(ASPIRIN)
        self.molecule.compute_aromaticity_if_needed()
        m = self.molecule
        self.atoms = [AtomInfo(m.atomic_number(i), m.ncon(i), m.hcount(i),
                              m.is_aromatic(i), m.ring_bond_count(i),
                              m.fused_system_size(i)) for i in range(m.natoms())]
        self.adjacency = [[] for _ in self.atoms]
        for bond in m.bonds():
            # Aromaticity must be checked before the underlying Kekule bond type.
            kind = ':' if bond.is_aromatic() else '=' if bond.is_double_bond() else '-'
            self.adjacency[bond.a1()].append((bond.a2(), kind))
            self.adjacency[bond.a2()].append((bond.a1(), kind))
        for neighbors in self.adjacency:
            neighbors.sort()
        self.templates = {(atom, context): self._template(atom, context)
                          for atom in range(m.natoms()) for context in CONTEXTS}

    def _template(self, root, context):
        """Render a small rooted source subgraph, leaving atoms as integer tokens.

        This traversal is at most 13 atoms deep. Large generated expressions are
        represented separately and never traversed recursively in Python.
        """
        radius = {'atom': 0, 'bond': 1, 'radius2': 2, 'molecule': 13}[context]
        allowed = {i for i in range(len(self.atoms))
                   if self.molecule.bonds_between(root, i) <= radius}
        visited, edges, children, closures = {}, set(), {}, []

        def visit(atom):
            visited[atom] = len(visited)
            children[atom] = []
            for other, bond in self.adjacency[atom]:
                edge = tuple(sorted((atom, other)))
                if other not in allowed or edge in edges:
                    continue
                edges.add(edge)
                if other in visited:
                    closures.append((other, atom, bond))
                else:
                    children[atom].append((other, bond))
                    visit(other)

        visit(root)
        rings = {atom: [] for atom in visited}
        for number, (first, second, bond) in enumerate(closures):
            rings[first].append(RingToken(number, bond, True))
            rings[second].append(RingToken(number, bond, False))

        def emit(atom):
            result = [atom, *rings[atom]]
            for i, (other, bond) in enumerate(children[atom]):
                branch = i + 1 < len(children[atom])
                if branch:
                    result.append('(')
                result.append(bond)
                result.extend(emit(other))
                if branch:
                    result.append(')')
            return result

        return emit(root)


class Generator:
    def __init__(self, model=None, seed=2026, contexts=('bond', 'radius2', 'molecule'),
                 max_depth=32, alternatives=1, lillymol_extensions=False):
        if max_depth < 0 or alternatives < 1:
            raise ValueError('max_depth must be nonnegative and alternatives must be positive')
        if not contexts or any(context not in CONTEXTS for context in contexts):
            raise ValueError('contexts must be chosen from ' + ','.join(CONTEXTS))
        self.model = model if model is not None else AspirinModel()
        self.rng = random.Random(seed)
        self.contexts = tuple(contexts)
        self.max_depth = max_depth
        self.alternatives = alternatives
        self.extensions = lillymol_extensions
        self.parts, self.size, self.expandable = self._instantiate(0, 'molecule', 0)
        self.expansions = 0
        self.depth = 0

    def _instantiate(self, root, context, depth):
        parts, leaves, size, labels = [], [], 0, {}
        for token in self.model.templates[root, context]:
            if isinstance(token, int):
                node = AtomNode(token, depth,
                                self.model.atoms[token].predicate(self.rng, self.extensions))
                parts.append(node)
                size += len(node.predicate) + 2
                if depth < self.max_depth:
                    leaves.append(node)
            else:
                if isinstance(token, RingToken):
                    if token.number not in labels:
                        unused = [i for i in range(1, 100) if i not in labels.values()]
                        labels[token.number] = self.rng.choice(unused)
                    number = labels[token.number]
                    token = (token.bond if token.opening else '') + (
                        str(number) if number < 10 else f'%{number}')
                parts.append(token)
                size += len(token)
        return parts, size, leaves

    def expand(self):
        if not self.expandable:
            raise ValueError(f'No expandable atoms remain at --max_depth {self.max_depth}; '
                             'increase the limit or request fewer expansions/a smaller target')
        # Swap/pop avoids linear-time deletion from a growing list of leaves.
        index = self.rng.randrange(len(self.expandable))
        node = self.expandable[index]
        self.expandable[index] = self.expandable[-1]
        self.expandable.pop()
        node.alternatives = []
        for _ in range(self.alternatives):
            context = self.rng.choice(self.contexts)
            parts, size, leaves = self._instantiate(node.atom, context, node.depth + 1)
            node.alternatives.append(parts)
            self.expandable.extend(leaves)
            # Each alternative adds $(...), and a comma or the initial semicolon.
            self.size += size + 4
        self.depth = max(self.depth, node.depth + 1)
        self.expansions += 1

    def smarts(self):
        """Serialize iteratively, including expressions deeper than Python's limit."""
        stack = list(reversed(self.parts))
        chunks = []
        while stack:
            token = stack.pop()
            if isinstance(token, str):
                chunks.append(token)
            elif token.alternatives is None:
                chunks.append('[' + token.predicate + ']')
            else:
                chunks.append('[' + token.predicate + ';')
                expansion = []
                for i, parts in enumerate(token.alternatives):
                    if i:
                        expansion.append(',')
                    expansion.append('$(')
                    expansion.extend(parts)
                    expansion.append(')')
                stack.append(']')
                stack.extend(reversed(expansion))
        return ''.join(chunks)


def nonnegative(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError('must be nonnegative')
    return number


def positive(value):
    number = nonnegative(value)
    if number == 0:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    limit = parser.add_mutually_exclusive_group()
    limit.add_argument('--nsmarts', type=nonnegative, help='Number of expansions (default 100); 0 writes the baseline')
    limit.add_argument('--target_bytes', type=positive, help='Expand until SMARTS reaches this size; may overshoot by one expansion')
    parser.add_argument('--write_last', action='store_true', help='Write only the final SMARTS, recommended for large targets')
    parser.add_argument('--seed', type=int, default=2026)
    parser.add_argument('--max_depth', type=nonnegative, default=32, help='Maximum nested $() levels (default 32)')
    parser.add_argument('--alternatives', type=positive, default=1, help='ORed recursive contexts per expansion (default 1)')
    parser.add_argument('--contexts', default='bond,radius2,molecule', help='Comma-separated atom,bond,radius2,molecule contexts')
    extensions = parser.add_mutually_exclusive_group()
    extensions.add_argument('--lillymol_extensions', action='store_true', help='Include /IWfss and >/< relational predicates')
    extensions.add_argument('--no_lillymol_extensions', dest='lillymol_extensions', action='store_false', help='Use only portable SMARTS (default)')
    parser.set_defaults(lillymol_extensions=False)
    args = parser.parse_args(argv)
    contexts = tuple(dict.fromkeys(item.strip() for item in args.contexts.split(',')))
    if not contexts or any(item not in CONTEXTS for item in contexts):
        parser.error('--contexts must be a comma-separated selection from ' + ','.join(CONTEXTS))
    try:
        generator = Generator(seed=args.seed, contexts=contexts, max_depth=args.max_depth,
                              alternatives=args.alternatives,
                              lillymol_extensions=args.lillymol_extensions)
        count = args.nsmarts if args.nsmarts is not None else 100
        written = 0
        while (generator.size < args.target_bytes if args.target_bytes is not None
               else generator.expansions < count):
            generator.expand()
            if not args.write_last:
                print(f'{generator.smarts()} #{generator.expansions}')
                written += 1
        if args.write_last or written == 0:
            print(f'{generator.smarts()} #{generator.expansions}')
        mode = 'LillyMol extensions' if args.lillymol_extensions else 'portable SMARTS'
        print(f'{generator.expansions} expansions, {generator.size} SMARTS bytes, '
              f'maximum recursion depth {generator.depth}, {mode}', file=sys.stderr)
    except (ValueError, OSError) as error:
        print(f'mkaspirin: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
