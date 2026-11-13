#!/usr/bin/env python3
"""Exercise stored-neighbour reassignment for leader and Taylor–Butina."""
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tfrecord_utils import tfrecord

exe = str(Path(sys.argv[1]).resolve())

def parse(text):
    groups = []
    for block in text.split('|\n'):
        names, distances = [], {}
        for line in block.splitlines():
            if line.startswith('PCN<'):
                names.append(line[4:-1])
            elif line.startswith('DIST<'):
                distances[names[-1]] = float(line[5:-1])
            elif line.startswith('CSIZE<'):
                size = int(line[6:-1])
        if names:
            assert len(names) == size
            groups.append((names, distances))
    return groups

with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    def graph(names, edges):
        result = {name: {} for name in names}
        for a, b, d in edges:
            result[a][b] = d
            result[b][a] = d
        return result

    def write(data, serialized):
        path = root / ('input.rio' if serialized else 'input.nn')
        # Reverse index order among equal-distance entries deliberately: ties must
        # resolve by molecule input order, not neighbour-list order.
        neighbours = {a: sorted(b.items(), key=lambda x: (x[1], -list(data).index(x[0])))
                      for a, b in data.items()}
        if serialized:
            path.write_bytes(b''.join(tfrecord(a, b) for a, b in neighbours.items()))
        else:
            path.write_text(''.join(f'$SMI<C>\nPCN<{a}>\n' +
                ''.join(f'PCN<{b}>\nDIST<{d}>\n' for b, d in nbrs) + '|\n'
                for a, nbrs in neighbours.items()))
        return path

    def run(path, serialized, mode, *args):
        result = subprocess.run([exe, '-h', *(['-r'] if serialized else []), *mode,
                                 *args, str(path)], capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stderr
        return result

    def check(data, serialized, mode, cutoff, *args):
        path = write(data, serialized)
        baseline = run(path, serialized, mode, *args)
        initial = parse(baseline.stdout)
        original = {item: g[0] for g, _ in initial for item in g}
        sizes = {g[0]: len(g) for g, _ in initial}
        assignment = dict(original)
        for item in original:
            if sizes[original[item]] > cutoff:
                continue
            eligible = [other for other in data[item] if other in original and
                        sizes[original[other]] > cutoff]
            if eligible:
                closest = min(eligible, key=lambda other: (data[item][other], list(data).index(other)))
                assignment[item] = original[closest]
        expected = {}
        for group, _ in initial:
            old_root = group[0]
            members = [item for item in data if assignment.get(item) == old_root]
            if members:
                new_root = old_root if old_root in members else members[0]
                expected[new_root] = set(members)
        result = run(path, serialized, mode, *args, '-R', str(cutoff))
        actual = parse(result.stdout)
        assert {g[0]: set(g) for g, _ in actual} == expected, (actual, expected)
        for group, distances in actual:
            leader = group[0]
            for item in group[1:]:
                d = data[leader].get(item, data[item].get(leader))
                if d is None:
                    assert item not in distances
                else:
                    assert abs(distances[item] - d) < 0.0001
        return result, baseline

    data = graph(['s', 'a', 'b', 'c', 'd', 'e', 'f', 't'], [
        ('a', 'b', .1), ('c', 'd', .1), ('e', 'f', .1),
        ('s', 'b', .4), ('s', 'c', .45), ('s', 'a', .8),
        ('t', 'd', .4), ('t', 'b', .4)])
    split_pair = graph(['c', 'd', 'a', 'b', 'h', 'e', 'f', 'g', 's'], [
        ('c', 'd', .1), ('a', 'b', .1), ('a', 'h', .1),
        ('e', 'f', .1), ('e', 'g', .1), ('c', 'b', .35), ('c', 'e', .5),
        ('d', 'f', .3), ('d', 'a', .6)])
    partial = graph(['c', 'd', 'a', 'b', 'h'], [
        ('c', 'd', .1), ('a', 'b', .1), ('a', 'h', .1), ('c', 'b', .35)])
    for serialized in (False, True):
        for mode in ([], ['-T']):
            for cap in ([], ['-m', '2']):
                result, _ = check(data, serialized, mode, 1, '-t', '.2', *cap)
            check(split_pair, serialized, mode, 2, '-t', '.2')
            check(partial, serialized, mode, 2, '-t', '.2')
            asymmetric = {a: dict(b) for a, b in data.items()}
            del asymmetric['a']['s']  # New leader distance exists only in s's list.
            check(asymmetric, serialized, mode, 1, '-t', '.2')
            only_large = graph(['a', 'b', 'c', 'd'], [('a', 'b', .1), ('c', 'd', .1)])
            result, baseline = check(only_large, serialized, mode, 1, '-t', '.2')
            assert result.stdout == baseline.stdout
            result, baseline = check(data, serialized, mode, 100, '-t', '.2')
            assert result.stdout == baseline.stdout and 'retained' in result.stderr
            check(data, serialized, mode, 1, '-t', '.2', '-C', '1')
            check(data, serialized, mode, 1, '-t', '.2', '-s', '4')
            # Threshold workers must have private memberships and reassignment state.
            path = write(data, serialized)
            stem = root / 'clusters'
            run(path, serialized, mode, '-t', '.2,.5', '-R', '1', '-B', str(stem))
            for index, threshold in enumerate(('.2', '.5')):
                serial = run(path, serialized, mode, '-t', threshold, '-R', '1')
                assert (root / f'clusters{index}.ldr').read_text() == serial.stdout
            smiles = run(path, serialized, mode, '-t', '.2', '-R', '1', '-a').stdout
            assert sum(line.startswith('C ') for line in smiles.splitlines()) == len(data)

        # Leader-specific assertions: attraction by a nonleader and input-order tie.
        result, _ = check(data, serialized, [], 1, '-t', '.2')
        groups = {g[0]: (set(g), ds) for g, ds in parse(result.stdout)}
        assert {'s', 't'} <= groups['a'][0] and groups['a'][1]['s'] == .8
        assert 't' not in groups['a'][1] and 'omitted unavailable leader distances' in result.stderr
        result, _ = check(split_pair, serialized, [], 2, '-t', '.2')
        groups = {g[0]: set(g) for g, _ in parse(result.stdout)}
        assert 'c' in groups['a'] and 'd' in groups['e']
        result, _ = check(partial, serialized, [], 2, '-t', '.2')
        groups = {g[0]: set(g) for g, _ in parse(result.stdout)}
        assert groups['d'] == {'d'} and 'c' in groups['a']
        path = write(data, serialized)
        for options in (['-t', '.2', '-R', '0'], ['-J', '.,1', '-R', '1']):
            result = subprocess.run([exe, *(['-r'] if serialized else []),
                                     *options, str(path)], capture_output=True)
            assert result.returncode != 0
print('nn_leader_and_jp reassignment checks passed')
