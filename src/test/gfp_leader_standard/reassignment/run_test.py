#!/usr/bin/env python3
"""Compare both reassignment modes against independent distance searches."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())
source = Path(sys.argv[2]) / 'fingerprints.gfp'
records = [s + '|\n' for s in source.read_text().split('|\n') if s.strip()]

def parse(text):
    result = []
    for record in text.split('|\n'):
        ids = re.findall(r'^PCN<([^>]+)>', record, re.M)
        if ids:
            size = int(re.search(r'CSIZE<(\d+)>', record)[1])
            assert size == len(ids)
            result.append((ids, list(map(float, re.findall(r'DIST<([^>]+)>', record)))))
    return result

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    def run(path, *args):
        r = subprocess.run([exe, '-h', '1', *args, str(path)], text=True,
                           capture_output=True, timeout=30)
        assert r.returncode == 0, r.stderr
        return r

    # Distances are measured independently by putting each item first in turn.
    ids = [re.search(r'PCN<([^>]+)>', r)[1] for r in records]
    matrix = {}
    for i in range(len(records)):
        path = tmp / 'ordered.gfp'
        path.write_text(records[i] + ''.join(records[:i] + records[i+1:]))
        groups = parse(run(path, '-t', '0.999999').stdout)
        assert len(groups) == 1
        names, distances = groups[0]
        matrix.update({(ids[i], name): distance for name, distance in zip(names, [0] + distances)})

    # Singleton first, two larger groups, a pair, and more singletons.
    expanded = []
    for i, copies in enumerate([1, 3, 3, 2, 1, 1]):
        for j in range(copies):
            expanded.append(re.sub(r'PCN<[^>]+>', f'PCN<m{i}_{j}>', records[i]))
    path = tmp / 'pool.gfp'
    path.write_text(''.join(expanded))
    def distance(a, b):
        return matrix[ids[int(a[1:].split('_')[0])], ids[int(b[1:].split('_')[0])]]

    for threshold in ('0.00001', '0.35'):
        original = parse(run(path, '-t', threshold).stdout)
        for cutoff in (1, 2):
            for suffix in ('', ':any'):
                targets = [g for g, _ in original if len(g) > cutoff]
                expected = {g[0]: list(g) for g in targets}
                assert targets
                order = [re.search(r'PCN<([^>]+)>', r)[1] for r in expanded]
                target_members = [name for name in order if any(name in g for g in targets)]
                if not suffix:
                    target_members = [name for name in target_members if name in expected]
                for group, _ in original:
                    if len(group) > cutoff:
                        continue
                    for name in group:
                        closest = min(target_members, key=lambda other: distance(name, other))
                        destination = next(g[0] for g in targets if closest in g)
                        expected[destination].append(name)
                actual = parse(run(path, '-t', threshold, '-R', str(cutoff) + suffix).stdout)
                assert {g[0]: set(g) for g, _ in actual} == {k: set(v) for k, v in expected.items()}
                for group, distances in actual:
                    for name, d in zip(group[1:], distances):
                        assert abs(d - distance(group[0], name)) < 0.00002
                # Determinism across worker counts.
                assert run(path, '-t', threshold, '-R', str(cutoff) + suffix).stdout == run(
                    path, '-h', '4', '-t', threshold, '-R', str(cutoff) + suffix).stdout

    # Exercise the parallel branch with more than 32 reassignment candidates.
    large = tmp / 'parallel.gfp'
    large.write_text(''.join(re.sub(r'PCN<[^>]+>', f'PCN<p{i}_{j}>', records[i])
                             for i, count in enumerate((40, 60)) for j in range(count)))
    serial = run(large, '-t', '0.00001', '-R', '40').stdout
    parallel = run(large, '-h', '4', '-t', '0.00001', '-R', '40').stdout
    assert serial == parallel
    groups = parse(serial)
    assert len(groups) == 1 and groups[0][0][0] == 'p1_0' and len(groups[0][0]) == 100

    baseline = run(path, '-t', '0.00001')
    no_targets = run(path, '-t', '0.00001', '-R', '100')
    assert baseline.stdout == no_targets.stdout and 'retaining original' in no_targets.stderr
    assert run(path, '-t', '0.00001', '-n', '1', '-R', '1').stdout == run(
        path, '-t', '0.00001', '-n', '1').stdout
    previous = tmp / 'previous.gfp'
    previous.write_text(expanded[0])
    excluded = parse(run(path, '-t', '0.00001', '-A', str(previous), '-R', '1').stdout)
    assert all('m0_0' not in g for g, _ in excluded)
    proto = run(path, '-t', '0.00001', '-R', '1', '-k', '.').stdout
    assert 'PCN<' not in proto and 'name:' in proto
    for value in ('0', '0:any', '1:bad', '1:any:any', ':any'):
        bad = subprocess.run([exe, '-t', '0.2', '-R', value, str(path)], capture_output=True)
        assert bad.returncode != 0
print('gfp_leader_standard reassignment checks passed')
