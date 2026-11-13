#!/usr/bin/env python3
"""Check final reassignment against all original destination members."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())
fixture = Path(__file__).resolve().parents[2] / 'gfp_leader_standard/reassignment/in/fingerprints.gfp'
records = [r + '|\n' for r in fixture.read_text().split('|\n') if r.strip()]
ids = [re.search(r'PCN<([^>]+)>', r)[1] for r in records]

def parse(text):
    groups = []
    for record in text.split('|\n'):
        names = re.findall(r'^PCN<([^>]+)>', record, re.M)
        if names:
            assert int(re.search(r'CSIZE<(\d+)>', record)[1]) == len(names)
            groups.append((names, list(map(float, re.findall(r'^DIST<([^>]+)>', record, re.M)))))
    return groups

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    def run(path, *args):
        r = subprocess.run([exe, *args, str(path)], capture_output=True, text=True, timeout=30)
        assert r.returncode == 0, r.stderr
        return r

    expanded = [re.sub(r'PCN<[^>]+>', f'PCN<m{i}_{j}>', records[i])
                for i, copies in enumerate((1, 3, 3, 2, 1, 1)) for j in range(copies)]
    path = tmp / 'pool.gfp'
    path.write_text(''.join(expanded))
    for metric in ([], ['-V', 'a=0.3', '-V', 'b=1.7']):
        matrix = {}
        for i in range(6):
            ordered = tmp / 'ordered.gfp'
            ordered.write_text(records[i] + ''.join(records[:i] + records[i+1:6]))
            groups = parse(run(ordered, '-t', '1', *metric).stdout)
            assert len(groups) == 1
            names, distances = groups[0]
            matrix.update({(ids[i], name): d for name, d in zip(names, [0] + distances)})
        def distance(a, b):
            return matrix[ids[int(a[1:].split('_')[0])], ids[int(b[1:].split('_')[0])]]

        for args, cutoff in ((['-t', '0.00001'], 1), (['-t', '0.35'], 2),
                             (['-t', '0.00001', '-m', '2'], 1),
                             (['-t', '0.00001', '-S', 'col=2'], 1)):
            # Score ordering must not change input-order nearest-member tie breaking.
            if '-S' in args:
                scored = []
                for index, record in enumerate(expanded):
                    scored.append(re.sub(r'PCN<([^>]+)>', rf'PCN<\1 {index}>', record))
                path.write_text(''.join(scored))
            else:
                path.write_text(''.join(expanded))
            original = parse(run(path, *args, *metric).stdout)
            targets = [g for g, _ in original if len(g) > cutoff]
            assert targets
            expected = {g[0]: set(g) for g in targets}
            order = re.findall(r'^PCN<([^>]+)>', path.read_text(), re.M)
            eligible = [item for item in order if any(item in g for g in targets)]
            for group, _ in original:
                if len(group) > cutoff:
                    continue
                for item in group:
                    nearest = min(eligible, key=lambda target: distance(target, item))
                    dest = next(g[0] for g in targets if nearest in g)
                    expected[dest].add(item)
            actual = parse(run(path, *args, *metric, '-R', str(cutoff)).stdout)
            assert {g[0]: set(g) for g, _ in actual} == expected
            for group, distances in actual:
                for item, d in zip(group[1:], distances):
                    assert abs(d - distance(group[0], item)) < 0.0002
            assert sum(len(g) for g, _ in actual) == len(expanded)
            sorted_actual = parse(run(path, *args, *metric, '-R', str(cutoff), '-r').stdout)
            assert {g[0]: set(g) for g, _ in sorted_actual} == expected
            assert all(ds == sorted(ds) for _, ds in sorted_actual)

    path.write_text(''.join(expanded))
    baseline = run(path, '-t', '0.00001').stdout
    no_target = run(path, '-t', '0.00001', '-R', '100')
    assert no_target.stdout == baseline and 'retaining original' in no_target.stderr
    assert run(path, '-t', '0.00001', '-C', '1').stdout == run(
        path, '-t', '0.00001', '-C', '1', '-R', '1').stdout
    previous = tmp / 'previous.gfp'
    previous.write_text(expanded[0])
    excluded = run(path, '-t', '0.00001', '-A', str(previous), '-R', '1').stdout
    assert 'PCN<m0_0>' not in excluded
    # Full search ignores the initial early-abandon threshold.
    assert run(path, '-t', '0.00001', '-R', '1').stdout == run(
        path, '-t', '0.00001', '-X', '0.00001', '-R', '1').stdout
    echoed = run(path, '-t', '0.00001', '-R', '1', '-E', 'ALL').stdout
    assert echoed.count('FPIW<') == len(expanded)
    proto = run(path, '-t', '0.00001', '-R', '1', '-k', '.').stdout
    assert 'PCN<' not in proto and 'name:' in proto
    sorted_output = parse(run(path, '-t', '0.00001', '-R', '1', '-r').stdout)
    assert all(ds == sorted(ds) for _, ds in sorted_output)
    for args in (['-R', '0'], ['-R', '1', '-D', 'sbc=0.1']):
        r = subprocess.run([exe, '-t', '0.2', *args, str(path)], capture_output=True)
        assert r.returncode != 0
print('gfp_leader reassignment checks passed')
