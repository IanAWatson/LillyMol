#!/usr/bin/env python3
"""Compare compact candidates with a simple full-pool greedy reference."""
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())
fixture = Path(__file__).resolve().parents[1] / 'reassignment/in/fingerprints.gfp'
records = [r + '|\n' for r in fixture.read_text().split('|\n') if r.strip()]
ids = [re.search(r'PCN<([^>]+)>', r)[1] for r in records]

def parse(text):
    groups = []
    for record in text.split('|\n'):
        names = re.findall(r'^PCN<([^>]+)>', record, re.M)
        if names:
            assert len(names) == int(re.search(r'CSIZE<(\d+)>', record)[1])
            groups.append(names)
    return groups

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    def run(path, *args):
        r = subprocess.run([exe, '-h', '4', *args, str(path)], text=True,
                           capture_output=True, timeout=30)
        assert r.returncode == 0, r.stderr
        return r.stdout

    # Gather pair distances independently of candidate filtering: each query
    # forms one sphere containing the entire small fixture.
    matrix = {}
    for i in range(len(records)):
        path = tmp / 'distances.gfp'
        path.write_text(records[i] + ''.join(records[:i] + records[i+1:]))
        text = run(path, '-t', '.999999')
        assert len(parse(text)) == 1
        names = re.findall(r'^PCN<([^>]+)>', text, re.M)
        distances = [0] + list(map(float, re.findall(r'^DIST<([^>]+)>', text, re.M)))
        matrix.update({(i, ids.index(name)): d for name, d in zip(names, distances)})

    interleaved = list(range(12)) * 60
    blocks = sorted(interleaved)
    shuffled = list(interleaved)
    random.Random(7391).shuffle(shuffled)
    for ordering in (interleaved, blocks, shuffled, interleaved[:501], interleaved[:502]):
        path = tmp / 'pool.gfp'
        path.write_text(''.join(re.sub(r'PCN<[^>]+>', f'PCN<m{i}>', records[k])
                                for i, k in enumerate(ordering)))
        for threshold in (.00001, .35, .7):
            for exclusion in (False, True):
                initial = {i for i, k in enumerate(ordering)
                           if exclusion and matrix[0, k] <= threshold}
                reference = []
                selected = set(initial)
                while len(selected) < len(ordering):
                    centre = next(i for i in range(len(ordering)) if i not in selected)
                    members = [i for i, k in enumerate(ordering) if i not in selected
                               and matrix[ordering[centre], k] <= threshold]
                    assert members[0] == centre
                    reference.append([f'm{i}' for i in members])
                    selected.update(members)
                previous = tmp / 'previous.gfp'
                previous.write_text(records[0])
                args = ['-A', str(previous)] if exclusion else []
                for batch in (1, 2):
                    for limit in (1, 3, 7, 1000):
                        actual = parse(run(path, '-t', str(threshold), '-p', str(batch),
                                           '-n', str(limit), *args))
                        assert actual == reference[:limit], (threshold, exclusion, batch, limit)
print('gfp_leader_standard compaction checks passed')
