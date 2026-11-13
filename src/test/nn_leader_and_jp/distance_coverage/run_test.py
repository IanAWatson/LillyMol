#!/usr/bin/env python3
"""Coverage warnings are advisory, per-file, and based on unpruned input."""
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tfrecord_utils import tfrecord

exe = str(Path(sys.argv[1]).resolve())
prefix = "Warning: distance coverage in '"
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    for serialized in (False, True):
        def write(name, records):
            path = root / name
            if serialized:
                path.write_bytes(b''.join(tfrecord(owner, nbrs) for owner, nbrs in records))
            else:
                path.write_text(''.join('$SMI<C>\nPCN<' + owner + '>\n' +
                    ''.join(f'PCN<{nbr}>\nDIST<{dist}>\n' for nbr, dist in nbrs) +
                    '|\n' for owner, nbrs in records))
            return path

        def run(args, files):
            result = subprocess.run([exe, *(['-r'] if serialized else []), *args,
                                     *map(str, files)], capture_output=True, text=True,
                                    timeout=15)
            assert result.returncode == 0, result.stderr
            assert all(line.startswith(prefix) for line in result.stderr.splitlines()), result.stderr
            return result

        first = write('first', [('a', [('b', 0.25)]), ('b', [('a', 0.25)])])
        quiet = run(['-t', '0.25'], [first])
        assert not quiet.stderr
        warned = run(['-t', '0.5'], [first])
        assert warned.stderr.count(prefix) == 1
        assert 'thresholds 0.5 exceed maximum stored distance 0.25' in warned.stderr
        second = write('second', [('c', [('d', 0.75)]), ('d', [('c', 0.75)])])
        result = run(['-s', '4', '-t', '0.5'], [first, second])
        assert result.stderr.count(prefix) == 1 and str(first) in result.stderr
        result = run(['-t', '0.25,0.5,1', '-B', str(root / 'clusters')], [first])
        assert 'thresholds 0.5 1 exceed' in result.stderr
        repeated = run(['-t', '0.25', '-t', '0.5', '-t', '1', '-B', str(root / 'repeated')], [first])
        assert repeated.stderr == result.stderr
        # The long edge is discarded by -s, but still counts towards observed coverage.
        assert not run(['-s', '1', '-t', '0.25'], [first]).stderr
        # Do not inspect an unread suffix or later file.
        tail = write('tail', [('a', []), ('b', [('a', 0.75)])])
        result = run(['-s', '1', '-t', '0.5'], [tail, root / 'missing'])
        assert result.stderr.count(prefix) == 1 and 'no stored neighbour distances' in result.stderr
        zero = write('zero', [('a', [('a', 0.0)])])
        assert not run(['-h', '-t', '0'], [zero]).stderr
        assert 'maximum stored distance 0' in run(['-h', '-t', '0.5'], [zero]).stderr
        assert not run(['-J', '.,1'], [tail]).stderr
print('nn_leader_and_jp distance coverage checks passed')
