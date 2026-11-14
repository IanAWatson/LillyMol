#!/usr/bin/env python3
"""Paired leader scans must produce exactly the ordinary leader output."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile

exe = str(Path(sys.argv[1]).resolve())
fixture = Path(__file__).resolve().parents[1] / 'reassignment/in/fingerprints.gfp'
records = [r + '|\n' for r in fixture.read_text().split('|\n') if r.strip()]
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    path = tmp / 'pool.gfp'
    # More than 500 candidates exercises both parallel scans. Copies also test
    # fallback when the next two prospective leaders belong to the same sphere.
    path.write_text(''.join(re.sub(r'PCN<([^>]+)>', rf'PCN<\1_{copy}>', r)
                            for copy in range(51) for r in records))
    previous = tmp / 'previous.gfp'
    previous.write_text(records[0])
    def run(*args):
        result = subprocess.run([exe, *args, str(path)], capture_output=True,
                                text=True, timeout=30)
        assert result.returncode == 0, result.stderr
        return result.stdout
    for threshold in ('0.00001', '0.35', '0.7'):
        for extra in ([], ['-n', '1'], ['-n', '2'], ['-n', '3'], ['-n', '7'],
                      ['-A', str(previous)], ['-R', '2'], ['-R', '2:any'], ['-k', '.']):
            reference = run('-t', threshold, '-h', '1', *extra)
            for threads in ('1', '4'):
                for batch in ('1', '2'):
                    actual = run('-t', threshold, '-h', threads, '-p', batch, *extra)
                    assert actual == reference, (threshold, extra, threads, batch)
print('gfp_leader_standard paired leader checks passed')
