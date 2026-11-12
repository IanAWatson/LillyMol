#!/usr/bin/env python3
"""Check external query loading and failures detected before randomisation."""
from pathlib import Path
import os
import subprocess
import sys

exe = sys.argv[1]


def run(name, options, smiles, success=True, env=None):
    source = Path(name + '.smi')
    source.write_text(smiles + '\n')
    result = subprocess.run(
        [exe, '-v', '-n', '10', '-F', name + '.fail.smi', *options, str(source)],
        capture_output=True, text=True, timeout=30, env=env)
    assert result.returncode >= 0, result.stderr
    assert (result.returncode == 0) == success, result.stderr
    return result


query_dir = Path('external_queries')
query_dir.mkdir(exist_ok=True)
query = query_dir / 'amide.textproto'
query.write_text('# Query file starts with comments\n\nquery { smarts: "[ND1H1]=[CD3]-[OD1H]" }\n')
external = query_dir / 'standardise.textproto'
external.write_text('query_file: "amide.textproto" smiles: "N-C=O" name: "amide"\n')
result = run('relative_query', ['-g', 'EXT:' + str(external)], 'CC(=[NH])[OH]')
assert '1 molecules changed' in result.stderr, result.stderr
assert '0 molecules with failures' in result.stderr, result.stderr
assert Path('relative_query.fail.smi').read_text() == '', result.stderr

# Legacy MSI queries use the same initial SMARTS atom ordering.
legacy = query_dir / 'amide.qry'
legacy.write_text('# Legacy query\n\n(0 Query\n (A C smarts "[ND1H1]=[CD3]-[OD1H]")\n)\n')
external.write_text('query_file: "amide.qry" smiles: "N-C=O"\n')
result = run('legacy_query', ['-g', 'EXT:' + str(external)], 'CC(=N)O')
assert '1 molecules changed' in result.stderr, result.stderr
assert Path('legacy_query.fail.smi').read_text() == '', result.stderr

# Expansion occurs before resolving a relative filename against the proto file.
env = dict(os.environ, STANDARDISE_QUERY_FILE='amide.textproto')
external.write_text('query_file: "${STANDARDISE_QUERY_FILE}" smiles: "N-C=O"\n')
result = run('expanded_query', ['-g', 'EXT:' + str(external)], 'CC(=N)O', env=env)
assert '1 molecules changed' in result.stderr, result.stderr
assert Path('expanded_query.fail.smi').read_text() == '', result.stderr

# A query can exclude a template atom from its embeddings. Such a mapping
# must be rejected before any changes, even though the query has both atoms.
query.write_text('query { query_atom { id: 0 atom_properties { atomic_number: 6 } } '
                 'query_atom { id: 1 atom_properties { atomic_number: 8 } '
                 'single_bond: 0 include_in_embedding: false } }\n')
external.write_text('query_file: "amide.textproto" smiles: "C=O"\n')
result = run('omitted_atom', ['-g', 'EXT:' + str(external)], 'CO')
assert '0 molecules changed' in result.stderr, result.stderr
assert ('invalid template atom mapping' in result.stderr or
        'embedding shorter than template' in result.stderr), result.stderr

external.write_text('query_file: "missing.textproto" smiles: "N-C=O"\n')
run('missing_query', ['-g', 'EXT:' + str(external)], 'CC(=N)O', success=False)
external.write_text('smarts: "C" smiles: "CC"\n')
run('short_query', ['-g', 'EXT:' + str(external)], 'CC', success=False)

# The replacement fragment is valid, but imposing its triple bond on this
# branched target produces an invalid carbon valence. The tester must report it.
external.write_text('smarts: "[CD1]-[CD3]" smiles: "C#C"\n')
result = run('introduced_valence', ['-g', 'EXT:' + str(external)], 'CC(C)C input')
assert '1 molecules with failures' in result.stderr, result.stderr
assert '1 valence errors introduced' in result.stderr, result.stderr
assert 'input valence' in Path('introduced_valence.fail.smi').read_text(), result.stderr
run('abort_on_valence', ['-f', '-g', 'EXT:' + str(external)], 'CC(C)C input', success=False)

result = run('invalid_start', ['-g', 'amidine'], '[CH5] invalid')
assert 'Abnormal starting valence' in result.stderr, result.stderr
assert '0 valence errors introduced' in result.stderr, result.stderr
assert Path('invalid_start.fail.smi').read_text() == '', result.stderr

# Invalid inputs may still be repaired; measuring their initial valence must
# not prevent the standardiser from running.
external.write_text('smarts: "[#6]" smiles: "[CH3-]"\n')
result = run('repaired_start', ['-g', 'EXT:' + str(external)], '[CH5] repair')
assert '1 molecules changed' in result.stderr, result.stderr
assert 'Abnormal starting valence' not in result.stderr, result.stderr
assert '0 valence errors introduced' in result.stderr, result.stderr
assert Path('repaired_start.fail.smi').read_text() == '', result.stderr

print('Standardisation correctness tests passed')

# Opposite hydrazone preferences must be resolved before reading molecules.
for name, options in [
    ('hydrazone_conflict_all', ['-g', 'all', '-g', 'beta-enehydrazine']),
    ('hydrazone_conflict_explicit', ['-g', 'beta-enehydrazine', '-g', 'beta-hydrazone']),
]:
    result = run(name, options, 'CC(=O)CC(C)=NNC', success=False)
    assert 'mutually exclusive' in result.stderr, result.stderr
    assert '-g all -g -beta-hydrazone -g beta-enehydrazine' in result.stderr, result.stderr
