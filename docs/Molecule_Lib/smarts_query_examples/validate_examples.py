#!/usr/bin/env python3
"""Validate the documented examples using LillyMol and, optionally, RDKit."""
import argparse
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import tempfile

from generate_smarts import patterns


def check_names(actual, expected, caption):
    expected = set(expected)
    if actual != expected:
        raise ValueError(f'{caption}: missing {sorted(expected - actual)}, '
                         f'unexpected {sorted(actual - expected)}')


def main():
    root = Path(__file__).resolve().parent
    repo = Path(os.environ.get('LILLYMOL_HOME', root.parents[2]))
    local_exe = repo/'bin'/platform.system()/'tsubstructure'
    default_exe = str(local_exe) if local_exe.is_file() else shutil.which('tsubstructure')
    if default_exe is None:
        default_exe = 'tsubstructure'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tsubstructure', default=default_exe)
    parser.add_argument('--rdkit', action='store_true', help='Also validate the standard SMARTS with RDKit')
    args = parser.parse_args()
    # These calls deliberately use the CLI, not LillyMol Python bindings.
    exe = str(Path(shutil.which(args.tsubstructure) or args.tsubstructure).resolve())
    print(f'LillyMol executable: {exe}')
    if args.rdkit:
        from rdkit import Chem, __version__
        print(f'RDKit {__version__}')
    cases = json.loads((root/'cases.json').read_text())
    generated = patterns()
    checks = 0
    with tempfile.TemporaryDirectory(prefix='smarts query examples ') as tmp:
        for case in cases:
            name = case['name']
            source = root/case['input']
            command = [exe, '-q', 'PROTO:'+str(root/(name+'.textproto')), '-m', '-']
            if 'labels' in case:
                command.extend(['-j', '1', '-j', 'same'])
            result = subprocess.run([*command, str(source)], cwd=tmp,
                                    capture_output=True, text=True, timeout=60)
            if result.returncode:
                raise ValueError(f'{name}: query failed: {result.stderr}')
            records = [line.split(maxsplit=1) for line in result.stdout.splitlines()]
            check_names({fields[1] for fields in records}, case['expected'], name+' proto')
            for smiles, identifier in records:
                if 'labels' in case:
                    labels = len(re.findall(r'\[1(?=[A-Za-z*])', smiles))
                    if labels != case['labels'][identifier]:
                        raise ValueError(f'{name}/{identifier}: expected {case["labels"][identifier]} '
                                         f'labelled atoms, got {labels}: {smiles}')
            checks += 1
            if case.get('standard'):
                path = root/case['standard']
                smarts = path.read_text().split()[0]
                if smarts != generated[path.stem]:
                    raise ValueError(f'{path}: differs from generate_smarts.py')
                result = subprocess.run([exe, '-q', 'S:'+str(path), '-m', '-', str(source)],
                                        cwd=tmp, capture_output=True, text=True, timeout=60)
                if result.returncode:
                    raise ValueError(f'{name}: standard SMARTS failed: {result.stderr}')
                check_names({line.split(maxsplit=1)[1] for line in result.stdout.splitlines()},
                            case['expected'], name+' standard/LillyMol')
                checks += 1
                if args.rdkit:
                    query = Chem.MolFromSmarts(smarts)
                    if query is None:
                        raise ValueError(f'{name}: RDKit rejected SMARTS')
                    matched = set()
                    for line in source.read_text().splitlines():
                        smiles, identifier = line.split(maxsplit=1)
                        molecule = Chem.MolFromSmiles(smiles)
                        if molecule is None:
                            raise ValueError(f'{name}/{identifier}: RDKit rejected SMILES')
                        if molecule.HasSubstructMatch(query):
                            matched.add(identifier)
                    check_names(matched, case['expected'], name+' standard/RDKit')
                    checks += 1
            print(f'{name}: PASS ({len(case["expected"])} matching molecules)')
    print(f'{checks} checks passed; region and matched-atom isotope counts also verified.')


if __name__ == '__main__':
    main()
