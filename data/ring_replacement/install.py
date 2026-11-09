#!/usr/bin/env python3
"""Expand the distributed ring collection using the host's filename convention."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import tempfile


SAFE_NAME = re.compile(r'rings_(?:[0-9]+(?:Al|Ar)?)+\.smi')
GENERATED_NAME = re.compile(r'rings_(?:[0-9]+(?:Al|Ar|A|a)?)+\.smi')


def install(directory: Path, system: str):
    if system not in ('Linux', 'Darwin'):
        raise ValueError(f'Unsupported platform {system!r}; specify Linux or Darwin')
    manifest = json.loads((directory / 'manifest.json').read_text())
    if not manifest:
        raise ValueError('Empty ring collection manifest')
    outputs = set()
    # Validate and stage the complete collection before replacing any installed
    # files. Temporary files are on the same filesystem for atomic renames.
    with tempfile.TemporaryDirectory(prefix='.ring-install-', dir=directory) as tmp:
        staging = Path(tmp)
        for name, expected in manifest.items():
            if not SAFE_NAME.fullmatch(name):
                raise ValueError(f'Invalid ring filename {name!r}')
            output = name if system == 'Darwin' else name.replace('Al', 'A').replace('Ar', 'a')
            if output.casefold() in {n.casefold() for n in outputs} and system == 'Darwin':
                raise ValueError(f'Case-insensitive filename collision: {output}')
            outputs.add(output)
            data = gzip.decompress((directory / (name + '.gz')).read_bytes())
            if (len(data) != expected['bytes'] or
                    hashlib.sha256(data).hexdigest() != expected['sha256'] or
                    len(data.splitlines()) != expected['records']):
                raise ValueError(f'Ring collection checksum or size mismatch: {name}')
            (staging / output).write_bytes(data)
        for name in outputs:
            os.replace(staging / name, directory / name)
        # Remove stale generated rings, including the other platform's naming
        # scheme, without touching the compressed sources or unrelated files.
        for path in directory.glob('rings_*.smi'):
            if GENERATED_NAME.fullmatch(path.name) and path.name not in outputs:
                path.unlink()
    return len(outputs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--platform', choices=('Linux', 'Darwin'), default=platform.system())
    args = parser.parse_args()
    try:
        count = install(args.directory, args.platform)
    except (OSError, ValueError, EOFError) as exc:
        parser.exit(1, f'Cannot install ring replacements: {exc}\n')
    print(f'Installed {count} ring replacement files for {args.platform}')


if __name__ == '__main__':
    main()
