#!/usr/bin/env python3
"""Build a size-limited fragment distribution from dicer_to_topological_types output."""

import argparse
import io
from pathlib import Path
import re
import tarfile
import tempfile


FILENAME = re.compile(r'DICER_([12])_([0-9]+)\.textproto')
NAT = re.compile(rb'\bnat:\s*(\d+)\b')
COUNT = re.compile(rb'\bn:\s*(\d+)\b')


def build(source: Path, output: Path, max_atoms: int = 11):
    if max_atoms < 1:
        raise ValueError('Maximum atom count must be positive')
    inputs = []
    for path in source.glob('DICER_*.textproto'):
        match = FILENAME.fullmatch(path.name)
        if match:
            inputs.append((int(match[1]), int(match[2]), path))
    if not inputs:
        raise ValueError(f'No DICER_1_* or DICER_2_* textprotos in {source}')
    summary = []
    # Replace the archive only after every selected input has been processed.
    with tempfile.TemporaryDirectory(prefix='.dicer-distribution-', dir=output.parent) as tmp:
        archive = Path(tmp) / output.name
        with tarfile.open(archive, 'w:xz', preset=6) as tar:
            for attachments, number, path in sorted(inputs):
                if attachments == 1 and number > max_atoms:
                    continue
                selected = []
                with path.open('rb') as inp:
                    for line_number, line in enumerate(inp, 1):
                        nat, count = NAT.search(line), COUNT.search(line)
                        if not nat or not count:
                            raise ValueError(f'{path}:{line_number}: missing nat or n')
                        size = int(nat[1])
                        if attachments == 1 and size != number:
                            raise ValueError(f'{path}:{line_number}: nat disagrees with filename')
                        if size <= max_atoms:
                            selected.append((int(count[1]), line.rstrip(b'\r\n') + b'\n'))
                if not selected:
                    continue
                selected.sort(key=lambda item: (-item[0], item[1]))
                data = b''.join(line for _, line in selected)
                name = f'FRAG_{attachments}_{number}.textproto'
                member = tarfile.TarInfo(name)
                member.size = len(data)
                member.mode = 0o644
                tar.addfile(member, io.BytesIO(data))
                summary.append((name, len(selected), len(data)))
        if not summary:
            raise ValueError('No fragments satisfy the atom-count limit')
        archive.replace(output)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parent / 'dicer_fragments.tar.xz')
    parser.add_argument('--max_atoms', type=int, default=11)
    args = parser.parse_args()
    try:
        summary = build(args.source, args.output, args.max_atoms)
    except (OSError, ValueError, tarfile.TarError) as exc:
        parser.exit(1, f'Cannot build fragment distribution: {exc}\n')
    for name, count, size in summary:
        print(f'{name}: {count} fragments, {size} bytes')
    print(f'Total: {sum(count for _, count, _ in summary)} fragments in '
          f'{len(summary)} files; archive {args.output.stat().st_size} bytes')


if __name__ == '__main__':
    main()
