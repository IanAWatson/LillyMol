#!/usr/bin/env python3
"""Front-end for LillyMol ring_replacement.

This wrapper validates collection/ring-system names and assembles a
ring_replacement command line with the appropriate -R arguments.

Wrapper options accept both hyphenated and underscored spellings, for example:
  --collections-dir and --collections_dir
  --list-rings and --list_rings

Options not known to this wrapper are passed through unchanged to the underlying
ring_replacement executable.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
from pathlib import Path
from typing import Sequence


DEFAULT_RING_REPLACEMENT = "ring_replacement"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Front-end for LillyMol ring_replacement.",
        allow_abbrev=False,
    )

    parser.add_argument(
        "--collections-dir",
        "--collections_dir",
        required=True,
        type=Path,
        dest="collections_dir",
        help="Directory containing collection subdirectories.",
    )

    parser.add_argument(
        "--collection",
        help=(
            "Collection name. Must correspond to a directory under "
            "--collections-dir. Required unless --list-collections is used."
        ),
    )

    parser.add_argument(
        "--ring",
        action="append",
        dest="rings",
        help=(
            "Ring system identifier, e.g. 3a4A, 6a, 6A, 5a6a, or 6a6a6a. "
            "Can be specified multiple times."
        ),
    )

    list_group = parser.add_mutually_exclusive_group()

    list_group.add_argument(
        "--list-collections",
        "--list_collections",
        action="store_true",
        dest="list_collections",
        help="List available collections and exit.",
    )

    list_group.add_argument(
        "--list-rings",
        "--list_rings",
        action="store_true",
        dest="list_rings",
        help="List available ring systems for --collection and exit.",
    )

    parser.add_argument(
        "--ring-replacement",
        "--ring_replacement",
        default=DEFAULT_RING_REPLACEMENT,
        dest="ring_replacement",
        help=(
            "Path to the ring_replacement executable. "
            f"Defaults to {DEFAULT_RING_REPLACEMENT!r}."
        ),
    )

    parser.add_argument(
        "--dry-run",
        "--dry_run",
        action="store_true",
        dest="dry_run",
        help="Print the ring_replacement command but do not run it.",
    )

    parser.add_argument(
        "--inexact",
        action="store_true",
        dest="inexact",
        help="Allow substituents at any open position in a replacement ring. More structures, higher risk",
    )

    return parser


def validate_collections_dir(collections_dir: Path) -> None:
    if not collections_dir.is_dir():
        raise FileNotFoundError(
            f"--collections-dir is not a directory: {collections_dir}"
        )


def validate_collection(collections_dir: Path, collection: str) -> Path:
    collection_dir = collections_dir / collection

    if not collection_dir.is_dir():
        raise FileNotFoundError(
            f"Collection {collection!r} is invalid: "
            f"directory does not exist: {collection_dir}\n"
            "Use --list-collections to see valid collections."
        )

    return collection_dir


def ring_replacement_directory(collection_dir: Path, collection: str) -> Path:
    rr_dir = collection_dir / f"ring_replacement_{collection}"

    if not rr_dir.is_dir():
        raise FileNotFoundError(
            f"Collection {collection!r} is invalid: "
            f"ring replacement directory does not exist: {rr_dir}"
        )

    return rr_dir


def list_collections(collections_dir: Path) -> list[str]:
    """Return names of directories that look like valid collections."""
    collections: list[str] = []

    for path in collections_dir.iterdir():
        if not path.is_dir():
            continue

        collection = path.name
        rr_dir = path / f"ring_replacement_{collection}"

        # Only show directories that look like valid collections.
        if rr_dir.is_dir():
            collections.append(collection)

    return sorted(collections)


def list_ring_systems(ring_dir: Path, collection: str) -> list[str]:
    """Return ring-system names inferred from collection_ring.smi files."""
    prefix = f"{collection}_"
    suffix = ".smi"

    ring_systems: list[str] = []

    for path in ring_dir.iterdir():
        if not path.is_file():
            continue

        name = path.name

        if not name.startswith(prefix):
            continue

        if not name.endswith(suffix):
            continue

        ring_system = name[len(prefix) : -len(suffix)]
        if ring_system:
            ring_systems.append(ring_system)

    return sorted(ring_systems)


def validate_ring_files(
    ring_dir: Path,
    collection: str,
    rings: Sequence[str],
) -> list[Path]:
    ring_files: list[Path] = []

    for ring_system in rings:
        ring_file = ring_dir / f"{collection}_{ring_system}.smi"

        if not ring_file.is_file():
            raise FileNotFoundError(
                f"Ring system {ring_system!r} is invalid for collection {collection!r}: "
                f"file does not exist: {ring_file}\n"
                f"Use --collection {collection} --list-rings to see valid ring systems."
            )

        ring_files.append(ring_file)

    return ring_files


def build_command(
    ring_replacement: str,
    ring_files: Sequence[Path],
    passthrough_args: Sequence[str],
) -> list[str]:
    command = [ring_replacement]

    for ring_file in ring_files:
        command.extend(["-R", str(ring_file)])

    command.extend(passthrough_args)
    return command


def print_items(items: Sequence[str]) -> None:
    for item in items:
        print(item)


def main() -> int:
    parser = build_parser()
    args, passthrough_args = parser.parse_known_args()

    try:
        validate_collections_dir(args.collections_dir)

        if args.list_collections:
            print_items(list_collections(args.collections_dir))
            return 0

        if not args.collection:
            parser.error("--collection is required unless --list-collections is used")

        collection_dir = validate_collection(args.collections_dir, args.collection)
        ring_dir = ring_replacement_directory(collection_dir, args.collection)

        if args.list_rings:
            print_items(list_ring_systems(ring_dir, args.collection))
            return 0

        if not args.rings:
            parser.error(
                "At least one --ring option is required unless --list-rings is used"
            )

        ring_files = validate_ring_files(ring_dir, args.collection, args.rings)

    except FileNotFoundError as err:
        parser.error(str(err))

    exe = args.ring_replacement
    if args.inexact:
      exe = 'ring_replacement_inexact'
    command = build_command(
        ring_replacement=exe,
        ring_files=ring_files,
        passthrough_args=passthrough_args,
    )

    print(shlex.join(command))

    if args.dry_run:
        return 0

    completed = subprocess.run(command, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
