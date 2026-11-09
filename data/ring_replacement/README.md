# Ring replacement collection

This directory distributes replacement rings extracted from ChEMBL. Install the
collection before using it:

```shell
python3 "$LILLYMOL_HOME/data/ring_replacement/install.py"
```

The normal build invokes this through `src/uncompress_and_install.sh`.
Installation requires Python 3 and its standard library only. The `to_linux.sh`
and `to_mac.sh` scripts also install the collection and can be called from any
working directory.

## Portable distribution and installed names

Only gzip-compressed, case-safe source files are tracked. For example,
`rings_5Al.smi.gz` contains aliphatic five-membered rings and
`rings_5Ar.smi.gz` contains aromatic five-membered rings. Hidden filenames would
not avoid case collisions on case-insensitive filesystems.

Installation retains the compressed sources and expands them to:

| Source | Linux | Mac |
| --- | --- | --- |
| `rings_5Al.smi.gz` | `rings_5A.smi` | `rings_5Al.smi` |
| `rings_5Ar.smi.gz` | `rings_5a.smi` | `rings_5Ar.smi` |
| `rings_56.smi.gz` | `rings_56.smi` | `rings_56.smi` |

Names without aromaticity qualifiers contain aromaticity-independent queries.
The naming conversion applies to each ring in a fused system, for example
`rings_5Al6Ar.smi.gz` becomes `rings_5A6a.smi` on Linux.

`manifest.json` records the uncompressed byte count, SHA-256 checksum and record
count for every source. The installer validates and stages the entire collection
before replacing any outputs, then replaces individual files with atomic
renames. It refreshes existing outputs, restores missing files, and removes stale
`rings_*.smi` files from previous collections or naming schemes. Keep custom ring
collections elsewhere. Generated files are ignored by Git; do not commit them.

For testing or installing into a separate directory containing the compressed
sources and manifest:

```shell
python3 install.py --directory /path/to/collection --platform Darwin
```

`common_rings.txt` is an informational subset list with current record counts.
Its filenames use the Linux convention; substitute `Al` for `A` and `Ar` for
`a` when using the corresponding installed files on Mac.

## Collection provenance and format

The current 183 files were supplied from
`/d1/ian/CHEMBL/ring_replacement_chembl` on 2026-11-09, using Chembl 37, replacing the previous
ChEMBL 36 collection. The supplied filenames were renamed from `chembl_` to
`rings_`; record contents are unchanged.

The collection includes 132 aromaticity-specific files and 51
aromaticity-independent files. Uncompressed size is 34,498,161 bytes; gzip sources
occupy about 4.6 MB.

Although installed filenames end in `.smi`, each line is a replacement-ring
textproto. Most LillyMol tools can read these with `-i textproto`.
`ring_replacement` also supports `-D` and `-F` for subsetting replacement rings.
