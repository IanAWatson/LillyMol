#!/usr/bin/env python3
"""Render dicer_fragment_activity JSON; no molecule toolkit or protobuf runtime needed."""
import argparse
import json
from pathlib import Path
import sys

from activity_report import SORT_FIELDS, build_html, ranked_rows


def nonnegative(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError('must be non-negative')
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='Activity report JSON, or - for stdin')
    parser.add_argument('-o', '--output', required=True, help='HTML output, or - for stdout')
    parser.add_argument('--top_n', type=nonnegative, default=100, help='Maximum HTML rows (default 100; 0 keeps all)')
    parser.add_argument('--sort', choices=SORT_FIELDS, default='abs_median', help='Rank by effect or support (default abs_median)')
    parser.add_argument('--ascending', action='store_true', help='Sort low to high, e.g. most negative median differences first')
    parser.add_argument('--min_matches', type=nonnegative, default=1, help='Minimum matched molecules (default 1)')
    parser.add_argument('--min_unmatched', type=nonnegative, default=1, help='Minimum unmatched molecules (default 1)')
    parser.add_argument('--title', default='Dicer Fragment Activity Exploration')
    args = parser.parse_args(argv)
    if args.input != '-' and args.output != '-':
        source, destination = Path(args.input), Path(args.output)
        if (source.resolve() == destination.resolve() or
                (source.exists() and destination.exists() and source.samefile(destination))):
            parser.error('HTML output would overwrite the JSON input')
    try:
        if args.input == '-':
            data = json.load(sys.stdin)
        else:
            data = json.loads(Path(args.input).read_text())
        rows = ranked_rows(data, args.sort, args.ascending, args.top_n,
                           args.min_matches, args.min_unmatched)
        document = build_html(data, args.title, rows=rows)
        if args.output == '-':
            sys.stdout.write(document)
        else:
            Path(args.output).write_text(document)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == '__main__':
    sys.exit(main())
