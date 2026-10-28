#!/usr/bin/env python3
"""Explore activity distributions associated with substructure motifs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from activity_report import SORT_FIELDS, build_html, ranked_rows

from google.protobuf import text_format

from substructure_exploration_pb2 import SubstructureExploration

from lillymol import Molecule, Standardise
from lillymol_io import MolReaderContext
from lillymol_query import SubstructureQuery


@dataclass
class MoleculeRecord:
    smiles: str
    identifier: str
    activity: float
    mol: Molecule


@dataclass
class Example:
    smiles: str
    labelled_smiles: str
    identifier: str
    activity: float
    nhits: int


@dataclass
class QuerySpec:
    name: str
    source: str
    source_type: str
    query: SubstructureQuery
    ncon: list[int] = field(default_factory=list)
    charge_assigner: bool = False


@dataclass
class QueryAccumulator:
    spec: QuerySpec
    matched: list[float] = field(default_factory=list)
    unmatched: list[float] = field(default_factory=list)
    count_values: dict[str, list[float]] = field(default_factory=dict)
    examples: list[Example] = field(default_factory=list)
    molecules_with_multiple_matches: int = 0
    max_matches: int = 0


def die(message: str) -> None:
    print(message, file=sys.stderr)
    sys.exit(1)


def separator_for_file(path: str) -> str | None:
    if path.lower().endswith('.csv'):
        return ','
    return None


def read_config(path: str) -> SubstructureExploration:
    message = SubstructureExploration()
    try:
        text = Path(path).read_text()
    except OSError as exc:
        die(f"Cannot read config '{path}': {exc}")

    try:
        text_format.Parse(text, message)
    except text_format.ParseError as exc:
        die(f"Cannot parse config '{path}': {exc}")

    return message


def build_queries(config: SubstructureExploration) -> list[QuerySpec]:
    result = []
    for ndx, proto in enumerate(config.query, 1):
        qry = SubstructureQuery()
        source_type = proto.WhichOneof('query_spec')
        if source_type == 'smarts':
            source = proto.smarts
            if not qry.build_from_smarts(source):
                die(f"Query {ndx}: invalid SMARTS '{source}'")
        elif source_type == 'query':
            source = text_format.MessageToString(proto.query, as_one_line=True)
            if not qry.construct_from_proto(proto.query):
                die(f"Query {ndx}: cannot construct inline SubstructureQuery")
        elif source_type == 'fname':
            source = proto.fname
            if source.endswith('.qry'):
                ok = qry.read_msi(source)
            else:
                ok = qry.read_proto(source) or qry.read_msi(source)
            if not ok:
                die(f"Query {ndx}: cannot read query file '{source}'")
        else:
            die(f"Query {ndx}: no smarts, fname or inline query specified")

        if not proto.HasField('name') or not proto.name:
            die(f"Query {ndx}: every query must have a name")

        qry.set_find_unique_embeddings_only(True)
        result.append(QuerySpec(
            name=proto.name,
            source=source,
            source_type=source_type,
            query=qry,
            ncon=list(proto.ncon),
            charge_assigner=proto.charge_assigner if proto.HasField('charge_assigner') else False,
        ))

    names = [q.name for q in result]
    duplicate_names = sorted(name for name in set(names) if names.count(name) > 1)
    if duplicate_names:
        die(f"Duplicate query names: {', '.join(duplicate_names)}")

    return result


def read_activity(path: str) -> dict[str, float]:
    result: dict[str, float] = {}
    sep = separator_for_file(path)
    try:
        with open(path, newline='') as reader:
            if sep == ',':
                records = csv.reader(reader)
                header = next(records, None)
                if header is None:
                    die(f"{path}: empty activity file")
                for line_number, row in enumerate(records, 2):
                    if not row:
                        continue
                    if len(row) < 2:
                        die(f"{path}:{line_number}: activity records need at least two columns")
                    add_activity_record(result, row[0], row[1], path, line_number)
            else:
                header = reader.readline()
                if not header:
                    die(f"{path}: empty activity file")
                for line_number, line in enumerate(reader, 2):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    fields = line.split()
                    if len(fields) < 2:
                        die(f"{path}:{line_number}: activity records need at least two columns")
                    add_activity_record(result, fields[0], fields[1], path, line_number)
    except OSError as exc:
        die(f"Cannot read activity file '{path}': {exc}")

    if not result:
        die(f"{path}: no activity values read")

    return result


def add_activity_record(activity: dict[str, float], identifier: str, value: str,
                        fname: str, line_number: int) -> None:
    identifier = identifier.strip()
    if not identifier:
        die(f"{fname}:{line_number}: empty identifier")
    if identifier in activity:
        die(f"{fname}:{line_number}: duplicate activity for '{identifier}'")
    try:
        activity[identifier] = float(value)
    except ValueError:
        die(f"{fname}:{line_number}: invalid activity '{value}'")


def read_molecules(path: str, activity: dict[str, float], okmissing_activity: bool, standardize: bool = True) -> tuple[list[MoleculeRecord], int]:
    result: list[MoleculeRecord] = []
    missing_activity = 0
    seen: set[str] = set()

    standardiser = Standardise()
    if standardize:
        standardiser.activate_all()

    with MolReaderContext(path) as reader:
        for mol in reader:
            identifier = str(mol.name()).split()[0]
            if not identifier:
                die(f"{path}: molecule with no identifier")
            if identifier in seen:
                die(f"{path}: duplicate molecule identifier '{identifier}'")
            seen.add(identifier)

            if identifier not in activity:
                if not okmissing_activity:
                    die(f"{path}: missing activity for '{identifier}'")
                missing_activity += 1
                continue

            if standardize:
                standardiser.process(mol)

            result.append(MoleculeRecord(
                smiles=str(mol.smiles()),
                identifier=identifier,
                activity=activity[identifier],
                mol=mol,
            ))

    if not result:
        die(f"{path}: no molecules with activity")

    unused_activity = sorted(set(activity) - seen)
    if unused_activity:
        print(f"Warning: {len(unused_activity)} activity records not found in smiles file", file=sys.stderr)

    return result, missing_activity


def percentile(sorted_values: list[float], pct: float) -> float | None:
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * pct / 100.0
    lower = math.floor(pos)
    upper = math.ceil(pos)
    if lower == upper:
        return sorted_values[int(pos)]
    fraction = pos - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def summary(values: Iterable[float]) -> dict[str, float | int | None]:
    vals = sorted(values)
    if not vals:
        return {
            'n': 0,
            'min': None,
            'q1': None,
            'median': None,
            'mean': None,
            'q3': None,
            'max': None,
        }

    return {
        'n': len(vals),
        'min': vals[0],
        'q1': percentile(vals, 25),
        'median': percentile(vals, 50),
        'mean': statistics.fmean(vals),
        'q3': percentile(vals, 75),
        'max': vals[-1],
    }


def common_bins(values: list[float], nbins: int) -> dict[str, object]:
    minval = min(values)
    maxval = max(values)
    if minval == maxval:
        edges = [minval, maxval]
        return {'min': minval, 'max': maxval, 'edges': edges}

    nbins = max(1, nbins)
    width = (maxval - minval) / nbins
    edges = [minval + i * width for i in range(nbins)] + [maxval]
    return {'min': minval, 'max': maxval, 'edges': edges}


def histogram(values: Iterable[float], edges: list[float]) -> list[int]:
    counts = [0] * (len(edges) - 1)
    if not counts:
        return counts

    minval = edges[0]
    maxval = edges[-1]
    if minval == maxval:
        counts[0] = sum(1 for _ in values)
        return counts

    width = (maxval - minval) / len(counts)
    for value in values:
        bucket = int((value - minval) / width)
        if bucket < 0:
            bucket = 0
        elif bucket >= len(counts):
            bucket = len(counts) - 1
        counts[bucket] += 1
    return counts


def bucket_for_count(nhits: int, buckets: list[str]) -> str:
    for bucket in buckets:
        if bucket.endswith('+'):
            try:
                threshold = int(bucket[:-1])
            except ValueError:
                continue
            if nhits >= threshold:
                return bucket
        else:
            try:
                if nhits == int(bucket):
                    return bucket
            except ValueError:
                continue
    return buckets[-1]


def parse_count_buckets(text: str) -> list[str]:
    buckets = [token.strip() for token in text.split(',') if token.strip()]
    if not buckets:
        die('--count-buckets must specify at least one bucket')
    if buckets[0] != '0':
        die('--count-buckets must include 0 as the first bucket')
    return buckets


def labelled_example(record: MoleculeRecord, matches, nhits: int) -> Example:
    mol = Molecule(record.mol)
    for match in matches or []:
        mol.set_isotopes(match, 1)
    return Example(
        smiles=record.smiles,
        labelled_smiles=str(mol.aromatic_smiles()),
        identifier=record.identifier,
        activity=record.activity,
        nhits=nhits,
    )


def analyse(args: argparse.Namespace) -> dict[str, object]:
    config = read_config(args.config)
    queries = build_queries(config)
    if not queries:
        die(f"{args.config}: no queries")

    unsupported = [q.name for q in queries if q.charge_assigner or q.ncon]
    if unsupported:
        print(
            'Warning: charge_assigner and ncon query metadata are recorded but not applied by this script: ' +
            ', '.join(unsupported),
            file=sys.stderr,
        )

    activity = read_activity(args.activity)
    molecules, missing_activity = read_molecules(args.smiles, activity, args.okmissing_activity, args.standardize)
    bins = common_bins([m.activity for m in molecules], args.bins)
    edges = bins['edges']
    count_buckets = parse_count_buckets(args.count_buckets)

    accumulators = [QueryAccumulator(q) for q in queries]
    for acc in accumulators:
        acc.count_values = {bucket: [] for bucket in count_buckets}

    for record in molecules:
        for acc in accumulators:
            matches = acc.spec.query.substructure_search_matches(record.mol)
            nhits = len(matches) if matches is not None else 0
            acc.max_matches = max(acc.max_matches, nhits)
            if nhits > 1:
                acc.molecules_with_multiple_matches += 1

            bucket = bucket_for_count(nhits, count_buckets)
            acc.count_values[bucket].append(record.activity)

            if nhits > 0:
                acc.matched.append(record.activity)
                if len(acc.examples) < args.max_examples:
                    acc.examples.append(labelled_example(record, matches, nhits))
            else:
                acc.unmatched.append(record.activity)

    results = []
    for ndx, acc in enumerate(accumulators, 1):
        if len(acc.matched) < args.min_matches:
            continue
        matched_summary = summary(acc.matched)
        unmatched_summary = summary(acc.unmatched)
        delta_mean = None
        if matched_summary['mean'] is not None and unmatched_summary['mean'] is not None:
            delta_mean = matched_summary['mean'] - unmatched_summary['mean']
        delta_median = None
        if matched_summary['median'] is not None and unmatched_summary['median'] is not None:
            delta_median = matched_summary['median'] - unmatched_summary['median']

        results.append({
            'index': ndx,
            'name': acc.spec.name,
            'source_type': acc.spec.source_type,
            'source': acc.spec.source,
            'ncon': acc.spec.ncon,
            'charge_assigner': acc.spec.charge_assigner,
            'matched': {
                'summary': matched_summary,
                'histogram': histogram(acc.matched, edges),
            },
            'unmatched': {
                'summary': unmatched_summary,
                'histogram': histogram(acc.unmatched, edges),
            },
            'delta_mean': delta_mean,
            'delta_median': delta_median,
            'multiplicity': {
                'molecules_with_multiple_matches': acc.molecules_with_multiple_matches,
                'max_matches': acc.max_matches,
                'profiles': {
                    bucket: {
                        'summary': summary(values),
                        'histogram': histogram(values, edges),
                    }
                    for bucket, values in acc.count_values.items()
                },
            },
            'examples': [example.__dict__ for example in acc.examples],
        })

    return {
        'metadata': {
            'smiles': args.smiles,
            'activity': args.activity,
            'config': args.config,
            'standardize': args.standardize,
            'molecules_read_with_activity': len(molecules),
            'molecules_missing_activity': missing_activity,
            'queries': len(queries),
            'queries_reported': len(results),
            'min_matches': args.min_matches,
            'histogram_bins': args.bins,
            'multiplicity_mode': args.multiplicity,
            'count_buckets': count_buckets,
        },
        'global': {
            'summary': summary(m.activity for m in molecules),
            'histogram': histogram((m.activity for m in molecules), edges),
            'bins': bins,
        },
        'queries': results,
    }


def write_json(data: dict[str, object], path: str | None) -> None:
    text = json.dumps(data, indent=2, sort_keys=True) + '\n'
    if path is None or path == '-':
        sys.stdout.write(text)
    else:
        Path(path).write_text(text)


def command_analyse(args: argparse.Namespace) -> int:
    data = analyse(args)
    write_json(data, args.output)
    if args.html:
        Path(args.html).write_text(build_html(data, args.title))
    return 0


def command_html(args: argparse.Namespace) -> int:
    try:
        data = json.loads(Path(args.input).read_text())
    except OSError as exc:
        die(f"Cannot read '{args.input}': {exc}")
    except json.JSONDecodeError as exc:
        die(f"Cannot parse JSON '{args.input}': {exc}")
    if args.sort is not None or args.top_n:
        rows = ranked_rows(data, args.sort or 'abs_median', args.ascending,
                           args.top_n, args.min_matches, args.min_unmatched)
    else:
        rows = [row for row in data['queries']
                if row['matched']['summary']['n'] >= args.min_matches and
                row['unmatched']['summary']['n'] >= args.min_unmatched]
    Path(args.output).write_text(build_html(data, args.title, rows=rows))
    return 0


def nonnegative_integer(value: str) -> int:
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError('must be a non-negative integer')
    return result


def add_analyse_args(parser: argparse.ArgumentParser) -> None:
    default_config = Path(__file__).resolve().parents[1] / 'data' / 'substructure_exploration.textproto'
    parser.add_argument('--config', default=str(default_config), help='SubstructureExploration textproto; default curated catalogue')
    parser.add_argument('--standardize', action=argparse.BooleanOptionalAction, default=True,
                        help='Apply all LillyMol chemical standardisations before matching (default enabled)')
    parser.add_argument('-A', '--activity', required=True, help='Activity file with header, id/activity columns')
    parser.add_argument('-o', '--output', default='-', help='Output JSON file, or - for stdout')
    parser.add_argument('--html', help='Optional HTML report written from the same analysis')
    parser.add_argument('--title', default='Substructure Activity Exploration', help='HTML report title')
    parser.add_argument('--bins', type=int, default=20, help='Number of common histogram bins')
    parser.add_argument('--min_matches', type=nonnegative_integer, default=1,
                        help='Minimum number of matched molecules for a query to appear in JSON/HTML (default 1; 0 keeps all queries)')
    parser.add_argument('--max_examples', type=int, default=3, help='Matched example molecules retained per query')
    parser.add_argument('--okmissing_activity', action='store_true', help='Allow molecules with no activity value; they are skipped')
    parser.add_argument('--multiplicity', choices=['presence', 'count'], default='presence', help='Primary multiplicity interpretation')
    parser.add_argument('--count-buckets', default='0,1,2,3+', help='Hit count buckets stored in JSON')
    parser.add_argument('smiles', help='Input SMILES file')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='command')

    analyse_parser = subparsers.add_parser('analyse', help='Analyse SMILES/activity/query data')
    add_analyse_args(analyse_parser)
    analyse_parser.set_defaults(func=command_analyse)

    html_parser = subparsers.add_parser('html', help='Convert analysis JSON to a static HTML report')
    html_parser.add_argument('input', help='Analysis JSON')
    html_parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    html_parser.add_argument('--title', default='Substructure Activity Exploration', help='HTML report title')
    html_parser.add_argument('--top_n', type=nonnegative_integer, default=0,
                             help='Limit HTML rows after ranking (default 0: all)')
    html_parser.add_argument('--sort', choices=SORT_FIELDS,
                             help='Rank by effect or support; --top_n defaults to abs_median')
    html_parser.add_argument('--ascending', action='store_true', help='Rank low to high')
    html_parser.add_argument('--min_matches', type=nonnegative_integer, default=0,
                             help='Minimum matched molecules when rendering existing JSON')
    html_parser.add_argument('--min_unmatched', type=nonnegative_integer, default=0,
                             help='Minimum unmatched molecules when rendering existing JSON')
    html_parser.set_defaults(func=command_html)

    args = parser.parse_args(argv)
    if not hasattr(args, 'func'):
        parser.print_help(sys.stderr)
        return 1
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
