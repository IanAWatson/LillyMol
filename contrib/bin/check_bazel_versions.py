#!/usr/bin/env python3
"""Compare literal bazel_dep declarations with local or online BCR metadata.

Informational only: no files are changed and no Git or Bazel commands are run.
"""

import argparse
import ast
from dataclasses import dataclass
from datetime import date
import json
from pathlib import Path
import re
import sys
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


# Bazel's relaxed SemVer order, including numeric identifier spelling as a
# tie-breaker. Build metadata does not affect ordering. See:
# https://github.com/bazelbuild/bazel/blob/master/src/main/java/com/google/devtools/build/lib/bazel/bzlmod/Version.java
VERSION = re.compile(r'([A-Za-z0-9.]+)(?:-([A-Za-z0-9.-]+))?(?:\+[A-Za-z0-9.-]+)?')
DATED_RELEASE = re.compile(r'(\d{4}-\d{2}-\d{2})(?:\.bcr\.\d+)?')
OVERRIDES = {'local_path_override', 'git_override', 'archive_override',
             'single_version_override', 'multiple_version_override'}


class CheckError(Exception):
    """An input or registry record cannot be checked reliably."""


def identifiers(value):
    result = []
    for part in value.split('.'):
        if not part:
            raise CheckError(f'Empty identifier in version {value!r}')
        numeric = part.isascii() and part.isdigit()
        if numeric and int(part) > 2**64 - 1:
            raise CheckError(f'Numeric version identifier too large: {part}')
        result.append((0, int(part), part) if numeric else (1, 0, part))
    return tuple(result)


def version_key(value):
    """Return a sort key matching Bazel, rather than Python package versions."""
    if value == '':
        # Bazel uses the empty version for non-registry overrides and sorts it last.
        return (1, (), 0, ())
    match = VERSION.fullmatch(value)
    if match is None:
        raise CheckError(f'Invalid Bazel version: {value!r}')
    release, prerelease = match.groups()
    return (0, identifiers(release), 1 if prerelease is None else 0,
            identifiers(prerelease) if prerelease is not None else ())


def is_prerelease(value):
    """Use Bazel syntax, except that ISO-date releases are not preview builds."""
    normalized = value.split('+', 1)[0]
    dated = DATED_RELEASE.fullmatch(normalized)
    if dated:
        try:
            date.fromisoformat(dated[1])
            return False
        except ValueError:
            pass
    return VERSION.fullmatch(value)[2] is not None


def literal_keyword(call, name, default=None):
    for keyword in call.keywords:
        if keyword.arg == name:
            if isinstance(keyword.value, ast.Constant) and isinstance(keyword.value.value, str):
                return keyword.value.value
            raise CheckError(f'Line {call.lineno}: {name} must be a literal string')
    return default


def read_module(path):
    """Read declarations without executing Starlark or matching commented text.

    MODULE.bazel's ordinary call syntax is also valid Python. Only top-level
    calls with literal names/versions are supported; dynamic declarations and
    includes fail explicitly instead of giving an incomplete report.
    """
    try:
        tree = ast.parse(path.read_text(), filename=str(path))
    except (OSError, SyntaxError, UnicodeError) as exc:
        raise CheckError(f'Cannot parse {path}: {exc}') from exc
    declarations = {}
    overrides = {}
    direct_calls = {id(node.value) for node in tree.body
                    if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)}
    for call in ast.walk(tree):
        if not isinstance(call, ast.Call) or not isinstance(call.func, ast.Name):
            continue
        function = call.func.id
        if function not in OVERRIDES | {'bazel_dep', 'include'}:
            continue
        if function == 'include':
            raise CheckError(f'Line {call.lineno}: include() is not supported; check included files separately')
        if id(call) not in direct_calls or call.args or any(k.arg is None for k in call.keywords):
            raise CheckError(f'Line {call.lineno}: {function} must be a top-level call with explicit keywords')
        name = literal_keyword(call, 'name' if function == 'bazel_dep' else 'module_name')
        if not name or name in {'.', '..'} or not re.fullmatch(r'[A-Za-z0-9_.-]+', name):
            raise CheckError(f'Line {call.lineno}: missing or invalid module name')
        if function == 'bazel_dep':
            if name in declarations:
                raise CheckError(f'Duplicate bazel_dep for {name}')
            declarations[name] = literal_keyword(call, 'version', '')
        else:
            overrides[name] = function
    if not declarations:
        raise CheckError(f'{path}: no bazel_dep declarations found')
    return declarations, overrides


class Registry:
    def __init__(self, directory=None, timeout=15):
        self.directory = directory
        self.timeout = timeout
        if directory is not None and not (directory / 'modules').is_dir():
            raise CheckError(f'Not a registry root (missing modules directory): {directory}')

    def read(self, name):
        if self.directory is not None:
            path = self.directory / 'modules' / name / 'metadata.json'
            try:
                text = path.read_text()
            except FileNotFoundError as exc:
                raise CheckError('Module metadata not found in registry clone') from exc
            except (OSError, UnicodeError) as exc:
                raise CheckError(str(exc)) from exc
        else:
            url = f'https://bcr.bazel.build/modules/{quote(name, safe="")}/metadata.json'
            try:
                request = Request(url, headers={'User-Agent': 'LillyMol-check-bazel-versions'})
                with urlopen(request, timeout=self.timeout) as response:
                    text = response.read().decode('utf-8')
            except HTTPError as exc:
                raise CheckError(f'Registry HTTP error {exc.code}') from exc
            except (URLError, OSError, UnicodeError) as exc:
                raise CheckError(f'Cannot read registry: {exc}') from exc
        try:
            metadata = json.loads(text)
        except ValueError as exc:
            raise CheckError(f'Invalid metadata JSON: {exc}') from exc
        if not isinstance(metadata, dict):
            raise CheckError('Metadata must be an object')
        versions = metadata.get('versions')
        yanked = metadata.get('yanked_versions', {})
        if not isinstance(versions, list) or any(not isinstance(v, str) or not v for v in versions):
            raise CheckError('Metadata needs a list of nonempty version strings')
        if not isinstance(yanked, dict) or any(not isinstance(k, str) or not isinstance(v, str)
                                               for k, v in yanked.items()):
            raise CheckError('Invalid yanked_versions mapping')
        return metadata


@dataclass
class Result:
    name: str
    declared: str
    latest: str
    status: str
    notes: str = ''


def check_dependency(name, declared, metadata, include_prereleases=False):
    current = declared.split('+', 1)[0]
    current_key = version_key(declared)
    notes = []
    versions = metadata['versions']
    yanked = metadata.get('yanked_versions', {})
    if current not in versions:
        notes.append('Declared version is not listed in this registry')
    if current in yanked:
        notes.append('Declared version is yanked: ' + yanked[current])
    if metadata.get('deprecated'):
        notes.append('Module deprecated: ' + str(metadata['deprecated']))
    # Validate all versions, including excluded ones. Never silently pick an
    # older version because the newest record has an unfamiliar format.
    ordered = sorted(versions, key=version_key)
    eligible = [v for v in ordered if v not in yanked and
                (include_prereleases or not is_prerelease(v))]
    if not eligible:
        return Result(name, declared, '-', 'NO_ELIGIBLE_VERSION', '; '.join(notes))
    latest = eligible[-1]
    latest_key = version_key(latest)
    status = ('UPDATE_AVAILABLE' if current_key < latest_key else
              'DECLARED_NEWER' if current_key > latest_key else 'CURRENT')
    if current in yanked:
        notes.insert(0, 'Comparison: ' + status)
        status = 'YANKED'
    return Result(name, declared, latest, status, '; '.join(notes))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module', nargs='?', type=Path, default=Path('MODULE.bazel'),
                        help='MODULE.bazel to inspect (default: current directory)')
    parser.add_argument('--registry', type=Path,
                        help='Local registry root; otherwise query https://bcr.bazel.build')
    parser.add_argument('--include-prereleases', action='store_true',
                        help='Also consider prerelease versions (ISO-date releases are always considered)')
    parser.add_argument('--timeout', type=float, default=15, help='HTTPS timeout in seconds (default 15)')
    args = parser.parse_args(argv)
    if not 0 < args.timeout < float('inf'):
        parser.error('--timeout must be positive and finite')
    try:
        declarations, overrides = read_module(args.module)
        registry = Registry(args.registry, args.timeout)
    except CheckError as exc:
        print(f'Error: {exc}', file=sys.stderr)
        return 1
    source = str(args.registry.resolve()) if args.registry else 'https://bcr.bazel.build'
    print(f'Module: {args.module.resolve()}')
    print(f'Registry: {source}')
    if args.registry:
        print('Local snapshot only; the clone is not fetched or updated.')
    print('Comparing declared versions, not Bazel-resolved transitive versions.')
    rows = []
    errors = 0
    for name, declared in declarations.items():
        if name in overrides:
            rows.append(Result(name, declared or '-', '-', 'SKIPPED', overrides[name]))
        elif not declared:
            rows.append(Result(name, '-', '-', 'SKIPPED', 'No declared version'))
        else:
            try:
                rows.append(check_dependency(name, declared, registry.read(name), args.include_prereleases))
            except CheckError as exc:
                errors += 1
                rows.append(Result(name, declared, '-', 'ERROR', str(exc)))
    columns = [('MODULE', 'DECLARED', 'LATEST', 'STATUS')]
    columns += [(r.name, r.declared, r.latest, r.status) for r in rows]
    widths = [max(len(row[i]) for row in columns) for i in range(4)]
    for index, row in enumerate(columns):
        print('  '.join(value.ljust(width) for value, width in zip(row, widths)))
        if index and rows[index - 1].notes:
            print('  ' + rows[index - 1].notes)
    print(f'{len(rows)} dependencies; {sum(r.status == "UPDATE_AVAILABLE" for r in rows)} updates available; '
          f'{sum(r.status == "SKIPPED" for r in rows)} skipped; {errors} errors.')
    return 1 if errors else 0


if __name__ == '__main__':
    sys.exit(main())
