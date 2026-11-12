#!/usr/bin/env python3
"""Compare descriptor tables exactly, allowing named columns an absolute tolerance."""

import argparse
from decimal import Decimal, InvalidOperation
from itertools import zip_longest
import sys


def compare_lines(expected, actual, tolerances):
    expected, actual = iter(expected), iter(actual)
    header = next(expected, '').split()
    if not header or header != next(actual, '').split():
        raise ValueError('Descriptor headers differ or are missing')
    if len(header) != len(set(header)):
        raise ValueError('Duplicate descriptor columns')
    if header[0] in tolerances or not tolerances.keys() <= set(header[1:]):
        raise ValueError('Tolerances must name existing descriptor columns')
    for line_number, (left, right) in enumerate(zip_longest(expected, actual), 2):
        if left is None or right is None:
            raise ValueError(f'Line {line_number}: row counts differ')
        left, right = left.split(), right.split()
        if len(left) != len(header) or len(right) != len(header):
            raise ValueError(f'Line {line_number}: column count differs from header')
        for column, x, y in zip(header, left, right):
            if x == y:
                continue
            if column in tolerances:
                try:
                    a, b = Decimal(x), Decimal(y)
                    if a.is_finite() and b.is_finite() and abs(a - b) <= tolerances[column]:
                        continue
                except InvalidOperation:
                    pass
            raise ValueError(f'Line {line_number}, {column}: expected {x}, got {y}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tolerance', action='append', default=[], metavar='COLUMN=VALUE')
    parser.add_argument('expected')
    parser.add_argument('actual')
    args = parser.parse_args()
    tolerances = {}
    try:
        for specification in args.tolerance:
            column, value = specification.split('=', 1)
            tolerance = Decimal(value)
            if not tolerance.is_finite() or tolerance < 0 or column in tolerances:
                raise ValueError(f'Invalid tolerance: {specification}')
            tolerances[column] = tolerance
        with open(args.expected) as expected, open(args.actual) as actual:
            compare_lines(expected, actual, tolerances)
    except (OSError, ValueError, InvalidOperation) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
