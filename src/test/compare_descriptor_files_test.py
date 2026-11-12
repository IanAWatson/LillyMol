"""Ensure numeric tolerance cannot hide changes to a descriptor table's structure."""
import unittest
from decimal import Decimal

from compare_descriptor_files import compare_lines


class CompareDescriptorFilesTest(unittest.TestCase):
    def compare(self, expected, actual, tolerances=None):
        compare_lines(expected.splitlines(), actual.splitlines(),
                      tolerances if tolerances is not None else {'value': Decimal('0.001')})

    def test_one_printed_unit_is_allowed(self):
        self.compare('Id value count\na 5.235 2\n', 'Id value count\na 5.236 2\n')
        self.compare('Id value\na -162.763\n', 'Id value\na -162.762\n')

    def test_larger_difference_is_rejected(self):
        with self.assertRaises(ValueError):
            self.compare('Id value\na 5.235\n', 'Id value\na 5.237\n')

    def test_untolerated_columns_and_identifiers_are_exact(self):
        for actual in ['Id value count\na 1 3', 'Id value count\nb 1 2',
                       'Id value count\na 1 2.0']:
            with self.assertRaises(ValueError):
                self.compare('Id value count\na 1 2', actual)

    def test_headers_row_counts_and_widths_are_exact(self):
        for actual in ['Id other\na 1', 'Id value', 'Id value\na 1\nb 2',
                       'Id value\na', 'Id value\na 1 extra']:
            with self.assertRaises(ValueError):
                self.compare('Id value\na 1', actual)

    def test_missing_and_nonfinite_values_cannot_use_tolerance(self):
        for value in ['.', 'NaN', 'Infinity']:
            with self.assertRaises(ValueError):
                self.compare('Id value\na 1', f'Id value\na {value}')
        self.compare('Id value\na .', 'Id value\na .')

    def test_tolerance_column_must_be_a_descriptor(self):
        for column in ['Id', 'missing']:
            with self.assertRaises(ValueError):
                self.compare('Id value\na 1', 'Id value\na 1', {column: Decimal('0.001')})


if __name__ == '__main__':
    unittest.main()
