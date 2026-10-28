"""Executable regression tests; all fixtures are molecule-independent proto data."""
import json
from pathlib import Path
import statistics
import struct
import subprocess
import sys
import tempfile
import unittest

EXE = Path(sys.argv.pop(1)).resolve()


def crc32c(data):
    crc = 0xffffffff
    for byte in data:
        crc ^= byte
        for _ in range(8):
            crc = (crc >> 1) ^ (0x82f63b78 if crc & 1 else 0)
    return crc ^ 0xffffffff


def masked(data):
    crc = crc32c(data)
    return (((crc >> 15) | (crc << 17)) + 0xa282ead8) & 0xffffffff


def tfrecord(payload):
    length = struct.pack('<Q', len(payload))
    return length + struct.pack('<I', masked(length)) + payload + struct.pack('<I', masked(payload))


def varint(value):
    result = bytearray()
    while value >= 128:
        result.append((value & 127) | 128)
        value >>= 7
    result.append(value)
    return bytes(result)


def string_field(number, value):
    value = value.encode() if isinstance(value, str) else value
    return varint(number * 8 + 2) + varint(len(value)) + value


def molecule_proto(name, smiles, fragments):
    result = string_field(1, name) + string_field(2, smiles)
    for fragment in fragments:
        encoded = string_field(2, fragment.get('smi', ''))
        encoded += b'\x20' + varint(fragment.get('nat', 1))
        encoded += b'\x28' + varint(fragment.get('n', 1))
        if 'usmi' in fragment:
            encoded += string_field(9, fragment['usmi'])
        result += string_field(5, encoded)
    return result


def molecule_text(name, smiles, fragments):
    result = f'name:{json.dumps(name)} smiles:{json.dumps(smiles)}'
    for fragment in fragments:
        result += ' fragment {' + ' '.join(
            f'{key}:{json.dumps(value)}' for key, value in fragment.items()) + '}'
    return result + '\n'


class FragmentActivityTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.activity = self.root / 'activity.txt'
        self.activity.write_text('ID activity\nA 0\nB 2\nC 4\nD 6\nUNSEEN 100\n')
        self.rows = [
            ('A', 'CC', [{'smi': 'C', 'nat': 1, 'n': 99}, {'smi': 'C', 'nat': 1}]),
            ('B', 'CN', [{'smi': 'C', 'nat': 1}, {'smi': 'N', 'nat': 1}]),
            ('C', 'N', [{'smi': 'N', 'nat': 1}]),
            ('D', 'O', []),
        ]
        self.text = self.root / 'diced.textproto'
        self.write_rows()

    def write_rows(self):
        self.text.write_text(''.join(molecule_text(*row) for row in self.rows))

    def run_tool(self, *options, inputs=None, success=True):
        if inputs is None:
            inputs = [self.text]
        result = subprocess.run([str(EXE), '-A', str(self.activity), *map(str, options),
                                 *map(str, inputs)], capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def report(self, *options, **kwargs):
        return json.loads(self.run_tool('-json', '-', *options, **kwargs).stdout)

    def test_population_support_quantiles_and_histograms(self):
        report = self.report('-bins', '3')
        self.assertEqual(report['global']['summary'], dict(
            n=4, min=0, q1=1.5, median=3, mean=3, q3=4.5, max=6))
        self.assertEqual(report['global']['histogram'], [1, 1, 2])
        self.assertEqual(report['metadata']['activities_unused'], 1)
        self.assertEqual([row['source'] for row in report['queries']], ['C', 'N'])
        carbon = report['queries'][0]
        self.assertEqual(carbon['matched']['summary'], dict(
            n=2, min=0, q1=.5, median=1, mean=1, q3=1.5, max=2))
        self.assertEqual(carbon['unmatched']['summary'], dict(
            n=2, min=4, q1=4.5, median=5, mean=5, q3=5.5, max=6))
        self.assertEqual(carbon['matched']['histogram'], [1, 1, 0])
        self.assertEqual(carbon['unmatched']['histogram'], [0, 0, 2])
        self.assertEqual(carbon['delta_mean'], -4)
        self.assertEqual(carbon['delta_median'], -4)
        self.assertEqual([e['identifier'] for e in carbon['examples']], ['A', 'B'])
        rows = [line.split() for line in self.run_tool().stdout.splitlines()]
        self.assertEqual(rows[0], 'Frag Parent Natoms N Fraction Minval Ave Maxval Diff'.split())
        self.assertEqual(rows[1][:4], ['C', 'A', '1', '2'])
        self.assertEqual(list(map(float, rows[1][4:])), [.5, 0, 1, 2, -2])

    def test_binary_matches_text(self):
        binary = self.root / 'diced.tfdata'
        binary.write_bytes(b''.join(tfrecord(molecule_proto(*row)) for row in self.rows))
        self.assertEqual(self.report(), self.report('-tfdata', inputs=[binary]))

    def test_corrupt_binary_is_an_error(self):
        good = tfrecord(molecule_proto(*self.rows[0]))
        bad_crc = bytearray(good)
        bad_crc[-1] ^= 1
        length = struct.pack('<Q', 2**63)
        oversized = length + struct.pack('<I', masked(length))
        for content in [good + tfrecord(b'\xff'), good + b'\x01', good[:-1],
                        bytes(bad_crc), oversized, tfrecord(b'')]:
            with self.subTest(content=content):
                binary = self.root / 'bad.tfdata'
                binary.write_bytes(content)
                result = self.run_tool('-tfdata', inputs=[binary], success=False)
                self.assertEqual(result.stdout, '')

    def test_duplicate_and_missing_identifiers(self):
        self.run_tool(inputs=[self.text, self.text], success=False)
        self.rows[1] = ('A additional name', 'CN', [])
        self.write_rows()
        self.assertIn('Duplicate molecule', self.run_tool(success=False).stderr)
        self.rows[1] = ('UNKNOWN', 'CN', [])
        self.write_rows()
        self.assertIn('No activity', self.run_tool(success=False).stderr)
        self.activity.write_text('ID value\nA 0\nA 1\n')
        self.assertIn('Duplicate activity', self.run_tool(success=False).stderr)

    def test_support_limits_examples_and_constant_activity(self):
        self.assertEqual(self.report('-support', '3')['queries'], [])
        self.assertEqual(self.report('-unmatched_support', '3')['queries'], [])
        self.assertTrue(all(not row['examples'] for row in self.report('-examples', '0')['queries']))
        self.rows = [('A', 'CC', [{'smi': 'CC', 'nat': 2}]), ('B', 'N', [])]
        self.write_rows()
        self.assertEqual(self.report('-C', '1')['queries'], [])
        self.activity.write_text('ID value\nA 5\nB 5\n')
        report = self.report()
        self.assertEqual(report['global']['bins']['edges'], [5, 5])
        self.assertEqual(report['global']['histogram'], [2])
        self.assertEqual(report['queries'][0]['delta_median'], 0)
        self.rows = [('A', 'C', []), ('B', 'N', [])]
        self.write_rows()
        self.assertEqual(self.report()['queries'], [])

    def test_coordinate_free_identity_and_validation(self):
        self.rows = [('A', 'CC', [{'smi': 'C{{1,2,3}}', 'usmi': 'C', 'nat': 1}]),
                     ('B', 'CN', [{'smi': 'C{{4,5,6}}', 'usmi': 'C', 'nat': 1}]),
                     ('C', 'N', []), ('D', 'O', [])]
        self.write_rows()
        self.assertEqual(self.report()['queries'][0]['matched']['summary']['n'], 2)
        self.assertEqual(self.report()['queries'][0]['source'], 'C')
        self.rows[1][2][0]['nat'] = 2
        self.write_rows()
        self.assertIn('Inconsistent atom counts', self.run_tool(success=False).stderr)

    def test_empty_and_invalid_inputs_and_clobbering(self):
        self.run_tool(inputs=[], success=False)
        self.text.write_text('')
        self.run_tool(success=False)
        self.write_rows()
        for filename in [self.activity, self.text]:
            before = filename.read_bytes()
            self.run_tool('-json', filename, success=False)
            self.assertEqual(filename.read_bytes(), before)
        alias = self.root / 'activity-hardlink.txt'
        alias.hardlink_to(self.activity)
        before = self.activity.read_bytes()
        self.run_tool('-json', alias, success=False)
        self.assertEqual(self.activity.read_bytes(), before)
        directory = self.root / 'not-a-file'
        directory.mkdir()
        self.run_tool(inputs=[self.text, directory], success=False)
        for flag, value in [('-bins', '0'), ('-examples', '-1'), ('-support', '-1')]:
            self.run_tool(flag, value, success=False)
        self.activity.write_text('ID value\nA nan\n')
        self.run_tool(success=False)

    def test_complement_statistics_against_explicit_groups(self):
        values = [(i * 7) % 11 - 5 for i in range(23)]
        memberships = {('C' * (j + 1)): {i for i in range(23) if (i * 3 + j * 7) % 13 < j % 10 + 1}
                       for j in range(14)}
        self.activity.write_text('ID value\n' + ''.join(f'M{i} {value}\n' for i, value in enumerate(values)))
        self.rows = [(f'M{i}', 'C' * 14, [{'smi': key, 'nat': len(key)}
                       for key, members in memberships.items() if i in members]) for i in range(23)]
        self.write_rows()
        report = self.report('-bins', '7')
        for row in report['queries']:
            members = memberships[row['source']]
            for field, group in [('matched', [values[i] for i in members]),
                                 ('unmatched', [v for i, v in enumerate(values) if i not in members])]:
                with self.subTest(fragment=row['source'], group=field):
                    actual = row[field]['summary']
                    self.assertEqual(actual['n'], len(group))
                    self.assertAlmostEqual(actual['mean'], statistics.fmean(group))
                    self.assertEqual(actual['median'], statistics.median(group))
                    if len(group) > 1:
                        quantiles = statistics.quantiles(group, method='inclusive')
                        self.assertAlmostEqual(actual['q1'], quantiles[0])
                        self.assertAlmostEqual(actual['q3'], quantiles[2])
                    self.assertEqual(sum(row[field]['histogram']), len(group))


if __name__ == '__main__':
    unittest.main()
