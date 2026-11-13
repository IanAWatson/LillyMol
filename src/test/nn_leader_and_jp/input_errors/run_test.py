#!/usr/bin/env python3
"""Distinguish input corruption from clean EOF, preserving explicit -s limits."""
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tfrecord_utils import frame, masked_crc, string_field, tfrecord

EXE = str(Path(sys.argv[1]).resolve())
MODES = (['-h', '-t', '0.2'], ['-T', '-h', '-t', '0.2'], ['-J', '.,1'])
VALID = '$SMI<C>\nPCN<m0>\n|\n'


class InputErrors(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)

    def write(self, name, content):
        path = self.directory / name
        if isinstance(content, bytes):
            path.write_bytes(content)
        else:
            path.write_text(content)
        return path

    def run_tool(self, files, mode, serialized=False, limit=None, success=False):
        args = [EXE, *mode]
        if serialized:
            args.append('-r')
        if limit is not None:
            args.extend(['-s', str(limit)])
        result = subprocess.run([*args, *map(str, files)], capture_output=True,
                                text=True, timeout=15)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertGreater(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, '', 'An input error must not produce partial clusters')
        self.assertNotIn('AddressSanitizer', result.stderr)
        self.assertNotIn('runtime error:', result.stderr)
        return result

    def corrupt_frames(self):
        good = tfrecord('m1')
        bad_length_crc = bytearray(good)
        bad_length_crc[8] ^= 1
        bad_data_crc = bytearray(good)
        bad_data_crc[-1] ^= 1
        enormous_length = struct.pack('<Q', 1 << 63)
        return [good[:1], good[:11], good[:13], good[:-1],
                bytes(bad_length_crc), bytes(bad_data_crc),
                frame(b'\x12\x05m'),  # Invalid protobuf with valid frame CRCs.
                frame(b''), frame(string_field(1, 'C')),  # Missing identifiers.
                frame(string_field(2, '')), frame(string_field(2, ' ')),
                enormous_length + masked_crc(enormous_length)]

    def test_serialized_errors_after_valid_prefix(self):
        for mode in MODES:
            for i, bad in enumerate(self.corrupt_frames()):
                with self.subTest(mode=mode, bad=i):
                    path = self.write('bad.rio', tfrecord('m0') + bad)
                    result = self.run_tool([path], mode, serialized=True)
                    self.assertIn('Cannot extract identifiers', result.stderr)

    def test_serialized_error_in_middle(self):
        path = self.write('bad.rio', tfrecord('m0') + frame(b'\x12\x05m') + tfrecord('m1'))
        for mode in MODES:
            self.run_tool([path], mode, serialized=True)

    def test_serialized_unread_tail_ignored_only_after_limit(self):
        path = self.write('bad.rio', tfrecord('m0') + frame(b'\x12\x05m'))
        for mode in MODES:
            result = self.run_tool([path], mode, serialized=True, limit=1, success=True)
            self.assertEqual(re.findall(r'^PCN<([^>]+)>', result.stdout, re.MULTILINE), ['m0'])
            self.run_tool([path], mode, serialized=True, limit=2)

    def test_tdt_malformed_records_after_valid_prefix(self):
        bad_records = [
            '$SMI<C>\nPCN<>\n|\n', '$SMI<C>\nPCN< >\n|\n',
            '$SMI<C>\nPCN<\tm1>\n|\n', '$SMI<C>\nPCN<m1\n|\n',
            '$SMI<C>\nPCN<m1>\nPCN<>\nDIST<0>\n|\n',
            '$SMI<C>\n|\n', 'PCN<m1>\n|\n', '|\n',
            '$SMI<C>\nPCN<m1>\n', 'unrecognised trailing data\n',
            '$SMI<C\nPCN<m1>\n|\n',
        ]
        for mode in MODES:
            for i, bad in enumerate(bad_records):
                with self.subTest(mode=mode, bad=i):
                    path = self.write('bad.nn', VALID + bad)
                    self.run_tool([path], mode)

    def test_empty_or_malformed_numeric_values(self):
        for distance in ('DIST<>', 'DIST<0.1', 'DIST<text>'):
            path = self.write('bad.nn', VALID + '$SMI<C>\nPCN<m1>\nPCN<m0>\n' + distance + '\n|\n')
            for mode in MODES:
                self.run_tool([path], mode)

    def test_tdt_unread_tail_ignored_only_after_limit(self):
        path = self.write('bad.nn', VALID + '$SMI<C>\nPCN<>\n|\n')
        for mode in MODES:
            result = self.run_tool([path], mode, limit=1, success=True)
            self.assertEqual(re.findall(r'^PCN<([^>]+)>', result.stdout, re.MULTILINE), ['m0'])
            self.run_tool([path], mode, limit=2)

    def test_bad_later_file_does_not_accept_earlier_files(self):
        for serialized in (False, True):
            first = self.write('first', tfrecord('m0') if serialized else VALID)
            bad = self.write('second', b'\x00' if serialized else '$SMI<C>\nPCN<>\n|\n')
            for mode in MODES:
                self.run_tool([first, bad], mode, serialized, limit=99)
                self.run_tool([first, bad], mode, serialized, limit=1, success=True)

    def test_read_error_is_not_eof(self):
        # Opening a directory may succeed, but reading its contents must fail.
        for serialized in (False, True):
            first = self.write('first', tfrecord('m0') if serialized else VALID)
            for mode in MODES:
                self.run_tool([first, self.directory], mode, serialized, limit=99)

    def test_clean_eof_empty_later_file_and_empty_smiles(self):
        for serialized in (False, True):
            if serialized:
                content = tfrecord('m0') + frame(string_field(2, 'm1'))
                empty = b''
            else:
                content = VALID + '\n$SMI<>\nPCN<m1>\n|\n\n'
                empty = ''
            first = self.write('first', content)
            second = self.write('second', empty)
            for mode in MODES:
                result = self.run_tool([first, second], mode, serialized, limit=99, success=True)
                self.assertEqual(re.findall(r'^PCN<([^>]+)>', result.stdout, re.MULTILINE),
                                 ['m0', 'm1'])
                self.assertEqual(result.stderr, '')


if __name__ == '__main__':
    result = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(InputErrors))
    if not result.wasSuccessful():
        sys.exit(1)
    print('nn_leader_and_jp input error checks passed')
