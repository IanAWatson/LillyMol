"""Correctness checks for tcount's record selection, token modes and file handling."""

import pathlib
import re
import subprocess
import sys
import tempfile
import unittest

_EXECUTABLE = str(pathlib.Path(sys.argv.pop(1)).resolve())


class TcountTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._directory = pathlib.Path(self._tmp.name)
        self._input = self._directory / "input.txt"
        self._input.write_text("a b\nc d\ne f\n")

    def run_tool(self, *options, inputs=("input.txt",), success=True, stdin=""):
        result = subprocess.run([_EXECUTABLE, *options, *inputs], cwd=self._directory,
                                capture_output=True, text=True, input=stdin, timeout=20)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        if not success:
            self.assertGreater(result.returncode, 0, result.stderr)
        return result

    def counts(self, result):
        return {int(columns): int(records.replace(",", "")) for records, columns in
                re.findall(r"^([\d,]+) records had (\d+) tokens", result.stderr, re.M)}

    def test_record_limit_is_exact(self):
        result = self.run_tool("-n", "1", "-C", "counts.txt")
        self.assertEqual(self.counts(result), {2: 1})
        self.assertIn("1 records read", result.stderr)
        self.assertEqual((self._directory / "counts.txt").read_text(), "Line Tcount\n1 2\n")

    def test_record_limit_applies_per_file(self):
        (self._directory / "second.txt").write_text("a b\nc d\n")
        result = self.run_tool("-n", "1", "-v", inputs=("input.txt", "second.txt"))
        self.assertEqual(result.stderr.count("1 records read"), 2)
        self.assertIn("Read 2 lines from 2 files", result.stderr)

    def test_default_and_empty_field_modes(self):
        self._input.write_text(",a,,b,\n")
        self.assertEqual(self.counts(self.run_tool("-d", "comma")), {2: 1})
        self.assertEqual(self.counts(self.run_tool("-d", "comma", "-u")), {5: 1})

    def test_quoted_fields_include_leading_and_trailing_empty_fields(self):
        self._input.write_text(',"a,b",\n')
        for options in (("-q",), ("-q", "-u")):
            with self.subTest(options=options):
                self.assertEqual(self.counts(self.run_tool("-d", "comma", *options)), {3: 1})

    def test_doubled_quotes(self):
        self._input.write_text('"a""b,c",d\n')
        self.assertEqual(self.counts(self.run_tool("-q", "-d", "comma")), {2: 1})

    def test_empty_record_has_no_tokens_in_all_modes(self):
        self._input.write_text("\n")
        for options in ((), ("-u",), ("-q",), ("-k",), ("-q", "-u", "-k")):
            with self.subTest(options=options):
                self.assertEqual(self.counts(self.run_tool(*options)), {0: 1})

    def test_backslash_mode_counts_tokens_not_separators(self):
        for text, count in (("abc", 1), ("a,b", 2), (r"a\,b,c", 2),
                            (r"a\,,b", 2), (r"a\\,b", 2), (",a,,b,", 2)):
            with self.subTest(text=text):
                self._input.write_text(text + "\n")
                self.assertEqual(self.counts(self.run_tool("-k", "-d", "comma")), {count: 1})

    def test_backslash_mode_combines_with_empty_fields(self):
        self._input.write_text(r"a\,,,b" + "\n")
        self.assertEqual(self.counts(self.run_tool("-k", "-u", "-d", "comma")), {3: 1})

    def test_backslash_escaped_quotes(self):
        self._input.write_text(r',"a\",b",c' + "\n")
        self.assertEqual(self.counts(self.run_tool("-q", "-k", "-d", "comma")), {3: 1})

    def test_unmatched_quote_warns_without_failing(self):
        self._input.write_text('a,"b,c\n')
        result = self.run_tool("-q", "-d", "comma")
        self.assertEqual(self.counts(result), {2: 1})
        self.assertIn("unmatched quote", result.stderr)

    def test_translation_preserves_quoted_and_escaped_separators(self):
        cases = (((',"a,b",c\n'), ("-q",), '|"a,b"|c\n'),
                 ((r"a\,b,c" + "\n"), ("-k",), r"a\,b|c" + "\n"))
        for text, options, expected in cases:
            with self.subTest(options=options):
                self._input.write_text(text)
                columns = 3 if "-q" in options else 2
                self.run_tool(*options, "-d", "comma", "-t", str(columns),
                              "-T", "selected.txt", "-o", "vbar")
                self.assertEqual((self._directory / "selected.txt").read_text(), expected)

    def test_uppercase_literal_token_and_record_delimiters(self):
        self._input.write_text("aXb\n")
        self.assertEqual(self.counts(self.run_tool("-d", "X")), {2: 1})
        self._input.write_text("a bXc dX")
        self.assertEqual(self.counts(self.run_tool("-c", "X")), {2: 2})

    def test_named_record_delimiters(self):
        for name in ("DOS", "windows", "CRLF", "newline", "linux", "unix", "LF"):
            for newline in (b"\n", b"\r\n"):
                with self.subTest(name=name, newline=newline):
                    self._input.write_bytes(newline.join((b"a b", b"c d", b"")))
                    result = self.run_tool("-c", name, "-v")
                    self.assertEqual(self.counts(result), {2: 2})
                    length = 3 if newline == b"\n" or name.lower() in ("dos", "windows", "crlf") else 4
                    self.assertIn(f"records between {length} and {length}", result.stderr)

    def test_windows_mode_consumes_crlf_while_lf_mode_preserves_cr(self):
        self._input.write_bytes(b"a b\r\nc d\r\n")
        for mode in ("DOS", "windows", "CRLF", "newline", "linux", "LF"):
            with self.subTest(mode=mode):
                self.run_tool("-c", mode, "-t", "2", "-T", "selected.txt")
                expected = b"a b\nc d\n" if mode.lower() in ("dos", "windows", "crlf") else self._input.read_bytes()
                self.assertEqual((self._directory / "selected.txt").read_bytes(), expected)

    def test_empty_windows_record_differs_from_literal_cr_in_lf_record(self):
        self._input.write_bytes(b"\r\n")
        self.assertEqual(self.counts(self.run_tool("-c", "DOS")), {0: 1})
        self.assertEqual(self.counts(self.run_tool("-c", "linux")), {1: 1})
        self.assertEqual(self.counts(self.run_tool()), {0: 1})  # Historical default.

    def test_windows_records_with_eof_offset(self):
        self._input.write_bytes(b"a b\r\nc d\r\n")
        result = self.run_tool("-c", "DOS", "-e", "5", "-v")
        self.assertEqual(self.counts(result), {2: 1})
        self.assertIn("records between 3 and 3", result.stderr)

    def test_named_custom_record_delimiter(self):
        self._input.write_text("a b;c d;")
        self.assertEqual(self.counts(self.run_tool("-c", "semic")), {2: 2})
        self._input.write_bytes(b"a b\rc d\r")
        for name in ("CR", "^m"):
            self.assertEqual(self.counts(self.run_tool("-c", name)), {2: 2})

    def test_empty_delimiters_are_rejected(self):
        for flag in ("-c", "-d", "-o"):
            for value in ("", "empty"):
                with self.subTest(flag=flag, value=value):
                    self.run_tool(flag, value, success=False)

    def test_record_and_token_delimiter_help(self):
        for flag, text in (("-c", "DOS"), ("-d", "comma")):
            with self.subTest(flag=flag):
                result = self.run_tool(flag, "help", inputs=())
                self.assertIn(text, result.stderr)

    def test_removed_length_option_is_rejected(self):
        result = self.run_tool("-L", "2", success=False)
        self.assertIn("unrecognised option '-L'", result.stderr)

    def test_eof_offset_preserves_aligned_records(self):
        self._input.write_text("a b\nc d\n")
        for offset, records in ((8, 2), (4, 1), (5, 1), (3, 0), (1, 0)):
            with self.subTest(offset=offset):
                result = self.run_tool("-e", str(offset))
                expected = {2: records} if records else {}
                self.assertEqual(self.counts(result), expected)
                self.assertIn(f"{records} records read", result.stderr)

    def test_eof_offset_combines_with_record_limit(self):
        self._input.write_text("a b\nc d\ne f\n")
        result = self.run_tool("-e", "8", "-n", "1", "-C", "counts.txt")
        self.assertEqual(self.counts(result), {2: 1})
        self.assertEqual((self._directory / "counts.txt").read_text(), "Line Tcount\n1 2\n")

    def test_eof_offset_respects_custom_record_delimiter(self):
        self._input.write_text("a b;c d;e f;")
        self.assertEqual(self.counts(self.run_tool("-c", "semic", "-e", "8")), {2: 2})

    def test_length_statistics_are_per_file(self):
        self._input.write_text("x\n")
        (self._directory / "long.txt").write_text("abcdefgh\n")
        result = self.run_tool("-v", inputs=("input.txt", "long.txt"))
        self.assertIn("'input.txt' records between 1 and 1", result.stderr)
        self.assertIn("'long.txt' records between 8 and 8", result.stderr)
        self.assertNotIn("ave 4.5", result.stderr)

    def test_abort_counts_records_and_continues_other_files(self):
        self._input.write_text("a b\nc\ne f\n")
        (self._directory / "second.txt").write_text("a b\nc d\n")
        result = self.run_tool("-b", "-g", "-v", inputs=("input.txt", "second.txt"))
        self.assertIn("2 records read from 'input.txt'", result.stderr)
        self.assertIn("Read 4 lines from 2 files", result.stderr)
        self.assertIn("1 non rectangular files found", result.stderr)

    def test_reported_and_unreported_records(self):
        self._input.write_text("a b\nc\nx y\nz z\n")
        self.run_tool("-r", "-R", "reported.txt", "-U", "unreported.txt")
        self.assertEqual((self._directory / "reported.txt").read_text(), "c\nx y\n")
        self.assertEqual((self._directory / "unreported.txt").read_text(), "a b\nz z\n")

    def test_duplicate_token_selections_write_each_destination(self):
        self.run_tool("-t", "2", "-T", "first.txt", "-t", "2", "-T", "second.txt")
        self.assertEqual((self._directory / "first.txt").read_text(), self._input.read_text())
        self.assertEqual((self._directory / "second.txt").read_text(), self._input.read_text())

    def test_input_output_conflicts_are_rejected_before_truncation(self):
        original = self._input.read_text()
        for options in (("-r", "-R", "input.txt"), ("-r", "-U", "./input.txt"),
                        ("-C", "input.txt"), ("-t", "2", "-T", "input.txt")):
            with self.subTest(options=options):
                self.run_tool(*options, success=False)
                self.assertEqual(self._input.read_text(), original)

    def test_output_cannot_replace_later_input(self):
        path = self._directory / "second.txt"
        path.write_text("original\n")
        self.run_tool("-r", "-R", "second.txt", inputs=("input.txt", "second.txt"), success=False)
        self.assertEqual(path.read_text(), "original\n")

    def test_hard_link_and_symlink_input_aliases_are_rejected(self):
        original = self._input.read_text()
        alias = self._directory / "alias.txt"
        for hardlink in (False, True):
            with self.subTest(hardlink=hardlink):
                if hardlink:
                    alias.hardlink_to(self._input)
                else:
                    alias.symlink_to(self._input)
                self.run_tool("-C", "alias.txt", success=False)
                self.assertEqual(self._input.read_text(), original)
                alias.unlink()

    def test_output_destinations_must_differ(self):
        self.run_tool("-r", "-R", "out.txt", "-C", "./out.txt", success=False)
        self.assertFalse((self._directory / "out.txt").exists())

    def test_output_comparison_error_is_not_ignored(self):
        (self._directory / "loop").symlink_to("loop")
        self.run_tool("-C", "loop/counts.txt", success=False)
        self.assertEqual(self._input.read_text(), "a b\nc d\ne f\n")

    def test_buffered_output_failures_return_nonzero(self):
        if not pathlib.Path("/dev/full").exists():
            self.skipTest("Requires /dev/full")
        self._input.write_text("a b\nc\n")
        for options in (("-C", "/dev/full"), ("-r", "-R", "/dev/full"),
                        ("-r", "-U", "/dev/full"), ("-t", "2", "-T", "/dev/full")):
            with self.subTest(options=options):
                result = self.run_tool(*options, success=False)
                self.assertIn("output", result.stderr)

    def test_stdin_and_empty_input(self):
        result = self.run_tool(inputs=("-",), stdin="a b\nc d\n")
        self.assertEqual(self.counts(result), {2: 2})
        self._input.write_text("")
        result = self.run_tool("-v")
        self.assertIn("0 records read", result.stderr)
        self.assertNotIn("records between", result.stderr)

    def test_human_readable_counts(self):
        self._input.write_text("a b\n" * 1000)
        result = self.run_tool("-h")
        self.assertEqual(self.counts(result), {2: 1000})
        self.assertIn("1,000 records read", result.stderr)


if __name__ == "__main__":
    unittest.main()
