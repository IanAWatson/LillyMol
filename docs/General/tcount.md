# tcount
`tcount` is like `wc` except that it also counts the number of tokens
on each line. This can be particularly helpful when dealing with tabular
files which should be tabular. Debugging problems in such files can
be difficult.

Usage 
```
Usage: tcount <options> <fname1> <fname2> ...
  -n <number>       process only the first <number> records from each file
  -e <offset>       start processing <offset> bytes before EOF
  -r                report different token count from previous record
  -f                report different token count from first record
  -R <name>         write reported   records to <name>
  -U <name>         write unreported records to <name>
  -t <count>        write records with <count> tokens to -T <fname>
  -T <fname>        file name for corresponding -t
  -o <char>         when writing -T files, use <char> as output delimiter
  -C <fname>        write a descriptor file of 'line.number tcount'
  -b                abort on different token count from previous record
  -w                write reported records to stderr
  -d <delim>        delimiter is any character, or a token word(e.g. 'tab'). see '-d help' for a list (default ' ')
  -c <delim>        record delimiter: newline/linux/DOS or a character (see -c help)
  -u                consecutive delimiters mean empty words in between
  -g                report non rectangular files
  -q                data is quoted tokens
  -k                data may contain backslashed token separators
  -v                verbose output
```
The `-r` option can be particularly useful, since rather than just accumulating
statics on how many rows have a how many tokens, it will report any occurrence
of differing token counts from one line to the next.

For reasons that are not clear, tcount can also be faster than wc on large files.
I think it might be related to shell variable LANG, but not sure.

## Delimiters and quoted fields

Prefer a named record delimiter: `-c newline` or `-c linux` for LF records,
and `-c DOS` for Windows CRLF records. Names are case-insensitive.
`DOS`, `windows`, and `CRLF` enable the reader's Windows mode: LF ends the
record and the preceding CR is removed. `newline`, `linux`, `unix`, and `LF`
select LF-only handling, preserving any preceding CR as record data.
Without `-c`, the historical default accepts both LF and CRLF and strips CR.
For another record separator, supply a literal character or a character name,
such as `-c semic`. `-c CR` and `-c '^m'` select a literal carriage return;
use `-c DOS` for CRLF files. `-c help` lists the choices.

Token delimiters use `-d`, for example `-d tab`, `-d comma`, or `-d X`.
Literal characters retain their case. `-d help` lists character names.
By default repeated token delimiters are collapsed and leading/trailing
delimiters do not add tokens. `-u` counts empty fields between delimiters,
including leading and trailing empty fields. An empty record has zero tokens.

`-q` ignores delimiters within double quotes and retains empty fields. It
handles doubled quotes inside quoted fields. `-k` makes backslash escape the
next character, including a token separator, quote, or another backslash.
These options combine: `-q -u` keeps quote handling, and `-q -k` also allows
backslash-escaped quotes. Unmatched quotes produce a warning; multiline
quoted fields are not supported.

```shell
tcount -c DOS -d comma -q data.csv
tcount -c linux -d tab -u data.tsv
```

When `-o` changes the token delimiter in `-T` output, quoted and escaped
separators remain unchanged. Output records end with LF. Repeated `-t` values
can send matching records to more than one `-T` destination.

## Record selection and statistics

`-n N` processes at most N records from each input. `-e N` starts N bytes
before EOF: a record beginning exactly at that offset is retained; a partial
initial record is discarded. Counts, `-C` record numbers, and the `-n` limit
apply to the complete records processed after that point. `-C` numbering
starts at one for each input file.

Verbose record-length statistics are computed separately for each input file.
With multiple inputs, `-v` also reports the total number of processed records.
`-b` stops the current file at the first change in token count, includes that
record in its statistics, and continues to subsequent input files.

The unfinished `-L` per-column length option has been removed.


## Exit status

`tcount` returns zero after examining its inputs, even when a file is
non-rectangular. Options such as `-g`, `-r`, and `-b` do not turn column-count
differences into a failing exit status.

Missing input files, fatal read errors, failed end-of-file offset operations
(`-e`), and failures to open, write, or close requested output files return a nonzero
status. Output files must differ from all input files and from one another;
checks include relative/absolute aliases, symlinks and hard links.
When processing multiple inputs, an input-file failure does not prevent later
files from being examined, but the final exit status remains nonzero.
