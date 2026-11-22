#include <cstdlib>

#include <fstream>
#include <limits>
#include <memory>
#include <vector>
#include <iostream>

#include <cassert>

using std::cerr;

static int verbose = 0;

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/string_with_commas.h"

#include "tokens_in_quoted_string.h"

static uint64_t records_to_process = 0;

struct FileStatistics {
  extending_resizable_array<uint64_t> _token_counts;
  Accumulator_Int<uint64_t> _record_length;
  uint64_t _records = 0;
  int _first_token_count = -1;
  int _previous_token_count = -1;
  bool _rectangular = true;
};

static DelimitedTextOptions token_options;

//  Sometimes it is the end of the file which is of interest

static off_t offset_from_eof = 0;

static int abort_on_different_token_count = 0;

static int report_different_token_count = 0;

static int report_different_from_first_record = 0;

static char token_delimiter = ' ';

static char output_token_delimiter = ' ';

static int consecutive_delimiters_mean_empty_words_between = 0;

static char record_delimiter = '\n';
// Retain the historical tolerant default; an explicit -c selects the mode.
static bool dos_record_terminators = true;

/*
  When reporting a record different, do we just report the record
  number of do we also write the record
*/

static int write_reported_records_to_cerr = 0;

static IWString new_fname_unreported;
static IWString new_fname_reported;

static std::ofstream new_fstream_unreported;
static std::ofstream new_fstream_reported;

/*
  When dealing with multiple files, it is sometimes convenient to get
  a count of the number of non-rectangular files
*/

static int report_non_rectangular_files = 0;

static uint64_t non_rectangular_files_found = 0;

static int data_is_quoted_tokens = 0;

static std::ofstream stream_for_count;

static int contains_backslashed_separators = 0;

// We can optionally write large numbers with commas.
static int human_readable_output = 0;

// Accumulated across all files.
static uint64_t total_lines_read = 0;

/*
  When writing different column counts to different files, we need a
  means of specifying that
*/

class Stream_for_Columns : public std::ofstream {
 private:
  const int _tokens;
  const IWString _fname;

 public:
  Stream_for_Columns(int, const char *);

  const IWString& fname() const { return _fname; }

  int tokens() const {
    return _tokens;
  }
};

Stream_for_Columns::Stream_for_Columns(int nc, const char *fname) : _tokens(nc), _fname(fname) {
  open(fname, std::ios::out);

  if (!good()) {
    cerr << "Cannot open '" << fname << "' for " << _tokens << " stream\n";
  } else if (verbose) {
    cerr << "Stream '" << fname << "' initialised for " << _tokens << " tokens\n";
  }

  return;
}

static resizable_array_p<Stream_for_Columns> cstreams;

static bool
process_to_streams(const const_IWSubstring& buffer, int ntokens) {
  for (Stream_for_Columns* stream : cstreams) {
    if (ntokens != stream->tokens()) {
      continue;
    }
    if (token_delimiter == output_token_delimiter) {
      *stream << buffer << '\n';
    } else {
      IWString translated;
      CountDelimitedTokens(buffer, token_options, 0, &translated, output_token_delimiter);
      *stream << translated << '\n';
    }
    if (!stream->good()) {
      cerr << "Error writing '" << stream->fname() << "'\n";
      return false;
    }
  }
  return true;  // A record need not match any requested column count.
}

template class resizable_array_p<Stream_for_Columns>;
template class resizable_array_base<Stream_for_Columns *>;

static void
DisplayDelimiterOptions(int rc) {
  char_name_to_char_usage("-d");
  std::exit(rc);
}

static void
DisplayRecordDelimiterOptions(int rc) {
  cerr << "Record delimiters: newline/linux/unix/LF for LF records,\n"
       << "DOS/windows/CRLF for Windows records (LF separator plus CR stripping).\n"
       << "Explicit LF mode preserves CR as data; the default accepts LF or CRLF.\n"
       << "Use a literal character or a character name for another separator.\n"
       << "Use ^m or CR for a literal carriage-return separator.\n";
  char_name_to_char_usage("-c");
  std::exit(rc);
}

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << "Counts tokens per line. Designed for identifying tabular and non-tabular files\n";
  cerr << "  -n <number>       process only the first <number> records from each file\n";
  cerr << "  -e <offset>       start processing <offset> bytes before EOF\n";
  cerr << "  -r                report different token count from previous record\n";
  cerr << "  -f                report different token count from first record\n";
  cerr << "  -R <name>         write reported   records to <name>\n";
  cerr << "  -U <name>         write unreported records to <name>\n";
  cerr << "  -t <count>        write records with <count> tokens to -T <fname>\n";
  cerr << "  -T <fname>        file name for corresponding -t\n";
  cerr << "  -o <char>         when writing -T files, use <char> as output delimiter\n";
  cerr << "  -C <fname>        write a descriptor file of 'line.number tcount'\n";
  cerr << "  -b                abort on different token count from previous record\n";
  cerr << "  -w                write reported records to stderr\n";
  cerr << "  -c <delim>        record delimiter: newline/linux/DOS or a character (see -c help)\n";
  cerr << "  -u                consecutive delimiters mean empty words in between\n";
  cerr << "  -g                report non rectangular files\n";
  cerr << "  -q                quoted fields; delimiters inside quotes are data\n";
  cerr << "  -k                backslash escapes the next character (combines with -q/-u)\n";
  cerr << "  -h                human readable output, 1000000 written as 1,000,000\n";
  cerr << "  -v                verbose output\n";
// clang-format on

  exit(rc);
}

/*
  An R file might have one fewer token on the first record than all other records
*/

static int
header_record_has_only_different_token_count(
    const extending_resizable_array<uint64_t> &token_counts,
    const int tokens_on_first_record) {
  resizable_array<int> non_zero_token_counts;

  int nt = token_counts.number_elements();

  for (int i = 1; i < nt; i++) {
    if (0 == token_counts[i]) {
      continue;
    }

    non_zero_token_counts.add(i);

    if (1 == non_zero_token_counts.number_elements())  // first non zero token count
    {
      if (1 != token_counts[i]) {  // if an R file, there must be a single header record
                                   // with that token count
        return 0;
      }

      if (tokens_on_first_record != i) {
        return 0;
      }
    } else if (non_zero_token_counts.number_elements() > 2) {
      return 0;
    }
  }

  return 2 == non_zero_token_counts.number_elements();
}

static int
CountTokens(const const_IWSubstring& buffer, uint64_t line_number) {
  // Keep the ordinary counting paths as simple as before.
  if (data_is_quoted_tokens || contains_backslashed_separators) {
    return CountDelimitedTokens(buffer, token_options, line_number);
  }
  if (consecutive_delimiters_mean_empty_words_between) {
    return buffer.nwords_single_delimiter(token_delimiter);
  }
  return buffer.nwords(token_delimiter);
}

static bool
ProcessInput(iwstring_data_source& input, FileStatistics& stats) {
  const_IWSubstring buffer;
  while ((records_to_process == 0 || stats._records < records_to_process) &&
         input.next_record(buffer)) {
    ++stats._records;
    const int tokens = CountTokens(buffer, stats._records);
    ++stats._token_counts[tokens];
    if (verbose) {
      stats._record_length.extra(buffer.length());
    }
    if (verbose > 2) {
      cerr << "Record " << stats._records << " had " << tokens << " words and "
           << buffer.length() << " characters\n";
    }
    if (!process_to_streams(buffer, tokens)) {
      return false;
    }
    if (stream_for_count.is_open()) {
      stream_for_count << stats._records << ' ' << tokens << '\n';
      if (!stream_for_count.good()) {
        cerr << "Error writing per-record token counts\n";
        return false;
      }
    }
    if (stats._first_token_count < 0) {
      stats._first_token_count = tokens;
    }
    bool reported = false;
    if (stats._previous_token_count >= 0 && stats._previous_token_count != tokens) {
      stats._rectangular = false;
      if (abort_on_different_token_count) {
        cerr << "Line " << stats._records << " has " << tokens << " tokens\n"
             << "Previous line had " << stats._previous_token_count << '\n';
        break;
      }
      if (report_different_token_count) {
        cerr << "Token count mismatch line " << stats._records << ". " << tokens
             << " vs " << stats._previous_token_count << '\n';
        reported = true;
      }
    }
    if (report_different_from_first_record && tokens != stats._first_token_count) {
      cerr << "Token count mismatch, line " << stats._records << ". First record had "
           << stats._first_token_count << '\n';
      reported = true;
    }
    stats._previous_token_count = tokens;
    if (reported) {
      if (write_reported_records_to_cerr) {
        cerr << buffer << '\n';
      }
      if (new_fstream_reported.is_open()) {
        new_fstream_reported << buffer << '\n';
        if (!new_fstream_reported.good()) {
          cerr << "Error writing reported records\n";
          return false;
        }
      }
    } else if (new_fstream_unreported.is_open()) {
      new_fstream_unreported << buffer << '\n';
      if (!new_fstream_unreported.good()) {
        cerr << "Error writing unreported records\n";
        return false;
      }
    }
  }
  return input.good();
}

static void
ReportFile(const char* fname, const FileStatistics& stats) {
  const bool looks_like_R = header_record_has_only_different_token_count(
      stats._token_counts, stats._first_token_count);
  for (int i = 0; i < stats._token_counts.number_elements(); ++i) {
    const uint64_t count = stats._token_counts[i];
    if (count == 0) {
      continue;
    }
    if (human_readable_output) {
      cerr << iwstring::with_commas(count);
    } else {
      cerr << count;
    }
    cerr << " records had " << i << " tokens";
    if (i == stats._first_token_count && looks_like_R) {
      cerr << " HEADER";
    }
    cerr << '\n';
  }
  if (human_readable_output) {
    cerr << iwstring::with_commas(stats._records);
  } else {
    cerr << stats._records;
  }
  cerr << " records read from '" << fname << "'";
  if (verbose && stats._records != 0) {
    cerr << " records between " << stats._record_length.minval() << " and "
         << stats._record_length.maxval();
    if (stats._record_length.n() > 1) {
      cerr << " ave " << stats._record_length.average();
    }
  }
  cerr << '\n';
}

// Skip a partial initial record, but retain a complete record at the offset.
static bool
PositionNearEnd(iwstring_data_source& input) {
  const auto size = input.file_size();
  if (size < 0 || offset_from_eof > size) {
    cerr << "Requested offset from EOF " << offset_from_eof << " larger than file "
         << size << '\n';
    return false;
  }
  const off_t start = size - offset_from_eof;
  char previous = record_delimiter;
  if (start > 0 && (!input.seekg(start - 1) || input.read_bytes(&previous, 1) != 1)) {
    cerr << "Cannot inspect record boundary at offset " << start << '\n';
    return false;
  }
  if (!input.seekg(start)) {
    cerr << "Cannot seek to offset " << start << '\n';
    return false;
  }
  if (previous != record_delimiter) {
    const_IWSubstring discarded;
    input.next_record(discarded);
  }
  return input.good();
}

static int
tcount(const char* fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 1;
  }
  if (verbose) {
    cerr << fname << '\n';
  }
  input.set_record_delimiter(record_delimiter);
  input.set_dos(dos_record_terminators);
  if (offset_from_eof && !PositionNearEnd(input)) {
    return 1;
  }
  FileStatistics stats;
  const bool ok = ProcessInput(input, stats);
  total_lines_read += stats._records;
  if (!stats._rectangular) {
    ++non_rectangular_files_found;
  }
  ReportFile(fname, stats);
  if (!ok) {
    cerr << "Error processing '" << fname << "'\n";
    return 1;
  }
  return 0;
}

// Validate all destinations before opening any of them, including aliases.
static bool
CheckOutputNames(const Command_Line& cl) {
  std::vector<IWString> outputs;
  for (char flag : {'U', 'R', 'C', 'T'}) {
    if (flag != 'T' && cl.option_count(flag) > 1) {
      cerr << "Specify -" << flag << " only once\n";
      return false;
    }
    IWString name;
    for (int i = 0; cl.value(flag, name, i); ++i) {
      if (name.empty()) {
        cerr << "Empty output filename for -" << flag << '\n';
        return false;
      }
      for (const char* input : cl) {
        if (const_IWSubstring(input) == "-") {
          continue;
        }
        std::error_code error;
        const auto same = iwmisc::SameFile(name, IWString(input), error);
        if (!same) {
          cerr << "Cannot compare output '" << name << "' with input '" << input
               << "': " << error.message() << '\n';
          return false;
        }
        if (*same) {
          cerr << "Output '" << name << "' would overwrite input '" << input << "'\n";
          return false;
        }
      }
      for (const IWString& previous : outputs) {
        std::error_code error;
        const auto same = iwmisc::SameFile(name, previous, error);
        if (!same || *same) {
          cerr << "Output destinations '" << name << "' and '" << previous
               << "' conflict";
          if (error) {
            cerr << ": " << error.message();
          }
          cerr << '\n';
          return false;
        }
      }
      outputs.push_back(name);
    }
  }
  return true;
}

static bool
CloseOutput(std::ofstream& stream, const IWString& name) {
  if (!stream.is_open()) {
    return true;
  }
  stream.close();  // close also flushes buffered data; retain any failure status.
  if (!stream.good()) {
    cerr << "Error closing output '" << name << "'\n";
    return false;
  }
  return true;
}

static bool
FinishOutputs(const Command_Line& cl) {
  bool ok = true;
  if (!CloseOutput(new_fstream_reported, new_fname_reported)) {
    ok = false;
  }
  if (!CloseOutput(new_fstream_unreported, new_fname_unreported)) {
    ok = false;
  }
  if (stream_for_count.is_open() && !CloseOutput(stream_for_count, cl.string_value('C'))) {
    ok = false;
  }
  for (Stream_for_Columns* stream : cstreams) {
    if (!CloseOutput(*stream, stream->fname())) {
      ok = false;
    }
  }
  return ok;
}

static int
tcount(int argc, char **argv) {
  Command_Line cl(argc, argv, "c:d:wt:T:U:R:fbrvn:e:guqC:o:kh");

  if (cl.unrecognised_options_encountered()) {
    usage(2);
  }

  DisplayDashHelpIfRequested(cl, 'd', DisplayDelimiterOptions);
  DisplayDashHelpIfRequested(cl, 'c', DisplayRecordDelimiterOptions);
  verbose = cl.option_count('v');

  if (cl.option_present('d')) {
    IWString tmp;
    cl.value('d', tmp);
    if (tmp == "help") {
      DisplayDelimiterOptions(0);
    }
    if (!char_name_to_char(tmp) || tmp.length() != 1) {
      cerr << "Sorry, the delimiter character can be one character only, or one of a "
              "list of tokens.  see -d help for details\n";
      usage(13);
    }

    token_delimiter = tmp[0];
  }

  if (cl.option_present('c')) {
    IWString value = cl.string_value('c');
    if (value.length() > 1) {
      value.to_lowercase();
    }
    dos_record_terminators = false;
    if (value == "dos" || value == "windows" || value == "crlf") {
      // The reader's Windows mode consumes CRLF as one record terminator.
      value = '\n';
      dos_record_terminators = true;
    } else if (value == "newline" || value == "linux" || value == "unix" || value == "lf") {
      value = '\n';
    } else if (value == "^m" || value == "cr") {
      value = '\r';
    } else if (!char_name_to_char(value) || value.length() != 1) {
      cerr << "Invalid record delimiter; use '-c help' for choices\n";
      return 1;
    }
    record_delimiter = value[0];
  }

  if (cl.option_present('u')) {
    consecutive_delimiters_mean_empty_words_between = 1;
    if (verbose) {
      cerr << "Consecutive delimiters mean empty words in between adjacent delimiters\n";
    }
  }
  if (cl.option_present('q')) {
    data_is_quoted_tokens = 1;
    if (verbose) {
      cerr << "Data is quoted tokens\n";
    }
  }

  if (cl.option_present('h')) {
    human_readable_output = 1;
    if (verbose) {
      cerr << "Numbers of records written in human readable form\n";
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', records_to_process) || records_to_process <= 0) {
      cerr << "The -n option requires a positive whole number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will read only the first " << records_to_process
           << " records from each file\n";
    }
  }

  if (cl.option_present('e')) {
    uint64_t tmp;
    if (!cl.value('e', tmp) || tmp == 0 ||
        tmp > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
      cerr << "Must specify a value offset from eof\n";
      usage(8);
    }

    offset_from_eof = tmp;

    if (verbose) {
      cerr << "Will start processing " << offset_from_eof << " bytes from eof\n";
    }
  }

  if (cl.option_present('b')) {
    abort_on_different_token_count = 1;
    if (verbose) {
      cerr << "Will stop processing on a differing token count\n";
    }
  }

  if (cl.option_present('r')) {
    report_different_token_count = 1;
    if (verbose) {
      cerr << "Will report all instances where token count is different from previous "
              "record\n";
    }
  }

  if (cl.option_present('f')) {
    report_different_from_first_record = 1;

    if (verbose) {
      cerr << "Will report instances where token count is different from first record\n";
    }
  }

  if (cl.option_present('g')) {
    report_non_rectangular_files = 1;
    if (verbose) {
      cerr << "Will report non-rectangular files\n";
    }
  }

  if (cl.option_present('k')) {
    contains_backslashed_separators = 1;

    if (verbose) {
      cerr << "Input may contain backslashed token separators\n";
    }
  }

  if (cl.empty()) {
    usage(1);
  }
  if (!CheckOutputNames(cl)) {
    return 1;
  }
  token_options = {token_delimiter, data_is_quoted_tokens != 0,
                   contains_backslashed_separators != 0,
                   consecutive_delimiters_mean_empty_words_between != 0};

  if (cl.option_present('U')) {
    if (cl.option_present('b')) {
      cerr << "The -U and -b options are inconsistent\n";
      usage(18);
    }

    if (!cl.option_present('r') && !cl.option_present('f')) {
      cerr << "The -U option only makes sense with either the -r or -f options\n";
      usage(73);
    }

    cl.value('U', new_fname_unreported);

    new_fstream_unreported.open(new_fname_unreported.null_terminated_chars(),
                                std::ios::out);
    if (!new_fstream_unreported.good()) {
      cerr << "Cannot open '" << new_fname_unreported << "' for output\n";
      return 8;
    }

    if (verbose) {
      cerr << "Unreported records will be written to '" << new_fname_unreported << "'\n";
    }
  }

  if (cl.option_present('R')) {
    if (cl.option_present('b')) {
      cerr << "The -R and -b options are inconsistent\n";
      usage(13);
    }

    if (!cl.option_present('r') && !cl.option_present('f')) {
      cerr << "The -R option only makes sense with either the -r or -f options\n";
      usage(73);
    }

    cl.value('R', new_fname_reported);

    new_fstream_reported.open(new_fname_reported.null_terminated_chars(), std::ios::out);
    if (!new_fstream_reported.good()) {
      cerr << "Cannot open reported stream '" << new_fname_reported << "'\n";
      return 9;
    }

    if (verbose) {
      cerr << "Reported records will be written to '" << new_fname_reported << "'\n";
    }
  }

  if (cl.option_present('o')) {
    IWString o = cl.string_value('o');
    if (!char_name_to_char(o) || o.length() != 1) {
      cerr << "Unrecognised output token directive '" << o << "'\n";
      usage(1);
    }
    output_token_delimiter = o[0];
  }

  int nt = cl.option_count('t');
  if (nt != cl.option_count('T')) {
    cerr << "There must be the same number of -t and -T options\n";
    usage(34);
  }

  for (int i = 0; i < nt; i++) {
    int t;
    if (!cl.value('t', t, i) || t < 0) {
      cerr << "The -t option requires a whole non-negative number\n";
      usage(5);
    }

    IWString fname;
    (void)cl.value('T', fname, i);

    auto stream = std::make_unique<Stream_for_Columns>(t, fname.null_terminated_chars());
    if (!stream->good()) {
      return 9;
    }
    cstreams.add(stream.release());
  }

  if (cl.option_present('w')) {
    write_reported_records_to_cerr = 1;
    if (verbose) {
      cerr << "Will echo reported records to stderr\n";
    }
  }

  if (cl.option_present('C')) {
    const char *c = cl.option_value('C');

    stream_for_count.open(c, std::ios::out);

    if (!stream_for_count.good()) {
      cerr << "Cannot open stream for tcount '" << c << "'\n";
      return 1;
    }

    if (verbose) {
      cerr << "Data on column counts written to '" << c << "'\n";
    }

    stream_for_count << "Line Tcount\n";
  }

  int rc = 0;
  for (const char* fname: cl) {
    if (tcount(fname) != 0) {
      rc = 1;  // Keep processing other files, but never lose an earlier error.
    }
  }

  if (report_non_rectangular_files) {
    cerr << non_rectangular_files_found << " non rectangular files found\n";
  }

  if (verbose && cl.size() > 1) {
    cerr << "Read " << iwstring::with_commas(total_lines_read) << " lines from "
         << cl.size() << " files\n";
  }

  if (!FinishOutputs(cl)) {
    rc = 1;
  }
  return rc;
}

int
main(int argc, char **argv) {
  return tcount(argc, argv);
}
