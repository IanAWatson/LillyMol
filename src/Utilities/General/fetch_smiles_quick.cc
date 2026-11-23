/*
  Fetches records from a smiles file based on identifiers in another file
*/

#include <stdlib.h>

#include <fstream>
#include <limits>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int identifier_column_in_identifier_file = 0;

static int identifier_column_in_smiles_file = 1;

static std::unique_ptr<RE2> identifier_regexp_in_smiles_file;

static int ignore_duplicate_identifiers_in_identifier_file = 0;

static int duplicate_identifiers_in_identifier_file = 0;

static int ignore_column_not_present = 0;

static std::ofstream stream_for_not_in_identifier_file;
static std::ofstream stream_for_not_in_smiles_file;

static int write_stream_for_not_in_smiles_as_smiles = 0;

static int items_written = 0;
static int not_in_identifier_file = 0;

static int include_identifier_file_info = 1;

static IWString string_before_identifier_file_info(' ');

static int strip_leading_zeros_from_identifiers = 0;

static int erase_identifiers_written = 1;

static int block_size = 8192;

static int invert_fetching_operation = 0;

static int identifier_file_is_descriptor_file = 0;

static char identifier_file_column_separator = ' ';
static char smiles_file_column_separator = ' ';


// Useful for fetching things out of textproto files.
static int tokens_are_quoted = 0;

struct OutputFiles {
  IWString _selected;  // Empty in normal mode, where selected records use stdout.
  IWString _unrequested;
  IWString _missing;
};

static bool
FlushIfNeeded(IWString_and_File_Descriptor& output) {
  const bool ok = output.write_if_buffer_holds_more_than(block_size) != 0;
  if (!ok) {
    cerr << "Error writing selected records\n";
  }
  return ok;
}

static bool
FinishBufferedOutput(IWString_and_File_Descriptor& output, const char* name) {
  if (!output.is_open()) {
    return true;
  }
  if (!output.close()) {
    cerr << "Error closing output '" << name << "'\n";
    return false;
  }
  return true;
}

static bool
FinishSideOutput(std::ofstream& output, const char* label) {
  if (!output.is_open()) {
    return true;
  }
  output.close();
  if (!output.good()) {
    cerr << "Error closing " << label << " output\n";
    return false;
  }
  return true;
}

static bool
FinishOutputs(IWString_and_File_Descriptor& output) {
  bool ok = output.flush() != 0;  // stdout is borrowed, so do not close it here.
  if (!ok) {
    cerr << "Error flushing selected output\n";
  }
  if (!FinishSideOutput(stream_for_not_in_identifier_file, "unrequested (-X)")) {
    ok = false;
  }
  if (!FinishSideOutput(stream_for_not_in_smiles_file, "missing (-Y)")) {
    ok = false;
  }
  return ok;
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
  cerr << "Fetches records from one file based on identifiers in one or more other file(s)\n";
  cerr << prog_name << " identifier_file smiles_file > newfile\n";
  cerr << " -c <col>       identifier column in identifier file\n";
  cerr << " -C <col>       identifier column in smiles file\n";
  cerr << " -C RX=<rx>     identifier is whichever column(s) match <rx>\n";
  cerr << " -K <id>        one or more identifiers to fetch (multiple -K options)\n";
  cerr << " -d             ignore duplicate identifiers in identifier file\n";
  cerr << " -q             quietly ignore duplicate identifiers in identifier file\n";
  cerr << " -a             write all instances of identifiers in smiles file\n";
  cerr << "                by default, only the first is written\n";
  cerr << " -X <fname>     write smiles records not in <identifier_file> to <fname>\n";
  cerr << " -Y <fname>     write identifiers not in <smiles_file> to <fname>\n";
  cerr << " -w             write the -Y file as a smiles file (swap columns)\n";
  cerr << " -k             suppress addition of info from identifier file\n";
  cerr << " -n <string>    string to insert between record and info from identifier file\n";
  cerr << " -x             invert behaviour, selection becomes deselection\n";
  cerr << " -z             strip leading zero's from identifiers\n";
  cerr << " -j             identifier file is descriptor file, skip header record\n";
  cerr << " -b             stop processing identifier file on error, but continue processing (dangerous)\n";
  cerr << " -f             do not fail if a specified column is not present in either the identifier file or smiles file\n";
  cerr << " -i <char>      column separator in identifier file\n";
  cerr << " -I <char>      column separator in smiles file\n";
  cerr << " -S <stem>      first files are identifier files, last is haystack. Create many subsets\n";
  cerr << " -u <suffix>    suffix for -S files created\n";
  cerr << " -g <ndx>       start number for files created (-g 1 for dopattern)\n";
  cerr << " -e             tokens may be quoted\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
separate_into_id_and_data_2(const const_IWSubstring& buffer, int col, const char sep,
                            IWString& id, IWString& zdata) {
  if (0 == col) {  // the most common case
    //  cerr << "Extracting column " << col << " from '" << buffer << "'\n";
    int ispace = buffer.index(sep);

    if (ispace < 0) {  // just one token on the line
      id = buffer;
      return 1;
    }

    const char* s = buffer.rawchars();
    id.set(s, ispace);

    zdata.set(s + ispace + 1, buffer.length() - ispace - 1);

    //  cerr << "Split into '" << id << "' and '" << zdata << "'\n";

    return 1;
  }

  if (col < 0) {
    int nw;

    if (' ' == sep) {
      nw = buffer.nwords();
    } else {
      nw = buffer.nwords(sep);
    }

    if (-col > nw) {
      cerr << "Record contains only " << nw << " columns, cannot fetch column " << col
           << '\n';
      return 0;
    }

    col = nw + col;
  }

  int j = 0;
  const_IWSubstring token;

  for (int c = 0; buffer.NextWord(token, j, sep); ++c) {
    if (c == col) {
      id = token;
    } else if (include_identifier_file_info) {
      zdata.append_with_spacer(token);
    }
  }

  if (0 == id.length()) {
    if (ignore_column_not_present) {
      return 0;
    }

    cerr << "Cannot extract column " << col << " from '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

static int
separate_into_id_and_data(const const_IWSubstring& buffer, int col, const char sep,
                          IWString& id, IWString& zdata) {
  if (!separate_into_id_and_data_2(buffer, col, sep, id, zdata)) {
    return 0;
  }

  if (strip_leading_zeros_from_identifiers) {
    id.remove_leading_chars('0');
  }

  return 1;
}

static int
read_identifiers_to_fetch(const const_IWSubstring& buffer, int record_number,
                          IW_STL_Hash_Map_String& identifiers_to_fetch) {
  IWString id, zdata;
  if (separate_into_id_and_data(buffer, identifier_column_in_identifier_file,
                                identifier_file_column_separator, id, zdata)) {
    ;
  } else if (ignore_column_not_present) {
    return 1;
  } else {
    return 0;
  }

  IW_STL_Hash_Map_String::const_iterator f = identifiers_to_fetch.find(id);
  if (f == identifiers_to_fetch.end()) {
    identifiers_to_fetch[id] = zdata;
    return 1;
  }

  if (ignore_duplicate_identifiers_in_identifier_file) {
    if (1 == ignore_duplicate_identifiers_in_identifier_file) {
      cerr << "Ignoring duplicate identifer '" << id << "', line " << record_number
           << "\n";
    }
    duplicate_identifiers_in_identifier_file++;
    return 1;
  }

  cerr << "read_identifiers_to_fetch:duplicate identifier '" << id << "'\n";
  return 0;
}

static int
read_identifiers_to_fetch(iwstring_data_source& input,
                          IW_STL_Hash_Map_String& identifiers_to_fetch) {
  const_IWSubstring buffer;

  if (identifier_file_is_descriptor_file) {
    (void)input.next_record(buffer);
  }

  while (input.next_record(buffer)) {
    if (!read_identifiers_to_fetch(buffer, input.lines_read(), identifiers_to_fetch)) {
      cerr << "Fatal error processing identifier file '" << buffer << "', line "
           << input.lines_read() << '\n';
      return 0;
    }
  }

  if (!input.good()) {
    cerr << "Error reading identifier file '" << input.fname() << "'\n";
    return 0;
  }
  return !identifiers_to_fetch.empty();
}

// `translate_tabs` should be set when the caller has not specified a column
// separator for this file. Tabs are then converted to spaces, which makes a
// tab separated file work without any extra options. When a separator has been
// specified with -i the tabs in the file are significant and must be preserved.
static int
read_identifiers_to_fetch(const char* fname,
                          IW_STL_Hash_Map_String& identifiers_to_fetch,
                          int translate_tabs) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_translate_tabs(translate_tabs);

  return read_identifiers_to_fetch(input, identifiers_to_fetch);
}

static int
handle_record_not_in_identifier_file(const const_IWSubstring& buffer,
                                     const IWString& id) {
  if (verbose > 2) {
    cerr << "Smiles file identifier '" << id << "' not requested\n";
  }

  not_in_identifier_file++;

  if (stream_for_not_in_identifier_file.is_open()) {
    stream_for_not_in_identifier_file << buffer << '\n';
    if (!stream_for_not_in_identifier_file.good()) {
      cerr << "Error writing unrequested output (-X)\n";
      return 0;
    }
  }

  return 1;
}

template <typename T> void
Unquote(T& token) {
  token.remove_leading_chars(1);
  token.chop();
}

static int
fetch_smiles_quick(const const_IWSubstring& buffer, const IWString& id,
                   const IWString& zdata, IW_STL_Hash_Map_String& identifiers_to_fetch,
                   IWString& output_buffer) {
  IW_STL_Hash_Map_String::const_iterator f = identifiers_to_fetch.find(id);

  bool is_match = false;

  if (f != identifiers_to_fetch.end()) {
    is_match = true;
  } else if (tokens_are_quoted && id[0] == '"') {
    IWString tmp(id);
    Unquote(tmp);
    f = identifiers_to_fetch.find(tmp);
    if (f != identifiers_to_fetch.end()) {
      is_match = true;
    }
  }

  // if (is_match)
  //   cerr << "Got match for '" << id << "'\n";

  if (invert_fetching_operation) {
    if (is_match) {
      is_match = false;
      //    cerr << "Switched to non-match for '" << id << "'\n";
    } else {
      is_match = true;
    }
  }

  if (!is_match) {
    return handle_record_not_in_identifier_file(buffer, id);
  }

  output_buffer << buffer;

  if (invert_fetching_operation) {  // the iterator may be bad
    ;
  } else if (include_identifier_file_info) {
    if (!f->second.empty()) {
      output_buffer << string_before_identifier_file_info << (*f).second;
    }
  }

  output_buffer << '\n';

  items_written++;

  if (erase_identifiers_written) {
    identifiers_to_fetch.erase(id);
  }

  return 1;
}

static int
fetch_smiles_quick(const const_IWSubstring& buffer,
                   IW_STL_Hash_Map_String& identifiers_to_fetch,
                   IWString& output_buffer) {
  IWString id, zdata;

  if (identifier_regexp_in_smiles_file) {
    int i = 0;
    const_IWSubstring token;
    int col = 0;
    while (buffer.nextword(token, i)) {
      if (!iwre2::RE2PartialMatch(token, *identifier_regexp_in_smiles_file)) {
        continue;
      }

      if (identifiers_to_fetch.contains(token)) {
        ;  // good
      } else if (tokens_are_quoted && token[0] == '"') {
        Unquote(token);
        cerr << "unquoted to '" << token << "'\n";
      } else {
        continue;
      }

      id = token;
      zdata = buffer;
      zdata.remove_word(col);

      return fetch_smiles_quick(buffer, id, zdata, identifiers_to_fetch, output_buffer);
    }

    return handle_record_not_in_identifier_file(buffer, "");
  } else if (!separate_into_id_and_data(buffer, identifier_column_in_smiles_file,
                                        smiles_file_column_separator, id, zdata)) {
    if (ignore_column_not_present) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return fetch_smiles_quick(buffer, id, zdata, identifiers_to_fetch, output_buffer);
  }
}

static int
fetch_smiles_quick(iwstring_data_source& input,
                   IW_STL_Hash_Map_String& identifiers_to_fetch,
                   IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!fetch_smiles_quick(buffer, identifiers_to_fetch, output)) {
      cerr << "Fatal error processing '" << buffer << "', line " << input.lines_read()
           << '\n';
      return 0;
    }

    if (!FlushIfNeeded(output)) {
      return 0;
    }

    if (0 == identifiers_to_fetch.size()) {
      return 1;
    }
  }

  if (!input.good()) {
    cerr << "Error reading smiles file '" << input.fname() << "'\n";
    return 0;
  }
  return 1;
}

static size_t
fetch_smiles_quick_multiple_files(iwstring_data_source& input,
                                  const IW_STL_Hash_Map_String& file_contents,
                                  IWString_and_File_Descriptor& output,
                                  IWString_and_File_Descriptor& xout,
                                  IWString_and_File_Descriptor& yout) {
  IW_STL_Hash_Set written;

  bool fill_written_hash = false;

  if (xout.is_open()) {
    fill_written_hash = true;
  }

  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    IWString id, zdata;
    if (!separate_into_id_and_data(buffer, identifier_column_in_identifier_file,
                                   identifier_file_column_separator, id, zdata)) {
      return 0;
    }

    const auto f = file_contents.find(id);
    if (f != file_contents.end()) {
      output << (*f).second;
      if (!FlushIfNeeded(output)) {
        return 0;
      }
      if (fill_written_hash) {
        written.insert(id);
      }
      continue;
    }

    //  The identifier was not in the haystack

    if (yout.is_open()) {
      yout << id << '\n';
      if (!yout.write_if_buffer_holds_more_than(block_size)) {
        cerr << "Error writing subset missing output\n";
        return 0;
      }
    }
  }

  if (!input.good()) {
    cerr << "Error reading subset identifier file '" << input.fname() << "'\n";
    return 0;
  }
  if (xout.is_open()) {
    for (const auto& c : file_contents) {
      if (written.contains(c.first)) {
        continue;
      }

      xout << c.second;
      if (!xout.write_if_buffer_holds_more_than(block_size)) {
        cerr << "Error writing subset unrequested output\n";
        return 0;
      }
    }
  }

  return 1;
}

static size_t
fetch_smiles_quick_multiple_files(const char* ifile,
                                  const IW_STL_Hash_Map_String& file_contents,
                                  const char* ofile, const char* xfile,
                                  const char* yfile) {
  iwstring_data_source input(ifile);

  if (!input.good()) {
    cerr << "fetch_smiles_quick_multiple_files:cannot open input file '" << ifile
         << "'\n";
    return 0;
  }

  IWString_and_File_Descriptor output;

  if (!output.open(ofile)) {
    cerr << "fetch_smiles_quick_multiple_files:cannot open output file '" << ofile
         << "'\n";
    return 0;
  }

  IWString_and_File_Descriptor xout, yout;

  if (0 == ::strlen(xfile)) {
    ;
  } else if (xout.open(xfile)) {
    ;
  } else {
    cerr << "fetch_smiles_quick_multiple_files:cannot open xfile '" << xfile << "'\n";
    return 0;
  }

  if (0 == ::strlen(yfile)) {
    ;
  } else if (yout.open(yfile)) {
    ;
  } else {
    cerr << "fetch_smiles_quick_multiple_files:cannot open yfile '" << yfile << "'\n";
    return 0;
  }

  bool ok = fetch_smiles_quick_multiple_files(input, file_contents, output, xout, yout);
  if (!FinishBufferedOutput(output, ofile)) {
    ok = false;
  }
  if (!FinishBufferedOutput(xout, xfile)) {
    ok = false;
  }
  if (!FinishBufferedOutput(yout, yfile)) {
    ok = false;
  }
  return ok;
}

static int
read_file_contents(iwstring_data_source& input, IW_STL_Hash_Map_String& file_contents) {
  IWString buffer;
  while (input.next_record(buffer)) {
    IWString id;
    int i = 0;
    if (!buffer.nextword(id, i) || !buffer.nextword(id, i)) {
      cerr << "read_file_contents:cannot read identifier from '" << buffer << "'\n";
      return 0;
    }

    if (!file_contents.contains(id)) {
      buffer << '\n';
      file_contents[id] = buffer;
    } else if (ignore_duplicate_identifiers_in_identifier_file) {
      if (1 == ignore_duplicate_identifiers_in_identifier_file) {
        cerr << "Ignoring duplicate identifer '" << id << "', line " << input.lines_read()
             << "\n";
      }
      duplicate_identifiers_in_identifier_file++;
      continue;
    } else {
      cerr << "read_file_contents:duplicate identifier '" << buffer << "'\n";
      return 0;
    }
  }

  if (!input.good()) {
    cerr << "Error reading subset haystack '" << input.fname() << "'\n";
    return 0;
  }
  return !file_contents.empty();
}

static int
read_file_contents(const char* fname, IW_STL_Hash_Map_String& file_contents) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "read_file_contents:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_file_contents(input, file_contents);
}

static int
fetch_smiles_quick_multiple_files(const Command_Line& cl,
                                  const std::vector<OutputFiles>& outputs) {
  IW_STL_Hash_Map_String file_contents;
  if (!read_file_contents(cl.last_item(), file_contents)) {
    cerr << "fetch_smiles_quick_multiple_files:cannot read file contents '"
         << cl.last_item() << "'\n";
    return 0;
  }
  for (int i = 0; i < cl.number_elements() - 1; ++i) {
    const OutputFiles& files = outputs[i];
    if (!fetch_smiles_quick_multiple_files(cl[i], file_contents,
            files._selected.AsString().c_str(), files._unrequested.AsString().c_str(),
            files._missing.AsString().c_str())) {
      cerr << "Cannot process fetching from '" << cl[i] << "'\n";
      return 0;
    }
  }
  return 1;
}

// See the comment on read_identifiers_to_fetch for the meaning of
// `translate_tabs`. Note that because this was previously unconditional, a tab
// separated smiles file has always worked without needing -I.
static int
fetch_smiles_quick(const char* fname, IW_STL_Hash_Map_String& identifiers_to_fetch,
                   IWString_and_File_Descriptor& output, int translate_tabs) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  input.set_dos(1);
  input.set_translate_tabs(translate_tabs);

  if (verbose > 0) {
    cerr << "Processing '" << fname << "'\n";
  }

  return fetch_smiles_quick(input, identifiers_to_fetch, output);
}

static int
determine_inter_column_separator(const Command_Line& cl, const char flag, char& sep) {
  IWString i = cl.string_value(flag);

  if (!char_name_to_char(i)) {
    cerr << "Unrecognised column separator '" << i << "'\n";
    return 0;
  }

  sep = i[0];

  return 1;
}

// Check all input and output identities before opening any output. -F is an
// input too, even though it is not one of the positional arguments.
static bool
CheckOutputNames(const Command_Line& cl, const std::vector<OutputFiles>& outputs) {
  std::vector<IWString> inputs;
  for (const char* input : cl) {
    if (const_IWSubstring(input) != "-") {
      inputs.emplace_back(input);
    }
  }
  if (cl.option_present('F') && cl.string_value('F') != "-") {
    inputs.push_back(cl.string_value('F'));
  }
  // Text readers expand environment variables; gzip readers use literal names.
  for (IWString& input : inputs) {
    if (!input.ends_with(".gz")) {
      if (auto expanded = input.ExpandEnvironmentVariables()) {
        input = *expanded;
      }
    }
  }
  std::vector<IWString> previous;
  for (const OutputFiles& files : outputs) {
    for (const IWString* name : {&files._selected, &files._unrequested, &files._missing}) {
      if (name->empty()) {  // An optional output was not requested.
        continue;
      }
      IWString path(*name);
      // Match buffered text output filename handling. Gzip outputs and
      // normal-mode ofstreams use literal names instead.
      if (!files._selected.empty() && !path.ends_with(".gz")) {
        if (path.length() > 2 && path.starts_with(">>")) {
          path.remove_leading_chars(2);
        }
        if (auto expanded = path.ExpandEnvironmentVariables()) {
          path = *expanded;
        }
      }
      for (const IWString& input : inputs) {
        std::error_code error;
        const auto same = iwmisc::SameFile(path, input, error);
        if (!same) {
          cerr << "Cannot compare output '" << *name << "' with input '" << input
               << "': " << error.message() << '\n';
          return false;
        }
        if (*same) {
          cerr << "Cannot overwrite input '" << input << "' with output '" << *name << "'\n";
          return false;
        }
      }
      for (const IWString& other : previous) {
        std::error_code error;
        const auto same = iwmisc::SameFile(path, other, error);
        if (!same || *same) {
          cerr << "Conflicting output files '" << *name << "' and '" << other << "'";
          if (error) {
            cerr << ": " << error.message();
          }
          cerr << '\n';
          return false;
        }
      }
      previous.push_back(path);
    }
  }
  return true;
}

static int
fetch_smiles_quick(int argc, char** argv) {
  Command_Line cl(argc, argv, "vc:C:X:Y:wkdzB:an:xF:qbjK:S:u:g:fi:I:e");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.option_present('a')) {
    erase_identifiers_written = 0;

    if (cl.option_present('Y')) {
      cerr << "Sorry, the -a and -Y options do not work together\n";
      return 6;
    }

    if (verbose) {
      cerr << "Will write all instances of identifiers in the smiles file\n";
    }
  }

  // Tabs are translated to spaces only when the column separator for that file
  // has not been specified. Specifying -i or -I means the tabs in that file are
  // significant, and translating them would defeat the option.

  int translate_tabs_in_identifier_file = 1;
  int translate_tabs_in_smiles_file = 1;

  if (cl.option_present('i')) {
    if (!determine_inter_column_separator(cl, 'i', identifier_file_column_separator)) {
      return 1;
    }
    translate_tabs_in_identifier_file = 0;
  }

  if (cl.option_present('I')) {
    if (!determine_inter_column_separator(cl, 'I', smiles_file_column_separator)) {
      return 1;
    }
    translate_tabs_in_smiles_file = 0;
  }

  if (cl.option_present('e')) {
    tokens_are_quoted = 1;
    if (verbose) {
      cerr << "Will strip quotes from tokens\n";
    }
  }

  if (cl.option_present('n')) {
    if (cl.option_present('k')) {
      cerr << "The -n and -k options are mutually exclusive\n";
      usage(5);
    }

    cl.value('n', string_before_identifier_file_info);

    if (verbose) {
      cerr << "Will put '" << string_before_identifier_file_info
           << " before identifier file info\n";
    }
  }

  if (cl.option_present('B')) {
    if (!cl.value('B', block_size) || block_size < 1) {
      cerr << "INvalid block size (-B option)\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Will use " << block_size << " as the output block size\n";
    }
  }

  if (cl.option_present('d')) {
    ignore_duplicate_identifiers_in_identifier_file = 1;

    if (verbose) {
      cerr << "Will ignore duplicates in identifier file\n";
    }
  }

  if (cl.option_present('q')) {
    ignore_duplicate_identifiers_in_identifier_file = 2;

    if (verbose) {
      cerr << "Will quietly ignore duplicates in identifier file\n";
    }
  }

  if (cl.option_present('f')) {
    ignore_column_not_present = cl.option_count('f');

    if (verbose) {
      cerr << "Will ignore records (in both identifier file and target file) where not "
              "enough columns\n";
    }
  }

  if (cl.option_present('z')) {
    strip_leading_zeros_from_identifiers = 1;

    if (verbose) {
      cerr << "Will strip leading zero's from identifiers\n";
    }
  }

  if (cl.option_present('c')) {
    if (!cl.value('c', identifier_column_in_identifier_file) ||
        0 == identifier_column_in_identifier_file) {
      cerr << "The column in identifier file option (-c) must be a valid column number\n";
      usage(5);
    }

    if (verbose) {
      cerr << "Identifiers in identifier file in column "
           << identifier_column_in_identifier_file << '\n';
    }

    if (identifier_column_in_identifier_file > 0) {
      identifier_column_in_identifier_file--;
    }
  }

  if (cl.option_present('C')) {
    const_IWSubstring c = cl.string_value('C');

    if (c.starts_with("RX=")) {
      c.remove_leading_chars(3);
      if (!iwre2::RE2Reset(identifier_regexp_in_smiles_file, c)) {
        cerr << "Invalid smiles file identifier regexp '" << c << "'\n";
        return 5;
      }

      if (verbose) {
        cerr << "Identifiers in smiles file must match '"
             << identifier_regexp_in_smiles_file->pattern() << "'\n";
      }

      identifier_column_in_smiles_file = -1;
    }
    //  else if (! cl.value('C', identifier_column_in_smiles_file) ||
    //  identifier_column_in_smiles_file < 1)
    else if (!cl.value('C', identifier_column_in_smiles_file)) {
      cerr << "The column in smiles file option (-C) must be a whole +ve number\n";
      usage(5);
    } else if (verbose) {
      cerr << "Identifiers in smiles file in column " << identifier_column_in_smiles_file
           << '\n';
    }

    if (identifier_column_in_smiles_file > 0) {
      identifier_column_in_smiles_file--;
    }
  }

  if (cl.option_present('k')) {
    include_identifier_file_info = 0;

    if (verbose) {
      cerr << "Suppress addition of identifier file info\n";
    }
  }

  int fail_if_read_identifiers_to_fetch_fails = 1;

  if (cl.option_present('b')) {
    fail_if_read_identifiers_to_fetch_fails = 0;

    if (verbose) {
      cerr << "Will stop processing identifier file on error, and continue with "
              "identifiers found so far\n";
    }
  }

  if (cl.option_present('j')) {
    identifier_file_is_descriptor_file = 1;

    if (verbose) {
      cerr << "Will skip header record in identifier file\n";
    }
  }

  IW_STL_Hash_Map_String identifiers_to_fetch;

  identifiers_to_fetch.rehash(1000000);

  // Where do we start scanning the command line

  int clstart = 1;

  if (cl.option_present('F')) {  // nobody uses this any more
    const char* f = cl.option_value('F');
    if (!read_identifiers_to_fetch(f, identifiers_to_fetch,
                                   translate_tabs_in_identifier_file)) {
      cerr << "Cannot read identifiers to fetch from '" << f << "'\n";
      return 6;
    }

    clstart = 0;
  } else if (cl.option_present('K')) {
    IWString k;
    for (int i = 0; cl.value('K', k, i); ++i) {
      int j = 0;
      IWString token;
      while (k.nextword(token, j, '.')) {
        identifiers_to_fetch[token] = "";
      }
    }
    clstart = 0;
  } else if (cl.number_elements() < 2) {
    cerr << "Must specify at least two files\n";
    usage(2);
  } else if (cl.option_present('S')) {
    ;
  } else if (read_identifiers_to_fetch(cl[0], identifiers_to_fetch,
                                     translate_tabs_in_identifier_file)) {
    ;
  } else if (fail_if_read_identifiers_to_fetch_fails ||
             0 == identifiers_to_fetch.size()) {
    cerr << "Cannot read identifiers to fetch from '" << cl[0] << "'\n";
    return 6;
  } else {
    cerr << "Warning, possibly incomplete data read from '" << cl[0]
         << "', allowed by -b option\n";
  }

  if (verbose) {
    cerr << "Will fetch " << identifiers_to_fetch.size() << " identifiers\n";
  }

  if (cl.option_present('x')) {
    invert_fetching_operation = 1;

    if (verbose) {
      cerr << "Will discard rather than fetch identifiers\n";
    }
  }

  if (cl.option_present('w')) {
    write_stream_for_not_in_smiles_as_smiles = 1;

    if (verbose) {
      cerr << "Will write stream for not in smiles file as smiles\n";
    }
  }

  for (char flag : {'X', 'Y', 'S'}) {
    if (cl.option_count(flag) > 1 ||
        (cl.option_present(flag) && cl.string_value(flag).empty())) {
      cerr << "Specify one non-empty filename or stem for -" << flag << '\n';
      return 1;
    }
  }
  if (cl.option_present('S')) {
    if (cl.size() < 2) {
      cerr << "Subset mode requires an identifier file and a haystack\n";
      return 1;
    }
    const IWString stem = cl.string_value('S');
    const IWString xstem = cl.string_value('X');
    const IWString ystem = cl.string_value('Y');
    const IWString suffix = cl.string_value('u');
    int seq_begin = 0;
    if (cl.option_present('g') && (!cl.value('g', seq_begin) || seq_begin < 0)) {
      cerr << "The starting subset index (-g) must be a non-negative integer\n";
      return 1;
    }
    if (cl.number_elements() - 2 > std::numeric_limits<int>::max() - seq_begin) {
      cerr << "Subset file indices exceed the supported range\n";
      return 1;
    }
    std::vector<OutputFiles> outputs(cl.number_elements() - 1);
    for (int i = 0; i < cl.number_elements() - 1; ++i) {
      outputs[i]._selected << stem << (seq_begin + i) << suffix;
      if (!xstem.empty()) {
        outputs[i]._unrequested << xstem << (seq_begin + i);
      }
      if (!ystem.empty()) {
        outputs[i]._missing << ystem << (seq_begin + i);
      }
    }
    if (!CheckOutputNames(cl, outputs)) {
      return 1;
    }
    return fetch_smiles_quick_multiple_files(cl, outputs) ? 0 : 2;
  }

  OutputFiles files;
  if (cl.option_present('X')) {
    files._unrequested = cl.string_value('X');
  }
  if (cl.option_present('Y')) {
    files._missing = cl.string_value('Y');
    if (write_stream_for_not_in_smiles_as_smiles) {
      files._missing.EnsureEndsWith(".smi");
    }
  }
  if (!CheckOutputNames(cl, {files})) {
    return 1;
  }
  if (!files._unrequested.empty()) {
    stream_for_not_in_identifier_file.open(files._unrequested.AsString());
    if (!stream_for_not_in_identifier_file.good()) {
      cerr << "Cannot open unrequested output '" << files._unrequested << "'\n";
      return 4;
    }
  }
  if (!files._missing.empty()) {
    stream_for_not_in_smiles_file.open(files._missing.AsString());
    if (!stream_for_not_in_smiles_file.good()) {
      cerr << "Cannot open missing output '" << files._missing << "'\n";
      return 4;
    }
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;
  for (int i = clstart; i < cl.number_elements(); i++) {
    if (!fetch_smiles_quick(cl[i], identifiers_to_fetch, output,
                            translate_tabs_in_smiles_file)) {
      rc = 1;
      break;
    }
  }

  if (0 == items_written) {
    cerr << "Warning, nothing written!\n";
  }

  if (verbose) {
    cerr << items_written << " items written\n";
    if (not_in_identifier_file) {
      cerr << not_in_identifier_file << " items not in identifier file\n";
    }
    if (duplicate_identifiers_in_identifier_file) {
      cerr << "ignored " << duplicate_identifiers_in_identifier_file
           << " duplicate identifiers in identifier file\n";
    }
  }

  if (stream_for_not_in_smiles_file.is_open() && identifiers_to_fetch.size()) {
    if (verbose) {
      cerr << identifiers_to_fetch.size() << " identifiers not in smiles file\n";
    }

    for (IW_STL_Hash_Map_String::const_iterator i = identifiers_to_fetch.begin();
         i != identifiers_to_fetch.end(); ++i) {
      if (write_stream_for_not_in_smiles_as_smiles) {
        if ((*i).second.length()) {
          stream_for_not_in_smiles_file << (*i).second << ' ' << (*i).first;
        } else {
          stream_for_not_in_smiles_file << (*i).first;
        }
      } else {
        stream_for_not_in_smiles_file << (*i).first;
        if ((*i).second.length()) {
          stream_for_not_in_smiles_file << ' ' << (*i).second;
        }
      }
      stream_for_not_in_smiles_file << '\n';
      if (!stream_for_not_in_smiles_file.good()) {
        cerr << "Error writing missing output (-Y)\n";
        rc = 1;
        break;
      }
    }
  }

  if (!FinishOutputs(output)) {
    rc = 1;
  }
  return rc;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = fetch_smiles_quick(argc, argv);

  return rc;
}
