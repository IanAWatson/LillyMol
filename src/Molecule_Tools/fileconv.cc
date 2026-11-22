/*
  File conversion utility.
*/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <memory>
#include <limits>
#include <assert.h>


using std::cerr;

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/output.h"

#include "fileconv_opts.h"

static const char *prog_name;

namespace fileconv {

// When fileconv is processing, there can be a number of output streams.
struct LocalOptions {
  int verbose = 0;
  uint64_t molecules_read = 0;
  uint64_t molecules_changed = 0;
  uint64_t molecules_written = 0;
  int audit_input = 0;
  int debug_print_each_molecule = 0;

  // If a -S option is specified.
  IWString output_file_stem;

  // The main output stream - the -S option, or a per file stream.
  Molecule_Output_Object output;

  // Log messages about rejections.
  IWString_and_File_Descriptor reject_log;

  // Rejected molecules
  Molecule_Output_Object rejections_output;

  // Smiles before any filtering.
  // As part of the smiles infrastructure, I wanted to simultaneously
  // produce llyg.smi and llygO.smi. This allows that...
  IWString_and_File_Descriptor stream_for_smiles_before_filters;

  int connection_table_errors_allowed = 0;

  // Optionally we can write connection table errors to a destination.
  IWString connection_table_error_file;

  Report_Progress report_progress;

  // functions
  private:
    std::vector<IWString> _input_files;
    std::vector<IWString> _molecule_files;
    std::vector<IWString> _output_files;
    IWString _rejection_stem;
    IWString _rejection_log_name;
    IWString _prefiltered_smiles_file;
    uint64_t _main_files_created = 0;
    uint64_t _rejection_files_created = 0;

    int GatherInputFiles(const Command_Line& cl);
    int ReserveOutputFile(const IWString& fname, bool remember = true);
    int ReserveMoleculeOutputs(const IWString& stem,
                              const Molecule_Output_Object& output, int keep_suffix, bool remember = true);
    int PrepareDynamicOutput(Molecule& m, Molecule_Output_Object& output, uint64_t& count);
    int OpenAuxiliaryOutputs();

    // If requested, write `smiles` to stream_for_smiles_before_filters
    int MaybeWritePostTransformSmiles(const IWString& smiles);

    // If there is a rejection stream open, write `m`.
    int MaybeWriteRejection(Molecule& m,
                             const IWString& rejection_reason);

    // Parse the -L option.
    int SetupRejectionsFile(Command_Line& cl, char flag);

    // Parse the -Y B4F= option.
    int OpenPreFilteredSmilesFile(Command_Line& cl, char flag);

    void SetDebugPrintEachMolecule(Command_Line& cl, char flag);

    int ParseBadHandlingoptions(Command_Line& cl, char flag);

    int SetupReporting(Command_Line& cl, char flag);

  public:
    // Initialise settings from command line options.
    int Build(Command_Line& cl);

    // Set verbosity, connection table errors... on a newly opened
    // input file.
    void SetupInputFile(data_source_and_type<Molecule>& input) const;

    // Does it seem likely that we will be able to open an output stream
    // for processing input from `fname`.
    int OutputFileSeemsOk(const char * fname) const;

    // If needed, open the output stream needed for processing `fname`.
    int MaybeOpenNewOutputFIle(const char * fname, int first_call);

    // A successful processing. If active, write it.
    int MaybeWrite(Molecule& m);

    // The most commonly used method.
    int Fileconv(Molecule& m, fileconv::FileconvConfig& config);

    int ReportResults(const Command_Line& cl, std::ostream& output) const;
    int FinishOutputs();
};

int
LocalOptions::GatherInputFiles(const Command_Line& cl) {
  for (const char* fname : cl) {
    _input_files.emplace_back(fname);
    const_IWSubstring arg(fname);
    if (!arg.starts_with("F:")) {
      _molecule_files.emplace_back(fname);
      continue;
    }
    arg.remove_leading_chars(2);
    _input_files.emplace_back(arg);
    iwstring_data_source input(arg);
    if (!input.good()) {
      cerr << "Cannot open input file list '" << arg << "'\n";
      return 0;
    }
    input.set_strip_trailing_blanks(1);
    input.set_strip_leading_blanks(1);
    input.set_skip_blank_lines(1);
    IWString line;
    while (input.next_record(line)) {
      if (line.starts_with('#')) {
        continue;
      }
      _input_files.push_back(line);
      _molecule_files.push_back(line);
    }
  }
  return 1;
}

int
LocalOptions::ReserveOutputFile(const IWString& fname, bool remember) {
  if (fname == "-") {
    return 1;
  }
  if (fname.empty()) {
    cerr << "Empty output filename\n";
    return 0;
  }
  for (const IWString& input : _input_files) {
    if (input != "-") {
      const auto same = iwmisc::SameFile(fname, input);
      if (!same) {
        cerr << "Cannot compare file paths\n";
        return 0;
      }
      if (*same) {
        cerr << "Output '" << fname << "' would overwrite input '" << input << "'\n";
        return 0;
      }
    }
  }
  for (const IWString& output : _output_files) {
    const auto same = iwmisc::SameFile(fname, output);
    if (!same) {
      cerr << "Cannot compare file paths\n";
      return 0;
    }
    if (*same) {
      cerr << "Output destinations '" << fname << "' and '" << output << "' collide\n";
      return 0;
    }
  }
  if (remember) {
    _output_files.push_back(fname);
  }
  return 1;
}

int
LocalOptions::ReserveMoleculeOutputs(const IWString& stem,
                                   const Molecule_Output_Object& output, int keep_suffix, bool remember) {
  if (stem == "-") {
    return 1;
  }
  for (FileType type : output.output_types()) {
    IWString fname;
    if (!create_file_with_appropriate_name(stem, fname, type, keep_suffix) ||
        !ReserveOutputFile(fname, remember)) {
      return 0;
    }
  }
  return 1;
}

int
LocalOptions::OpenAuxiliaryOutputs() {
  if (!_rejection_stem.empty() && !rejections_output.new_stem(_rejection_stem, 1)) {
    cerr << "Cannot open rejected molecule output '" << _rejection_stem << "'\n";
    return 0;
  }
  if (!_rejection_log_name.empty() && !reject_log.open(_rejection_log_name.null_terminated_chars())) {
    cerr << "Cannot open rejection log '" << _rejection_log_name << "'\n";
    return 0;
  }
  if (!_prefiltered_smiles_file.empty() &&
      !stream_for_smiles_before_filters.open(_prefiltered_smiles_file.null_terminated_chars())) {
    cerr << "Cannot open pre-filtered smiles '" << _prefiltered_smiles_file << "'\n";
    return 0;
  }
  return 1;
}

// Flush before streams are destroyed so failures affect the exit status.
static int
FlushMolecules(Molecule_Output_Object& output) {
  int ok = 1;
  for (auto* stream : output) {
    stream->flush();
    if (!stream->good()) {
      ok = 0;
    }
  }
  return ok;
}

int
LocalOptions::PrepareDynamicOutput(Molecule& m, Molecule_Output_Object& output,
                                   uint64_t& count) {
  IWString stem;
  if (output.name_token_for_file_name() >= 0) {
    if (!m.name().word(output.name_token_for_file_name(), stem)) {
      cerr << "Cannot obtain output filename from '" << m.name() << "'\n";
      return 0;
    }
  } else {
    if (output.molecules_written() % output.molecules_per_file() != 0) {
      return 1;
    }
    stem = output.stem();
    stem << count;
  }
  // Generated names are not accumulated: checking them must not become quadratic
  // in the number of molecules. Check against inputs and named output streams.
  if (!ReserveMoleculeOutputs(stem, output, 0, false) || !FlushMolecules(output)) {
    return 0;
  }
  ++count;
  return 1;
}

int
LocalOptions::FinishOutputs() {
  int ok = FlushMolecules(output);
  if (!FlushMolecules(rejections_output)) {
    ok = 0;
  }
  if (reject_log.is_open() && !reject_log.close()) {
    ok = 0;
  }
  if (stream_for_smiles_before_filters.is_open() && !stream_for_smiles_before_filters.close()) {
    ok = 0;
  }
  std::cout.flush();
  if (!std::cout.good()) {
    ok = 0;
  }
  if (!ok) {
    cerr << "fileconv: error finalising output\n";
  }
  return ok;
}

int
LocalOptions::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (!GatherInputFiles(cl)) {
    return 0;
  }

  if (! SetupRejectionsFile(cl, 'L')) {
    cerr << "Cannot initialise rejection files (-L)\n";
    return 0;
  }

  if (! OpenPreFilteredSmilesFile(cl, 'Y')) {
    return 0;
  }

  SetDebugPrintEachMolecule(cl, 'Y');

  if (!SetupReporting(cl, 'Y')) {
    return 0;
  }

  audit_input = cl.option_present('a');

  if (! ParseBadHandlingoptions(cl, 'B')) {
    return 0;
  }

  if (!connection_table_error_file.empty() && !ReserveOutputFile(connection_table_error_file)) {
    return 0;
  }

  // If just auditing the input, no output to be set up.
  if (audit_input) {
    return OpenAuxiliaryOutputs();
  }

  // Setup output options.
  if (cl.option_present('S')) {
    if (1 != cl.option_count('S')) {
      cerr << "Sorry, only one -S option possible\n";
      return 0;
    }

    cl.value('S', output_file_stem);
    if (output_file_stem.empty()) {
      cerr << "The -S output stem cannot be empty\n";
      return 0;
    }
    if (verbose)
      cerr << "New files will be created with stem '" << output_file_stem << "'\n";
  }

  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);

    if (verbose)
      cerr << "Default output type smiles\n";
  } else if (! output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output types\n";
    return 0;
  }

  if (cl.option_present('G')) {
    output.set_molecules_per_file(1);
    if (verbose)
      cerr << "Each molecule written to its own file (-G option)\n";
  }

  if (output.name_token_for_file_name() >= 0 || output.molecules_per_file() > 0) {
    // Names are checked just before writing each generated file.
  } else if (!output_file_stem.empty()) {
    if (!ReserveMoleculeOutputs(output_file_stem, output, 0)) {
      return 0;
    }
  } else if (output.name_token_for_file_name() < 0) {
    for (const IWString& input : _molecule_files) {
      if (!ReserveMoleculeOutputs(input, output, 0)) {
        return 0;
      }
    }
  }
  return OpenAuxiliaryOutputs();
}

void
DisplayDashBOptions(int rc) {
  std::cerr << " -B <nn>        ignore as many as <nn> otherwise fatal input errors\n";
  std::cerr << " -B log=<fname> echo (unchanged) rejected input records to <fname>\n";

  std::exit(rc);
}

int
LocalOptions::ParseBadHandlingoptions(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  const_IWSubstring b;
  for (int i = 0; cl.value(flag, b, i); ++i) {
    if (b.starts_with("log=")) {
      b.remove_leading_chars(4);
      connection_table_error_file = b;

      if (verbose)
        cerr << "Will log connection table errors to '" << connection_table_error_file << "'\n";
    } else if ("help" == b) {
      DisplayDashBOptions(0);
    } else if (b.numeric_value(connection_table_errors_allowed) &&
               connection_table_errors_allowed >= 0) {
      if (verbose)
        cerr << connection_table_errors_allowed << " connection table errors allowed\n";
    } else {
      cerr << "Invalid -B qualifier '" << b << "'\n";
      return 0;
    }
  }

  return 1;
}

void
LocalOptions::SetDebugPrintEachMolecule(Command_Line& cl, char flag) {
  if (! cl.option_present(flag)) {
    return;
  }
  const_IWSubstring y;
  for (int i = 0; cl.value(flag, y, i); ++i) {
    if (y == "dbg") {
      debug_print_each_molecule = 1;
      return;
    }
  }
}

int
LocalOptions::SetupRejectionsFile(Command_Line& cl, char flag) {
  if (!cl.option_present(flag)) {
    return 1;
  }

  if (!cl.option_present('o'))
    rejections_output.add_output_type(FILE_TYPE_SMI);
  else if (!rejections_output.determine_output_types(cl)) {
    cerr << "Cannot discern output types for rejections file\n";
    return 0;
  }

  _rejection_stem = cl.option_value(flag);
  if (_rejection_stem.empty()) {
    cerr << "The rejection output stem cannot be empty\n";
    return 0;
  }
  if (!ReserveMoleculeOutputs(_rejection_stem, rejections_output, 1)) {
    return 0;
  }
  // Keep the historical log name unless it is also a molecule output filename.
  // With -L -, only rejected molecules are written to stdout.
  if (_rejection_stem != "-") {
    _rejection_log_name = _rejection_stem;
    for (FileType type : rejections_output.output_types()) {
      IWString fname;
      create_file_with_appropriate_name(_rejection_stem, fname, type, 1);
      const auto same = iwmisc::SameFile(fname, _rejection_log_name);
      if (!same) {
        cerr << "Cannot compare file paths\n";
        return 0;
      }
      if (*same) {
        _rejection_log_name << ".log";
        break;
      }
    }
    if (!ReserveOutputFile(_rejection_log_name)) {
      return 0;
    }
  }

  return 1;
}

// Parse '-Y B4F=<fname>' and open stream_for_smiles_before_filters.
int
LocalOptions::OpenPreFilteredSmilesFile(Command_Line& cl, char flag) {
  if (! cl.option_present(flag)) {
    return 1;
  }
  IWString y;
  for (int i = 0; cl.value(flag, y, i); ++i) {
    if (! y.starts_with("B4F=")) {
      continue;
    }

    y.remove_leading_chars(4);
    IWString fname(y);
    if ('-' == y)
      ;
    else if (!fname.ends_with(".smi"))
      fname << ".smi";

    if (!ReserveOutputFile(fname)) {
      return 0;
    }
    _prefiltered_smiles_file = fname;

    return 1;
  }

  return 1;
}

int
LocalOptions::SetupReporting(Command_Line& cl, char flag) {
  const_IWSubstring y;
  for (int i = 0; cl.value(flag, y, i); ++i) {
    if (! y.starts_with("rpt=")) {
      continue;
    }
    y.remove_leading_chars(4);
    uint32_t rpt;
    if (! y.numeric_value(rpt) || rpt == 0) {
      cerr << "Invalid -Y rpt=... directive '" << y << "'\n";
      return 0;
    }
    report_progress.set_report_every(rpt);
    return 1;
  }

  return 1;
}

void
LocalOptions::SetupInputFile(data_source_and_type<Molecule>& input) const {
  if (verbose > 1) {
    input.set_verbose(1);
  }

  if (connection_table_errors_allowed) {
    input.set_connection_table_errors_allowed(connection_table_errors_allowed);
    if (connection_table_error_file.length())
      input.set_connection_table_error_file(connection_table_error_file);
  }
}

int
LocalOptions::OutputFileSeemsOk(const char * fname) const {
  if (audit_input) {
    return 1;
  }

  // output file names generated by the output object
  if (output.name_token_for_file_name() >= 0)  {
    return 1;
  }

  // Single output file for each input.
  if (output_file_stem.length()) {
    return 1;
  }

  if (output.would_use_name(fname)) {
    cerr << "LocalOptions::OutputFileSeemsOk: input '" << fname << "' and output must be distinct\n";
    return 0;
  }

  return 1;
}


int
LocalOptions::MaybeOpenNewOutputFIle(const char * fname, int first_call) {
  if (audit_input) {  // Nothing to do.
    return 1;
  }

  // One file for all input files, open once.
  if (output_file_stem.length()) {
    if (first_call)
      return output.new_stem(output_file_stem);
    else
      return output.ok();
  }

  // New output file for each molecule.
  if (output.name_token_for_file_name() >= 0) {
    return 1;
  }

  // New output file for this input.
  if (!FlushMolecules(output) || !output.new_stem(fname)) {
    return 0;
  }

  return 1;
}

int
LocalOptions::MaybeWritePostTransformSmiles(const IWString& smiles) {
  if (smiles.empty()) {
    return 1;
  }
  if (! stream_for_smiles_before_filters.is_open()) {
    return 1;
  }

  // organic subset
  stream_for_smiles_before_filters << smiles << '\n';
  return stream_for_smiles_before_filters.write_if_buffer_holds_more_than(32768);
}

int
LocalOptions::MaybeWriteRejection(Molecule& m,
                             const IWString& rejection_reason) {
  if (reject_log.is_open()) {
    reject_log << molecules_read << " '" << m.name() << "' REASON " << rejection_reason << '\n';
    if (!reject_log.write_if_buffer_holds_more_than(4096)) {
      return 0;
    }
  }

  if (rejections_output.active()) {
    m << " REJ " << rejection_reason;
    if ((rejections_output.name_token_for_file_name() >= 0 ||
         rejections_output.molecules_per_file() > 0) &&
        !PrepareDynamicOutput(m, rejections_output, _rejection_files_created)) {
      return 0;
    }
    if (!rejections_output.write(m)) {
      return 0;
    }
  }

  return 1;
}

int
LocalOptions::ReportResults(const Command_Line& cl, std::ostream& output) const {
  if (cl.number_elements() > 1)
    output << molecules_read << " molecules read, ";
  if (!audit_input)
    output << molecules_written << " molecules written";
  if (molecules_changed)
    output << ' ' << molecules_changed << " molecules changed";
  output << '\n';

  return output.good();
}

int
LocalOptions::Fileconv(Molecule& m,
                       fileconv::FileconvConfig& config) {

  if (report_progress()) {
    cerr << "Read " << molecules_read << " molecules";
    if (! audit_input) {
      cerr << ", wrote " << molecules_written;
      if (molecules_changed) {
        cerr << molecules_changed << " molecules changed";
      }
    }
    cerr << '\n';
  }

  if (debug_print_each_molecule) {
    m.compute_aromaticity_if_needed();
    m.debug_print(cerr);
  }
    
  FileconvResult result = config.Process(m);
  if (result.error) {
    cerr << "Fatal error processing " << m.smiles() << ' ' << m.name() << '\n';
    return 0;
  }

  if (!MaybeWritePostTransformSmiles(result.smiles_before_filters)) {
    return 0;
  }

  if (result.rejected) {
    return MaybeWriteRejection(m, result.rejection_reason);
  }

  if (result.molecule_changed) {
    molecules_changed++;
  }

  if (audit_input) {  // Nothing more to do.
    return 1;
  }

  if (output.name_token_for_file_name() >= 0 || output.molecules_per_file() > 0) [[unlikely]] {
    if (!PrepareDynamicOutput(m, output, _main_files_created)) {
      return 0;
    }
  }

  if (! output.write(m)) {
    return 0;
  }

  molecules_written++;

  return 1;
}

int
Fileconv(data_source_and_type<Molecule> & input,
         fileconv::FileconvConfig& config,
         LocalOptions& local_options) {
  Molecule * m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    local_options.molecules_read++;

    if (! local_options.Fileconv(*m, config)) {
      return 0;
    }
  }

  if (input.stopped_because_of_error())
    return 0;

  return 1;
}

int
Fileconv(const char *fname,
         FileType input_type,
         fileconv::FileconvConfig& config,
         LocalOptions & local_options) {
  assert(nullptr != fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    if (FILE_TYPE_INVALID == input_type) {
      cerr << "Cannot determine input type for '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << prog_name << ": cannot open '" << fname << "'\n";
    return 0;
  }

  local_options.SetupInputFile(input);

  if (! local_options.OutputFileSeemsOk(fname)) {
    return 0;
  }
// Set up the output object for this file stem
// If there is a new stem for output files, make sure we get it.
// Make sure we deal properly with the case of multiple input files
// and a single output file (via the -S option)

  static bool first_call = true;

  if (! local_options.MaybeOpenNewOutputFIle(fname, first_call)) {
    cerr << "Cannot establish output for '" << fname << "'\n";
    first_call = false;
    return 0;
  }

  first_call = false;

  return Fileconv(input, config, local_options);
}

int
FileconvListOfFiles(iwstring_data_source & input,
                    FileType input_type,
                    fileconv::FileconvConfig& config,
                    LocalOptions& local_options) {
  input.set_strip_trailing_blanks(1);
  input.set_skip_blank_lines(1);
  input.set_strip_leading_blanks(1);

  IWString buffer;

  while (input.next_record(buffer)) {
    if (buffer.starts_with('#'))
      continue;

    if (local_options.verbose > 1)
      cerr << "Processing '" << buffer << "'\n";

    if (! Fileconv(buffer.null_terminated_chars(), input_type, config, local_options)) {
      cerr << "Fatal error processing file '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
FileconvListOfFiles(const char * fname,
                    FileType input_type,
                    fileconv::FileconvConfig& config,
                    LocalOptions& local_options) {
  const char * tmp = fname + 2;

  iwstring_data_source input(tmp);

  if (! input.good()) {
    cerr << "Cannot open list of files to process '" << fname << "'\n";
    return 0;
  }

  return FileconvListOfFiles(input, input_type, config, local_options);
}

int
AllFilesRecognisedByType(const Command_Line& cl) {
  for (const char * fname: cl) {
    const_IWSubstring tmp(fname);
    if (tmp.starts_with("F:")) {
      continue;
    }
    if (!discern_file_type_from_name(fname)) {
      cerr << "Cannot determine file type from '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

static int
Fileconv(int argc, char ** argv) {
  Command_Line cl(argc, argv, "PN:T:eB:I:s:g:h:H:t:n:L:S:GA:K:w:W:m:X:c:C:Q:O:E:f:F:vVi:ao:r:R:p:J:Y:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    fileconv::Usage(1);
  }

  DisplayDashHelpIfRequested(cl, 'J', Structure_Fixing::DisplayDashJOptions);

  DisplayDashHelpIfRequested(cl, 'h', fileconv::DisplayDashhOptions);
  DisplayDashHelpIfRequested(cl, 'p', fileconv::DisplayPDirectives);
  DisplayDashHelpIfRequested(cl, 'Y', fileconv::DisplayDashYOptions);
  DisplayDashHelpIfRequested(cl, 'I', fileconv::DisplayDashIOptions);
  DisplayDashHelpIfRequested(cl, 's', fileconv::DisplayChiralityOptions);
  DisplayDashHelpIfRequested(cl, 'f', fileconv::DisplayfOptions);
  DisplayDashHelpIfRequested(cl, 'F', fileconv::DisplayFOptions);
  DisplayDashHelpIfRequested(cl, 'O', fileconv::DisplayOrganicOptions);
  DisplayDashHelpIfRequested(cl, 'Q', fileconv::DisplayChargeOptions);

  DisplayDashHelpIfRequested(cl, 'B', DisplayDashBOptions);

  fileconv::FileconvConfig config;
  if (! config.Build(cl)) {
    cerr << "Cannot initialise fileconv\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;
  LocalOptions local_options;
  if (! local_options.Build(cl)) {
    cerr << "Cannot initialise local options\n";
    return 1;
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      fileconv::Usage(6);
    }
  }

  if (FILE_TYPE_INVALID != input_type) {  // great, explicitly specified
  } else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {  // reading from a pipe, assume smiles input
    if (local_options.verbose)
      cerr << "Assuming smiles input from pipe read\n";
    input_type = FILE_TYPE_SMI;
  } else if (AllFilesRecognisedByType(cl)) {
  } else {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    fileconv::Usage(1);
  }

  int rc = 0;
  for (const char * fname: cl) {
    const const_IWSubstring as_string(fname);

    if (local_options.verbose)
      cerr << "Processing '" << fname << "'\n";

    int myrc;
    if (as_string.length() > 2 && as_string.starts_with("F:"))
      myrc = FileconvListOfFiles(fname, input_type, config, local_options);
    else
      myrc = Fileconv(fname, input_type, config, local_options);

    if (!myrc) {
      rc = 1;
      break;
    }
  }

  if (!local_options.FinishOutputs()) {
    rc = 1;
  }

  if (local_options.verbose) {
    local_options.ReportResults(cl, cerr);
    config.ReportResults(cl, std::cerr);
  }

  return rc;
}

}  // namespace fileconv
//#define DLL_FLAG=1
#ifndef DLL_FLAG

int
main(int argc, char **argv) {
  prog_name = argv[0];

  int rc = fileconv::Fileconv(argc, argv);

  return rc;
}

#endif

#ifdef HZ_CSHARP

extern "C" int
fileconv_csharp(int argc, char **argv) {
  return fileconv::Fileconv(argc, argv);
}

#endif

