// Use SAFE smiles representations to do molecular transformations.

#include <iostream>

// Deliberate choice of flat_hash_map because we need key pointer stability.
#include "absl/container/flat_hash_map.h"

#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"
#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace safe_replacement {

using std::cerr;

// By convention the Usage function tells how to use the tool.
void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  cerr << R"(
)";
  // clang-format on
  ::exit(rc);
}

// A class that holds all the information needed for the
// application. The idea is that this should be entirely
// self contained. If it were moved to a separate header
// file, then unit tests could be written and tested
// separately.
class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    int _molecules_read = 0;

    resizable_array_p<Substructure_Query> _fragment_must_contain;
    resizable_array_p<Substructure_Query> _fragment_must_not_contain;

    absl::flat_hash_map<IWString, dicer_data::DicerFragment> _fragment;

    // private functions

    int ReadFragments(IWString& fname);
    int ReadFragments(iwstring_data_source& input);
    int ReadFragment(const const_IWSubstring& line);

  public:
    Options();

    // Get user specified command line directives.
    int Initialise(Command_Line& cl);

    int verbose() const {
      return _verbose;
    }

    // After each molecule is read, but before any processing
    // is attempted, do any preprocessing transformations.
    int Preprocess(Molecule& m);

    // Helpful when the -i option is not given.
    int MaybeDiscernInputType(const char * fname);

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove all chirality\n";
    }
  }

  if (! cl.option_present('F')) {
    cerr << "Must specify fragments to read via the -F option\n";
    return 0;
  }

  if (cl.option_present('F')) {
    IWString fname;
    for (int i = 0; cl.value('F', fname, i); ++i) {
      if (! ReadFragments(fname)) {
        cerr << "Cannot read fragments from '" << fname << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      cerr << "Read " << _fragment.size() << " fragments\n";
    }
  }

  return 1;
}

int
Options::ReadFragments(IWString& fname) {
  iwstring_data_source input;
  if (! input.open(fname.null_terminated_chars())) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragments(input);
}

int
Options::ReadFragments(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadFragment(buffer)) {
      cerr << "Invalid textproto '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadFragment(const const_IWSubstring& line) {
  google::protobuf::io::ArrayInputStream zero_copy_array(line.data(), line.nchars());
  dicer_data::DicerFragment proto;
  if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
    cerr << "Options:ReadFragment:cannot parse proto " << line << '\n';
    return 0;
  }

  Molecule m;
  if (! m.build_from_smiles(proto.smi()) {
  }
  IWString key(proto.smi());
  _fragment.emplace(key, std::move(proto));

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  // Other information about what has happened.

  return 1;
}

// If the input type is known, return it.
// Otherwise examine the file name's suffix to 
// determine the type.
int
Options::MaybeDiscernInputType(const char * fname) {

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  // Do the actual work, probably write something to `output`.

  // The IWString_and_File_Descriptor object needs to be flushed.
  output.write_if_buffer_holds_more_than(4092);

  return 1;
}

int
SafeReplacement2(Options& options,
                 Molecule& m,
                 IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
SafeReplacement2(Options& options,
                 const_IWSubstring& buffer,
                 IWString_and_File_Descriptor& output) {
  Molecule m;
  if (! m.build_from_smiles(buffer)) {
    cerr << "SafeReplacement2:invalid smiles '" << buffer << "'\n";
    return 0;
  }

  return SafeReplacement2(options, m, output);
}

int
SafeReplacement(Options& options,
                iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! SafeReplacement2(options, buffer, output)) {
      cerr << "SafeReplacement2::fatal error processing '" << buffer << "'\n";
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
SafeReplacement(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input;
  if (! input.open(fname)) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return SafeReplacement(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:H:N:T:A:lcg:i:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }


  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! SafeReplacement(options, fname, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace safe_replacement

int
main(int argc, char ** argv) {

  int rc = safe_replacement::Main(argc, argv);

  return rc;
}
