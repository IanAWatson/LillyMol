// Generate a molecular formula fingerprint

#include <cstdint>
#include <iostream>
#include <sstream>

#include "Foundational/cmdline/cmdline.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/mformula.h"
#include "Utilities/GFP_Tools/fixed_size_counted_fingerprint.h"

namespace formula_fingperprint {

using std::cerr;
using mformula::MFormula;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(Fingerprint based on the molecular formula.
 -J <tag>       tag to use. The tag prefix determines the kind of fingerprint written.
                FC...  fixed size counted fingerprint (default FCFML), as the Python GFP.formula()
                NC...  sparse counted fingerprint
 -f             function as a TDT filter.
 -L s,d         log scaling of counts. Emphasises rarer elements, see header file
                for discussion. The default is -L 200,41, as in Python. Use -L none for raw counts.
 -H             do not count atoms bearing hydrogens. For approximate substructure searching a
                molecule must have at least the counts of the query, but a query has hydrogens on
                atoms that are substituted in a molecule that contains it. The default tag
                becomes FCFMLNH, as in Python GFP.formula(hydrogens=False), so these cannot be
                compared with fingerprints that include the hydrogen counts.
 -i ...         standard input       control options.
 -E ....        standard Element     control options.
 -A ....        standard Aromaticity control options.
 -v             verbose output.
)";
  // clang-format off

  ::exit(rc);
}

// The log scaling used by the Python GFP.formula().
constexpr float kPythonLogS = 200.0f;
constexpr float kPythonLogD = 41.0f;

class Options {
  private:
    int _verbose;

    int _function_as_tdt_filter;

    // Determined by the tag prefix. FC fixed size counted, otherwise sparse (NC).
    int _fixed_counted = 1;

    // The -H option turns this off.
    int _include_hydrogens = 1;

    // The -L option
    // The default is the same as the Python GFP.formula().
    float _log_s = kPythonLogS;
    float _log_d = kPythonLogD;

    IWString _smiles_tag;
    IWString _identifier_tag;
    IWString _tag;

  public:
    Options();

    int Initialise(Command_Line& cl);

    const IWString& smiles_tag() const {
      return _smiles_tag;
    }

    int WriteFingerprint(Molecule& m, IWString_and_File_Descriptor& output);

    int WriteFixedCounted(const MFormula& mf, Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _function_as_tdt_filter = 0;
  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
  _tag = "FCFML<";
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('f')) {
    _function_as_tdt_filter = 1;
    if (_verbose) {
      cerr << "Will function as a tdt filter\n";
    }
  }

  if (cl.option_present('H')) {
    _include_hydrogens = 0;
    _tag = "FCFMLNH<";
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    _tag.EnsureEndsWith('<');
  }

  // The tags used by the Python GFP.formula() must mean the same thing here, or
  // fingerprints that cannot be compared would carry the same tag.
  if (_tag == "FCFML<" && ! _include_hydrogens) {
    cerr << "Options::Initialise:tag FCFML is for fingerprints that include the hydrogen counts, "
            "use FCFMLNH with -H\n";
    return 0;
  }
  if (_tag == "FCFMLNH<" && _include_hydrogens) {
    cerr << "Options::Initialise:tag FCFMLNH is for fingerprints without the hydrogen counts, "
            "specify -H\n";
    return 0;
  }

  // The tag prefix determines what kind of fingerprint is written, so the tag
  // and the data can never disagree.
  if (_tag.starts_with("FC")) {
    _fixed_counted = 1;
  } else if (_tag.starts_with("NC")) {
    _fixed_counted = 0;
  } else {
    cerr << "Options::Initialise:invalid tag '" << _tag
         << "', must start with FC (fixed size counted) or NC (sparse counted)\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Writing " << (_fixed_counted ? "fixed size counted" : "sparse counted")
         << " fingerprints with tag " << _tag << '\n';
  }

  if (cl.option_present('L')) {
    const_IWSubstring l = cl.string_value('L');

    if (l == "none") {
      _log_s = 0.0f;
      _log_d = 0.0f;
    } else {
      const_IWSubstring s, d;
      if (! l.split(s, ',', d) || s.empty() || d.empty()) {
        cerr << "Options::Initialise:invalid -L specifier '" << l << "' must be <float>,<float> or none\n";
        return 0;
      }
      if (! s.numeric_value(_log_s) || _log_s <= 0.0f) {
        cerr << "Options::Initialise:invalid log S value '" << l << "'\n";
        return 0;
      }
      if (! d.numeric_value(_log_d) || _log_d <= 0.0f) {
        cerr << "Options::Initialise:invalid log D value '" << l << "'\n";
        return 0;
      }
    }
  }

  if (_verbose) {
    if (_log_s > 0.0f) {
      cerr << "Log scaling enabled with s=" << _log_s << " d=" << _log_d << '\n';
    } else {
      cerr << "No log scaling, raw counts\n";
    }
  }

  // The tag does not record the scaling, so a different scaling gives a
  // fingerprint with the same tag, but which cannot be compared with Python.
  if (_fixed_counted && (_log_s != kPythonLogS || _log_d != kPythonLogD)) {
    cerr << "Warning: the Python GFP.formula() uses log scaling " << kPythonLogS << ','
         << kPythonLogD << ", these fingerprints will not be comparable with it\n";
  }

  return 1;
}

int
Options::WriteFingerprint(Molecule& m,
                          IWString_and_File_Descriptor& output) {
  
  MFormula mf;
  mf.set_include_hydrogens(_include_hydrogens);
  if (_log_s > 0.0f) {
    mf.set_log_scaling_factors(_log_s, _log_d);
  }

  mf.Build(m);

  IWString fp;
  if (_fixed_counted) {
    return WriteFixedCounted(mf, m, output);
  }

  mf.ToSparseFingerprint(fp);

  if (_function_as_tdt_filter) {
    output << _tag << fp << ">\n";
    return 1;
  }

  output << _smiles_tag << m.smiles() << ">\n";
  output << _identifier_tag << m.name() << ">\n";
  output << _tag << fp << ">\n";
  output << "|\n";

  return 1;
}

int
Options::WriteFixedCounted(const MFormula& mf, Molecule& m,
                           IWString_and_File_Descriptor& output) {
  int counts[mformula::kMFOther + 1];
  if (! mf.ToFixedCountedFingerprint(counts, mformula::kMFOther + 1)) {
    cerr << "Options::WriteFixedCounted:cannot generate fingerprint " << m.name() << '\n';
    return 0;
  }

  Fixed_Size_Counted_Fingerprint_uchar fp;
  if (! fp.construct_from_array_of_ints(counts, mformula::kMFOther + 1)) {
    cerr << "Options::WriteFixedCounted:cannot build fingerprint " << m.name() << '\n';
    return 0;
  }

  std::ostringstream tmp;
  if (! fp.write_daylight_ascii_representation(tmp, _tag, 1)) {
    cerr << "Options::WriteFixedCounted:cannot write fingerprint " << m.name() << '\n';
    return 0;
  }

  if (! _function_as_tdt_filter) {
    output << _smiles_tag << m.smiles() << ">\n";
    output << _identifier_tag << m.name() << ">\n";
  }

  output << tmp.str().c_str();

  if (! _function_as_tdt_filter) {
    output << "|\n";
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  return 1;
}

int
FormulaFingerprintFilterRecord(const_IWSubstring line,  // note local copy
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  line.remove_leading_chars(options.smiles_tag().length());
  line.chop();

  Molecule m;
  if (! m.build_from_smiles(line)) {
    cerr << "FormulaFingerprintFilterRecord:invalid smiles\n";
    return 0;
  }

  return options.WriteFingerprint(m, output);
}

int
FormulaFingerprintFilter(iwstring_data_source& input,
                         Options& options,
                         IWString_and_File_Descriptor& output) {

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
    if (! buffer.starts_with(options.smiles_tag())) {
      continue;
    }

    if (! FormulaFingerprintFilterRecord(buffer, options, output)) {
      cerr << "FormulaFingerprintFilter:error processing\n";
      cerr << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
FormulaFingerprintFilter(const char* fname,
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "FormulaFingerprintFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return FormulaFingerprintFilter(input, options, output);
}

int
FormulaFingerprint(Molecule& m,
                   Options& options,
                   IWString_and_File_Descriptor& output) {
  return options.WriteFingerprint(m, output);
}

int
FormulaFingerprint(data_source_and_type<Molecule>& input,
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    FormulaFingerprint(*m, options, output);
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
FormulaFingerprint(const char* fname,
                   FileType input_type,
                         Options& options,
                         IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);
                          
  if (! input.good()) {
    cerr << "FormulaFingerprintFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return FormulaFingerprint(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:EJ:fi:L:H");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  Options options;

  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  const int function_as_tdt_filter = cl.option_present('f');

  if (function_as_tdt_filter && cl.option_present('i')) {
    cerr << "Cannot combine TDT filter (-f) and input file (-i) options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficent arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    int rc;
    if (function_as_tdt_filter) {
      rc = FormulaFingerprintFilter(fname, options, output);
    } else {
      rc = FormulaFingerprint(fname, input_type, options, output);
    }

    if (rc == 0) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace formula_fingperprint

int 
main(int argc, char ** argv) {
  int rc = formula_fingperprint::Main(argc, argv);

  return rc;
}
