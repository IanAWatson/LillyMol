// We have a list of de-functionalised reagents.
// Given a set of molecule, for each reagent, perform a substructure seach
// on that molecule and identify matched atoms.

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>

#include "absl/container/flat_hash_set.h"

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/mformula.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/reagent_substructure_search.pb.h"
#else
#include "reagent_substructure_search.pb.h"
#endif

namespace reagent_substructure_search {

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
// clang-format on
// clang-format off
  cerr << R"(Reagent substructure search.
Consumed
 -R <fname>  recordio file of reagent_substructure_search::ReagentData serialized protos,
                such as might come from Enamine2Reagents.
 -m          minimum number of atoms in reagents.
 -M          maximum number of atoms in reagents.
 -r <ncon>   allowed values for number of connections.
 -O <fname>  write labelled molecules to <fname>.
 -c          remove chirality.
 -v          verbose output.
)";
// clang-format on

  ::exit(rc);
}

#ifdef NOT_NEEDED_JJAA
class MoleculeAndFormula {
  private:
    Molecule* _m = nullptr;
    mformula::MFormula _formula;

    Molecule_to_Match _target;

  public:
    int Initialise(const std::string& smiles);

    const Molecule& mol() const {
      return *_m;
    }
    Molecule& mol() {
      return *_m;
    }

    const mformula::MFormula& formula() const {
      return _formula;
    }
    mformula::MFormula& formula() {
      return _formula;
    }

};

int
MoleculeAndFormula::Initialise(const std::string& smiles) {
  if (! _m->build_from_smiles(smiles)) {
    cerr << "MoleculeAndFormula::Initialise:invalid smiles '" << smiles << "'\n";
    return 0;
  }

  _formula.Build(*_m);

  _target.initialise_molecule(_m);
}
#endif

// As we identify matched atoms in a target molecule, we need to keep
// track of which atoms are hit by the queries. For each of the _natoms
// atoms in the molecule, keep track of which Reagent matches.
class Coverage {
  private:
    const int _natoms;

    absl::flat_hash_set<uint32_t>* _by_atom;

  public:
    Coverage(int n);
    ~Coverage();

    // An indication that reagent number `reagent` has matched `atom`.
    void Extra(int atom, int reagent);

    float FractionAtomsCovered() const;

    uint32_t Count(int atom) const {
      return _by_atom[atom].size();
    }
};

Coverage::Coverage(int n) : _natoms(n) {
  _by_atom = new absl::flat_hash_set<uint32_t>[_natoms];
}

Coverage::~Coverage() {
  delete [] _by_atom;
}

void
Coverage::Extra(int atom, int reagent) {
  _by_atom[atom].insert(reagent);
}

float
Coverage::FractionAtomsCovered() const {
  int atoms_matched = 0;
  for (int i = 0; i < _natoms; ++i) {
    if (! _by_atom[i].empty()) {
      ++atoms_matched;
    }
  }

  return iwmisc::Fraction<float>(atoms_matched, _natoms);
}

class Reagent {
  private:
    ReagentData _proto;

    // Derived form the smiles, will always be present.
    mformula::MFormula _formula_smiles;

    // will come into existence if we need the aromatic formula.
    std::unique_ptr<Molecule> _m;

    std::unique_ptr<mformula::MFormula> _formula;

    // If we decide we need to do searches, this will come into existence.
    std::unique_ptr<Substructure_Query> _query;

    // Our position in some external array.
    uint32_t _ndx;

  // private functions.

  public:
    Reagent();

    // `proto` is destroyed.
    bool BuildFromProto(ReagentData& proto,
                uint32_t min_natoms,
                uint32_t max_natoms,
                const bool* _ok_ncon);

    void set_ndx(int s) {
      _ndx = s;
    }

    int natoms() const {
      return _proto.natoms();
    }

    int SubstructureSearch(Molecule_to_Match& target,
                           const mformula::MFormula& rhs_mf,
                           const mformula::MFormula& rhs_mf_arom,
                           Coverage& coverage);
};

Reagent::Reagent() {
  _ndx = 0;
}

bool
Reagent::BuildFromProto(ReagentData& proto,
                uint32_t min_natoms,
                uint32_t max_natoms,
                const bool* _ok_ncon) {
  if (proto.natoms() < min_natoms) {
    return false;
  }
  if (proto.natoms() > max_natoms) {
    return false;
  }

  if (! _ok_ncon[proto.ncon()]) {
    return false;
  }

  if (! _formula_smiles.Build(proto.smiles())) {
    cerr << "Reagent::BuildFromProto:invalid smiles '" << proto.smiles() << "'\n";
    return 0;
  }

  _proto = std::move(proto);

  return true;
}

int
Reagent::SubstructureSearch(Molecule_to_Match& target,
                const mformula::MFormula& rhs_mf,
                const mformula::MFormula& rhs_mf_arom,
                Coverage& coverage) {
  if (static_cast<int32_t>(_proto.natoms()) > target.natoms()) {
    return 0;
  }

  if (! _formula_smiles.IsSubset(rhs_mf)) {
    return 0;
  }

  if (! _formula) {
    _m = std::make_unique<Molecule>(_proto.natoms());
    if (! _m->build_from_smiles(_proto.smiles())) {
      cerr << "Reagent::SubstructureSearch:invalid smiles '" << _proto.smiles() << '\n';
      return 0;
    }

    _formula = std::make_unique<mformula::MFormula>();
    _formula->Build(*_m);
  }

  if (! _formula->IsSubset(rhs_mf_arom)) {
    return 0;
  }

  if (! _query) {
    static Molecule_to_Query_Specifications mqs;
    mqs.set_substituents_only_at_isotopic_atoms(1);
    mqs.set_must_have_substituent_at_every_isotopic_atom(1);

    _query = std::make_unique<Substructure_Query>();
    if (! _query->create_from_molecule(*_m, mqs)) {
      cerr << "Reagent::SubstructureSearch:cannot create query\n";
      return 0;
    }
  }

  Substructure_Results sresults;
  if (! _query->substructure_search(target, sresults)) {
    return 0;
  }

  for (const Set_of_Atoms* e : sresults.embeddings()) {
    for (const atom_number_t a : *e) {
      coverage.Extra(a, _ndx);
    }
  }

  return 1;
}

constexpr int kMaxNcon = 3;

class Options {
  private:
    int _verbose = 0;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    // Not a part of all applications, just an example...
    Element_Transformations _element_transformations;

    uint32_t _min_natoms;
    uint32_t _max_natoms;

    // We only read Reagents that have their ncon() value set in
    // this array.
    bool _ok_ncon[kMaxNcon + 1];

    // Read from the -R file.
    resizable_array_p<Reagent> _reagent;

    uint64_t _molecules_read = 0;

    extending_resizable_array<uint32_t> _matches;

    Accumulator<double> _acc_coverage;

    IWString_and_File_Descriptor _stream_for_labelled_molecules;

    // private functions

    int ReadReagents(IWString& fname);
    int ReadReagents(iw_tf_data_record::TFDataReader& reader);
    int AddReagent(reagent_substructure_search::ReagentData& proto);

    int HandleNoMatches(Molecule& m);

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

  _min_natoms = 0;
  _max_natoms = std::numeric_limits<uint32_t>::max();

  std::fill_n(_ok_ncon, 4, true);
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

  if (cl.option_present('m')) {
    if (! cl.value('m', _min_natoms)) {
      cerr << "Options::Initialise:invalid minimum atom count -m\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only consider reagents with at least " << _min_natoms << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_natoms)) {
      cerr << "Options::Initialise:invalid maximum atom count -M\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only consider reagents with at most " << _max_natoms << " atoms\n";
    }
  }

  if (cl.option_present('r')) {
    std::fill_n(_ok_ncon, kMaxNcon + 1, true);
    uint32_t r;
    for (int i = 0; cl.value('r', r, i); ++i) {
      if (r >= kMaxNcon) {
        cerr << "Invalid number of connections " << r << '\n';
        return 0;
      }
    }
    _ok_ncon[r] = true;
  }

  if (! cl.option_present('R')) {
    cerr << "Must specify one or more TFDataReader binary files of\n";
    cerr << "reagent_substructure_search::ReagentData protos via the -R option\n";
    return 0;
  } else {
    IWString fname;
    for (int i = 0; cl.value('R', fname, i); ++i) {
      if (! ReadReagents(fname)) {
        cerr << "Options::Initialise:cannot read reagents '" << fname << "'\n";
        return 0;
      }
    }
  }

  if (_verbose) {
    cerr << "Read " << _reagent.size() << " reagents\n";
  }

  if (cl.option_present('O')) {
    IWString fname = cl.string_value('O');
    fname.EnsureEndsWith(".smi");
    if (! _stream_for_labelled_molecules.open(fname)) {
      cerr << "Options::Initialise:cannot open '" << fname << "' -O\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Labelled molecules written to '" << fname << "'\n";
    }
  }

  return 1;
}

int
Options::ReadReagents(IWString& fname) {
  iw_tf_data_record::TFDataReader reader(fname);
  if (! reader.good()) {
    cerr << "Options::ReadReagents:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadReagents(reader);
}

int
Options::ReadReagents(iw_tf_data_record::TFDataReader& reader) {
  while (1) {
    std::optional<reagent_substructure_search::ReagentData> maybe_proto =
      reader.ReadProto<reagent_substructure_search::ReagentData>();
    if (maybe_proto) {
    } else if (reader.eof()) {
      return 1;
    } else {
      cerr << "Options::ReadReagents:error reading proto\n";
      return 0;
    }

    if (! AddReagent(*maybe_proto)) {
      cerr << "Options::ReadReagent:invalid proto\n";
      // Will fail if it has been destroyed
      cerr << maybe_proto->ShortDebugString() << '\n';
      return 0;
    }
  }

  return _reagent.size();
}

int
Options::AddReagent(reagent_substructure_search::ReagentData& proto) {
  std::unique_ptr<Reagent> rgnt = std::make_unique<Reagent>();

  if (! rgnt->BuildFromProto(proto, _min_natoms, _max_natoms, _ok_ncon)) {
    return 0;
  }

  rgnt->set_ndx(_reagent.size());

  _reagent << rgnt.release();

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";

  output << "Coverage btw " << _acc_coverage.minval() << " and " <<
            _acc_coverage.maxval() << " ave " << _acc_coverage.average() << '\n';

  for (int i = 0; i < _matches.number_elements(); ++i) {
    if (_matches[i]) {
      output << _matches[i] << " molecules matched " << i << " reagents\n";
    }
  }

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
    m.revert_all_directional_bonds_to_non_directional();
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

  const int matoms = m.natoms();

  Molecule_to_Match target(&m);

  Coverage coverage(matoms);

  mformula::MFormula formula;
  formula.set_consider_aromatic(0);
  formula.Build(m);

  mformula::MFormula formula_arom;
  formula_arom.set_consider_aromatic(1);
  formula_arom.Build(m);

  uint32_t matches_found = 0;
  for (Reagent* r : _reagent) {
    if (! r->SubstructureSearch(target, formula, formula_arom, coverage)) {
      continue;
    }
    ++matches_found;
  }

  if (_verbose) {
    ++_matches[matches_found];
  }

  if (matches_found == 0) {
    return HandleNoMatches(m);
  }

  float f = coverage.FractionAtomsCovered();

  _acc_coverage.extra(f);

  if (_stream_for_labelled_molecules.is_open()) {
    // Make a copy - do we really need to????
    Molecule tmp(m);
    for (int i = 0; i < matoms; ++i) {
      isotope_t iso = coverage.Count(i);
      if (iso > 0) {
        tmp.set_isotope(i, iso);
      }
    }

    _stream_for_labelled_molecules << tmp.aromatic_smiles() << ' ' << m.name() << '\n';

    _stream_for_labelled_molecules.write_if_buffer_holds_more_than(32768);
  }

  output.write_if_buffer_holds_more_than(4092);

  return 1;
}

int
Options::HandleNoMatches(Molecule& m) {
  return 1;
}

int
ApplicationName(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
ApplicationName(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! ApplicationName(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
ApplicationName(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "ApplicationName:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return ApplicationName(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:H:N:T:A:lcg:i:m:M:r:O:");

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
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! ApplicationName(options, fname, input_type, output)) {
      cerr << "ApplicationName::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace reagent_substructure_search

int
main(int argc, char ** argv) {

  int rc = reagent_substructure_search::Main(argc, argv);

  return rc;
}
