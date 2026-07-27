#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_preprocessing.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/substructure.h"

#include "Molecule_Tools/chemotypes.h"
#include "Molecule_Tools/dicer_fragments.pb.h"

#include "absl/container/flat_hash_map.h"
#include "google/protobuf/text_format.h"

namespace chemotypes_main {

using std::cerr;
using molecule_processing::MoleculePreprocessing;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on

  cerr << R"(Extracts query-defined molecular chemotypes.
The first query that matches, and the first embedding of that query, define the seed ring system.

Options:
 -q <query>     LillyMol query file defining the chemotype seed.
 -s <smarts>    SMARTS defining the chemotype seed.
 -n <n>         include <n> directly adjacent ring systems, default 0.
 -r <n>         minimum number of rings required for processing, default 0.
 -D <dist>      maximum bond distance to an adjacent ring system.
 -u             include one-hop atoms attached to retained ring atoms.
 -x             ignore singly connected attached atoms with -u.
 -t             with -n, include all adjacent ring systems tied at the cutoff distance.
 -P <atype>     atom typing specification; non-terminal attachment atoms become labelled.
 -I <iso>       label retained ring exit-point atoms with isotope <iso>; incompatible with -P.
 -p <text>      write the parent before each chemotype; append <text> to parent name.
                Use -p . or -p def for no parent annotation.
 -F <fname>      write accumulated dicer_data::DicerFragment textproto summary.
 -z i|ignore    ignore molecules not matching any query.
 -z f|first     use the first embedding if a query matches multiple ring systems.
 -S <stem>      output file stem; default stdout.
 -i <type>      input type.
 -o <type>      output type, default smiles.
 -g ...         chemical standardisation.
 -l             reduce to largest fragment.
 -c             remove chirality.
 -A ...         aromaticity options.
 -E ...         element options.
 -v             verbose output.
)";

  ::exit(rc);
}

class Options {
 private:
  int _verbose = 0;

  MoleculePreprocessing _preprocessing;

  int _min_rings = 0;
  int _ignore_molecules_not_matching_query = 0;
  int _write_parent = 0;
  IWString _parent_annotation;
  IWString _summary_fname;

  resizable_array_p<Substructure_Query> _queries;
  chemotypes::ChemotypeOptions _chemotype_options;
  chemotypes::ChemotypeScratch _scratch;

  Atom_Typing_Specification _atom_typing;
  Atom_Typing_Specification* _atom_typing_ptr = nullptr;

  absl::flat_hash_map<IWString, dicer_data::DicerFragment> _chemotype;

  uint64_t _molecules_read = 0;
  uint64_t _molecules_written = 0;
  uint64_t _parent_molecules_written = 0;
  uint64_t _unique_chemotypes = 0;
  uint64_t _duplicate_chemotypes = 0;
  uint64_t _molecules_not_matching_queries = 0;
  uint64_t _molecules_below_min_rings = 0;
  uint64_t _molecules_matched_query_without_ring_atom = 0;
  uint64_t _molecules_with_ambiguous_query_matches = 0;
  uint64_t _atom_typing_failures = 0;

 public:
  int Initialise(Command_Line& cl);
  int Preprocess(Molecule& m);
  int Process(Molecule& m, Molecule_Output_Object& output);
  int AccumulateChemotype(Molecule& m, const IWString& parent_name);
  int WriteSummary();
  int Report(std::ostream& output) const;

  int verbose() const {
    return _verbose;
  }
};

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! _preprocessing.Initialise(cl)) {
    cerr << "Cannot initialise molecule preprocessing\n";
    return 0;
  }

  if (! cl.option_present('q') && ! cl.option_present('s')) {
    cerr << "Must specify one or more queries via -q or -s\n";
    return 0;
  }

  if (cl.option_present('q')) {
    if (! process_queries(cl, _queries, _verbose, 'q')) {
      cerr << "Cannot read queries (-q)\n";
      return 0;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Query> query = std::make_unique<Substructure_Query>();
      if (! query->create_from_smarts(smarts)) {
        cerr << "Invalid smarts '" << smarts << "'\n";
        return 0;
      }
      _queries << query.release();
    }
  }

  for (Substructure_Query* query : _queries) {
    query->set_find_unique_embeddings_only(1);
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _chemotype_options.adjacent_ring_systems_to_include) ||
        _chemotype_options.adjacent_ring_systems_to_include < 0) {
      cerr << "The adjacent ring system count (-n) must be a non-negative integer\n";
      return 0;
    }
  }


  if (cl.option_present('r')) {
    if (! cl.value('r', _min_rings) || _min_rings < 0) {
      cerr << "The minimum ring count (-r) must be a non-negative integer\n";
      return 0;
    }
  }

  if (cl.option_present('D')) {
    if (! cl.value('D', _chemotype_options.max_distance) ||
        _chemotype_options.max_distance <= 0) {
      cerr << "The max distance (-D) must be a positive integer\n";
      return 0;
    }
  }

  if (cl.option_present('u')) {
    _chemotype_options.include_attached_atoms = 1;
  }

  if (cl.option_present('p')) {
    _write_parent = 1;
    const_IWSubstring p = cl.string_value('p');
    if (p != '.' && p != "def" && p != "default") {
      _parent_annotation << ' ' << p;
    }
  }

  if (cl.option_present('F')) {
    _summary_fname = cl.string_value('F');
  }

  if (cl.option_present('x')) {
    _chemotype_options.ignore_singly_connected_attached_atoms = 1;
  }

  if (cl.option_present('t')) {
    _chemotype_options.include_tied_adjacent_ring_systems = 1;
  }


  if (cl.option_present('I')) {
    uint64_t iso;
    if (! cl.value('I', iso) || iso == 0) {
      cerr << "The exit point isotope (-I) must be a positive integer\n";
      return 0;
    }
    _chemotype_options.isotope_for_exit_points = static_cast<isotope_t>(iso);
  }

  if (cl.option_present('P')) {
    const_IWSubstring p = cl.string_value('P');
    if (! _atom_typing.build(p)) {
      cerr << "Invalid atom typing specification (-P) '" << p << "'\n";
      return 0;
    }
    _atom_typing_ptr = &_atom_typing;
  }


  if (_chemotype_options.isotope_for_exit_points != 0 && _atom_typing_ptr != nullptr) {
    cerr << "The fixed exit point isotope option (-I) cannot be used with atom typing (-P)\n";
    return 0;
  }

  if (cl.option_present('z')) {
    IWString z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i' || z == "ignore") {
        _ignore_molecules_not_matching_query = 1;
      } else if (z == 'f' || z == "first") {
        _chemotype_options.choose_first_embedding = true;
      } else {
        cerr << "Unrecognised -z qualifier '" << z << "'\n";
        return 0;
      }
    }
  }

  if (_verbose) {
    cerr << "Read " << _queries.number_elements() << " queries\n";
    cerr << "Will include " << _chemotype_options.adjacent_ring_systems_to_include
         << " adjacent ring systems\n";
    if (_min_rings > 0) {
      cerr << "Will only process molecules with at least " << _min_rings << " rings\n";
    }
    if (_chemotype_options.max_distance) {
      cerr << "Adjacent ring systems limited to distance "
           << _chemotype_options.max_distance << '\n';
    }
    if (_chemotype_options.include_attached_atoms) {
      cerr << "Will include one-hop atoms attached to retained ring atoms\n";
    }
    if (_chemotype_options.ignore_singly_connected_attached_atoms) {
      cerr << "Will ignore singly connected attached atoms\n";
    }
    if (_write_parent) {
      cerr << "Will write parent molecules before chemotypes";
      if (! _parent_annotation.empty()) {
        cerr << " with parent annotation '" << _parent_annotation << "'";
      }
      cerr << '\n';
    }
    if (! _summary_fname.empty()) {
      cerr << "Will write chemotype summary to '" << _summary_fname << "'\n";
    }
    if (_chemotype_options.include_tied_adjacent_ring_systems) {
      cerr << "Will include adjacent ring systems tied at the -n cutoff distance\n";
    }
    if (_chemotype_options.choose_first_embedding) {
      cerr << "Will use the first query embedding if multiple ring systems match\n";
    }
    if (_atom_typing_ptr != nullptr) {
      cerr << "Will label non-terminal attachment atoms with atom types\n";
    }
    if (_chemotype_options.isotope_for_exit_points != 0) {
      cerr << "Will label retained ring exit point atoms with isotope "
           << _chemotype_options.isotope_for_exit_points << '\n';
    }
  }

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_preprocessing.active()) {
    _preprocessing.Process(m);
  }

  return ! m.empty();
}

int
Options::Process(Molecule& m, Molecule_Output_Object& output) {
  ++_molecules_read;
  m.compute_aromaticity_if_needed();

  if (_min_rings > 0 && m.nrings() < _min_rings) {
    ++_molecules_below_min_rings;
    return 1;
  }

  const IWString parent_name = m.name();

  std::unique_ptr<Molecule> parent;
  if (_write_parent) {
    parent = std::make_unique<Molecule>(m);
  }

  chemotypes::ChemotypeQueryMatch match;
  const chemotypes::ChemotypeQueryMatchStatus status = chemotypes::ReduceToChemotype(
      m, _queries, _chemotype_options, _scratch, match, _atom_typing_ptr);

  switch (status) {
    case chemotypes::ChemotypeQueryMatchStatus::kMatched:
      if (parent) {
        if (! _parent_annotation.empty()) {
          parent->append_to_name(_parent_annotation);
        }
        if (! output.write(*parent)) {
          return 0;
        }
        ++_parent_molecules_written;
      }
      if (! AccumulateChemotype(m, parent_name)) {
        return 0;
      }
      ++_molecules_written;
      return output.write(m);

    case chemotypes::ChemotypeQueryMatchStatus::kNoQueryMatch:
      ++_molecules_not_matching_queries;
      if (_ignore_molecules_not_matching_query) {
        return 1;
      }
      cerr << "No query matched '" << m.name() << "'\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kMatchedQueryNoRingAtom:
      ++_molecules_matched_query_without_ring_atom;
      cerr << "Query matched no ring atom in '" << m.name() << "'\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kAmbiguousQueryMatches:
      ++_molecules_with_ambiguous_query_matches;
      cerr << "Query matched multiple ring systems in '" << m.name()
           << "', use -z first to use the first embedding\n";
      return 0;

    case chemotypes::ChemotypeQueryMatchStatus::kAtomTypingFailed:
      ++_atom_typing_failures;
      cerr << "Atom typing failed for '" << m.name() << "'\n";
      return 0;
  }

  return 0;
}

dicer_data::Isotope
DicerIsotope(const chemotypes::ChemotypeOptions& options,
             const Atom_Typing_Specification* atom_typing) {
  if (atom_typing != nullptr) {
    return dicer_data::ATYPE;
  }

  if (options.isotope_for_exit_points != 0) {
    return dicer_data::ATT;
  }

  return dicer_data::NONE;
}

int
Options::AccumulateChemotype(Molecule& m, const IWString& parent_name) {
  if (_summary_fname.empty()) {
    return 1;
  }

  const IWString& usmi = m.unique_smiles();
  auto iter = _chemotype.find(usmi);
  if (iter != _chemotype.end()) {
    iter->second.set_n(iter->second.n() + 1);
    ++_duplicate_chemotypes;
    return 1;
  }

  dicer_data::DicerFragment proto;
  proto.set_smi(usmi.data(), usmi.length());
  proto.set_par(parent_name.data(), parent_name.length());
  proto.set_nat(m.natoms());
  proto.set_n(1);

  const dicer_data::Isotope iso = DicerIsotope(_chemotype_options, _atom_typing_ptr);
  if (iso != dicer_data::NONE) {
    proto.set_iso(iso);
  }

  _chemotype.emplace(usmi, std::move(proto));
  ++_unique_chemotypes;
  return 1;
}

int
Options::WriteSummary() {
  if (_summary_fname.empty()) {
    return 1;
  }

  IWString_and_File_Descriptor output;
  if (! output.open(_summary_fname.null_terminated_chars())) {
    cerr << "Options::WriteSummary:cannot open '" << _summary_fname << "'\n";
    return 0;
  }

  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  for (const auto& [usmi, proto] : _chemotype) {
    if (! printer.PrintToString(proto, &buffer)) {
      cerr << "Options::WriteSummary:cannot write '" << proto.ShortDebugString() << "'\n";
      return 0;
    }
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  output.flush();
  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Wrote " << _molecules_written << " chemotypes\n";
  if (_parent_molecules_written) {
    output << "Wrote " << _parent_molecules_written << " parent molecules\n";
  }
  if (! _summary_fname.empty()) {
    output << "Accumulated " << _unique_chemotypes << " unique chemotypes and "
           << _duplicate_chemotypes << " duplicate chemotypes\n";
  }
  if (_molecules_not_matching_queries) {
    output << _molecules_not_matching_queries << " molecules did not match any query\n";
  }
  if (_molecules_below_min_rings) {
    output << _molecules_below_min_rings << " molecules had fewer than " << _min_rings
           << " rings\n";
  }
  if (_molecules_matched_query_without_ring_atom) {
    output << _molecules_matched_query_without_ring_atom
           << " molecules matched a query with no ring atom\n";
  }
  if (_molecules_with_ambiguous_query_matches) {
    output << _molecules_with_ambiguous_query_matches
           << " molecules matched queries in multiple ring systems\n";
  }
  if (_atom_typing_failures) {
    output << _atom_typing_failures << " atom typing failures\n";
  }

  return 1;
}

int
Chemotypes(Options& options, data_source_and_type<Molecule>& input,
           Molecule_Output_Object& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! options.Process(*m, output)) {
      return 0;
    }
  }

  return 1;
}

int
Chemotypes(Options& options, const char* fname, FileType input_type,
           Molecule_Output_Object& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Chemotypes:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return Chemotypes(options, input, output);
}

int
Chemotypes(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:g:cli:o:S:F:q:s:n:r:D:P:I:up:xtz:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    Usage(1);
  }
  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process element options\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (cl.number_elements() == 1 && 0 == strcmp(cl[0], "-")) {
    input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  Molecule_Output_Object output;
  if (! cl.option_present('o')) {
    output.add_output_type(FILE_TYPE_SMI);
  } else if (! output.determine_output_types(cl, 'o')) {
    cerr << "Cannot determine output type(s)\n";
    return 1;
  }

  if (cl.option_present('S')) {
    IWString stem = cl.string_value('S');
    if (output.would_overwrite_input_files(cl, stem)) {
      cerr << "Cannot overwrite input file(s) with stem '" << stem << "'\n";
      return 1;
    }
    if (! output.new_stem(stem)) {
      cerr << "Cannot open output stem '" << stem << "'\n";
      return 1;
    }
  } else {
    output.new_stem("-");
  }

  for (const char* fname : cl) {
    if (! Chemotypes(options, fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (! options.WriteSummary()) {
    return 1;
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace chemotypes_main

int
main(int argc, char** argv) {
  return chemotypes_main::Chemotypes(argc, argv);
}
