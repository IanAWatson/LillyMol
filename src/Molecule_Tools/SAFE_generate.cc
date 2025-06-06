// Consume the output of mol2SAFE and generate random variants.

#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <tuple>

#include "google/protobuf/io/zero_copy_stream_impl_lite.h"
#include "google/protobuf/text_format.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/matcher.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/molecule_filter_lib.h"
#include "Molecule_Tools/safe_generate_lib.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#include "Molecule_Tools/safe_generate.pb.h"
#else
#include "dicer_fragments.pb.h"
#include "safe_generate.pb.h"
#endif

namespace safe_generate {

using std::cerr;

constexpr char kPercent = '%';

constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';

using iwmatcher::Matcher;
using molecule_filter_lib::MoleculeFilter;

// We ignore any fragments with more than this many connections.
int max_ncon = 10;

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
  cerr << R"(Denovo generation of molecules from SAFE smiles.
Consumes the output from mol2SAFE, and the -L file is also from mol2SAFE.

mol2SAFE ... -I 1 -c -M 15 -z -S library.textproto file.smi > file.safe.smi
safe_generate ... -L library.textproto -C generate.config file.safe.smi > new_molecules.smi

 -C <fname>             safe_generate textproto configuration file.
 -a <ncon>              maximum value for number of connections processed - applied to library.
 -L <fname>             fragment library of SAFE fragments
 -Y <query>             queries for atoms that are allowed to change
 -N <query>             queries for atoms that are NOT allowed to change
 -z i                   ignore molecules not matching queries
 -x extra=<n>           the number of extra atoms in a fragment being substitued
 -x fewer=<n>           the number of fewer atoms in a fragment being substitued
 -n <n>                 number of molecules to generate by fragment swap from library
 -b <n>                 number of molecules to generate by breeding
 -e <n>                 for each molecule, exhaustively scan the library(s) making all possibilities.
                        For each molecule, a max of <n> variants are made. Use '-e 0' for no limit.
 -F <fname>             molecule filter textproto file - products must pass
 -p                     write the parent molecule before variants
 -s                     write SAFE smiles
 -k <number>            abandon if no new molecules formed for <number> steps
 -v                     verbose output
)";
// clang-format on

  ::exit(rc);
}

class SafedMolecule {
  private:
    Molecule _m;

    // the SAFE smiles from which we are built.
    IWString _smiles;

    // the indices where each fragment starts in _smiles
    resizable_array<int> _frag_start_smiles;
    // the indices where each fragment starts in _m;
    resizable_array<int> _frag_start_atom;

    resizable_array_p<SafeFragment> _frag;

    // When random items are selected, we limit the number of tries.
    int _max_attempts;

    std::mt19937 _rng;
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _dist;

    // We may have atom types.
    atom_type_t* _atype;

    // Private functions
    int GetFragment(atom_number_t zatom) const;

  public:
    SafedMolecule();
    ~SafedMolecule();

    int Build(const const_IWSubstring& buffer);

    int DebugPrint(std::ostream& output) const;

    Molecule& mol() {
      return _m;
    }

    const Molecule& mol() const {
      return _m;
    }

    const IWString& name() const {
      return _m.name();
    }

    const IWString& smiles() const {
      return _smiles;
    }

    int number_fragments() const {
      return _frag.number_elements();
    }

    const SafeFragment* fragment(int ndx) const {
      return _frag[ndx];
    }

    int SetupRng();

    int AssignAtomTypes(Atom_Typing_Specification& ats);

    const atom_type_t* AtomTypes() {
      return _atype;
    }
    atom_type_t AtomType(atom_number_t zatom) const {
      return _atype[zatom];
    }

    // Given a set of queries, perform the searches. For every match examine
    // the matched atoms, and any fragment that contains a matched atom is
    // marked for being OK to change.
    int IdentifyChangingFragments(resizable_array_p<Substructure_Query>& queries);
    int IdentifyUnChangingFragments(resizable_array_p<Substructure_Query>& queries);

    std::optional<int> ChooseFragment();

    // Return the index of a random fragment.
    std::optional<int> RandomFragment();

    // The index of the fragment closest in size to the atoms in `target`.
    std::optional<int> MatchingFragment(const SafeFragment* target) const;

    int NewSmiles(int f1_ndx, const SafeFragment& f2, IWString& new_smiles) const;
};

SafedMolecule::SafedMolecule() {
  std::random_device rd;
  _rng.seed(rd());
  _max_attempts = 10;  // arbitrary number.
  _atype = nullptr;
}

SafedMolecule::~SafedMolecule() {
  if (_atype != nullptr) {
    delete [] _atype;
  }
}

int
SafedMolecule::SetupRng() {
  if (_frag.empty()) {
    cerr << "SafedMolecule::SetupRng:empty\n";
    return 0;
  }

  _dist = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _frag.size() - 1);

  return 1;
}

int
SafedMolecule::DebugPrint(std::ostream& output) const {
  output << "SafedMolecule with " << _frag.size() << " fragments\n";
  output << _smiles << ' ' << _m.name() << '\n';

  for (int i = 0; i < _frag.number_elements(); ++i) {
    output << i << ' ';
    _frag[i]->DebugPrint(output);
  }
  return 1;
}

int
SafedMolecule::AssignAtomTypes(Atom_Typing_Specification& ats) {
  if (_atype == nullptr) {
    _atype = new atom_type_t[_m.natoms()];
  }

  return ats.assign_atom_types(_m, _atype);
}

int
SafedMolecule::Build(const const_IWSubstring& buffer) {
  const_IWSubstring smiles, id;
  if (! buffer.split(smiles, ' ', id) ||
        smiles.empty() || id.empty()) {
    cerr << "SafedMolecule::Build:cannot split to smiles id '" << buffer << "'\n";
    return 0;
  }

  if (! _m.build_from_smiles(smiles)) {
    cerr << "SafedMolecule::build_from_smiles:invalid smiles '" << smiles << "' from '" <<
             buffer << "'\n";
    return 0;
  }

  _smiles = smiles;

  _m.set_name(id);

  _frag_start_smiles << 0;
  _frag_start_atom << 0;

  const_IWSubstring token;
  int i = 0;
  int prev_i = i;
  int atom_count = 0;
  while (buffer.nextword(token, i, '.')) {
    std::unique_ptr<SafeFragment> f = std::make_unique<SafeFragment>();
    if (! f->Build(token)) {
      cerr << "SafedMolecule::Build:invalid fragment '" << token << "' from '" <<
              buffer << "'\n";
      return 0;
    }

    atom_count += f->natoms();
    f->set_name(id);
    _frag << f.release();
    _frag_start_smiles << prev_i;
    _frag_start_atom << atom_count;
    prev_i = i;
  }

#ifdef DEBUG_BUILD_SAFE_MOLECULE
  cerr << "From\n" << buffer << '\n';
  for (int s : _frag_start_smiles) {
    cerr << ' ' << s;
  }
  cerr << '\n';
#endif

  return 1;
}

#ifdef NOT_USED
int
SafedMolecule::Build(const dicer_data::DicerFragment& proto) {
  IWString smi = proto.smi();

  if (smi.empty()) {
    cerr << "SafedMolecule::Build:no smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (! _m.build_from_smiles(smi)) {
    cerr << "SafedMolecule::Build:invalid smiles '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  return 1;
}
#endif

// Return the fragment number where `zatom` is found.
// This could be made more efficient with an atom_to_fragment
// array, but this is not expected to be a bottleneck.
int
SafedMolecule::GetFragment(atom_number_t zatom) const {
  for (int i = _frag_start_atom.number_elements() - 1; i >= 0; --i) {
    if (zatom >= _frag_start_atom[i]) {
      return i;
    }
  }

  cerr << "SafedMolecule::GetFragment:cannot find fragment for atom " << zatom << '\n';
  cerr << "Have " << _m.natoms() << " atoms\n";
  cerr << "Fragment starts ";
  for (int f : _frag_start_atom) {
    cerr << ' ' << f;
  }
  cerr << '\n';
  cerr << "SafedMolecule::GetFragment:this should not happen\n";
  return -1;
}

int
SafedMolecule::IdentifyChangingFragments(resizable_array_p<Substructure_Query>& queries) {
  for (SafeFragment* f : _frag) {
    f->set_ok_to_select(0);
  }

  int rc = 0;

  Molecule_to_Match target(&_m);

  for (Substructure_Query* q : queries) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    for (const Set_of_Atoms* c : sresults.embeddings()) {
      for (atom_number_t a : *c) {
        int f = GetFragment(a);
        _frag[f]->set_ok_to_select(1);
        rc = 1;
      }
    }
  }

  return rc;
}

int
SafedMolecule::IdentifyUnChangingFragments(resizable_array_p<Substructure_Query>& queries) {
  int rc = 0;

  Molecule_to_Match target(&_m);

  for (Substructure_Query* q : queries) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    for (const Set_of_Atoms* c : sresults.embeddings()) {
      for (atom_number_t a : *c) {
        int f = GetFragment(a);
        _frag[f]->set_ok_to_select(0);
        rc = 1;
      }
    }
  }

  return rc;
}


std::optional<int>
SafedMolecule::ChooseFragment() {
  if (_frag.size() == 1) {
    return 0;
  }

  for (int i = 0; i < _max_attempts; ++i) {
    const int ndx = (*_dist)(_rng);
    if (! _frag[ndx]->ok_to_select()) {
      continue;
    }

    return ndx;
  }

  return std::nullopt;
}

std::optional<int>
SafedMolecule::RandomFragment() {
  for (int i = 0; i < _max_attempts; ++i) {
    const int ndx = (*_dist)(_rng);
    SafeFragment* f = _frag[ndx];
    if (f->ok_to_select()) {
      return ndx;
    }
  }

  return std::nullopt;
}

int
SafedMolecule::NewSmiles(int f1_ndx, const SafeFragment& f2, IWString& new_smiles) const {
  const SafeFragment* f1 = _frag[f1_ndx];
  f1->SameNumbers(f2, new_smiles);

  return 1;
}

std::optional<int>
SafedMolecule::MatchingFragment(const SafeFragment* target) const {
  int min_diff = std::numeric_limits<int>::max();
  int result = -1;

  const int natoms = _m.natoms();

  const int n = _frag.number_elements();
  for (int i = 0; i < n; ++i) {
    const SafeFragment* f = _frag[i];
    if (! f->ok_to_select()) {
      continue;
    }

    if (f->ncon() != target->ncon()) {
      continue;
    }

    int d = std::abs(natoms - f->natoms());
    if (d < min_diff) {
      min_diff = d;
      result = i;
    }
  }

  if (result >= 0) {
    return result;
  }

  return std::nullopt;
}

class LibraryIterator;

class Library {
  friend class LibraryIterator;

  private:
    // A mapping from <ncon, natoms> to the molecules fitting that description.
    absl::flat_hash_map<std::tuple<uint32_t, uint32_t>, resizable_array_p<SafeFragment>> _frag;

    // As we read the library, we keep track of the number of fragments
    // we discard because of too many connections
    int _excessive_connections;

    std::mt19937 _rng;

    int _max_attempts;

  public:
    Library();
    ~Library();

    int Build(IWString& fname);
    int Build(iwstring_data_source& input);
    int BuildMember(const const_IWSubstring& line);

    int SetupRng();

    // Change the isotopic label of each atom that has an isotopic label to `iso`.
    int SetIsotope(isotope_t iso);

    uint32_t size() const;

    const SafeFragment* GetFragment(int ncon, int natoms);
};

Library::Library() {
  _max_attempts = 10;  // arbitrary number.
  _excessive_connections = 0;
}

Library::~Library() {
}

int
Library::SetupRng() {
  std::random_device rd;
  _rng.seed(rd());

  return 1;
}

uint32_t
Library::size() const {
  uint32_t rc = 0;

  for (const auto& [k, v] : _frag) {
    rc += v.size();
  }

  return rc;
}

const SafeFragment*
Library::GetFragment(int ncon, int natoms) {
  std::tuple<uint32_t, uint32_t> key(ncon, natoms);
  auto iter = _frag.find(key);
  if (iter == _frag.end()) {
    return nullptr;
  }

  std::uniform_int_distribution<uint32_t> u(0, iter->second.size() - 1);
  uint32_t ndx = u(_rng);

  return iter->second.item(ndx);
}

int
Library::SetIsotope(isotope_t iso) {
  for (const auto& [_, v] : _frag) {
    for (SafeFragment* f : v) {
      Molecule& m = f->mol();
      const int matoms = m.natoms();
      for (int i = 0; i < matoms; ++i) {
        if (m.isotope(i)) {
          m.set_isotope(i, iso);
        }
      }
    }
  }

  return 1;
}

int
Library::Build(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Library::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  return Build(input);
}

int
Library::Build(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    // cerr << "Library building from '" << buffer << "'\n";
    if(! BuildMember(buffer)) {
      cerr << "Library::Build:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Library::BuildMember(const const_IWSubstring& line) {
  google::protobuf::io::ArrayInputStream zero_copy_array(line.data(), line.nchars());
  dicer_data::DicerFragment proto;
  if (!google::protobuf::TextFormat::Parse(&zero_copy_array, &proto)) {
    cerr << "Library:BuildMember:cannot parse proto " << line << '\n';
    return 0;
  }

  // cerr << "Proto built " << proto.ShortDebugString() << '\n';
  std::unique_ptr<SafeFragment> f = std::make_unique<SafeFragment>();
  if (! f->Build(proto)) {
    cerr << "line::BuildMember:cannot parse '" << line << "'\n";
    return 0;
  }

  // cerr << "Fragment contains " << f->ncon() << " connections, cmp " << max_ncon << '\n';
  if (f->ncon() > max_ncon) {
    ++_excessive_connections;
    return 1;
  }

  std::tuple<uint32_t, uint32_t> key(f->ncon(), f->natoms());
  auto iter = _frag.find(key);
  if (iter != _frag.end()) {
    iter->second << f.release();
    return 1;
  }

  resizable_array_p<SafeFragment> tmp;
  tmp << f.release();
  _frag.emplace(key, std::move(tmp));

  return 1;
}

class LibraryIterator {
  private:
    // The hash from a Library.
    const absl::flat_hash_map<std::tuple<uint32_t, uint32_t>, resizable_array_p<SafeFragment>>& _hash;

    // an iterator over _hash
    absl::flat_hash_map<std::tuple<uint32_t, uint32_t>, resizable_array_p<SafeFragment>>::const_iterator _iter;

    // within each value is a resizable array of SafeFragment*.
    SafeFragment* const * _frag_iter;

    // As we cycle over the items in _hash, this gets updated for each resizable_array_p.
    SafeFragment* const * _cend;

  public:
    LibraryIterator(const Library& lib);

    const SafeFragment* Next();
};

LibraryIterator::LibraryIterator(const Library& lib) : _hash(lib._frag), _iter(lib._frag.cbegin()) {
  _frag_iter = _iter->second.cbegin();
  _cend = _iter->second.cend();
}

const SafeFragment*
LibraryIterator::Next() {
  ++_frag_iter;
  if (_frag_iter == _iter->second.cend()) {
    ++_iter;
    if (_iter == _hash.cend()) {
      return nullptr;
    }
    _frag_iter = _iter->second.cbegin();
    _cend = _iter->second.cend();
  } 

  return *_frag_iter;
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

    safe_generate::Config _config;

    MoleculeFilter _filter;

    resizable_array_p<Library> _library;

    resizable_array_p<SafedMolecule> _mols;

    // The queries instantiated from the -Y and -N options.
    resizable_array_p<Substructure_Query> _can_change;
    resizable_array_p<Substructure_Query> _cannot_change;

    int _ignore_molecules_not_matching_queries;

    std::optional<uint32_t> _extra_atoms;
    std::optional<uint32_t> _fewer_atoms;

    std::optional<uint32_t> _extra_rings;
    std::optional<uint32_t> _fewer_rings;

    // Derived from the proto if set.
    std::optional<uint32_t> _max_formula_difference;

    std::optional<uint32_t> _max_distance_difference;

    Atom_Typing_Specification _atom_typing_specification;

    // The -X option.
    resizable_array_p<Substructure_Query> _discard_if_match;

    std::mt19937 _rng;

    // when asked to generate fragments, this is how many times we
    // try.
    int _max_attempts;

    // Random number generators for the _mols and _library arrays.
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _mols_dist;
    std::unique_ptr<std::uniform_int_distribution<uint32_t>> _libs_dist;

    // Do we write the starting molecule or not
    int _write_parent_molecule;

    // By default, we write aromatic smiles, but we can optionally write
    // SAFE smiles.
    int _write_safe_smiles;

    // we can impose limits on the size of fragments that are selected
    // for replacement.
    int _min_atoms_in_fragment;
    int _max_atoms_in_fragment;

    absl::flat_hash_set<IWString> _seen;

    int _remove_isotopes;

    // If we have not generated any new molecules for this many steps,
    // give up.
    int _abandon_for_futility;

    uint64_t _number_from_breeding;
    uint64_t _number_from_exhaustive;

    uint64_t _new_molecules_formed;
    uint64_t _rejected_by_bad_valence;
    uint64_t _rejected_by_discard_queries;
    uint64_t _rejected_by_seen_before;
    uint64_t _rejected_by_filter;
    uint64_t _rejected_by_adjacent_atoms;

  // Private functions.
    int ReadLibrary(IWString& fname);
    int ReadMolecules(iwstring_data_source& input);
    int ReadMolecule(const const_IWSubstring& line);

    int TransferFromConfig(const safe_generate::Config& proto);

    int IdentifyChangingFragments();
    int IdentifyChangingFragments(SafedMolecule& m);
    void IdentifyUnChangingFragments();
    void IdentifyUnChangingFragments(SafedMolecule& m);

    const SafeFragment* GetFragment(int ncon, int natoms);

    int Generate(SafedMolecule& m, Library& lib, IWString_and_File_Descriptor& output);
    int Generate(SafedMolecule& m,
                  const int f1_ndx, const SafeFragment& f2,
                  IWString_and_File_Descriptor& output);
    int AnyDiscardQueriesMatch(Molecule& m);
    int SeenBefore(Molecule& m);
    int OkAtomCount(const int natoms) const;
    int OkAtomCountDifference(int n1, int n2) const;
    int OkRingCountDifference(int n1, int n2) const;
    int OkDistanceDifference(int n1, int n2) const;
    int OkDifferences(const SafeFragment& f1, const SafeFragment& f2) const;
    int OkFormulaDifference(const SafeFragment& f1, const SafeFragment f2) const;
    int ProcessNewMolecule(Molecule& m, const IWString& name1,
                        const IWString& safe_smiles,
                        const SafeFragment& f2,
                        IWString_and_File_Descriptor& output);

    int Breed(IWString_and_File_Descriptor& output);
    int Breed(SafedMolecule& m1, SafedMolecule& m2, IWString_and_File_Descriptor& output);
    int SelectAtomCount(const SafeFragment& f);
    int SelectFragments(SafedMolecule& m1, const SafedMolecule& m2,
                int& f1, int& f2);
    int SelectFragmentsAtype(SafedMolecule& m1, const SafedMolecule& m,
                              int f1, int& f2);

    int MakeAllLibrary(const SafedMolecule& m, uint64_t max_make, IWString_and_File_Descriptor& output);
    int MakeAllLibrary(const SafedMolecule& m,
                        const Library& lib,
                        uint64_t max_make,
                        IWString_and_File_Descriptor& output);

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

    int ReadMolecules(const char* fname);

    const resizable_array_p<SafedMolecule>& mols() const {
      return _mols;
    }

    // The function that actually does the processing,
    // and may write to `output`.
    // You may instead want to use a Molecule_Output_Object if
    // Molecules are being written.
    // You may choose to use a std::ostream& instead of 
    // IWString_and_File_Descriptor.
    int Process(Molecule& mol, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;

    int SetupRng();

    int DoSubstructureSearches();

    int SetIsotope(isotope_t iso);

    int Generate(int ngenerate, IWString_and_File_Descriptor& output);
    int Generate(SafedMolecule& m, int ngenerate, IWString_and_File_Descriptor& output);
    int Breed(int nbreed, IWString_and_File_Descriptor& output);

    // For each molecule in the input, swap in each member of the library.
    // For each molecule, make a max of `max_make` variants.
    int MakeAllLibrary(uint64_t max_make, IWString_and_File_Descriptor& output);
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;

  _ignore_molecules_not_matching_queries = 0;

  _max_attempts = 100;
  _number_from_breeding = 0;
  _number_from_exhaustive = 0;

  std::random_device rd;
  _rng.seed(rd());

  _write_parent_molecule = 0;
  _write_safe_smiles = 0;

  _min_atoms_in_fragment = 0;
  _max_atoms_in_fragment = std::numeric_limits<int>::max();

  _remove_isotopes = 0;

  _abandon_for_futility = 0;

  _new_molecules_formed = 0;
  _rejected_by_bad_valence = 0;
  _rejected_by_discard_queries = 0;
  _rejected_by_seen_before = 0;
  _rejected_by_filter = 0;
  _rejected_by_adjacent_atoms = 0;
}

int
ReadQueries(Command_Line& cl, char flag, int verbose,
            resizable_array_p<Substructure_Query>& destination) {
  IWString q;
  for (int i = 0; cl.value(flag, q, i); ++i) {
    if (! process_cmdline_token(flag, q, destination, verbose)) {
      cerr << "ReadQueries:cannot process '" << q << "'\n";
      return 0;
    }
  }

  return destination.size();
}

int
Options::Initialise(Command_Line& cl) {

  _verbose = cl.option_count('v');

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
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

  if (cl.option_present('a')) {
    if (! cl.value('a', max_ncon) || max_ncon < 1) {
      cerr << "Options::Initialise:invalid max ncon value (-a)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only consider fragments with at most " << max_ncon <<
              " connections\n";
    }
  }

  if (cl.option_present('C')) {
    IWString fname = cl.string_value('C');
    std::optional<safe_generate::Config> maybe_config = 
                iwmisc::ReadTextProtoCommentsOK<safe_generate::Config>(fname);
    if (! maybe_config) {
      cerr << "Cannot initialise config (-C)\n";
      return 0;
    }

    _config = *maybe_config;
    TransferFromConfig(*maybe_config);
  }

  if (cl.option_present('T')) {
    if (!_element_transformations.construct_from_command_line(cl, _verbose, 'T'))
      Usage(8);
  }

  if (cl.option_present('Y')) {
    if (! ReadQueries(cl, 'Y', _verbose, _can_change)) {
      cerr << "Options::Initialise:cannot initialise -Y options\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _can_change.size() << " -Y queries\n";
    }
  }

  if (cl.option_present('N')) {
    if (! ReadQueries(cl, 'N', _verbose, _cannot_change)) {
      cerr << "Options::Initialise:cannot initialise -N options\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _cannot_change.size() << " -N queries\n";
    }
  }

  if (cl.option_present('X')) {
    if (! ReadQueries(cl, 'X', _verbose, _discard_if_match)) {
      cerr << "Options::Initialise:cannot initialise -X options\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _discard_if_match.size() << " -X queries\n";
    }
  }

  if (cl.option_present('z')) {
    const_IWSubstring z;
    for (int i = 0; cl.value('z', z, i); ++i) {
      if (z == 'i') {
        _ignore_molecules_not_matching_queries = 1;
        if (_verbose) {
          cerr << "Will ignore molecules not matching any of the can match queries\n";
        }
      } else  {
      }
    }
  }

  if (cl.option_present('L')) {
    IWString fname;
    for (int i = 0; cl.value('L', fname, i); ++i) {
      cerr << "Reading '" << fname << "'\n";
      if (! ReadLibrary(fname)) {
        cerr << "Options::Initialise:cannot read library file '" << fname << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      uint32_t nfrag = 0;
      for (const Library* lib : _library) {
        nfrag += lib->size();
      }
      cerr << "Read " << nfrag << " Library fragments from " <<
              _library.size() << " library files\n";
    }

    for (Library* lib : _library) {
      lib->SetupRng();
    }
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    if (! _filter.Build(fname)) {
      cerr << "Optionss::Initialise:cannot initialise filter '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Filter initialised '" << fname << "'\n";
    }
  } 

  if (cl.option_present('x')) {
    const_IWSubstring x;
    for (int i = 0; cl.value('x', x, i); ++i) {
      if (x.starts_with("extra=")) {
        x.remove_leading_chars(6);
        int tmp;
        if (! x.numeric_value(tmp)) {
          cerr << "Invalid 'extra=" << x << " directive\n";
          return 0;
        }
        _extra_rings = tmp;
      } else if (x .starts_with("fewer=")) {
        x.remove_leading_chars(6);
        int tmp;
        if (! x.numeric_value(tmp)) {
          cerr << "Invalid 'fewer=" << x << " directive\n";
          return 0;
        }
        _fewer_atoms = tmp;
      } else {
        cerr << "Unrecognised -x directive '" << x << "'\n";
        return 0;
      }
    }
  }

  if (cl.option_present('p')) {
    _write_parent_molecule = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('s')) {
    _write_safe_smiles = 1;
    if (_verbose) {
      cerr << "Will write SAFE smiles\n";
    }
  }

  if (cl.option_present('I')) {
    _remove_isotopes = 1;
    if (_verbose) {
      cerr << "Will remove isotopic labels from product molecules\n";
    }
  }

  if (cl.option_present('k')) {
    if (! cl.value('k', _abandon_for_futility)) {
      cerr << "Invalid abandon for futility (-k)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will abandon calculations unless a new molecule generated every " <<
                _abandon_for_futility << " steps\n";
    }
  }

  return 1;
}

int
Options::TransferFromConfig(const safe_generate::Config& proto) {
  if (proto.has_etrans()) {
    if (! _element_transformations.Build(proto.etrans())) {
      cerr << "Options::TransferFromConfig:invalid etrans " << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.has_ignore_molecules_not_matching_queries()) {
    _ignore_molecules_not_matching_queries = proto.ignore_molecules_not_matching_queries();
  }

  if (proto.has_min_atoms_in_fragment()) {
    _min_atoms_in_fragment = proto.min_atoms_in_fragment();
  }
  if (proto.has_max_atoms_in_fragment()) {
    _max_atoms_in_fragment = proto.max_atoms_in_fragment();
  }

  if (proto.has_extra_atoms()) {
    _extra_atoms = proto.extra_atoms();
  }
  if (proto.has_fewer_atoms()) {
    _fewer_atoms = proto.fewer_atoms();
  }
 
  for (const std::string& fname: proto.library()) {
    IWString tmp(fname);
    if (! ReadLibrary(tmp)) {
      cerr << "Options::TransferFromConfig::cannot read library '" << fname << "'\n";
      return 0;
    }
  }

  for (const std::string& q : proto.can_change()) {
    IWString tmp(q);
    constexpr char kFlag = 'Y';
    if (! process_cmdline_token(kFlag, tmp, _can_change, _verbose)) {
      cerr << "Options::TransferFromConfig:invalid can change " << q << '\n';
      return 0;
    }
  }

  for (const std::string& q : proto.cannot_change()) {
    IWString tmp(q);
    constexpr char kFlag = 'N';
    if (! process_cmdline_token(kFlag, tmp, _cannot_change, _verbose)) {
      cerr << "Options::TransferFromConfig:invalid cannot change " << q << '\n';
      return 0;
    }
  }

  for (const std::string& q :proto.discard_if_match()) {
    IWString tmp(q);
    constexpr char kFlag = 'X';
    if (! process_cmdline_token(kFlag, tmp, _discard_if_match, _verbose)) {
      cerr << "Options::TransferFromConfig:invalid discard if change " << q << '\n';
      return 0;
    }
  }

  if (proto.has_molecule_filter()) {
    if (! _filter.Build(proto.molecule_filter())) {
      cerr << "Options::TransferFromConfig:invalid molecule filter\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }
  }

  if (proto.has_max_formula_difference()) {
    _max_formula_difference = proto.max_formula_difference();
  }

  if (proto.has_max_distance_difference()) {
    _max_distance_difference = proto.max_distance_difference();
  }

  if (proto.has_atom_type()) {
    const_IWSubstring t = proto.atom_type();
    if (! _atom_typing_specification.build(t)) {
      cerr << "Options::TransferFromConfig:invalid atom typing specification " <<
        proto.atom_type() << '\n';
      return 0;
    }
    if (_verbose) {
      cerr << "Atom typing initialised " << proto.atom_type() << '\n';
    }
  }

  return 1;
}

int
Options::ReadLibrary(IWString& fname) {
  std::unique_ptr<Library> lib = std::make_unique<Library>();
  if (! lib->Build(fname)) {
    cerr << "Options::ReadLibrary:cannot read '" << fname << "'\n";
    return 0;
  }

  _library << lib.release();

  return _library.size();
}

int
Options::ReadMolecules(const char* fname) {
  if (_mols.empty()) {
    _mols.resize(100000);
  }

  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadMolecules:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadMolecules(input);
}

int
Options::ReadMolecules(iwstring_data_source& input) {
  const_IWSubstring line;
  while (input.next_record(line)) {
    if (! ReadMolecule(line)) {
      cerr << "Options::ReadMolecule:cannot process '" << line << "'\n";
      return 0;
    }
  }

  return _mols.number_elements();
}

int
Options::ReadMolecule(const const_IWSubstring& line) {
  std::unique_ptr<SafedMolecule> f = std::make_unique<SafedMolecule>();
  if (f == nullptr) {
    cerr << "Options::ReadMolecule:memory failure\n";
    return 0;
  }
  if (! f->Build(line)) {
    cerr << "Options::ReadMolecule:cannot process '" << line << "'\n";
    return 0;
  }

  _mols << f.release();

  return _mols.size();
}

int
Options::SetupRng() {
  if (_mols.empty()) {
    cerr << "Options::SetupRng: no molecules\n";
    return 0;
  }

  _mols_dist = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _mols.size() - 1);
  _libs_dist = std::make_unique<std::uniform_int_distribution<uint32_t>>(0, _library.size() - 1);

  for (SafedMolecule* f : _mols) {
    f->SetupRng();
  }

  return 1;
}

int
Options::SetIsotope(isotope_t iso) {
  for (Library* lib : _library) {
    lib->SetIsotope(iso);
  }

  return 1;
}

int
Options::DoSubstructureSearches() {
  if (_can_change.size() > 0) {
    if (! IdentifyChangingFragments()) {
      return 0;
    }
  }

  if (_cannot_change.size() > 0) {
    IdentifyUnChangingFragments();
  }

  return 1;
}

int
Options::AnyDiscardQueriesMatch(Molecule& m) {
  if (_discard_if_match.empty()) {
    return 0;
  }

  Molecule_to_Match target(&m);
  for (Substructure_Query* q : _discard_if_match) {
    if (q->substructure_search(target)) {
      return 1;
    }
  }

  return 0;
}

int
Options::IdentifyChangingFragments() {
  for (SafedMolecule* m : _mols) {
    if (! IdentifyChangingFragments(*m)) {
      cerr << "Options::IdentifyChangingFragments:no matches to " << m->name() << '\n';
      return _ignore_molecules_not_matching_queries;
    }
  }

  return 1;
}

int
Options::IdentifyChangingFragments(SafedMolecule& m) {
  return m.IdentifyChangingFragments(_can_change);
}

void
Options::IdentifyUnChangingFragments() {
  for (SafedMolecule* m : _mols) {
    m->IdentifyUnChangingFragments(_cannot_change);
  }
}

int
Options::Report(std::ostream& output) const {

  output << "Options: generated " << _new_molecules_formed << " molecules\n";
  output << _rejected_by_bad_valence << " _rejected_by_bad_valence\n";
  if (_discard_if_match.size() > 0) {
    output << _rejected_by_discard_queries << " _rejected_by_discard_queries\n";
  }
  output << _rejected_by_seen_before << " _rejected_by_seen_before\n";
  output << _rejected_by_filter << " _rejected_by_filter\n";
  output << _rejected_by_adjacent_atoms << " _rejected_by_adjacent_atoms\n";

  if (_number_from_breeding > 0) {
    output << "Breed " << _number_from_breeding << " molecules\n";
  }
  if (_number_from_exhaustive > 0) {
    output << "Exhaustive enumeration " << _number_from_exhaustive << '\n';
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
Options::Generate(int ngenerate, IWString_and_File_Descriptor& output) {
  for (SafedMolecule* m : _mols) {
    Generate(*m, ngenerate, output);

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Options::Generate(SafedMolecule& m, int ngenerate, IWString_and_File_Descriptor & output) {
  if (_write_parent_molecule) {
    output << m.smiles() << ' ' << m.name() << " parent\n";
  }

  int generated = 0;

  int max_attempts = ngenerate * 10;
  int last_successful = 0;

  for (int attempts = 0; attempts < max_attempts && generated < ngenerate; ++attempts) {
    int libindex = 0;
    if (_library.size() >= 1) {
      libindex = (*_libs_dist)(_rng);
    }
    if (Generate(m, *_library[libindex], output)) { 
      ++generated;
      last_successful = attempts;
    } else if (last_successful + _abandon_for_futility > attempts) {
      if (_verbose) {
        cerr << "Options::Breed:abandon after " << attempts << " attempts\n";
      }
      break;
    }
  }

  return 1;
}

// Fragment `f` is to be replaced. We need to decide on an atom count
// for the fragment that will replace it.
int
Options::SelectAtomCount(const SafeFragment& f) {
  // If nothing specified, return the same number of atoms.
  if (! _extra_atoms && ! _fewer_atoms) {
    return f.natoms();
  }

  // Need to select a value from a range. Establish the range.
  int aminval;
  // Make sure we have a valid minimum atom count.
  if (_fewer_atoms) {
    if (*_fewer_atoms >= static_cast<uint32_t>(f.natoms())) {
      aminval = 1;
    } else {
      aminval = f.natoms() - *_fewer_atoms;
    }
  } else {
    aminval = 1;
  }

  int amaxval;
  if (_extra_atoms) {
    amaxval = f.natoms() + *_extra_atoms;
  } else {
    amaxval = f.natoms() + 6;  // arbitrary number
  }

  std::uniform_int_distribution<uint32_t> u(aminval, amaxval);
  return u(_rng);
}

int
Options::Generate(SafedMolecule& m, Library& lib, IWString_and_File_Descriptor& output) {
//cerr << "Options;;Generate m ";
//m.DebugPrint(cerr);
  std::optional<int> f1_ndx = m.RandomFragment();
  if (f1_ndx == std::nullopt) {
    return 0;
  }

  const SafeFragment* f1 = m.fragment(*f1_ndx);
  // cerr << "f1_ndx " << *f1_ndx << " natoms " << f1->natoms() << " smiles " << f1->smiles() << '\n';
  // f1->DebugPrint(cerr);
  if (! OkAtomCount(f1->natoms())) {
    return 0;
  }

  const int natoms = SelectAtomCount(*f1);

  const SafeFragment* f2 = lib.GetFragment(f1->ncon(), natoms);
  if (f2 == nullptr) {
    return 0;
  }
#ifdef DEBUG_GENERATE
  cerr << f1->natoms() << " requested " << natoms << " atoms got " << f2->natoms() << " QQ " << (natoms == f2->natoms()) << '\n';
  cerr << f1->smiles() << ' ' << m.name() << '\n';
  Molecule mcopy(const_cast<SafeFragment*>(f1)->mol());
  cerr << mcopy.smiles() << " from molecule, nat " << mcopy.natoms() << " cmp " << f1->natoms() << '\n';
  cerr << "F2 smiles " << f2->smiles() << '\n';
#endif

  if (! OkDifferences(*f1, *f2)) {
    return 0;
  }

  if (! OkDifferences(*f1, *f2)) {
    return 0;
  }

  return Generate(m, *f1_ndx, *f2, output);
}

int
FormNewSmiles(const IWString& starting_smiles, int f1_ndx,
              const IWString& new_smiles,
              IWString& destination) {
  destination.reserve(starting_smiles.size() + 20);   // 20 is an arbitrary choice
  int frag_number = 0;
  int new_smiles_added = 0;
  for (int i = 0; i < starting_smiles.number_elements(); ++i) {
    const char c = starting_smiles[i];
    if (c == '.') {
      destination << c;
      ++frag_number;
    } else if (frag_number == f1_ndx) {
      if (! new_smiles_added) {
        destination << new_smiles;
        new_smiles_added = 1;
      }
    } else {
      destination << c;
    }
  }

  return 1;
}

// `f1` is a fragment number from within `m` and `f2` is a library fragment.
// Generate a molecule by replacing `f1` with `f2`.
int
Options::Generate(SafedMolecule& m,
                  const int f1_ndx, const SafeFragment& f2,
                  IWString_and_File_Descriptor& output) {
  const SafeFragment* f1 = m.fragment(f1_ndx);
  if (f1->ncon() != f2.ncon()) {
    cerr << "Options::Generate:ncon mismatch\n";
    return 0;
  }

  IWString new_smiles;
  if (! m.NewSmiles(f1_ndx, f2, new_smiles)) {
    return 0;
  }
  // cerr << "Fragment smiles " << f1->smiles() << '\n';
  // cerr << "from " << m.smiles() << " and " << f2.smiles() << " get " << new_smiles << '\n';

  IWString tmp;
  FormNewSmiles(m.smiles(), f1_ndx, new_smiles, tmp);

  Molecule newm;
  if (! newm.build_from_smiles(tmp)) {
    cerr << "Options::Generate:invalid smiles '" << tmp << "'\n";
    return 0;
  }

  return ProcessNewMolecule(newm, m.name(), tmp, f2, output);
}

int
Options::ProcessNewMolecule(Molecule& m, const IWString& name1,
                        const IWString& safe_smiles,
                        const SafeFragment& f2,
                        IWString_and_File_Descriptor& output) {
  output.write_if_buffer_holds_more_than(4096);

  ++_new_molecules_formed;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    m.unset_all_implicit_hydrogen_information(i);
    m.set_implicit_hydrogens_known(i, 0);
  }

  if (! m.valence_ok()) {
#ifdef DEBUG_GENERATE
    cerr << "Options::ProcessNewMolecule:invalid valence " << m.smiles() <<
            ' ' << name1 << '\n';
#endif
    ++_rejected_by_bad_valence;
    return 0;
  }

  if (AnyDiscardQueriesMatch(m)) {
#ifdef DEBUG_PROCESS_NEW_MOLECULE
    cerr << "Discard query\n";
#endif
    ++_rejected_by_discard_queries;
    return 0;
  }
  if (SeenBefore(m)) {
#ifdef DEBUG_PROCESS_NEW_MOLECULE
    cerr << "Seen before " << m.smiles() << '\n';
#endif
    ++_rejected_by_seen_before;
    return 0;
  }

  if (! _filter.Ok(m)) {
#ifdef DEBUG_PROCESS_NEW_MOLECULE
    cerr << "Filtered\n";
#endif
    ++_rejected_by_filter;
    return 0;
  }

  if (! BondsOk(m)) {
#ifdef DEBUG_PROCESS_NEW_MOLECULE
    cerr <<  m.smiles() << ' ' << m.name() << " adjacent bonds\n";
#endif
    ++_rejected_by_adjacent_atoms;
    return 0;
  }

  if (_remove_isotopes) {
    m.transform_to_non_isotopic_form();
  }

  // Form a new name.
  m << name1;

  m << " %% " << f2.name() << '.' << f2.ncon() << '.' << f2.natoms();

  if (_write_safe_smiles) {
    output << safe_smiles;
  } else {
    output << m.aromatic_smiles();
  }

  output << ' ' << m.name() << '\n';
  output << f2.smiles() << ' ' << f2.name() << ' ' << f2.nrings() << '\n';

  return 1;
}

int
Options::OkAtomCount(const int natoms) const {
  if (natoms < _min_atoms_in_fragment) {
    return 0;
  }

  if (natoms > _max_atoms_in_fragment) {
    return 0;
  }

  return 1;
}

int
Options::SeenBefore(Molecule& m) {
  auto f = _seen.find(m.unique_smiles());
  if (f != _seen.end()) {
    return 1;
  }

  _seen.insert(m.unique_smiles());
  return 0;
}

int
Options::Breed(int nbreed, IWString_and_File_Descriptor& output) {
  int generated = 0;
  int max_attempts = 50 * nbreed;
  int last_successful = 0;
  for (int attempts = 0; generated < nbreed && attempts < max_attempts; ++attempts) {
    int m1 = (*_mols_dist)(_rng);
    int m2 = 0;
    while (1) {
      m2 = (*_mols_dist)(_rng);
      if (m1 != m2) {
        break;
      }
    }
    // cerr << "Breeding with " << m1 << " and " << m2 << " generated " << generated << '\n';
    if (Breed(*_mols[m1], *_mols[m2], output)) {
      ++generated;
      last_successful = attempts;
    } else if (last_successful + _abandon_for_futility > attempts) {
      if (_verbose) {
        cerr << "Options::Breed:abandon after " << attempts << " attempts\n";
      }
      break;
    }
  }

  _number_from_breeding = generated;

  return 1;
}

// Identify a fragment number, `f1` from `m1` that is about the
// same size as a fragment `f2` from `m2`.
int
Options::SelectFragments(SafedMolecule& m1, const SafedMolecule& m2,
                int& f1, int& f2) {

  std::optional<int> maybef1 = m1.ChooseFragment();
  if (! maybef1) {
    return 0;
  }

  f1 = *maybef1;

  if (_atom_typing_specification.active()) {
    return SelectFragmentsAtype(m1, m2, f1, f2);
  }

  int atoms_in_f1 = m1.fragment(f1)->natoms();
  int ncon1 = m1.fragment(f1)->ncon();
  // cerr << "Finding frag with " << ncon1 << " connections, with " << atoms_in_f1 << " atoms\n";

  // Find the fragment in `m2` closest to atoms_in_f1.
  // Note that we do not respect the ok_to_select attribute?
  const int nf2 = m2.number_fragments();

  f2 = -1;
  int min_diff = std::numeric_limits<int>::max();
  for (int i = 0; i < nf2; ++i) {
    const SafeFragment* f = m2.fragment(i);
    if (ncon1 != f->ncon()) {
      continue;
    }

    int d = std::abs(f->natoms() - atoms_in_f1);
    if (d < min_diff) {
      min_diff = d;
      f2 = i;
    }
  }

  return f2 >= 0;
}

int
Options::SelectFragmentsAtype(SafedMolecule& m1, const SafedMolecule& m,
                              int f1, int& f2) {

  // TODO:ianwatson. Maybe implement this sometime, seems hard.
  cerr << "SelectFragmentsAtype:not implemented\n";
  std::exit(1);
#ifdef SELECTFRAGMENTSATYPE_DONE
  const SafeFragment* frag1 = m1.fragment(f1);

  const int atoms_in_f1 = frag1->natoms();

  // cerr << "Finding frag with " << ncon1 << " connections, with " << atoms_in_f1 << " atoms\n";

  // First determine the atom types that must be present in the second fragment.
  resizable_array<atom_type_t> types_needed;

  const Molecule& f1m = *frag1->mol();
  for (const Atom& a : f1m) {
    if (a.isotope() == 0) {
      continue;
    }
    types_needed << a.isotope();
  }

  int ncon1 = frag1->ncon();
  if (ncon1 != types_needed.number_elements()) {
    cerr << "Options::SelectFragmentsAtype:fragment has " << ncon << " connections but " <<
            "found " << types_needed << " isotopic atoms\n";
    return 0;
  }

  // Find the matching fragment in `m2` closest to atoms_in_f1.
  // Note that we do not respect the ok_to_select attribute?
  const int nf2 = m2.number_fragments();

  f2 = -1;
  int min_diff = std::numeric_limits<int>::max();
  for (int i = 0; i < nf2; ++i) {
    const SafeFragment* f = m2.fragment(i);
    if (ncon1 != f->ncon()) {
      continue;
    }

    int d = std::abs(f->natoms() - atoms_in_f1);
    if (d < min_diff) {
      min_diff = d;
      f2 = i;
    }
  }

  return f2 >= 0;
#endif
}

int
Options::OkDifferences(const SafeFragment& f1, const SafeFragment& f2) const {
  if (! OkAtomCountDifference(f1.natoms(), f2.natoms())) {
    return 0;
  }

  if (! OkRingCountDifference(f1.nrings(), f2.nrings())) {
    return 0;
  }

  if (! OkFormulaDifference(f1, f2)) {
    return 0;
  }

  if (f1.ncon() == 1) {
    // Singly connected fragments do not have a distance.
  } else if (! OkDistanceDifference(f1.distance(), f2.distance())) {
    return 0;
  }

  return 1;
}

// We are replacing a fragment with `n1` atoms with a fragment
// that contains `n2` atoms.
// Return true if that is OK wrt changes in atom count.
int
Options::OkAtomCountDifference(int n1, int n2) const {

  if (_extra_atoms == 0 && n2 > n1) {
    return 0;
  }

  if (_fewer_atoms == 0 && n2 < n1) {
    return 0;
  }

  if (n2 > n1 && (n2 - n1) > static_cast<int>(*_extra_atoms)) {
    return 0;
  }

  if (n2 < n1 && (n1 - n2) > static_cast<int>(*_fewer_atoms)) {
    return 0;
  }

  return 1;
}

int
Options::OkRingCountDifference(int n1, int n2) const {
  if (_extra_rings) {
    if (n1 < n2 && static_cast<uint32_t>(n2 - n1) > *_extra_rings) {
      return 0;
    }
  }

  if (_fewer_rings) {
    if (n1 > n2 && static_cast<uint32_t>(n1 - n2) < *_fewer_rings) {
      return 0;
    }
  }

  return 1;
}

int
Options::OkFormulaDifference(const SafeFragment& f1, const SafeFragment f2) const {
  if (! _max_formula_difference) {
    return 1;
  }

  return f1.FormulaDifference(f2) <= *_max_formula_difference;
}

int
Options::OkDistanceDifference(int d1, int d2) const {
  if (! _max_distance_difference) {
    return 1;
  }

  if (d1 == d2) {
    return 1;
  }

  if (d1 > d2) {
    return (d1 - d2) <= _max_distance_difference;
  } else {
    return (d2 - d1) <= _max_distance_difference;
  }
}

int
Options::Breed(SafedMolecule& m1, SafedMolecule& m2,
               IWString_and_File_Descriptor& output) {
  int f1_ndx, f2_ndx;
  if (! SelectFragments(m1, m2, f1_ndx, f2_ndx)) {
    return 0;
  }

  const SafeFragment* f1 = m1.fragment(f1_ndx);
  const SafeFragment* f2 = m2.fragment(f2_ndx);

  if (! OkDifferences(*f1, *f2)) {
    return 0;
  }

  IWString replacement_smiles;
  f1->SameNumbers(*f2, replacement_smiles);

#ifdef DEBUG_BREED
  cerr << "Calling FormNewSmiles with '" << m1.smiles() << "'\n";
  cerr << "Replacement " << replacement_smiles << '\n';
#endif
  IWString tmp;
  FormNewSmiles(m1.smiles(), f1_ndx, replacement_smiles, tmp);

  Molecule m;
  if (! m.build_from_smiles(tmp)) {
    cerr << "Options::Breed:invalid smiles '" << tmp << "'\n";
    cerr << "f1_ndx " << f1_ndx << ' ';
    m1.DebugPrint(cerr);
    cerr << "f2_ndx " << f2_ndx << ' ';
    m2.DebugPrint(cerr);
    f1->DebugPrint(cerr);
    f2->DebugPrint(cerr);
    cerr << "Replacement smiles " << replacement_smiles << '\n';
    return 0;
  }

  return ProcessNewMolecule(m, m1.name(), tmp, *f2, output);
}

int
Options::MakeAllLibrary(uint64_t max_make, IWString_and_File_Descriptor& output) {
  for (SafedMolecule* m : _mols) {
    MakeAllLibrary(*m, max_make, output);
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Options::MakeAllLibrary(const SafedMolecule& m,
                        uint64_t max_make, IWString_and_File_Descriptor& output) {
  for (const Library* lib : _library) {
    MakeAllLibrary(m, *lib, max_make, output);
    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Options::MakeAllLibrary(const SafedMolecule& m,
                        const Library& lib,
                        uint64_t max_make,
                        IWString_and_File_Descriptor& output) {
  LibraryIterator lib_iter(lib);
  const SafeFragment* f2;
  uint64_t rc = 0;
  while ((f2 = lib_iter.Next()) != nullptr) {
    std::optional<int> f1_ndx = m.MatchingFragment(f2);
    if (! f1_ndx) {
      continue;
    }
    if (! OkDifferences(*m.fragment(*f1_ndx), *f2)) {
      continue;
    }

    const SafeFragment* f1 = m.fragment(*f1_ndx);

    IWString replacement_smiles;
    f1->SameNumbers(*f2, replacement_smiles);

    IWString tmp;
    FormNewSmiles(m.smiles(), *f1_ndx, replacement_smiles, tmp);

    Molecule newmol;
    if (! newmol.build_from_smiles(tmp)) {
      cerr << "Options::MakeAllLibrary:invalid smiles formed " << tmp <<
              ' ' << m.name() << '\n';
      return 0;
    }

    rc += ProcessNewMolecule(newmol, m.name(), tmp, *f2, output);

    if (rc > max_make) {
      break;
    }
  }

  _number_from_exhaustive += rc;

  if (_verbose > 1) {
    cerr << "Enumerated " << rc << " fragments for " << m.name() << '\n';
  }

  return 1;
}

int
SafeGenerate(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:L:Y:N:a:C:z:X:F:x:n:pb:se:Ik:");

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

  for (const char * fname : cl) {
    if (verbose) {
      cerr << "Reading molecules from '" << fname << "'\n";
    }
    if (! options.ReadMolecules(fname)) {
      cerr << "SafeGenerate::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cerr << "Read " << options.mols().size() << " molecules\n";
  }

  options.SetupRng();

  options.DoSubstructureSearches();

  options.SetIsotope(2);

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('n')) {
    int ngenerate = 1;
    if (! cl.value('n', ngenerate) || ngenerate < 1) {
      cerr << "Invalid number to generate (-n)\n";
      return 1;
    }
    if (verbose) {
      cerr << "Will generate " << ngenerate << " variants\n";
    }

    options.Generate(ngenerate, output);
  }

  if (cl.option_present('b')) {
    int nbreed;
    if (! cl.value('b', nbreed) || nbreed < 1) {
      cerr << "The number to breed option (-b) must be a whole +ve number\n";
      Usage(1);
    }
    if (verbose) {
      cerr << "Will generate " << nbreed << " variants by breeding\n";
    }

    options.Breed(nbreed, output);
  }

  if (cl.option_present('e')) {
    uint64_t ngenerate;
    if (! cl.value('e', ngenerate)) {
      cerr << "The number to exhaustively make (-e) must be a whole +ve number\n";
      Usage(1);
    }
    if (verbose) {
      cerr << "Will generate " << ngenerate << " variants by exhaustive replacement\n";
    }

    options.MakeAllLibrary(ngenerate, output);
  }

  output.flush();


  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace safe_generate

int
main(int argc, char ** argv) {

  int rc = safe_generate::SafeGenerate(argc, argv);

  return rc;
}
