// Extract linkers from molecules

#include <iostream>
#include <limits>
#include <memory>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "Molecule_Tools.pb.h"
#endif

namespace smi2linker {

using std::cerr;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
  cerr << " -m <natoms>      discard fragments having fewer than <natoms> atoms\n";
  cerr << " -M <natoms>      discard fragments having more  than <natoms> atoms\n";
  cerr << " -P <atype>       standard atom typing specifications\n";
  cerr << " -I ...           fragments can be isotopically labelled, enter '-I help'\n";
  cerr << " -n               suppress normal output - for when the -F option is in use\n";
  cerr << " -F <fname>       write dicer_data::DicerFragment protos to <fname>\n";
  cerr << " -R <n>           report progress every <n> molecules processed\n";
  cerr << " -c               discard chirality\n";
  cerr << " -l               strip to largest fragment\n";
  cerr << " -v               verbose output\n";

  ::exit(rc);
}

// Different isotopic labels can be applied.
enum class Isotope {
  kNone,
  kConnection,
  kAtomicNumber,
};

class Options {
  private:
    int _verbose = 0;

    FileType _input_type = FILE_TYPE_INVALID;

    int _reduce_to_largest_fragment = 0;

    int _remove_chirality = 0;

    Chemical_Standardisation _chemical_standardisation;

    Element_Transformations _element_transformations;

    int _molecules_read = 0;

    int _molecules_with_no_rings = 0;

    int _fragments_generated;

    int _min_size;
    int _max_size;

    // We can choose to reject linkers that constitute a
    // large fraction of the atoms in the parent.
    float _max_fraction_parent_atoms;

    // It may be interesting to characterise zero length
    // linkers - bonds between two rings. The only way this
    // makes sense is to add the dummy atoms.
    int _process_zero_length_linkers;

    // One kind of annotation is to put an isotope on the
    // atoms that join to the rings.
    int _isotope_for_join_points;

    // Another annotation is to add Ar and Al atoms to the
    // fragment to designate whether it is an aromatic or
    // aromatic ring attached.
    int _alar_to_join_points;

    // As a modification of this, the Ar and Al atoms can have
    // an isotopic label indicating the ring size.
    int _label_alar_by_ring_size;

    // We can add a dummy atom at each join point, and the
    // isotope will be the atom type of the atom that used to
    // be there.
    int _atype_at_adjacent_atoms;

    // Token separator used in output files.
    IWString _sep;

    // By default, we write the fragments we generate to
    // the output stream. The alternative is to accumulate
    // data on fragments and write a summary proto.
    int _write_fragments;

    // If we are not writing fragments, we might be accumulating
    // fragment info and writing dicer_data::DicerFragment protos.
    int _accumulate_fragments;

    // If we are accumulating fragments, the destination.
    IWString_and_File_Descriptor _stream_for_fragment_summmary;

    // If atom typing is active.
    Atom_Typing_Specification _atom_typing;

    // To avoid passing around too many arguments, some variables
    // that are specific to the current molecule being processed.

    int* _xref;

    // The number of fragments generated by the current molecule.
    int _fragments_current_molecule;
    // Global statistics on how many fragments generated.
    extending_resizable_array<int> _fragments_per_molecule;

    // The assigned atom types of the current molecule.
    std::unique_ptr<uint32_t[]> _atype;

    // In order to create the global fingerprint we need a mapping
    // from smiles to data about that fragment.
    IW_STL_Hash_Map<IWString, dicer_data::DicerFragment> _fragments;

    Report_Progress _report_progress;

  // Private functions
    int InitialiseCurrentMoleculeData(Molecule& m);

    int OkSize(int natoms, int atoms_in_parent) const;

    int AddAlArToAdjacentAtoms(Molecule& m,
                       atom_number_t zatom,
                       const int* fragment_membership,
                       int fragment_number,
                       Molecule& frag);
    int AddDummyAtom(Molecule& m,
                       atom_number_t zatom,
                       const int* fragment_membership,
                       int fragment_number,
                       Molecule& frag) const;

    int ProcessFragment(Molecule& m,
                         const int* fragment_membership,
                         int fragment_number,
                         const int* join_point,
                         const int * adjacent_atoms,
                         IWString_and_File_Descriptor& output);
    int ProcessFragment(Molecule& frag,
                        dicer_data::Isotope iso_type,
                         IWString_and_File_Descriptor& output);

    int ProcessZeroLengthLinkers(Molecule& parent,
                        const int * is_ring,
                        IWString_and_File_Descriptor& output);
    int ProcessZeroLengthLinker(Molecule& parent,
                        atom_number_t a1,
                        atom_number_t a2,
                        IWString_and_File_Descriptor& output);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int MaybeDiscernInputType(const char * fname);

    FileType input_type() const {
      return _input_type;
    }

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;

    int verbose() const {
      return _verbose;
    }

    // This must be public since it is called after the last
    // molecule has been processed.
    // If accumulating fragments has not been requested, it is noop.
    int WriteFragmentStatisticsAsProto();
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _input_type = FILE_TYPE_INVALID;
  _molecules_read = 0;
  _molecules_with_no_rings = 0;

  _process_zero_length_linkers = 0;

  _min_size = 0;
  _max_size = std::numeric_limits<int>::max();
  // Initialised to something > 1.0
  _max_fraction_parent_atoms = 2.0f;

  _write_fragments = 1;
  _accumulate_fragments = 0;

  _fragments_generated = 0;

  _isotope_for_join_points = 0;
  _alar_to_join_points = 0;
  _label_alar_by_ring_size = 0;
  _atype_at_adjacent_atoms = 0;

  _sep = ' ';

  _xref = nullptr;
}

Options::~Options() {
  if (_xref != nullptr) {
    delete [] _xref;
  }
}

void
DisplayDashIQualifiers(std::ostream& output) {
  // clang-format off
  output << "The -I option controls isotopic labels applied to fragments\n";
  output << " -I join          isotope 1 applied to fragment atoms that join the rest of the molecule\n";
  output << " -I alar          Al (aliphatic) and Ar (arom) atms attached to the fragment\n";
  output << " -I alarsize      Al (aliphatic) and Ar (arom) atms attached to the fragment\n";
  output << "                  they are assigned an isotope based on the size of the ring\n";
  output << " -I atype         Dummy Ti atoms added with isotope being the atom type\n";
  output << "                  of the attached atom (non fragment) atom type\n";
  // clang-format on

  exit(0);
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
    if (! cl.value('m', _min_size) || _min_size < 0) {
      cerr << "The min fragment size option (-m) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only process fragments with at least " << _min_size << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_size) || _max_size < _min_size) {
      cerr << "The max fragment size option (-M) must be a whole +ve number >= " << _min_size << "\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will only process fragments with at most " << _max_size << " atoms\n";
    }
  }

  if (cl.option_present('f')) {
    if (! cl.value('f', _max_fraction_parent_atoms) ||
        _max_fraction_parent_atoms <= 0.0f ||
        _max_fraction_parent_atoms >= 1.0f) {
      cerr << "The max fraction parent atoms (-f) option must be a valid fraction\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will discard fragments that constitute > " << _max_fraction_parent_atoms << " of the parent atom count\n";
    }
  }

  int isotope_specifications = 0;
  if (cl.option_present('P')) {
    const IWString atype = cl.string_value('P');
    if (! _atom_typing.build(atype)) {
      cerr << "Invalid atom typeing '" << atype << "'\n";
      return 0;
    }
    _atype_at_adjacent_atoms = 1;
    ++isotope_specifications;
  }

  if (cl.option_present('z')) {
    _process_zero_length_linkers = 1;
    if (_verbose) {
      cerr << "Will process zero length linkers (biphenyl-like rings)\n";
    }
  }

  if (cl.option_present('R')) {
    if (! _report_progress.initialise(cl, 'R', _verbose)) {
      cerr << "Cannot initialise progress reporting (-R)\n";
      return 0;
    }
  }

  if (cl.option_present('I')) {
    IWString s;
    for (int i = 0; cl.value('I', s, i); ++i) {
      if (s == "join") {
        _isotope_for_join_points = 1;
        ++isotope_specifications;
      } else if (s == "alar") {
        _alar_to_join_points = 1;
        ++isotope_specifications;
      } else if (s == "alarsize") {
        _alar_to_join_points = 1;
        _label_alar_by_ring_size = 1;
        ++isotope_specifications;
      } else if (s == "atype") {
        _atype_at_adjacent_atoms = 1;
        if (! _atom_typing.active()) {
          cerr << "Atom typing not active, cannot do atom type isotope\n";
          return 0;
        }
        if (! cl.option_present('P')) {
          ++isotope_specifications;
        }
      } else if (s == "help") {
        DisplayDashIQualifiers(cerr);
        return 0;
      } else {
        cerr << "Unrecognised -I qualifier '" << s << "'\n";
        return 0;
      }
    }
  }

  if (isotope_specifications > 1) {
    cerr << "Multiple isotopic label specifications, cannot process\n";
    return 0;
  }

  if (cl.option_present('n')) {
    _write_fragments = 0;
    if (_verbose) {
      cerr << "Normal output suppressed\n";
    }
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    if (! _stream_for_fragment_summmary.open(fname.null_terminated_chars())) {
      cerr << "Options::Initialise:cannot open -F file '" << fname << "'\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will write fragment summary data to " << fname << '\n';
    }
    _accumulate_fragments = 1;
  }

  if (cl.option_present('i')) {
    if (! process_input_type(cl, _input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (1 == cl.number_elements() && 0 == strcmp(cl[0], "-")) {
    _input_type = FILE_TYPE_SMI;
  } else if (! all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  return 1;
}

int
Options::OkSize(int natoms, int atoms_in_parent) const {
  if (natoms < _min_size) {
    return 0;
  }
  if (natoms > _max_size) {
    return 0;
  }

  if (_max_fraction_parent_atoms >= 1.0f) {
    return 1;
  }

  float f = iwmisc::Fraction<float>(natoms, atoms_in_parent);
  if (f >= _max_fraction_parent_atoms) {
    return 0;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << _molecules_with_no_rings << " molecules had no rings\n";
  output << "Generated " << _fragments_generated << " fragments\n";
  for (int i = 0; i < _fragments_per_molecule.number_elements(); ++i) {
    if (_fragments_per_molecule[i]) {
      output << _fragments_per_molecule[i] << " molecules generated " << i << " fragments\n";
    }
  }
  return 1;
}

int
Options::MaybeDiscernInputType(const char * fname) {
  if (_input_type == FILE_TYPE_INVALID) {
    _input_type = discern_file_type_from_name(fname);
  }
  return 1;
}

int
Options::Preprocess(Molecule& m) {
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

// Recursively identify a linker region.
void
IdentifyFragment(const Molecule& m,
                 atom_number_t zatom,
                 const int* is_ring, 
                 int* already_done,
                 int fragment_number,
                 int* join_point,
                 int* adjacent_atoms) {
  const Atom* a = m.atomi(zatom);
  already_done[zatom] = fragment_number;
  // cerr << "Continue with atom " << zatom << '\n';

  for (const Bond * b : *a) {
    atom_number_t j = b->other(zatom);

    // cerr << "From atom " << zatom << " check " << j << " ring " << is_ring[j] << " already_done " << already_done[j] << '\n';
    if (is_ring[j]) {
      join_point[zatom] = fragment_number;
      adjacent_atoms[j] = fragment_number;
      continue;
    }

    if (already_done[j]) {
      continue;
    }

    IdentifyFragment(m, j, is_ring, already_done, fragment_number,
                     join_point, adjacent_atoms);
  }
}

std::unique_ptr<int[]>
IsRing(Molecule& m) {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> is_ring = std::make_unique<int[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    if (m.ring_bond_count(i) > 0) {
      is_ring[i] = 1;
    } else {
      is_ring[i] = 0;
    }
  }

  // Add doubly bonded atoms to rings.
  for (int i = 0; i < matoms; ++i) {
    if (is_ring[i]) {
      continue;
    }
    const Atom* a = m.atomi(i);
    if (a->ncon() != 1) {
      continue;
    }

    for (const Bond* b : *a) {
      if (! b->is_double_bond()) {
        continue;
      }

      atom_number_t j = b->other(i);
      if (is_ring[j]) {
        is_ring[i] = 1;
        break;
      }
    }
  }

  return is_ring;
}

int
AtomsInFragment(const int* already_done,
                int matoms,
                int fragment_number) {
  return std::count(already_done, already_done + matoms, fragment_number);
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;
  _fragments_current_molecule = 0;

  if (_report_progress()) {
    cerr << "Processed " << _molecules_read << " molecules, generated " << _fragments_generated << " fragments\n";
    if (_accumulate_fragments) {
    }
  }

  const int matoms = m.natoms();
  if (matoms == 0) {
    return 1;
  }

  if (m.nrings() == 0) {
    ++_molecules_with_no_rings;
    return 1;
  }

  // Write the parent, and a magic token to make this output
  // look like it came from dicer.
  if (_write_fragments) {
    output << m.smiles() << _sep << m.name() << " B=1\n";
  }

  InitialiseCurrentMoleculeData(m);

  std::unique_ptr<int[]> is_ring = IsRing(m);

  if (_process_zero_length_linkers && m.nrings() > 1) {
    ProcessZeroLengthLinkers(m, is_ring.get(), output);
  }

  std::unique_ptr<int[]> already_done(new_int(matoms));
  std::unique_ptr<int[]> join_point(new_int(matoms));
  std::unique_ptr<int[]> adjacent_atoms(new_int(matoms));

  int fragment_number = 1;
  for (int i = 0; i < matoms; ++i, ++fragment_number) {
    if (already_done[i] || is_ring[i]) {
      continue;
    }
    // cerr << "Begin fragment at " << i << '\n';
    IdentifyFragment(m, i, is_ring.get(), already_done.get(),
                     fragment_number, join_point.get(), adjacent_atoms.get());
    int atoms_in_fragment = AtomsInFragment(already_done.get(), matoms, fragment_number);
    if (! OkSize(atoms_in_fragment, m.natoms())) {
      continue;
    }
    ProcessFragment(m, already_done.get(), fragment_number,
                    join_point.get(), adjacent_atoms.get(),
                    output);
  }

  ++_fragments_per_molecule[_fragments_current_molecule];

  return 1;
}

int
Options::ProcessZeroLengthLinkers(Molecule& parent,
                        const int * is_ring,
                        IWString_and_File_Descriptor& output) {
  int rc = 0;
  for (const Bond* b : parent.bond_list()) {
    const atom_number_t a1 = b->a1();
    if (! is_ring[a1]) {
      continue;
    }
    const atom_number_t a2 = b->a2();
    if (! is_ring[a2]) {
      continue;
    }
    // We now have two atoms that are both rings.
    if (parent.in_same_ring(a1, a2)) {
      continue;
    }
    ProcessZeroLengthLinker(parent, a1, a2, output);
    ++rc;
  }

  return rc;
}

int
Options::ProcessZeroLengthLinker(Molecule& parent,
                        atom_number_t a1,
                        atom_number_t a2,
                        IWString_and_File_Descriptor& output) {
  static const Element* Al = get_element_from_atomic_number(13);
  static const Element* Ar = get_element_from_atomic_number(18);

  Molecule frag;
  frag.set_name(parent.name());
  if (_alar_to_join_points) {
    if (parent.is_aromatic(a1)) {
      frag.add(Ar);
    } else {
      frag.add(Al);
    }
    if (parent.is_aromatic(a2)) {
      frag.add(Ar);
    } else {
      frag.add(Al);
    }

    // We ignore the case where there might be a double bond between rings.
    frag.add_bond(0, 1, SINGLE_BOND);
    if (_label_alar_by_ring_size) {
      const Ring* r = parent.ring_containing_atom(a1);
      frag.set_isotope(0, r->size());
      r = parent.ring_containing_atom(a2);
      frag.set_isotope(1, r->size());
    }
    return ProcessFragment(frag, dicer_data::ALAR, output);
  }

  if (_atype_at_adjacent_atoms) {
    static const Element* Ti = get_element_from_atomic_number(22);
    frag.add(Ti);
    frag.add(Ti);
    frag.add_bond(0, 1, SINGLE_BOND);
    frag.set_isotope(0, _atype[a1]);
    frag.set_isotope(1, _atype[a1]);
    return ProcessFragment(frag, dicer_data::ATYPE, output);
  }

  cerr << "Not sure what typing to apply for zero length connectors\n";
  return 0;
}

int
Options::InitialiseCurrentMoleculeData(Molecule& m) {
  if (_xref != nullptr) {
    delete [] _xref;
  }
  _xref = new int[m.natoms()];

  if (_atype_at_adjacent_atoms) {
    _atype.reset(new uint32_t[m.natoms()]);
    _atom_typing.assign_atom_types(m, _atype.get());
  }

  return 1;
}

// Atom `zatom` in `m` is adjacent to one or more atoms that
// have been included in a fragment `fragment_membership[fragment_number]`
// Add Ar or Al atoms to `frag` at those points.
int
Options::AddAlArToAdjacentAtoms(Molecule& m,
                       atom_number_t zatom,
                       const int* fragment_membership,
                       int fragment_number,
                       Molecule& frag) {
  static const Element* Al = get_element_from_atomic_number(13);
  static const Element* Ar = get_element_from_atomic_number(18);

  assert(fragment_membership[zatom] != fragment_number);

  const Atom* a = m.atomi(zatom);

  int rc = 0;
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (fragment_membership[j] != fragment_number) {
      continue;
    }

    if (m.is_aromatic(zatom)) {
      frag.add(Ar);
    } else {
      frag.add(Al);
    }
    frag.add_bond(frag.natoms() - 1, _xref[j], SINGLE_BOND);
    if (_label_alar_by_ring_size) {
      const Ring* ri = m.ring_containing_atom(zatom);
      frag.set_isotope(frag.natoms() - 1, ri->size());
    }
    ++rc;
  }

  return rc;
}

// `zatom` has been marked as adjacent to a fragment atom.
// Identify the bonded atoms that are in the fragment, and
// add a dummy atom with zatom's atom type.
int
Options::AddDummyAtom(Molecule& m,
                       atom_number_t zatom,
                       const int* fragment_membership,
                       int fragment_number,
                       Molecule& frag) const {
  static const Element* Ti = get_element_from_atomic_number(22);

  assert(fragment_membership[zatom] != fragment_number);

  const Atom* a = m.atomi(zatom);

  int rc = 0;
  for (const Bond* b : *a) {
    atom_number_t j = b->other(zatom);
    if (fragment_membership[j] != fragment_number) {
      continue;
    }

    frag.add(Ti);
    frag.add_bond(frag.natoms() - 1, _xref[j], SINGLE_BOND);
    frag.set_isotope(frag.natoms() - 1, _atype[zatom]);
    ++rc;
  }

  return rc;
}

int
Options::ProcessFragment(Molecule& m,
                         const int* fragment_membership,
                         int fragment_number,
                         const int* join_point,
                         const int * adjacent_atoms,
                         IWString_and_File_Descriptor& output) {
  const int matoms = m.natoms();
  ++_fragments_current_molecule;

  Molecule frag;
  m.create_subset(frag, fragment_membership, fragment_number, _xref);
  frag.set_name(m.name());
  if (_isotope_for_join_points) {
    for (int i = 0; i < matoms; ++i) {
      if (join_point[i] == fragment_number) {
        frag.set_isotope(_xref[i], _isotope_for_join_points);
      }
    }
    return ProcessFragment(frag, dicer_data::ATT, output);
  }

  if (_alar_to_join_points) {
    for (int i = 0; i < matoms; ++i) {
      if (adjacent_atoms[i] != fragment_number) {
        continue;
      }
      AddAlArToAdjacentAtoms(m, i, fragment_membership,
                             fragment_number, frag);
    }
    return ProcessFragment(frag, dicer_data::ALAR, output);
  }

  if (_atype_at_adjacent_atoms) {
    for (int i = 0; i < matoms; ++i) {
      if (adjacent_atoms[i] != fragment_number) {
        continue;
      }
      AddDummyAtom(m, i, fragment_membership, fragment_number, frag);
    }
    return ProcessFragment(frag, dicer_data::ATYPE, output);
  }

  return ProcessFragment(frag, dicer_data::NONE, output);
}

int
Options::ProcessFragment(Molecule& frag,
                         dicer_data::Isotope iso_type,
                         IWString_and_File_Descriptor& output) {
  ++_fragments_generated;

  const IWString& usmi = frag.unique_smiles();

  if (_write_fragments) {
    output << usmi << _sep << frag.name() << _sep << frag.natoms() << '\n';
    output.write_if_buffer_holds_more_than(4096);
  }

  if (! _accumulate_fragments) {
    return 1;
  }

  auto iter = _fragments.find(usmi);
  if (iter != _fragments.end()) {
    auto current_count = iter->second.n();
    iter->second.set_n(current_count + 1);
    return 1;
  }

  auto [ins, _] = _fragments.emplace(std::make_pair(usmi, dicer_data::DicerFragment()));

  ins->second.set_smi(usmi.AsString());
  ins->second.set_par(frag.name().AsString());
  ins->second.set_n(1);
  ins->second.set_nat(frag.natoms());
  ins->second.set_iso(iso_type);

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::WriteFragmentStatisticsAsProto() {
  if (! _stream_for_fragment_summmary.active()) {
    return 1;
  }

  static google::protobuf::TextFormat::Printer printer;  
  printer.SetSingleLineMode(true);
  std::string buffer;

  for (const auto& [usmi, proto] : _fragments) {
    if (! printer.PrintToString(proto, &buffer)) {
      cerr << "WriteFragmentStatisticsAsProto:cannot write '" << proto.ShortDebugString() << "'\n";
      return 0;
    }
    // For loading a Berkeley db database, write the key first.
    _stream_for_fragment_summmary << usmi << ' ';

    _stream_for_fragment_summmary << buffer;
    _stream_for_fragment_summmary << '\n';
    _stream_for_fragment_summmary.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
Smi2Linker(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
Smi2Linker(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    if (! Smi2Linker(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
Smi2Linker(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  options.MaybeDiscernInputType(fname);

  data_source_and_type<Molecule> input(options.input_type(), fname);
  if (! input.good()) {
    cerr << "Smi2Linker:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return Smi2Linker(options, input, output);
}

int
Smi2Linker(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:clm:M:I:nF:P:zR:f:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
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
    if (! Smi2Linker(options, fname, output)) {
      cerr << "Smi2Linker::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  options.WriteFragmentStatisticsAsProto();

  return 0;
}

}  // namespace smi2linker

int
main(int argc, char ** argv) {

  int rc = smi2linker::Smi2Linker(argc, argv);

  return rc;
}
