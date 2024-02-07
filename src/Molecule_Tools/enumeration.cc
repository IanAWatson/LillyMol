// Systematic enumeration around a set of molecles.
// Takes a set of molecules and a fragment library,
// and adds all members of the fragment library to
// available sites in the molecules.

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <random>

#include "google/protobuf/text_format.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_set.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/etrans.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_fragments.pb.h"

namespace substituent_enumeration {

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
  cerr << R"(
Performs systematic enumeration of all available sites in a molecule.
Uses a library of fragments to add. These must be either smiles or
textproto represetation of DicerFragment protos.
Fragments added must be marked with an isotope. Only the first isotopic atom
in each fragment is used.
 -F <fname>             Read file of substituents as smiles.
 -F PROTO:<fname>       Read file of substituents as DicerFragment textproto.
 -F                     If any -F directive ends with ',arom' those fragments are only added to aroamtic sites.
 -m <natoms>            Only use substituents with at least <matoms>
 -M <natoms>            Only use substituents with at most <matoms>
 -y <query>             Specify the atom(s) in the starting molecule to which fragments can    be added.
 -n <query>             Specify the atom(s) in the starting molecule to which fragments cannot be added.
 -p                     write the parent molecule before any variants
 -Y ...                 other options, enter '-Y help' for info

 -v             verbose output
)";
// clang-format on

  ::exit(rc);
}

void
DisplayDashYOptions(std::ostream& output) {
  // clang-format off
  output << " -Y xiso           remove isotopic labels from product molecules\n";
  output << " -Y dmv=<n>        do not process any molecule that has more than <n> reactive sites\n";
  output << " -Y dns=<n>        if a molecule has more than <n> sites, randomly downsample to limit to approximately <n>\n";
  output << " -Y valence        check for OK valences in product molecules, and discard bad\n";
  output << " -Y stop=<n>       stop processing a molecule once it has generated <n> products\n";
  output << "                   note that this will be a function of the size of the fragment library as well as sites in the molecule\n";
  // clang-format on

  ::exit(0);
}

// Keep track of the sidechains being added.
// A molecule and the attachment point, and
// possibly with the number of times this is found in
// the originating collection.
class Fragment {
  private:
     Molecule* _frag = nullptr;

  // The atom by which this fragment is attached.
     atom_number_t _zatom = kInvalidAtomNumber;

  // What kind of atom is it, and is it a halogen
     atomic_number_t _atomic_number = kInvalidAtomicNumber;
     bool _is_halogen = false;

  // If this comes from a DicerFragment proto, the number
  // of instances.
     uint32_t _count = 0;

     // If set, this fragment will only attach to aromatic atoms in the
     // starting molecule.
     int _aromatic_only = 0;

  public:
    Fragment();
    ~Fragment();

    const Molecule* frag() const {
      return _frag;
    }

    int SetFragment(Molecule* f);

    int is_halogen() const {
      return _is_halogen;
    }

    atom_number_t zatom() const {
      return _zatom;
    }

    atomic_number_t atomic_number() const {
      return _atomic_number;
    }

    uint32_t count() const {
      return _count;
    }
    void set_count(uint32_t s) {
      _count = s;
    }

    int aromatic_only() const {
      return _aromatic_only;
    }
    void set_aromatic_only(int s) {
      _aromatic_only = s;
    }
};

Fragment::Fragment() {
}

bool
IsHalogem(atomic_number_t z) {
  if (z == 6) {
    return 0;
  }

  if (z == 7) {
    return 0;
  }

  if (z == 8) {
    return 0;
  }

  if (z == 9) {
    return 1;
  }
  if (z == 17) {
    return 1;
  }
  if (z == 35) {
    return 1;
  }
  if (z == 53) {
    return 1;
  }

  return 0;
}


// If we can identify an attachment point in `f` take ownership
// of `f`, otherwise return 0 and it will be up to the caller to
// dispose of `f`.
int
Fragment::SetFragment(Molecule* f) {
  const int matoms = f->natoms();

  for (int i = 0; i < matoms; ++i) {
    if (f->isotope(i) == 0) {
      continue;
    }

    _zatom = i;
    _atomic_number = f->atomic_number(i);
    _is_halogen = IsHalogem(f->atomic_number(i));
    _frag = f;
    return 1;
  }

  return 0;
}

Fragment::~Fragment() {
  if (_frag != nullptr) {
    delete _frag;
  }
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

    // A list of fragments that are to be added.
    resizable_array_p<Fragment> _fragment;

    int _min_atoms_in_fragment = 0;
    int _max_atoms_in_fragment = std::numeric_limits<int>::max();

    // We can specify atoms where substitution is OK or
    // atoms where substitution is not allowed.
    resizable_array_p<Substructure_Query> _ok_attach;
    resizable_array_p<Substructure_Query> _no_attach;

    int _write_parent = 0;

    // By default, we do not check the valence of the products.
    int _check_valence = 0;
    int _products_with_bad_valence = 0;

    // Even with discarding molecules with too many reactive sites and
    // downsampling, we can still end up with a molecule that generates
    // too many products.
    // Remember that this is counting the number of products formed,
    // not the number of sites processed, so it will need to be multiplied
    // by the size of the fragment library.
    int _stop_generating_if_more_than = 0;
    int _stopped_for_too_many_products = 0;

    // Do not generate duplicates.
    IW_STL_Hash_Set _seen;

    int _discard_if_too_many_sites = 0;
    int _discarded_for_too_many_sites = 0;

    // Another way of dealing with large molecules is to only randomly
    // do some of the sites. If there are too many sites, choose them 
    // randomly to get approximately this many.
    int _downsample_threshold = 0;

    std::default_random_engine _generator;

    int _molecules_read = 0;

    int _remove_isotopes = 0;

    // The number of sites in each starting molecule.
    extending_resizable_array<int> _nsites;
    // The number of molecules generated per starting molecule.
    // Generally this will be roughly nsites * number of fragments,
    // although that is an upper value.
    extending_resizable_array<int> _variants_generated;

    // Private functions

    int IdentifyAttachmentPoints(Molecule& m, int * attachment_point);
    int Seen(Molecule& m);
    int GenerateVariants(Molecule& m, atom_number_t zatom,
                          IWString_and_File_Descriptor& output);
    int OkProperties(Molecule& m);

    int ReadFragmentFromTextProto(const dicer_data::DicerFragment& proto, int aromatic_only);
    int ReadFragmentsFromTextProto(IWString& fname, int aromatic_only);
    int ReadFragmentsFromTextProto(iwstring_data_source& input, int aromatic_only);
    int ReadFragmentFromTextProto(const const_IWSubstring& buffer, int aromatic_only);

    int ReadFragmentsAsSmiles(IWString& fname, int aromatic_only);
    int ReadFragmentsAsSmiles(data_source_and_type<Molecule>& input, int aromatic_only);

    int OkFragment(Molecule& frag) const;

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
    if (! cl.value('m', _min_atoms_in_fragment)) {
      cerr << "Invalid min atoms in fragment (-c)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only use fragments with at least " << _min_atoms_in_fragment << " atoms\n";
    }
  }

  if (cl.option_present('M')) {
    if (! cl.value('M', _max_atoms_in_fragment)) {
      cerr << "Invalid max atoms in fragment (-c)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only use fragments with at most " << _max_atoms_in_fragment << " atoms\n";
    }
  }

  if (cl.option_present('F')) {
    IWString fname;
    for (int i = 0; cl.value('F', fname, i); ++i) {
      int only_aromatic = 0;
      if (fname.ends_with(",arom")) {
        fname.chop(5);
        only_aromatic = 1;
      }

      if (fname.starts_with("PROTO:")) {
        fname.remove_leading_chars(6);
        if (! ReadFragmentsFromTextProto(fname, only_aromatic)) {
          cerr << "Options::Initialise:cannot read testproto file '" << fname << "'\n";
          return 0;
        }
      } else {
        if (! ReadFragmentsAsSmiles(fname, only_aromatic)) {
          cerr << "Options::Initialise:cannot read smiles file '" << fname << "'\n";
          return 0;
        }
      }
    }
  }

  if (_fragment.empty()) {
    cerr << "Options::Initialise:no fragments - use the -F option to specify\n";
    return 0;
  }

  if (cl.option_present('y')) {
    if (! process_queries(cl, _ok_attach, _verbose, 'y')) {
      cerr << "Options::Initialise:cannot assemble ok attachment point queries (-y)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _ok_attach.size() << " ok attachment point queries\n";
    }
  }

  if (cl.option_present('n')) {
    if (! process_queries(cl, _no_attach, _verbose, 'n')) {
      cerr << "Options::Initialise:cannot assemble not ok attachment point queries (-n)\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Defined " << _no_attach.size() << " ok queries for excluding attachment points\n";
    }
  }

  if (cl.option_present('p')) {
    _write_parent = 1;
    if (_verbose) {
      cerr << "Will write the parent molecule\n";
    }
  }

  if (cl.option_present('Y')) {
    IWString y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "xiso") {
        _remove_isotopes = 1;
        if (_verbose) {
          cerr << "Will remove isotopic labels from products\n";
        }
      } else if (y.starts_with("dmv=")) {
        y.remove_leading_chars(4);
        if (! y.numeric_value(_discard_if_too_many_sites)) {
          cerr << "Invalid discard if too many sites directive 'dmv=" << y << "'\n";
          return 0;
        }
      } else if (y.starts_with("dns=")) {
        y.remove_leading_chars(4);
        if (! y.numeric_value(_downsample_threshold)) {
          cerr << "Invalid downsample threshold 'dns=" << y << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will downsample molecules with more than " << _downsample_threshold << " sites\n";
        }
      } else if (y == "valence") {
        _check_valence = 1;
        if (_verbose) {
          cerr << "Will discard products with bad valences\n";
        }
      } else if (y.starts_with("stop=")) {
        y.remove_leading_chars(5);
        if (! y.numeric_value(_stop_generating_if_more_than)) {
          cerr << "Invalid stop= directive '" << y << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Will abandon an input molecule if it generates more than " << _stop_generating_if_more_than << " products\n";
        }
      } else if (y == "help") {
        DisplayDashYOptions(cerr);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(cerr);
      }
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  if (_discarded_for_too_many_sites) {
    output << _discarded_for_too_many_sites << " input molecules discarded for more than " << _discard_if_too_many_sites << " sites\n";
  }
  if (_check_valence) {
    output << "Discarded " << _products_with_bad_valence << " products with bad valence\n";
  }
  if (_stop_generating_if_more_than) {
    output << _stopped_for_too_many_products << " enumerations stopped for generating " <<  _stop_generating_if_more_than << " products\n";
  }
  Accumulator_Int<int> acc;
  for (int i = 0; i < _nsites.number_elements(); ++i) {
    if (_nsites[i]) {
      output << _nsites[i] << " molecules had " << i << " sites for substitution\n";
      acc.extra(i, _nsites[i]);
    }
  }
  output << " ave " << acc.average() << " sites per starting molecule\n";
  for (int i = 0; i < _variants_generated.number_elements(); ++i) {
    if (_variants_generated[i]) {
      output << _variants_generated[i] << " molecules generated " << i << " variants\n";
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

  if (matoms == 0) {
    return 1;
  }

  if (_write_parent) {
    static constexpr char kSep = ' ';

    output << m.smiles() << kSep << m.name() << '\n';
  }

  std::unique_ptr<int[]> attachment_point = std::make_unique<int[]>(matoms);
  if (_ok_attach.size() || _no_attach.size()) {
    if (! IdentifyAttachmentPoints(m, attachment_point.get())) {
      cerr << "Options::Process:cannot identify attachment point atoms\n";
      return 0;
    }
  } else {
    std::fill_n(attachment_point.get(), matoms, 1);
  }

  for (int i = 0; i < matoms; ++i) {
    if (m.hcount(i) == 0) {
      attachment_point[i] = 0;
    }
  }

  const int nsites = std::count(attachment_point.get(), attachment_point.get() + matoms, 1);

  ++_nsites[nsites];

  if (_discard_if_too_many_sites) {
    if (nsites > _discard_if_too_many_sites) {
      ++_discarded_for_too_many_sites;
      return 0;
    }
  }
  // cerr << "Processing '" << m.name() << "' nsites " << nsites << "\n";

  std::unique_ptr<std::bernoulli_distribution> rng;
  if (nsites > _downsample_threshold) {
    float fraction = iwmisc::Fraction<double>(_downsample_threshold, nsites);
    rng = std::make_unique<std::bernoulli_distribution>(fraction);
  }

  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    if (! attachment_point[i]) {
      continue;
    }

    if (rng && ! (*rng)(_generator)) {
      continue;
    }

    rc += GenerateVariants(m, i, output);

    if(_stop_generating_if_more_than && rc > _stop_generating_if_more_than) {
      ++_stopped_for_too_many_products;
      break;
    }
  }

  output.write_if_buffer_holds_more_than(4092);

  ++_variants_generated[rc];

  return rc;
}

int
IdentifyAtoms(Molecule_to_Match& target,
              resizable_array_p<Substructure_Query>& queries,
              int* attachment_point,
              int flag) {
  int rc = 0;
  for (Substructure_Query* q : queries) {
    Substructure_Results sresults;
    if (q->substructure_search(target, sresults)) {
      sresults.each_embedding_set_vector(attachment_point, flag);
      ++rc;
    }
  }

  return rc;
}

int
Options::IdentifyAttachmentPoints(Molecule& m, int * attachment_point) {
  const int matoms = m.natoms();

  Molecule_to_Match target(&m);

  if (_ok_attach.size()) {
    std::fill_n(attachment_point, matoms, 0);
    return IdentifyAtoms(target, _ok_attach, attachment_point, 1);
  } else {
    std::fill_n(attachment_point, matoms, 1);
    return IdentifyAtoms(target, _no_attach, attachment_point, 0);
  }
}

// Return true if it is OK to join `zatom` in `m` to `frag`.
int
OkBondFormation(Molecule& m,
                atom_number_t zatom,
                const Fragment& frag) {
  if (frag.is_halogen() && ! m.is_aromatic(zatom)) {
    return 0;
  }

  if (frag.atomic_number() == 6) {
    return 1;
  }

  const atomic_number_t z = m.atomic_number(zatom);
  if (z == 6) {
    return 1;
  }

  // Must be two heteratoms, which we do not allow.
  return 0;
}

int
Options::GenerateVariants(Molecule& m, atom_number_t zatom,
                          IWString_and_File_Descriptor& output) {
  const int initial_matoms = m.natoms();

  static constexpr char kSep = ' ';

  static IWString percents(" %% ");

  // cerr << "For atom " << zatom << " in " << m.name() << " will generate " << _fragment.size() << " variants\n";
  int rc = 0;
  for (const Fragment* frag : _fragment) {
    if (! OkBondFormation(m, zatom, *frag)) {
      continue;
    }

    if (frag->aromatic_only() && ! m.is_aromatic(zatom)) {
      continue;
    }

    Molecule mcopy(m);
    mcopy += *frag->frag();
    mcopy.add_bond(zatom, initial_matoms + frag->zatom(), SINGLE_BOND);
    mcopy.unset_all_implicit_hydrogen_information(zatom);
    mcopy.unset_all_implicit_hydrogen_information(initial_matoms + frag->zatom());

    if (_remove_isotopes) {
      mcopy.unset_isotopes();
    }

    if (Seen(mcopy)) {
      continue;
    }

    if (! OkProperties(mcopy)) {
      continue;
    }

    mcopy.invalidate_smiles();

    if (_check_valence && ! mcopy.valence_ok()) {
      ++_products_with_bad_valence;
      cerr << "Bad valence " << mcopy.smiles() << ' ' << m.name() << '\n';
      continue;
    }

    output << mcopy.smiles() << kSep << m.name() << percents << frag->frag()->name();
    if (frag->count() > 0) {
      output << kSep << frag->count();
    }
    output << '\n';
    ++rc;
  }

  return rc;
}

int
Options::Seen(Molecule& m) {
  if (_seen.contains(m.unique_smiles())) {
    return 1;
  }

  _seen.insert(m.unique_smiles());

  return 0;
}

int
Options::ReadFragmentsFromTextProto(IWString& fname, int only_aromatic) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadFragmentsFromTextproto:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragmentsFromTextProto(input, only_aromatic);
}

int
Options::ReadFragmentsFromTextProto(iwstring_data_source& input, int only_aromatic) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadFragmentFromTextProto(buffer, only_aromatic)) {
      cerr << "Options::ReadFragmentsFromTextProto:cannot process " << buffer << '\n';
      return 0;
    }
  }

  return _fragment.size();
}

int
Options::ReadFragmentFromTextProto(const const_IWSubstring& buffer, int only_aromatic) {
  google::protobuf::io::ArrayInputStream input(buffer.data(), buffer.length());

  dicer_data::DicerFragment proto;
  if (! google::protobuf::TextFormat::Parse(&input, &proto)) {
    cerr << "DicerFragmentLookupImpl::Lookup:invalid db contents '";
    cerr << buffer << '\n';
    return 0;
  }

  return ReadFragmentFromTextProto(proto, only_aromatic);
}

int
Options::ReadFragmentFromTextProto(const dicer_data::DicerFragment& proto, int only_aromatic) {
  std::unique_ptr<Fragment> fragment = std::make_unique<Fragment>();

  Molecule* f = new Molecule();
  if (! f->build_from_smiles(proto.smi())) {
    cerr << "Options::ReadFragmentFromTextProto:invalid smiles\n";
    cerr << proto.ShortDebugString() << '\n';
    delete f;
    return 0;
  }

  if (! fragment->SetFragment(f)) {
    delete f;
    return 0;
  }

  if (! OkFragment(*f)) {
    return 0;
  }

  fragment->set_count(proto.n());
  fragment->set_aromatic_only(only_aromatic);

  _fragment << fragment.release();
  

  return 1;
}

int
Options::OkFragment(Molecule& frag) const {
  const int matoms = frag.natoms();
  if (matoms < _min_atoms_in_fragment) {
    return 0;
  }

  if (matoms > _max_atoms_in_fragment) {
    return 0;
  }

  return 1;
}

int
Options::ReadFragmentsAsSmiles(IWString& fname, int only_aromatic) {
  data_source_and_type<Molecule> input(FILE_TYPE_SMI, fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadFragmentsAsSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFragmentsAsSmiles(input, only_aromatic);
}

int
Options::ReadFragmentsAsSmiles(data_source_and_type<Molecule>& input, int only_aromatic) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Fragment> fragment = std::make_unique<Fragment>();
    if (! fragment->SetFragment(m)) {
      cerr << "Options::ReadFragmentsAsSmiles:cannot identify attachment in " <<
              m->smiles() << ' ' << m->name() << '\n';
      return 0;
    }

    if (! OkFragment(*m)) {
      return 0;
    }

    fragment->set_aromatic_only(only_aromatic);

    _fragment << fragment.release();
  }

  return _fragment.size();
}
int
Options::OkProperties(Molecule& m) {
  return 1;
}

int
SubstituentEnumeration(Options& options,
                Molecule& m,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
SubstituentEnumeration(Options& options,
                data_source_and_type<Molecule>& input,
                IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! options.Preprocess(*m)) {
      continue;
    }

    SubstituentEnumeration(options, *m, output);
  }

  return 1;
}

int
SubstituentEnumeration(Options& options,
             const char * fname,
             FileType input_type,
             IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "SubstituentEnumeration:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return SubstituentEnumeration(options, input, output);
}

int
SubstituentEnumeration(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:T:A:lcg:i:F:m:M:y:n:pY:");

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
    if (! SubstituentEnumeration(options, fname, input_type, output)) {
      cerr << "SubstituentEnumeration::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace substituent_enumeration

int
main(int argc, char ** argv) {

  int rc = substituent_enumeration::SubstituentEnumeration(argc, argv);

  return rc;
}
