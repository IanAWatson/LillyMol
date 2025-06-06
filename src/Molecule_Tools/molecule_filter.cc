// Molecule Filtering Tool.
// Designed as a first level filter for large collections.

#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/rotbond_common.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/molecule_filter_lib.h"
#include "Molecule_Tools/nvrtspsa.h"
#include "Molecule_Tools/xlogp.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/molecule_filter.pb.h"
#else
#include "Molecule_Tools.pb.h"
#endif

namespace molecule_filter {

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
// clang-format off
  cerr << R"(Filters molecules based on easy to compute molecular properties.
Designed for rapidly triaging large collections of molecules in the most
computationally efficient way possible.
 -F <fname>     textproto describing constraints on filter - required.
 -c             remove chirality
 -l             reduce to largest fragment
 -B <fname>     write rejected molecules to <fname>
 -v             verbose output
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

    alogp::ALogP _alogp;

    Chemical_Standardisation _chemical_standardisation;

    uint64_t _molecules_read = 0;
    uint64_t _molecules_passed = 0;

    MoleculeFilterData::Requirements _requirements;

    quick_rotbond::QuickRotatableBonds _rotbond;

    IWString_and_File_Descriptor _reject_stream;

    // Values accumulated based on rejections
    uint64_t _too_few_atoms = 0;
    uint64_t _too_many_atoms = 0;
    uint64_t _too_few_rings = 0;
    uint64_t _too_many_rings = 0;
    uint64_t _too_few_heteroatoms = 0;
    uint64_t _min_heteroatom_fraction = 0;
    uint64_t _max_heteroatom_fraction = 0;
    uint64_t _too_few_aromatic_rings = 0;
    uint64_t _too_many_aromatic_rings = 0;
    uint64_t _too_few_aliphatic_rings = 0;
    uint64_t _too_many_aliphatic_rings = 0;
    uint64_t _ring_system_too_large = 0;
    uint64_t _too_many_aromatic_rings_in_system = 0;
    uint64_t _ring_too_large = 0;
    uint64_t _non_organic = 0;
    uint64_t _isotope = 0;
    uint64_t _too_few_rotbond = 0;
    uint64_t _too_many_rotbond = 0;
    uint64_t _low_tpsa = 0;
    uint64_t _high_tpsa = 0;
    uint64_t _low_xlogp = 0;
    uint64_t _high_xlogp = 0;
    uint64_t _low_alogp = 0;
    uint64_t _high_alogp = 0;
    uint64_t _too_few_hba = 0;
    uint64_t _too_many_hba = 0;
    uint64_t _too_few_hbd = 0;
    uint64_t _too_many_hbd = 0;
    uint64_t _too_many_halogens = 0;
    uint64_t _too_long = 0;
    uint64_t _too_few_csp3 = 0;
    uint64_t _aromdens_too_high = 0;
    uint64_t _too_many_chiral = 0;
    uint64_t _too_many_fragments = 0;

    uint64_t _matches_exclusion_smarts = 0;
    uint64_t _no_match_required_smarts = 0;

  // Private functions
    int Process(Molecule& m);
    int Process(Molecule& m, int matoms, int nrings);
    int MaybeWriteToRejectStream(const const_IWSubstring& line);

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

    int Process(const const_IWSubstring& buffer, IWString_and_File_Descriptor& output);

    // After processing, report a summary of what has been done.
    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 0;
  _molecules_read = 0;
  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  nvrtspsa::set_display_psa_unclassified_atom_mesages(0);
  xlogp::SetIssueUnclassifiedAtomMessages(0);

  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_apply_zwitterion_correction(1);
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

  if (! cl.option_present('F')) {
    cerr << "Options::Initialise:must specify Requirements proto via the -F option\n";
    return 0;
  }

  if (cl.option_present('F')) {
    IWString fname = cl.string_value('F');
    std::optional<MoleculeFilterData::Requirements> maybe_proto =
        iwmisc::ReadTextProtoCommentsOK<MoleculeFilterData::Requirements>(fname);
    if (! maybe_proto) {
      cerr << "Options::Initialise:cannot read textproto '" << fname << "'\n";
      return 0;
    }

    _requirements = std::move(*maybe_proto);
  }

  if (_requirements.required_smarts_size() ||
      _requirements.must_not_have_smarts_size()) {
    cerr << "Options::Initialise:smarts not implemented\n";
    return 0;
  }

  if (_remove_chirality && _requirements.has_max_chiral()) {
    cerr << "Options::Initialise:removing chirality has been specified (-c)\n";
    cerr << "But the config file contains 'max_chiral'. Impossible\n";
    return 0;
  }

  if (_requirements.has_max_distance() && ! _reduce_to_largest_fragment) {
    cerr << "Options::Initialise:max distance specified, but largest fragment not selected (-l)\n";
    cerr << "Automatically enabling largest fragment selection\n";
    _reduce_to_largest_fragment = 1;
  }

  if (cl.option_present('B')) {
    IWString fname = cl.string_value('B');
    fname.EnsureEndsWith(".smi");
    if (! _reject_stream.open(fname.null_terminated_chars())) {
      cerr << "Options::Initialise:cannot open rejection file '" << fname << "'\n";
      return 0;
    } 

    if (_verbose) {
      cerr << "Rejected molecules written to " << fname << "'\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules, passed " << _molecules_passed <<
                ' ' << iwmisc::Fraction<float>(_molecules_passed, _molecules_read) << '\n';
  if (_requirements.has_exclude_non_organic()) {
    output << _non_organic << " non organic " << '\n';
  }
  if (_requirements.has_exclude_isotopes()) {
    output << _isotope << " isotope " << '\n';
  }
  if (_requirements.has_min_natoms()) {
    output << _too_few_atoms << " too few atoms " << _requirements.min_natoms() << '\n';
  }
  if (_requirements.has_max_natoms()) {
    output << _too_many_atoms << " too many atoms " << _requirements.max_natoms() << '\n';
  }
  if (_requirements.has_min_nrings()) {
    output << _too_few_rings << " too few rings " << _requirements.min_nrings() << '\n';
  }
  if (_requirements.has_max_nrings()) {
    output << _too_many_rings << " too many rings " << _requirements.max_nrings() << '\n';
  }

  if (_requirements.has_min_heteroatom_count()) {
    output << _too_few_heteroatoms << " too few heteroatoms " << _requirements.min_heteroatom_count() << '\n';
  }
  if (_requirements.has_min_heteroatom_fraction()) {
    output << _min_heteroatom_fraction << " min heteroatom fraction " << _requirements.min_heteroatom_fraction() << '\n';
  }
  if (_requirements.has_max_heteroatom_fraction()) {
    output << _max_heteroatom_fraction << " max heteroatom fraction " << _requirements.max_heteroatom_fraction() << '\n';
  }

  if (_requirements.has_min_aromatic_ring_count()) {
    output << _too_few_aromatic_rings << " too few aromatic rings " << _requirements.min_aromatic_ring_count() << '\n';
  }
  if (_requirements.has_max_aromatic_ring_count()) {
    output << _too_many_aromatic_rings << " too many aromatic rings " << _requirements.max_aromatic_ring_count() << '\n';
  }

  if (_requirements.has_min_aliphatic_ring_count()) {
    output << _too_few_aliphatic_rings << " too few aliphatic rings " << _requirements.min_aliphatic_ring_count() << '\n';
  }
  if (_requirements.has_max_aliphatic_ring_count()) {
    output << _too_many_aliphatic_rings << " too many aliphatic rings " << _requirements.max_aliphatic_ring_count() << '\n';
  }

  if (_requirements.has_max_ring_system_size()) {
    output << _ring_system_too_large << " ring systems too large " << _requirements.max_ring_system_size() << '\n';
  }
  if (_requirements.has_max_aromatic_rings_in_system()) {
    output << _too_many_aromatic_rings_in_system << " too many aromatic rings in system " << _requirements.max_aromatic_rings_in_system() << '\n';
  }

  if (_requirements.has_largest_ring_size()) {
    output << _ring_too_large << " ring too large " << _requirements.largest_ring_size() << '\n';
  }

  if (_requirements.has_min_rotatable_bonds()) {
    output << _too_few_rotbond << " too few rotatable bonds " << _requirements.min_rotatable_bonds() << '\n';
  }
  if (_requirements.has_max_rotatable_bonds()) {
    output << _too_many_rotbond << " too many rotatable bonds " << _requirements.max_rotatable_bonds() << '\n';
  }

  if (_requirements.has_min_tpsa()) {
    output << _low_tpsa << " low TPSA " << _requirements.min_tpsa() << '\n';
  }
  if (_requirements.has_max_tpsa()) {
    output << _high_tpsa << " high TPSA " << _requirements.max_tpsa() << '\n';
  }

  if (_requirements.has_min_alogp()) {
    output << _low_alogp << " low ALOGP " << _requirements.min_alogp() << '\n';
  }
  if (_requirements.has_max_alogp()) {
    output << _high_alogp << " high ALOGP " << _requirements.max_alogp() << '\n';
  }

  if (_requirements.has_min_xlogp()) {
    output << _low_xlogp << " low XLOGP " << _requirements.min_xlogp() << '\n';
  }
  if (_requirements.has_max_xlogp()) {
    output << _high_xlogp << " high XLOGP " << _requirements.max_xlogp() << '\n';
  }

  if (_requirements.has_min_hba()) {
    output << _too_few_hba << " too few HBA " << _requirements.min_hba() << '\n';
  }
  if (_requirements.has_max_hba()) {
    output << _too_many_hba << " too many HBA " << _requirements.max_hba() << '\n';
  }

  if (_requirements.has_min_hbd()) {
    output << _too_few_hbd << " too few HBD " << _requirements.min_hbd() << '\n';
  }
  if (_requirements.has_max_hbd()) {
    output << _too_many_hbd << " too many HBD " << _requirements.max_hbd() << '\n';
  }

  if (_requirements.has_max_halogen_count()) {
    output << _too_many_halogens << " too many halogens " << _requirements.max_halogen_count() << '\n';
  }

  if (_requirements.has_max_distance()) {
    output << _too_long << " molecules too long " << _requirements.max_distance() << '\n';
  }

  if (_requirements.has_min_sp3_carbon()) {
    output << _too_few_csp3 << " too few CSP3 " << _requirements.min_sp3_carbon() << '\n';
  }

  if (_requirements.has_max_aromatic_density()) {
    output << _aromdens_too_high << " aromatic density too high " << _requirements.max_aromatic_density() << '\n';
  }

  if (_requirements.has_max_chiral()) {
    output << _too_many_chiral << " too many chiral centres " << _requirements.max_chiral() << '\n';
  }

  if (_requirements.has_max_number_fragments()) {
    output << _too_many_fragments << " too many fragments " << _requirements.max_number_fragments() << '\n';
  }

  return 1;
}

int
CountHeteroatoms(const Molecule& m) {
  int rc = 0;
  m.each_atom_lambda([&rc](const Atom& a) {
    if (a.atomic_number() != 6) {
      ++rc;
    }
  });

  return rc;
}

int
AromaticRingCount(Molecule& m) {
  m.compute_aromaticity_if_needed();

  int rc = 0;
  for (const Ring* r : m.sssr_rings()) {
    if (r->is_aromatic()) {
      ++rc;
    }
  }

  return rc;
}

// Return true if we examine multiple fragments, which
// means that `largest_frag` will be different from `smiles`.
bool
LargestFragment(const const_IWSubstring& smiles,
                const_IWSubstring& largest_frag,
                int& natoms, int& nrings) {
  int max_atoms = 0;
  int i = 0;
  const_IWSubstring token;
  int fragments_examined = 0;
  while (smiles.nextword(token, i, '.')) {
    ++fragments_examined;
    int nri;
    int nat = count_atoms_in_smiles(token, nri);
    if (nat > max_atoms) {
      max_atoms = nat;
      nrings = nri;
      largest_frag = token;
    }
  }
  
  natoms = max_atoms;

  if (fragments_examined == 1) {
    return false;
  }
  return true;
}

// If chemical standardisation is in effect
int
Options::Process(const const_IWSubstring& line,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  const_IWSubstring smiles, id;
  line.split(smiles, ' ', id);

  bool smiles_changed = false;

  if (_requirements.has_max_number_fragments()) {
    const int nfrag = smiles.ccount('.') + 1;
    if (nfrag > _requirements.max_number_fragments()) {
      ++_too_many_fragments;
      MaybeWriteToRejectStream(line);
      return 0;
    }
  }

  const_IWSubstring largest_frag;
  int matoms = 0;
  int nrings = 0;
  if (_reduce_to_largest_fragment) {
    if (LargestFragment(smiles, largest_frag, matoms, nrings)) {
      smiles_changed = true;
    }
  } else {
    matoms = count_atoms_in_smiles(smiles, nrings);
    largest_frag = smiles;
  }

  if (matoms == 0) {
    return 0;
  }

  if (_requirements.has_min_natoms() && matoms < _requirements.min_natoms()) {
    ++_too_few_atoms;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  if (_requirements.has_max_natoms() && matoms > _requirements.max_natoms()) {
    ++_too_many_atoms;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  if (_requirements.has_min_nrings() && nrings < _requirements.min_nrings()) {
    ++_too_few_rings;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  if (_requirements.has_max_nrings() && nrings > _requirements.max_nrings()) {
    ++_too_many_rings;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  Molecule m;
  if (! m.build_from_smiles(largest_frag)) {
    cerr << "MoleculeFilterLine:invalid smiles '" << line << "'\n";
    return 0;
  }

  if (m.empty()) {
    cerr << "MoleculeFilterLine:no atoms '" << line << "'\n";
    return 0;
  }

  if (_requirements.has_max_chiral() &&
      m.chiral_centres() > _requirements.max_chiral()) {
    ++_too_many_chiral;
    MaybeWriteToRejectStream(line);
    return 0;
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
    smiles_changed = true;
  }

  if (_chemical_standardisation.active()) {
    if (_chemical_standardisation.process(m)) {
      smiles_changed = true;
      matoms = m.natoms();
    }
  }

  if (Process(m, matoms, nrings)) {
    ++_molecules_passed;

    // If the structure was altered, write the smiles, otherwise the original line
    if (smiles_changed) {
      output << m.smiles() << ' ' << m.name() << '\n';
    } else {
      output << line << '\n';
    }

    output.write_if_buffer_holds_more_than(4092);

    return 1;
  }

  MaybeWriteToRejectStream(line);

  return 0;
}

int
Options::MaybeWriteToRejectStream(const const_IWSubstring& line) {
  if (! _reject_stream.is_open()) {
    return 0;
  }

  _reject_stream << line << '\n';

  _reject_stream.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::Process(Molecule& m,
                 const int matoms,
                 const int nrings) {
  // cerr << "Processing " << m.smiles() << ' ' << m.name() << ' ' << matoms << '\n';
  if (m.natoms() != matoms) {
    cerr << "Atom count mismatch " << m.natoms() << " vs " << matoms << '\n';
  }

  if (_requirements.has_min_heteroatom_count() ||
      _requirements.has_min_heteroatom_fraction() ||
      _requirements.has_max_heteroatom_fraction()) {
    const int hac = CountHeteroatoms(m);

    if (_requirements.has_min_heteroatom_count() && hac < _requirements.min_heteroatom_count()) {
      ++_too_few_heteroatoms;
      return 0;
    }

    if (_requirements.has_min_heteroatom_fraction() ||
        _requirements.has_max_heteroatom_fraction()) {
      const float haf = iwmisc::Fraction<float>(hac, matoms);
      if (_requirements.has_min_heteroatom_fraction() && haf < _requirements.min_heteroatom_fraction()) {
        ++_min_heteroatom_fraction;
        return 0;
      }
      if (_requirements.has_max_heteroatom_fraction() && haf > _requirements.max_heteroatom_fraction()) {
        ++_max_heteroatom_fraction;
        return 0;
      }
    }
  }

  if (_requirements.has_exclude_non_organic() && ! m.organic_only()) {
    ++_non_organic;
    return 0;
  }

  if (_requirements.has_exclude_isotopes() && m.number_isotopic_atoms() > 0) {
    ++_isotope;
    return 0;
  }

  if (_requirements.has_max_chiral() &&
      m.chiral_centres() > _requirements.max_chiral()) {
    ++_too_many_chiral;
    return 0;
  }

  int arc = 0;
  int need_to_compute_aromatic_rings = 0;
  if (_requirements.has_min_aromatic_ring_count() ||
      _requirements.has_max_aromatic_ring_count() ||
      _requirements.has_max_aromatic_rings_in_system()) {
    need_to_compute_aromatic_rings = 1;
  }

  if (need_to_compute_aromatic_rings) {
    // Check nrings first. If not enough rings if all were aromatic...
    if (_requirements.has_min_aromatic_ring_count() && nrings < _requirements.min_aromatic_ring_count()) {
      ++_too_few_aromatic_rings;
      return 0;
    }

    arc = AromaticRingCount(m);
    if (_requirements.has_min_aromatic_ring_count() && arc < _requirements.min_aromatic_ring_count()) {
      ++_too_few_aromatic_rings;
      return 0;
    }
    if (_requirements.has_max_aromatic_ring_count() && arc > _requirements.max_aromatic_ring_count()) {
      ++_too_many_aromatic_rings;
      return 0;
    }
  }

  if (_requirements.has_min_aliphatic_ring_count() ||
      _requirements.has_max_aliphatic_ring_count()) {
    if (_requirements.has_min_aliphatic_ring_count() && nrings < _requirements.min_aliphatic_ring_count()) {
      ++_too_few_aliphatic_rings;
      return 0;
    }
    const int alring = nrings - AromaticRingCount(m);
    if (_requirements.has_min_aliphatic_ring_count() && alring < _requirements.min_aliphatic_ring_count()) {
      ++_too_few_aliphatic_rings;
      return 0;
    }
    if (_requirements.has_max_aliphatic_ring_count() && alring > _requirements.max_aliphatic_ring_count()) {
      ++_too_many_aliphatic_rings;
      return 0;
    }
  }

  if (_requirements.has_min_hba() || _requirements.has_max_hba() ||
      _requirements.has_min_hbd() || _requirements.has_max_hbd()) {
    int hba, hbd;
    molecule_filter_lib::RuleOfFive(m, hba, hbd);
    if (_requirements.has_min_hba() && hba < _requirements.min_hba()) {
      ++_too_few_hba;
      return 0;
    }
    if (_requirements.has_max_hba() && hba > _requirements.max_hba()) {
      ++_too_many_hba;
      return 0;
    }
    if (_requirements.has_min_hbd() && hbd < _requirements.min_hbd()) {
      ++_too_few_hbd;
      return 0;
    }
    if (_requirements.has_max_hbd() && hbd > _requirements.max_hbd()) {
      ++_too_many_hbd;
      return 0;
    }
  }

  if (_requirements.has_min_sp3_carbon()) {
    const int csp3 = molecule_filter_lib::Sp3Carbon(m);
    if (csp3 < _requirements.min_sp3_carbon()) {
      ++_too_few_csp3;
      return 0;
    }
  }

  if (_requirements.has_max_halogen_count()) {
    const int h = molecule_filter_lib::HalogenCount(m);
    if (h > _requirements.max_halogen_count()) {
      ++_too_many_halogens;
      return 0;
    }
  }

  if (_requirements.has_largest_ring_size()) {
    int rsze = m.ringi(nrings - 1)->number_elements();
    if (rsze > _requirements.largest_ring_size()) {
      ++_ring_too_large;
      return 0;
    }
  }

  if (_requirements.has_min_rotatable_bonds() || _requirements.has_max_rotatable_bonds()) {
    int rotb = _rotbond.Process(m);
    if (_requirements.has_min_rotatable_bonds() && rotb < _requirements.min_rotatable_bonds()) {
      ++_too_few_rotbond;
      return 0;
    }
    if (_requirements.has_max_rotatable_bonds() && rotb > _requirements.max_rotatable_bonds()) {
      ++_too_many_rotbond;
      return 0;
    }
  }

  if (_requirements.has_max_aromatic_density()) {
    const float aromdens = iwmisc::Fraction<float>(m.aromatic_atom_count(), matoms);
    if (aromdens > _requirements.max_aromatic_density()) {
      ++_aromdens_too_high;
      return 0;
    }
  }

  if (_requirements.has_max_distance() && matoms > _requirements.max_distance()) {
    const int d = m.longest_path();
    if (d > _requirements.max_distance()) {
      ++_too_long;
      return 0;
    }
  }

  if (_requirements.has_min_tpsa() || _requirements.has_max_tpsa()) {
    float tpsa = novartis_polar_surface_area(m);
    if (_requirements.has_min_tpsa() && tpsa < _requirements.min_tpsa()) {
      ++_low_tpsa;
      return 0;
    }
    if (_requirements.has_max_tpsa() && tpsa > _requirements.max_tpsa()) {
      ++_high_tpsa;
      return 0;
    }
  }

  // A temporary array that some external functions might need.
  std::unique_ptr<int[]> tmp;

  if (_requirements.has_max_ring_system_size() ||
      _requirements.has_max_aromatic_rings_in_system()) {
    if (nrings < _requirements.max_ring_system_size() &&
        nrings < _requirements.max_aromatic_rings_in_system()) {
      // no need to compute.
    }  else {
      const auto [max_ring_system_size, max_aromatic_rings_in_system] = 
          molecule_filter_lib::MaxRingSystemSize(m, tmp);
      if (_requirements.has_max_ring_system_size() && max_ring_system_size > _requirements.max_ring_system_size()) {
        ++_ring_system_too_large;
        return 0;
      }
      if (_requirements.has_max_aromatic_rings_in_system() &&
          max_aromatic_rings_in_system > _requirements.max_aromatic_rings_in_system()) {
        ++_too_many_aromatic_rings_in_system;
      }
    }
  }

  if (_requirements.has_min_xlogp() || _requirements.has_max_alogp()) {
    std::optional<double> x = _alogp.LogP(m);
    if (! x) {
    } else if (_requirements.has_min_alogp() && *x < _requirements.min_alogp()) {
      ++_low_alogp;
      return 0;
    }
    if (!x) {
    } else if (_requirements.has_max_alogp() && *x > _requirements.max_alogp()) {
      ++_high_alogp;
      return 0;
    }
  }

  if (_requirements.has_min_xlogp() || _requirements.has_max_xlogp()) {
    if (! tmp) {
      tmp.reset(new int[matoms]);
    }
    std::fill_n(tmp.get(), matoms, 0);
    std::optional<double> x = xlogp::XLogP(m, tmp.get());
    if (! x) {
    } else if (_requirements.has_min_xlogp() && *x < _requirements.min_xlogp()) {
      ++_low_xlogp;
      return 0;
    }
    if (!x) {
    } else if (_requirements.has_max_xlogp() && *x > _requirements.max_xlogp()) {
      ++_high_xlogp;
      return 0;
    }
  }


  return 1;
}

int
MoleculeFilterLine(Options& options,
                   const const_IWSubstring& line,
                   IWString_and_File_Descriptor& output) {

  options.Process(line, output);

  return 1;
}

int
MoleculeFilter(Options& options,
                iwstring_data_source& input,
                IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    MoleculeFilterLine(options, buffer, output);
  }

  return 1;
}

int
MoleculeFilter(Options& options,
             const char * fname,
             IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "MoleculeFilter:cannot open '" << fname << "'\n";
    return 0;
  }

  return MoleculeFilter(options, input, output);
}

int
MoleculeFilter(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:F:B:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

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
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  for (const char * fname : cl) {
    if (! MoleculeFilter(options, fname, output)) {
      cerr << "MoleculeFilter::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace molecule_filter

int
main(int argc, char ** argv) {

  int rc = molecule_filter::MoleculeFilter(argc, argv);

  return rc;
}
