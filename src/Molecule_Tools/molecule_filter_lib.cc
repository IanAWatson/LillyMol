#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <utility>

#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/proto_support.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/alogp.h"
#include "Molecule_Tools/xlogp.h"
#include "molecule_filter_lib.h"

namespace molecule_filter_lib {

using std::cerr;

namespace {

struct FeatureNameValue {
  std::string_view name;
  Feature feature;
};

constexpr std::array<FeatureNameValue, 27> kFeatureNameValues = {{
  {"natoms", Feature::kNatoms},
  {"nrings", Feature::kNrings},
  {"heteroatom_count", Feature::kHeteroatomCount},
  {"heteroatoms", Feature::kHeteroatomCount},
  {"heteroatom_fraction", Feature::kHeteroatomFraction},
  {"aromatic_ring_count", Feature::kAromaticRingCount},
  {"aromatic_rings", Feature::kAromaticRingCount},
  {"aliphatic_ring_count", Feature::kAliphaticRingCount},
  {"aliphatic_rings", Feature::kAliphaticRingCount},
  {"rotatable_bonds", Feature::kRotatableBonds},
  {"rotbond", Feature::kRotatableBonds},
  {"max_ring_system_size", Feature::kMaxRingSystemSize},
  {"ring_system_size", Feature::kMaxRingSystemSize},
  {"aromatic_rings_in_system", Feature::kAromaticRingsInSystem},
  {"tpsa", Feature::kTpsa},
  {"alogp", Feature::kAlogp},
  {"xlogp", Feature::kXlogp},
  {"hba", Feature::kHba},
  {"hbd", Feature::kHbd},
  {"largest_ring_size", Feature::kLargestRingSize},
  {"halogen_count", Feature::kHalogenCount},
  {"halogens", Feature::kHalogenCount},
  {"max_distance", Feature::kMaxDistance},
  {"longest_path", Feature::kMaxDistance},
  {"sp3_carbon", Feature::kSp3Carbon},
  {"aromatic_density", Feature::kAromaticDensity},
  {"chiral", Feature::kChiral},
}};

}  // namespace

std::optional<Feature>
FeatureFromName(std::string_view name) {
  if (name == "number_fragments" || name == "nfrag" || name == "fragments") {
    return Feature::kNumberFragments;
  }

  for (const FeatureNameValue& value : kFeatureNameValues) {
    if (name == value.name) {
      return value.feature;
    }
  }

  return std::nullopt;
}

std::string_view
FeatureName(Feature feature) {
  switch (feature) {
    case Feature::kNatoms:
      return "natoms";
    case Feature::kNrings:
      return "nrings";
    case Feature::kHeteroatomCount:
      return "heteroatom_count";
    case Feature::kHeteroatomFraction:
      return "heteroatom_fraction";
    case Feature::kAromaticRingCount:
      return "aromatic_ring_count";
    case Feature::kAliphaticRingCount:
      return "aliphatic_ring_count";
    case Feature::kRotatableBonds:
      return "rotatable_bonds";
    case Feature::kMaxRingSystemSize:
      return "max_ring_system_size";
    case Feature::kAromaticRingsInSystem:
      return "aromatic_rings_in_system";
    case Feature::kTpsa:
      return "tpsa";
    case Feature::kAlogp:
      return "alogp";
    case Feature::kXlogp:
      return "xlogp";
    case Feature::kHba:
      return "hba";
    case Feature::kHbd:
      return "hbd";
    case Feature::kLargestRingSize:
      return "largest_ring_size";
    case Feature::kHalogenCount:
      return "halogen_count";
    case Feature::kMaxDistance:
      return "max_distance";
    case Feature::kSp3Carbon:
      return "sp3_carbon";
    case Feature::kAromaticDensity:
      return "aromatic_density";
    case Feature::kChiral:
      return "chiral";
    case Feature::kNumberFragments:
      return "number_fragments";
  }

  return "unknown";
}


int
Utility::BuildFromProto(const molecule_filter_data::Utility& proto) {
  _name.clear();
  _weight = 1.0f;
  _x.clear();
  _y.clear();

  if (! proto.has_name() || proto.name().empty()) {
    cerr << "Utility::BuildFromProto:missing name\n";
    return 0;
  }
  _name = proto.name();
  std::optional<Feature> feature = FeatureFromName(_name);
  if (! feature) {
    cerr << "Utility::BuildFromProto:unrecognised feature name '" << _name << "'\n";
    return 0;
  }
  _feature = *feature;

  if (proto.has_weight()) {
    if (! std::isfinite(proto.weight()) || proto.weight() <= 0.0f) {
      cerr << "Utility::BuildFromProto:invalid weight " << proto.weight() << "\n";
      return 0;
    }
    _weight = proto.weight();
  }

  if (proto.point_size() < 2) {
    cerr << "Utility::BuildFromProto:" << _name << " must have at least two points\n";
    return 0;
  }

  std::vector<std::pair<double, double>> points;
  points.reserve(proto.point_size());
  for (const molecule_filter_data::Point& point : proto.point()) {
    const double x = point.x();
    const double y = point.y();
    if (! std::isfinite(x)) {
      cerr << "Utility::BuildFromProto:" << _name << " non-finite x value\n";
      return 0;
    }
    if (! std::isfinite(y) || y < 0.0 || y > 1.0) {
      cerr << "Utility::BuildFromProto:" << _name << " invalid y value " << y << "\n";
      return 0;
    }
    points.emplace_back(x, y);
  }

  std::sort(points.begin(), points.end(),
            [](const auto& p1, const auto& p2) { return p1.first < p2.first; });

  _x.reserve(points.size());
  _y.reserve(points.size());
  for (const auto& [x, y] : points) {
    if (! _x.empty() && x == _x.back()) {
      cerr << "Utility::BuildFromProto:" << _name << " duplicate x value " << x << "\n";
      _x.clear();
      _y.clear();
      return 0;
    }
    _x.push_back(x);
    _y.push_back(y);
  }

  return 1;
}

double
Utility::Value(double value) const {
  if (_x.empty()) {
    return 0.0;
  }

  if (value <= _x.front()) {
    return _y.front();
  }
  if (value >= _x.back()) {
    return _y.back();
  }

  for (int i = 1; i < static_cast<int>(_x.size()); ++i) {
    if (value > _x[i]) {
      continue;
    }

    const double fraction = (value - _x[i - 1]) / (_x[i] - _x[i - 1]);
    return _y[i - 1] + fraction * (_y[i] - _y[i - 1]);
  }

  return _y.back();
}

FeatureValues::FeatureValues(Molecule& m, const int matoms, const int nrings,
                             quick_rotbond::QuickRotatableBonds& rotbond,
                             alogp::ALogP& alogp,
                             xlogp::XLogPCalc& xlogp) :
    _m(m),
    _matoms(matoms),
    _nrings(nrings),
    _rotbond(rotbond),
    _alogp(alogp),
    _xlogp(xlogp) {
}

int
FeatureValues::HeteroatomCount() {
  if (! _heteroatom_count) {
    _heteroatom_count = molecule_filter_lib::CountHeteroatoms(_m);
  }
  return *_heteroatom_count;
}

int
FeatureValues::AromaticRingCount() {
  if (! _aromatic_ring_count) {
    _aromatic_ring_count = molecule_filter_lib::AromaticRingCount(_m);
  }
  return *_aromatic_ring_count;
}

int
FeatureValues::RotatableBonds() {
  if (! _rotatable_bonds) {
    _rotatable_bonds = _rotbond.Process(_m);
  }
  return *_rotatable_bonds;
}

void
FeatureValues::ComputeRingSystem() {
  if (_max_ring_system_size && _aromatic_rings_in_system) {
    return;
  }

  const auto [max_ring_system_size, aromatic_rings_in_system] = MaxRingSystemSize(_m, _tmp);
  _max_ring_system_size = max_ring_system_size;
  _aromatic_rings_in_system = aromatic_rings_in_system;
}

std::optional<double>
FeatureValues::ALogP() {
  if (! _alogp_value) {
    std::optional<double> value = _alogp.LogP(_m);
    if (! value) {
      return std::nullopt;
    }
    _alogp_value = *value;
  }

  return _alogp_value;
}

std::optional<double>
FeatureValues::XLogP() {
  if (! _xlogp_value) {
    if (! _tmp) {
      _tmp = std::make_unique<int[]>(_matoms);
    }
    std::fill_n(_tmp.get(), _matoms, 0);
    std::optional<double> value = _xlogp.LogP(_m, _tmp.get());
    if (! value) {
      return std::nullopt;
    }
    _xlogp_value = *value;
  }

  return _xlogp_value;
}

void
FeatureValues::ComputeHbaHbd() {
  if (_hba && _hbd) {
    return;
  }

  int hba;
  int hbd;
  RuleOfFive(_m, hba, hbd);
  _hba = hba;
  _hbd = hbd;
}

int
FeatureValues::HalogenCount() {
  if (! _halogen_count) {
    _halogen_count = molecule_filter_lib::HalogenCount(_m);
  }
  return *_halogen_count;
}

int
FeatureValues::MaxDistance() {
  if (! _max_distance) {
    _max_distance = _m.longest_path();
  }
  return *_max_distance;
}

int
FeatureValues::Sp3Carbon() {
  if (! _sp3_carbon) {
    _sp3_carbon = molecule_filter_lib::Sp3Carbon(_m);
  }
  return *_sp3_carbon;
}

int
FeatureValues::NumberFragments() {
  if (! _number_fragments) {
    _number_fragments = _m.number_fragments();
  }
  return *_number_fragments;
}

std::optional<double>
FeatureValues::Value(Feature feature) {
  switch (feature) {
    case Feature::kNatoms:
      return _matoms;
    case Feature::kNrings:
      return _nrings;
    case Feature::kHeteroatomCount:
      return HeteroatomCount();
    case Feature::kHeteroatomFraction:
      return iwmisc::Fraction<double>(HeteroatomCount(), _matoms);
    case Feature::kAromaticRingCount:
      return AromaticRingCount();
    case Feature::kAliphaticRingCount:
      return _nrings - AromaticRingCount();
    case Feature::kRotatableBonds:
      return RotatableBonds();
    case Feature::kMaxRingSystemSize:
      ComputeRingSystem();
      return *_max_ring_system_size;
    case Feature::kAromaticRingsInSystem:
      ComputeRingSystem();
      return *_aromatic_rings_in_system;
    case Feature::kTpsa:
      if (! _tpsa) {
        _tpsa = novartis_polar_surface_area(_m);
      }
      return *_tpsa;
    case Feature::kAlogp:
      return ALogP();
    case Feature::kXlogp:
      return XLogP();
    case Feature::kHba:
      ComputeHbaHbd();
      return *_hba;
    case Feature::kHbd:
      ComputeHbaHbd();
      return *_hbd;
    case Feature::kLargestRingSize:
      if (_nrings == 0) {
        return 0.0;
      }
      return _m.ringi(_nrings - 1)->number_elements();
    case Feature::kHalogenCount:
      return HalogenCount();
    case Feature::kMaxDistance:
      return MaxDistance();
    case Feature::kSp3Carbon:
      return Sp3Carbon();
    case Feature::kAromaticDensity:
      return iwmisc::Fraction<double>(_m.aromatic_atom_count(), _matoms);
    case Feature::kChiral:
      return _m.chiral_centres();
    case Feature::kNumberFragments:
      return NumberFragments();
  }

  return std::nullopt;
}

MoleculeFilter::MoleculeFilter() {
  _active = false;
}

int
MoleculeFilter::Build(IWString& fname) {
  std::optional<molecule_filter_data::Requirements> maybe_proto =
                iwmisc::ReadTextProto<molecule_filter_data::Requirements>(fname);
  if (! maybe_proto) {
    cerr << "MoleculeFilter::Build:cannot read '" << fname << "'\n";
    _active = false;
    _utilities.clear();
    return 0;
  }

  return Build(*maybe_proto);
}

int
MoleculeFilter::Build(const molecule_filter_data::Requirements& proto) {
  _active = false;
  _utilities.clear();

  _requirements = proto;

  InitialiseOptionalFeatures();

  if (! BuildUtilities()) {
    return 0;
  }

  _active = true;

  return 1;
}

int
MoleculeFilter::EvaluateUtilities(Molecule& m, const int matoms, const int nrings,
                                  std::vector<double>& per_feature_utility,
                                  double& overall_utility) {
  per_feature_utility.clear();
  overall_utility = 0.0;

  if (_utilities.empty()) {
    return 1;
  }

  FeatureValues feature_values(m, matoms, nrings, _rotbond, _alogp, _xlogp);
  per_feature_utility.reserve(_utilities.size());

  double weighted_sum = 0.0;
  double weight_sum = 0.0;
  double product = 1.0;
  double minval = 1.0;
  double maxval = 0.0;

  for (const Utility& utility : _utilities) {
    std::optional<double> raw_value = feature_values.Value(utility.feature());
    if (! raw_value) {
      cerr << "MoleculeFilter::EvaluateUtilities:cannot compute feature '" <<
              utility.name() << "'\n";
      per_feature_utility.clear();
      overall_utility = 0.0;
      return 0;
    }

    const double u = utility.Value(*raw_value);
    per_feature_utility.push_back(u);

    weighted_sum += utility.weight() * u;
    weight_sum += utility.weight();
    product *= u;
    if (u < minval) {
      minval = u;
    }
    if (u > maxval) {
      maxval = u;
    }
  }

  const molecule_filter_data::UtilityCombination combination =
      _requirements.has_utility_combination()
          ? _requirements.utility_combination()
          : molecule_filter_data::UTILITY_COMBINATION_WEIGHTED_AVERAGE;

  switch (combination) {
    case molecule_filter_data::UTILITY_COMBINATION_UNSPECIFIED:
    case molecule_filter_data::UTILITY_COMBINATION_WEIGHTED_AVERAGE:
      overall_utility = weighted_sum / weight_sum;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_WEIGHTED_SUM:
      overall_utility = weighted_sum;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_PRODUCT:
      overall_utility = product;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_MIN:
      overall_utility = minval;
      return 1;
    case molecule_filter_data::UTILITY_COMBINATION_MAX:
      overall_utility = maxval;
      return 1;
    default:
      cerr << "MoleculeFilter::EvaluateUtilities:unrecognised combination " <<
              static_cast<int>(combination) << "\n";
      return 0;
  }
}

int
MoleculeFilter::BuildUtilities() {
  const int nutility = _requirements.utility_size();
  if (nutility == 0) {
    return 1;
  }

  _utilities.reserve(nutility);
  for (const molecule_filter_data::Utility& proto : _requirements.utility()) {
    Utility utility;
    if (! utility.BuildFromProto(proto)) {
      _utilities.clear();
      return 0;
    }
    _utilities.push_back(std::move(utility));
  }

  return 1;
}

void
MoleculeFilter::InitialiseOptionalFeatures() {
  _rotbond.set_calculation_type(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  nvrtspsa::set_display_psa_unclassified_atom_mesages(0);
  _xlogp.SetIssueUnclassifiedAtomMessages(false);

  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_apply_zwitterion_correction(1);
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
    int nat = lillymol::count_atoms_in_smiles(token, nri);
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

std::tuple<int, int>
MaxRingSystemSize(Molecule& m, std::unique_ptr<int[]>& tmp) {
  const int matoms = m.natoms();

  m.compute_aromaticity_if_needed();

  if (! tmp) {
    tmp.reset(new int[matoms]);
  }
  std::fill_n(tmp.get(), matoms, 0);

  const int nrings = m.nrings();

  std::unique_ptr<int[]> ring_already_done = std::make_unique<int[]>(nrings);
  std::fill_n(ring_already_done.get(), nrings, 0);

  int max_system_size = 0;
  int max_aromatic_rings_in_system = 0;
  for (int i = 0; i < nrings; ++i) {
    if (ring_already_done[i]) {
      continue;
    }
    const Ring* ri = m.ringi(i);
    if (! ri->is_fused()) {
      continue;
    }

    int system_size = 1;
    int aromatic_rings_in_system;
    if (ri->is_aromatic()) {
      aromatic_rings_in_system = 1;
    } else {
      aromatic_rings_in_system = 0;
    }


    for (int j = i + 1; j < nrings; ++j) {
      if (ring_already_done[j]) {
        continue;
      }

      ring_already_done[j] = 1;
      const Ring* rj = m.ringi(j);
      if (ri->fused_system_identifier() == rj->fused_system_identifier()) {
        ++system_size;
        if (rj->is_aromatic()) {
          ++aromatic_rings_in_system;
        }
      }
    }
    if (system_size > max_system_size) {
      max_system_size = system_size;
    }
    if (aromatic_rings_in_system > max_aromatic_rings_in_system) {
      max_aromatic_rings_in_system = aromatic_rings_in_system;
    }
  }

  return std::make_tuple(max_system_size, max_aromatic_rings_in_system);
}

// Lifted from iwdescr.cc
void
RuleOfFive(Molecule & m, int& acceptor, int& donor) {
  acceptor = 0;
  donor = 0;

  const int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    atomic_number_t z = m.atomic_number(i);
    // Intercept the most common case.
    if (z == 6) {
      continue;
    }

    if (z == 7 || z == 8) {
    } else {
      continue;
    }

    ++acceptor;

    const int h = m.hcount(i);

    // acceptor
    if (0 == h) {
      continue;
    }

    if (7 == z && h > 1) {
      donor += 2;
    } else {
      donor += 1;
    }
  }
}

int
HalogenCount(const Molecule& m) {
  static std::vector<int> halogen = {
    0,  // 0
    0,  // 1
    0,  // 2
    0,  // 3
    0,  // 4
    0,  // 5
    0,  // 6
    0,  // 7
    0,  // 8
    1,  // 9
    0,  // 10
    0,  // 11
    0,  // 12
    0,  // 13
    0,  // 14
    0,  // 15
    0,  // 16
    1,  // 17
    0,  // 18
    0,  // 19
    0,  // 20
    0,  // 21
    0,  // 22
    0,  // 23
    0,  // 24
    0,  // 25
    0,  // 26
    0,  // 27
    0,  // 28
    0,  // 29
    0,  // 30
    0,  // 31
    0,  // 32
    0,  // 33
    0,  // 34
    1,  // 35
    0,  // 36
    0,  // 37
    0,  // 38
    0,  // 39
    0,  // 40
    0,  // 41
    0,  // 42
    0,  // 43
    0,  // 44
    0,  // 45
    0,  // 46
    0,  // 47
    0,  // 48
    0,  // 49
    0,  // 50
    0,  // 51
    0,  // 52
    1   // 53
  };

  int rc = 0;

  for (const Atom* a : m) {
    const uint32_t z = a->atomic_number();
    if (z < halogen.size()) {
      rc += halogen[z];
    }
  }

  return rc;
}

int 
Sp3Carbon(Molecule & m) {
  int rc = 0;

  for (const Atom* a : m) {
    if (a->atomic_number() != 6) {
      continue;
    }
    if (a->fully_saturated()) {
      ++rc;
    }
  }

  return rc;
}

int
Reject(RejectionReason& rejection_reason, RejectionReason reason) {
  rejection_reason = reason;
  return 0;
}

int
MoleculeFilter::Ok(Molecule& m) {
  RejectionReason rejection_reason;
  return Ok(m, rejection_reason);
}

int
MoleculeFilter::Ok(Molecule& m, RejectionReason& rejection_reason) {
  return Ok(m, m.natoms(), m.nrings(), rejection_reason);
}

int
MoleculeFilter::Ok(Molecule& m, const int matoms, const int nrings) {
  RejectionReason rejection_reason;
  return Ok(m, matoms, nrings, rejection_reason);
}

int
MoleculeFilter::Ok(Molecule& m, const int matoms, const int nrings,
                   RejectionReason& rejection_reason) {
  rejection_reason = RejectionReason::kPass;

  if (matoms == 0) {
    return Reject(rejection_reason, RejectionReason::kZeroAtoms);
  }

  if (_requirements.has_min_natoms() && matoms < _requirements.min_natoms()) {
    return Reject(rejection_reason, RejectionReason::kTooFewAtoms);
  }

  if (_requirements.has_max_natoms() && matoms > _requirements.max_natoms()) {
    return Reject(rejection_reason, RejectionReason::kTooManyAtoms);
  }

  if (_requirements.has_min_nrings() && nrings < _requirements.min_nrings()) {
    return Reject(rejection_reason, RejectionReason::kTooFewRings);
  }

  if (_requirements.has_max_nrings() && nrings > _requirements.max_nrings()) {
    return Reject(rejection_reason, RejectionReason::kTooManyRings);
  }

  if (_requirements.has_min_heteroatom_count() ||
      _requirements.has_min_heteroatom_fraction() ||
      _requirements.has_max_heteroatom_fraction()) {
    const int hac = CountHeteroatoms(m);

    if (_requirements.has_min_heteroatom_count() && hac < _requirements.min_heteroatom_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHeteroatoms);
    }

    if (_requirements.has_min_heteroatom_fraction() ||
        _requirements.has_max_heteroatom_fraction()) {
      const float haf = iwmisc::Fraction<float>(hac, matoms);
      if (_requirements.has_min_heteroatom_fraction() && haf < _requirements.min_heteroatom_fraction()) {
        return Reject(rejection_reason, RejectionReason::kMinHeteroatomFraction);
      }
      if (_requirements.has_max_heteroatom_fraction() && haf > _requirements.max_heteroatom_fraction()) {
        return Reject(rejection_reason, RejectionReason::kMaxHeteroatomFraction);
      }
    }
  }

  if (_requirements.has_exclude_non_organic() && ! m.organic_only()) {
    return Reject(rejection_reason, RejectionReason::kNonOrganic);
  }

  if (_requirements.has_exclude_isotopes() && m.number_isotopic_atoms() > 0) {
    return Reject(rejection_reason, RejectionReason::kIsotope);
  }

  if (_requirements.has_max_chiral() &&
      m.chiral_centres() > _requirements.max_chiral()) {
    return Reject(rejection_reason, RejectionReason::kTooManyChiral);
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
      return Reject(rejection_reason, RejectionReason::kTooFewAromaticRings);
    }

    arc = AromaticRingCount(m);
    if (_requirements.has_min_aromatic_ring_count() && arc < _requirements.min_aromatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewAromaticRings);
    }
    if (_requirements.has_max_aromatic_ring_count() && arc > _requirements.max_aromatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooManyAromaticRings);
    }
  }

  if (_requirements.has_min_aliphatic_ring_count() ||
      _requirements.has_max_aliphatic_ring_count()) {
    if (_requirements.has_min_aliphatic_ring_count() && nrings < _requirements.min_aliphatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewAliphaticRings);
    }
    const int alring = nrings - AromaticRingCount(m);
    if (_requirements.has_min_aliphatic_ring_count() && alring < _requirements.min_aliphatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooFewAliphaticRings);
    }
    if (_requirements.has_max_aliphatic_ring_count() && alring > _requirements.max_aliphatic_ring_count()) {
      return Reject(rejection_reason, RejectionReason::kTooManyAliphaticRings);
    }
  }

  if (_requirements.has_min_hba() || _requirements.has_max_hba() ||
      _requirements.has_min_hbd() || _requirements.has_max_hbd()) {
    int hba, hbd;
    RuleOfFive(m, hba, hbd);
    if (_requirements.has_min_hba() && hba < _requirements.min_hba()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHba);
    }
    if (_requirements.has_max_hba() && hba > _requirements.max_hba()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHba);
    }
    if (_requirements.has_min_hbd() && hbd < _requirements.min_hbd()) {
      return Reject(rejection_reason, RejectionReason::kTooFewHbd);
    }
    if (_requirements.has_max_hbd() && hbd > _requirements.max_hbd()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHbd);
    }
  }

  if (_requirements.has_min_sp3_carbon()) {
    const int csp3 = Sp3Carbon(m);
    if (csp3 < _requirements.min_sp3_carbon()) {
      return Reject(rejection_reason, RejectionReason::kTooFewSp3Carbon);
    }
  }

  if (_requirements.has_max_halogen_count()) {
    const int h = HalogenCount(m);
    if (h > _requirements.max_halogen_count()) {
      return Reject(rejection_reason, RejectionReason::kTooManyHalogen);
    }
  }

  if (_requirements.has_largest_ring_size() && nrings > 0) {
    int rsze = m.ringi(nrings - 1)->number_elements();
    if (rsze > _requirements.largest_ring_size()) {
      return Reject(rejection_reason, RejectionReason::kRingTooLarge);
    }
  }

  if (_requirements.has_min_rotatable_bonds() || _requirements.has_max_rotatable_bonds()) {
    int rotb = _rotbond.Process(m);
    if (_requirements.has_min_rotatable_bonds() && rotb < _requirements.min_rotatable_bonds()) {
      return Reject(rejection_reason, RejectionReason::kTooFewRotatableBonds);
    }
    if (_requirements.has_max_rotatable_bonds() && rotb > _requirements.max_rotatable_bonds()) {
      return Reject(rejection_reason, RejectionReason::kTooManyRotatableBonds);
    }
  }

  if (_requirements.has_max_aromatic_density()) {
    const float aromdens = iwmisc::Fraction<float>(m.aromatic_atom_count(), matoms);
    if (aromdens > _requirements.max_aromatic_density()) {
      return Reject(rejection_reason, RejectionReason::kAromaticDensityTooHigh);
    }
  }

  if (_requirements.has_max_distance() && matoms > _requirements.max_distance()) {
    const int d = m.longest_path();
    if (d > _requirements.max_distance()) {
      return Reject(rejection_reason, RejectionReason::kTooLong);
    }
  }

  if (_requirements.has_min_tpsa() || _requirements.has_max_tpsa()) {
    float tpsa = novartis_polar_surface_area(m);
    if (_requirements.has_min_tpsa() && tpsa < _requirements.min_tpsa()) {
      return Reject(rejection_reason, RejectionReason::kLowTpsa);
    }
    if (_requirements.has_max_tpsa() && tpsa > _requirements.max_tpsa()) {
      return Reject(rejection_reason, RejectionReason::kHighTpsa);
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
      const auto [max_ring_system_size, max_aromatic_rings_in_system] = MaxRingSystemSize(m, tmp);
      if (_requirements.has_max_ring_system_size() && max_ring_system_size > _requirements.max_ring_system_size()) {
        return Reject(rejection_reason, RejectionReason::kRingSystemTooLarge);
      }
      if (_requirements.has_max_aromatic_rings_in_system() &&
          max_aromatic_rings_in_system > _requirements.max_aromatic_rings_in_system()) {
        return Reject(rejection_reason, RejectionReason::kTooManyAromaticRingsInSystem);
      }
    }
  }

  if (_requirements.has_min_alogp() || _requirements.has_max_alogp()) {
    std::optional<double> x = _alogp.LogP(m);
    if (! x) {
    } else if (_requirements.has_min_alogp() && *x < _requirements.min_alogp()) {
      return Reject(rejection_reason, RejectionReason::kLowAlogp);
    }
    if (!x) {
    } else if (_requirements.has_max_alogp() && *x > _requirements.max_alogp()) {
      return Reject(rejection_reason, RejectionReason::kHighAlogp);
    }
  }

  if (_requirements.has_min_xlogp() || _requirements.has_max_xlogp()) {
    if (! tmp) {
      tmp.reset(new int[matoms]);
    }
    std::fill_n(tmp.get(), matoms, 0);
    std::optional<double> x = _xlogp.LogP(m, tmp.get());
    if (! x) {
    } else if (_requirements.has_min_xlogp() && *x < _requirements.min_xlogp()) {
      return Reject(rejection_reason, RejectionReason::kLowXlogp);
    }
    if (!x) {
    } else if (_requirements.has_max_xlogp() && *x > _requirements.max_xlogp()) {
      return Reject(rejection_reason, RejectionReason::kHighXlogp);
    }
  }

  return 1;
}


}  // namespace molecule_filter_lib
