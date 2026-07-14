#include "Utilities/GFP_Tools/gfp_context.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Molecule_Lib/iwmfingerprint.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/maccskeys_fn5.h"
#include "Molecule_Tools/mpr.h"
#include "Utilities/GFP_Tools/dyfp.h"

namespace gfp_context {
namespace {

constexpr char kSmilesTag[] = "$SMI<";
constexpr char kIdentifierTag[] = "PCN<";
constexpr char kMolecularPropertiesTag[] = "MPR<";

uint64_t
Fnv1a(uint64_t hash, const void* v, int n) {
  const unsigned char* data = reinterpret_cast<const unsigned char*>(v);
  for (int i = 0; i < n; ++i) {
    hash ^= data[i];
    hash *= 1099511628211ULL;
  }

  return hash;
}

uint64_t
Fnv1a(uint64_t hash, const IWString& s) {
  return Fnv1a(hash, s.rawchars(), s.length());
}

bool
TagFromDataitem(const const_IWSubstring& item, IWString& tag) {
  const int open_angle = item.index('<');
  if (open_angle < 0) {
    return false;
  }

  tag.set(item.rawchars(), open_angle + 1);
  return true;
}

int
CountTdts(iwstring_data_source& input, const char* fname) {
  const uint64_t result = input.count_records_starting_with("|");
  if (result > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
    std::cerr << "GFPList::ReadFile:too many TDT records in '" << fname << "'\n";
    return -1;
  }

  return static_cast<int>(result);
}

int
CompareTagsCanonical(const IWString& lhs, const IWString& rhs) {
  const int n1 = lhs.ends_with('<') ? lhs.length() - 1 : lhs.length();
  const int n2 = rhs.ends_with('<') ? rhs.length() - 1 : rhs.length();
  const int n = std::min(n1, n2);

  const int cmp = strncmp(lhs.rawchars(), rhs.rawchars(), n);
  if (cmp != 0) {
    return cmp;
  }
  if (n1 < n2) {
    return -1;
  }
  if (n1 > n2) {
    return 1;
  }
  return 0;
}

bool
NearestNeighbourBetter(const NearestNeighbour& lhs, const NearestNeighbour& rhs) {
  if (lhs.distance != rhs.distance) {
    return lhs.distance < rhs.distance;
  }
  return lhs.index < rhs.index;
}

struct NearestNeighbourMaxHeapCompare {
  bool
  operator()(const NearestNeighbour& lhs, const NearestNeighbour& rhs) const {
    // priority_queue places the highest priority item at the top. Return true
    // when `lhs` is better than `rhs`, so the worst retained neighbour is top.
    return NearestNeighbourBetter(lhs, rhs);
  }
};

using NearestNeighbourMaxHeap =
    std::priority_queue<NearestNeighbour, std::vector<NearestNeighbour>,
                        NearestNeighbourMaxHeapCompare>;

void
MaybeAddNearestNeighbour(NearestNeighbourMaxHeap& heap, int k,
                         NearestNeighbour candidate) {
  if (static_cast<int>(heap.size()) < k) {
    heap.push(candidate);
    return;
  }

  if (NearestNeighbourBetter(candidate, heap.top())) {
    heap.pop();
    heap.push(candidate);
  }
}

std::vector<NearestNeighbour>
SortedNearestNeighbours(NearestNeighbourMaxHeap& heap) {
  std::vector<NearestNeighbour> result;
  result.reserve(heap.size());
  while (!heap.empty()) {
    result.push_back(heap.top());
    heap.pop();
  }
  std::sort(result.begin(), result.end(), NearestNeighbourBetter);

  return result;
}

float
FixedBinaryTanimoto(const FixedBitVector& lhs, const FixedBitVector& rhs) {
  const int bic = lhs.bits().BitsInCommon(rhs.bits());
  const int denom = lhs.nset() + rhs.nset() - bic;
  if (denom == 0) {
    return 1.0f;
  }

  return static_cast<float>(bic) / static_cast<float>(denom);
}

}  // namespace

class StandardFingerprintGenerator {
 private:
  Molecular_Properties_Generator _mpr;
  MACCSKeys _mk;
  IWMFingerprintOptions _iw_options;
  Chemical_Standardisation _chemical_standardisation;
  int _preprocess = 1;
  int _tmp[2048];

  void Preprocess(Molecule& m);

 public:
  explicit StandardFingerprintGenerator(bool preprocess);

  int Generate(Molecule& m, GFPFingerprint& destination);
};

StandardFingerprintGenerator::StandardFingerprintGenerator(bool preprocess)
    : _preprocess(preprocess) {
  _chemical_standardisation.activate_all();
}

void
StandardFingerprintGenerator::Preprocess(Molecule& m) {
  m.reduce_to_largest_fragment_carefully();
  m.remove_all_chiral_centres();
  m.revert_all_directional_bonds_to_non_directional();
  _chemical_standardisation.process(m);
}

int
StandardFingerprintGenerator::Generate(Molecule& input, GFPFingerprint& destination) {
  std::unique_ptr<Molecule> copy;
  Molecule* m = &input;
  if (_preprocess) {
    copy = std::make_unique<Molecule>(input);
    Preprocess(*copy);
    m = copy.get();
  }

  int properties[8];
  _mpr(*m, properties);
  if (!destination.mutable_molecular_properties().BuildFromArray(properties, 8)) {
    return 0;
  }

  IWMFingerprint iwfp(_iw_options);
  iwfp.construct_fingerprint(*m);
  if (!destination.mutable_fixed_binary(0).BuildFromArray(
          iwfp.vector(), _iw_options.bits_per_fingerprint)) {
    return 0;
  }

  std::fill_n(_tmp, 2048, 0);
  _mk(*m, _tmp);
  if (!destination.mutable_fixed_binary(1).BuildFromArray(_tmp, _mk.nbits())) {
    return 0;
  }

  _mk.set_level_2_fingerprint(_tmp);
  if (!destination.mutable_fixed_binary(2).BuildFromArray(_tmp, _mk.nbits())) {
    return 0;
  }

  return 1;
}

MolecularProperties::~MolecularProperties() {
  delete[] _property;
}

MolecularProperties::MolecularProperties(const MolecularProperties& rhs) {
  *this = rhs;
}

MolecularProperties&
MolecularProperties::operator=(const MolecularProperties& rhs) {
  if (this == &rhs) {
    return *this;
  }

  delete[] _property;
  _nproperties = rhs._nproperties;
  if (_nproperties == 0) {
    _property = nullptr;
    return *this;
  }

  _property = new int[_nproperties];
  std::copy(rhs._property, rhs._property + _nproperties, _property);

  return *this;
}

MolecularProperties::MolecularProperties(MolecularProperties&& rhs) noexcept {
  *this = std::move(rhs);
}

MolecularProperties&
MolecularProperties::operator=(MolecularProperties&& rhs) noexcept {
  if (this == &rhs) {
    return *this;
  }

  delete[] _property;

  _nproperties = rhs._nproperties;
  _property = rhs._property;

  rhs._nproperties = 0;
  rhs._property = nullptr;

  return *this;
}

int
MolecularProperties::Build(const const_IWSubstring& buffer) {
  IWDYFP fp;
  if (!fp.construct_from_tdt_record(buffer)) {
    std::cerr << "MolecularProperties::Build:cannot parse '" << buffer << "'\n";
    return 0;
  }

  const int nproperties = fp.nbits() / 8;
  if (nproperties <= 0) {
    std::cerr << "MolecularProperties::Build:no properties in '" << buffer << "'\n";
    return 0;
  }

  if (_nproperties != nproperties) {
    delete[] _property;
    _nproperties = nproperties;
    _property = new int[_nproperties];
  }

  const unsigned char* bits = fp.bits();
  for (int i = 0; i < _nproperties; ++i) {
    _property[i] = bits[i];
  }

  return 1;
}

int
MolecularProperties::BuildFromArray(const int* properties, int nproperties) {
  if (nproperties <= 0) {
    return 0;
  }

  if (_nproperties != nproperties) {
    delete[] _property;
    _nproperties = nproperties;
    _property = new int[_nproperties];
  }

  std::copy(properties, properties + _nproperties, _property);
  return 1;
}

float
MolecularProperties::Similarity(const MolecularProperties& rhs) const {
  if (_nproperties != rhs._nproperties || _nproperties == 0) {
    return 0.0f;
  }

  float result = 0.0f;
  for (int i = 0; i < _nproperties; ++i) {
    const int lhs = _property[i];
    const int rhs_value = rhs._property[i];
    if (lhs == rhs_value) {
      result += 1.0f;
    } else if (lhs == 0 || rhs_value == 0) {
      result += 0.5f;
    } else if (lhs < rhs_value) {
      result += static_cast<float>(lhs) / static_cast<float>(rhs_value);
    } else {
      result += static_cast<float>(rhs_value) / static_cast<float>(lhs);
    }
  }

  return result / static_cast<float>(_nproperties);
}

int
FixedBitVector::Build(const const_IWSubstring& buffer) {
  if (!_bits.ConstructFromTdtRecordNset(buffer, _nset)) {
    std::cerr << "FixedBitVector::Build:cannot parse '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
FixedBitVector::BuildFromArray(const int* bits, int nbits) {
  if (nbits <= 0) {
    return 0;
  }

  if (!_bits.ConstructFromArrayIWBitsOrder(bits, nbits)) {
    return 0;
  }
  _nset = _bits.nset();

  return 1;
}

float
FixedBitVector::Tanimoto(const FixedBitVector& rhs) const {
  return FixedBinaryTanimoto(*this, rhs);
}

void
GFPFingerprint::FreeArrays() {
  delete[] _fixed_binary;
  delete[] _sparse;
  delete[] _fixed_counted;

  _fixed_binary = nullptr;
  _sparse = nullptr;
  _fixed_counted = nullptr;
}

GFPFingerprint::~GFPFingerprint() {
  FreeArrays();
}

GFPFingerprint::GFPFingerprint(GFPFingerprint&& rhs) noexcept {
  *this = std::move(rhs);
}

GFPFingerprint&
GFPFingerprint::operator=(GFPFingerprint&& rhs) noexcept {
  if (this == &rhs) {
    return *this;
  }

  FreeArrays();

  _context_hash = rhs._context_hash;
  _molecular_properties = std::move(rhs._molecular_properties);
  _fixed_binary = rhs._fixed_binary;
  _sparse = rhs._sparse;
  _fixed_counted = rhs._fixed_counted;

  rhs._context_hash = 0;
  rhs._fixed_binary = nullptr;
  rhs._sparse = nullptr;
  rhs._fixed_counted = nullptr;

  return *this;
}

int
GFPFingerprint::Allocate(const GFPContext& context) {
  FreeArrays();

  if (context.nfixed_binary()) {
    _fixed_binary = new FixedBitVector[context.nfixed_binary()];
  }
  if (context.nsparse()) {
    _sparse = new Sparse_Fingerprint[context.nsparse()];
  }
  if (context.nfixed_counted()) {
    _fixed_counted = new Fixed_Size_Counted_Fingerprint_uchar[context.nfixed_counted()];
  }

  _context_hash = context.context_hash();
  return 1;
}

int
GFPFingerprint::Build(const IW_TDT& tdt, const GFPContext& context) {
  if (_context_hash != context.context_hash()) {
    if (!Allocate(context)) {
      return 0;
    }
  }

  for (const Component& component : context.components()) {
    const_IWSubstring dataitem;
    if (!tdt.dataitem(component.tag, dataitem)) {
      std::cerr << "GFPFingerprint::Build:missing '" << component.tag << "'\n";
      return 0;
    }
    if (dataitem.ends_with('\n')) {
      dataitem.chop();
    }

    switch (component.kind) {
      case ComponentKind::kMolecularProperties:
        if (!_molecular_properties.Build(dataitem)) {
          return 0;
        }
        break;
      case ComponentKind::kFixedBinary:
        if (!_fixed_binary[component.index].Build(dataitem)) {
          return 0;
        }
        break;
      case ComponentKind::kSparse:
        if (!_sparse[component.index].construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPFingerprint::Build:cannot parse sparse '" << component.tag
                    << "'\n";
          return 0;
        }
        break;
      case ComponentKind::kFixedCounted:
        if (!_fixed_counted[component.index].construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPFingerprint::Build:cannot parse fixed counted '"
                    << component.tag << "'\n";
          return 0;
        }
        break;
    }
  }

  return 1;
}

void
GFPContext::ComputeHash() {
  uint64_t hash = 14695981039346656037ULL;
  for (const Component& component : _components) {
    const uint8_t kind = static_cast<uint8_t>(component.kind);
    hash = Fnv1a(hash, &kind, sizeof(kind));
    hash = Fnv1a(hash, &component.index, sizeof(component.index));
    hash = Fnv1a(hash, component.tag);
  }

  _context_hash = hash;
}

void
GFPContext::CanonicalizeComponents() {
  std::sort(_components.begin(), _components.end(),
            [](const Component& lhs, const Component& rhs) {
              return CompareTagsCanonical(lhs.tag, rhs.tag) < 0;
            });

  _nfixed_binary = 0;
  _nsparse = 0;
  _nfixed_counted = 0;
  _has_molecular_properties = false;

  for (Component& component : _components) {
    switch (component.kind) {
      case ComponentKind::kMolecularProperties:
        component.index = 0;
        _has_molecular_properties = true;
        break;
      case ComponentKind::kFixedBinary:
        component.index = _nfixed_binary;
        ++_nfixed_binary;
        break;
      case ComponentKind::kSparse:
        component.index = _nsparse;
        ++_nsparse;
        break;
      case ComponentKind::kFixedCounted:
        component.index = _nfixed_counted;
        ++_nfixed_counted;
        break;
    }
  }
}

void
GFPContext::BuildDefaultActiveList() {
  _active.clear();
  if (_components.empty()) {
    return;
  }

  const float weight = 1.0f / static_cast<float>(_components.size());
  for (Component& component : _components) {
    component.weight = weight;
    _active.push_back(ActiveComponent{component.kind, component.index, weight});
  }
}

GFPContext::GFPContext() = default;

GFPContext::~GFPContext() = default;

int
GFPContext::BuildFromTdt(const IW_TDT& tdt) {
  _components.clear();
  _active.clear();
  _nfixed_binary = 0;
  _nsparse = 0;
  _nfixed_counted = 0;
  _has_molecular_properties = false;
  _context_hash = 0;
  _standard_generator.reset();

  int i = 0;
  const_IWSubstring item;
  while (tdt.next_dataitem(item, i)) {
    IWString tag;
    if (!TagFromDataitem(item, tag)) {
      continue;
    }

    if (tag == kMolecularPropertiesTag) {
      _has_molecular_properties = true;
      _components.push_back(Component{ComponentKind::kMolecularProperties, tag, 0, 1.0f});
    } else if (tag.starts_with("FP")) {
      _components.push_back(
          Component{ComponentKind::kFixedBinary, tag, _nfixed_binary, 1.0f});
      ++_nfixed_binary;
    } else if (tag.starts_with("NC")) {
      _components.push_back(Component{ComponentKind::kSparse, tag, _nsparse, 1.0f});
      ++_nsparse;
    } else if (tag.starts_with("FC")) {
      _components.push_back(
          Component{ComponentKind::kFixedCounted, tag, _nfixed_counted, 1.0f});
      ++_nfixed_counted;
    }
  }

  if (_components.empty()) {
    std::cerr << "GFPContext::BuildFromTdt:no fingerprint components\n";
    return 0;
  }

  CanonicalizeComponents();
  ComputeHash();
  BuildDefaultActiveList();

  return 1;
}

int
GFPContext::BuildStandard(bool preprocess) {
  _components.clear();
  _active.clear();
  _nfixed_binary = 0;
  _nsparse = 0;
  _nfixed_counted = 0;
  _has_molecular_properties = true;
  _context_hash = 0;

  _components.push_back(
      Component{ComponentKind::kFixedBinary, IWString("FPIW<"), 0, 1.0f});
  _components.push_back(
      Component{ComponentKind::kFixedBinary, IWString("FPMK<"), 1, 1.0f});
  _components.push_back(
      Component{ComponentKind::kFixedBinary, IWString("FPMK2<"), 2, 1.0f});
  _components.push_back(
      Component{ComponentKind::kMolecularProperties, IWString("MPR<"), 0, 1.0f});
  _nfixed_binary = 3;

  CanonicalizeComponents();
  _standard_generator = std::make_unique<StandardFingerprintGenerator>(preprocess);
  ComputeHash();
  BuildDefaultActiveList();

  return 1;
}

std::vector<std::string>
GFPContext::Tags() const {
  std::vector<std::string> result;
  result.reserve(_components.size());
  for (const Component& component : _components) {
    result.emplace_back(component.tag.rawchars(), component.tag.length());
  }

  return result;
}

int
GFPContext::SetWeight(const const_IWSubstring& tag, float weight) {
  if (weight < 0.0f) {
    std::cerr << "GFPContext::SetWeight:invalid negative weight " << weight << '\n';
    return 0;
  }

  int found = 0;
  for (Component& component : _components) {
    if (component.tag == tag) {
      component.weight = weight;
      ++found;
    }
  }

  if (!found) {
    std::cerr << "GFPContext::SetWeight:no component '" << tag << "'\n";
    return 0;
  }

  _active.clear();
  float sum_weight = 0.0f;
  for (const Component& component : _components) {
    if (component.weight <= 0.0f) {
      continue;
    }
    _active.push_back(ActiveComponent{component.kind, component.index, component.weight});
    sum_weight += component.weight;
  }

  if (sum_weight == 0.0f) {
    std::cerr << "GFPContext::SetWeight:all weights are zero\n";
    return 0;
  }

  for (ActiveComponent& active : _active) {
    active.weight /= sum_weight;
  }

  return 1;
}

int
GFPContext::UseOnly(const std::vector<IWString>& tags) {
  _active.clear();

  for (const IWString& tag : tags) {
    bool found = false;
    for (const Component& component : _components) {
      if (component.tag != tag) {
        continue;
      }
      _active.push_back(
          ActiveComponent{component.kind, component.index, component.weight});
      found = true;
      break;
    }
    if (!found) {
      std::cerr << "GFPContext::UseOnly:no component '" << tag << "'\n";
      _active.clear();
      return 0;
    }
  }

  float sum_weight = 0.0f;
  for (const ActiveComponent& active : _active) {
    sum_weight += active.weight;
  }

  if (sum_weight == 0.0f) {
    std::cerr << "GFPContext::UseOnly:zero total weight\n";
    return 0;
  }

  for (ActiveComponent& active : _active) {
    active.weight /= sum_weight;
  }

  return 1;
}

void
GFPContext::UseAll() {
  _active.clear();

  float sum_weight = 0.0f;
  for (const Component& component : _components) {
    if (component.weight <= 0.0f) {
      continue;
    }
    _active.push_back(ActiveComponent{component.kind, component.index, component.weight});
    sum_weight += component.weight;
  }

  if (sum_weight == 0.0f) {
    BuildDefaultActiveList();
    return;
  }

  for (ActiveComponent& active : _active) {
    active.weight /= sum_weight;
  }
}

int
GFPContext::Fingerprint(Molecule& m, GFPFingerprint& result) {
  if (_standard_generator == nullptr) {
    std::cerr << "GFPContext::Fingerprint:context cannot generate fingerprints\n";
    return 0;
  }

  if (result.context_hash() != _context_hash && !result.Allocate(*this)) {
    return 0;
  }

  return _standard_generator->Generate(m, result);
}

float
GFPContext::Distance(const GFPFingerprint& lhs, const GFPFingerprint& rhs) const {
  if (lhs.context_hash() != _context_hash || rhs.context_hash() != _context_hash) {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float similarity = 0.0f;
  for (const ActiveComponent& active : _active) {
    switch (active.kind) {
      case ComponentKind::kMolecularProperties:
        similarity += active.weight *
                      lhs.molecular_properties().Similarity(rhs.molecular_properties());
        break;
      case ComponentKind::kFixedBinary:
        similarity +=
            active.weight *
            lhs.fixed_binary(active.index).Tanimoto(rhs.fixed_binary(active.index));
        break;
      case ComponentKind::kSparse:
        similarity +=
            active.weight * lhs.sparse(active.index).tanimoto(rhs.sparse(active.index));
        break;
      case ComponentKind::kFixedCounted:
        similarity +=
            active.weight *
            lhs.fixed_counted(active.index).tanimoto(rhs.fixed_counted(active.index));
        break;
    }
  }

  if (similarity > 1.0f && similarity < 1.0001f) {
    similarity = 1.0f;
  }

  return 1.0f - similarity;
}

GFPList::GFPList() : _context(std::make_shared<GFPContext>()) {
}

GFPList::GFPList(std::shared_ptr<GFPContext> context) : _context(std::move(context)) {
}

std::shared_ptr<GFPList>
GFPList::Standard(bool preprocess) {
  auto context = std::make_shared<GFPContext>();
  if (!context->BuildStandard(preprocess)) {
    return nullptr;
  }

  return std::make_shared<GFPList>(context);
}

std::shared_ptr<GFPList>
GFPList::StandardFromMolecules(const std::vector<Molecule*>& molecules, bool preprocess,
                               bool store_metadata) {
  auto result = GFPList::Standard(preprocess);
  if (result == nullptr) {
    return nullptr;
  }

  if (!result->AddMolecules(molecules, store_metadata)) {
    return nullptr;
  }

  return result;
}

int
GFPList::SetStoreMetadata(int store_metadata) {
  const int value = store_metadata ? 1 : 0;
  if (_store_metadata < 0) {
    _store_metadata = value;
    return 1;
  }

  if (_store_metadata == value) {
    return 1;
  }

  std::cerr << "GFPList::SetStoreMetadata:cannot mix metadata and no-metadata entries\n";
  return 0;
}

int
GFPList::Reserve(int size_hint) {
  if (size_hint <= 0) {
    return 1;
  }

  _molecular_properties.reserve(size_hint);
  _fixed_binary.reserve(size_hint * _context->nfixed_binary());
  _sparse.reserve(size_hint * _context->nsparse());
  _fixed_counted.reserve(size_hint * _context->nfixed_counted());
  if (_store_metadata != 0) {
    _smiles.reserve(size_hint);
    _ids.reserve(size_hint);
  }

  return 1;
}

int
GFPList::BuildContextIfNeeded(const IW_TDT& tdt) {
  if (_context->context_hash() != 0) {
    return 1;
  }

  return _context->BuildFromTdt(tdt);
}

int
GFPList::Add(const IW_TDT& tdt) {
  if (!BuildContextIfNeeded(tdt)) {
    return 0;
  }
  if (!SetStoreMetadata(1)) {
    return 0;
  }

  IWString smiles;
  if (!tdt.dataitem_value(kSmilesTag, smiles)) {
    std::cerr << "GFPList::Add:missing smiles\n";
    return 0;
  }
  IWString id;
  if (!tdt.dataitem_value(kIdentifierTag, id)) {
    std::cerr << "GFPList::Add:missing identifier\n";
    return 0;
  }

  for (const Component& component : _context->components()) {
    const_IWSubstring dataitem;
    if (!tdt.dataitem(component.tag, dataitem)) {
      std::cerr << "GFPList::Add:missing component '" << component.tag << "'\n";
      return 0;
    }
    if (dataitem.ends_with('\n')) {
      dataitem.chop();
    }

    switch (component.kind) {
      case ComponentKind::kMolecularProperties: {
        _molecular_properties.emplace_back();
        if (!_molecular_properties.back().Build(dataitem)) {
          return 0;
        }
        break;
      }
      case ComponentKind::kFixedBinary: {
        _fixed_binary.emplace_back();
        if (!_fixed_binary.back().Build(dataitem)) {
          return 0;
        }
        break;
      }
      case ComponentKind::kSparse: {
        _sparse.emplace_back();
        if (!_sparse.back().construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPList::Add:cannot parse sparse '" << component.tag << "'\n";
          return 0;
        }
        break;
      }
      case ComponentKind::kFixedCounted: {
        _fixed_counted.emplace_back();
        if (!_fixed_counted.back().construct_from_tdt_record(dataitem)) {
          std::cerr << "GFPList::Add:cannot parse fixed counted '" << component.tag
                    << "'\n";
          return 0;
        }
        break;
      }
    }
  }

  _smiles.push_back(smiles);
  _ids.push_back(id);
  ++_size;

  return 1;
}

int
GFPList::AddFingerprintComponents(const GFPFingerprint& fp) {
  if (fp.context_hash() != _context->context_hash()) {
    std::cerr << "GFPList::Add:incompatible fingerprint context\n";
    return 0;
  }

  for (const Component& component : _context->components()) {
    switch (component.kind) {
      case ComponentKind::kMolecularProperties:
        _molecular_properties.push_back(fp.molecular_properties());
        break;
      case ComponentKind::kFixedBinary:
        _fixed_binary.push_back(fp.fixed_binary(component.index));
        break;
      case ComponentKind::kSparse:
        _sparse.push_back(fp.sparse(component.index));
        break;
      case ComponentKind::kFixedCounted:
        _fixed_counted.push_back(fp.fixed_counted(component.index));
        break;
    }
  }

  ++_size;
  return 1;
}

int
GFPList::Add(const GFPFingerprint& fp) {
  if (!SetStoreMetadata(0)) {
    return 0;
  }

  return AddFingerprintComponents(fp);
}

int
GFPList::Add(const GFPFingerprint& fp, const IWString& smiles, const IWString& id) {
  if (!SetStoreMetadata(1)) {
    return 0;
  }

  if (!AddFingerprintComponents(fp)) {
    return 0;
  }

  _smiles.push_back(smiles);
  _ids.push_back(id);

  return 1;
}

int
GFPList::Add(Molecule& m) {
  return Add(m, true);
}

int
GFPList::Add(Molecule& m, bool store_metadata) {
  GFPFingerprint fp;
  if (!_context->Fingerprint(m, fp)) {
    return 0;
  }

  if (!store_metadata) {
    return Add(fp);
  }

  IWString smiles(m.smiles());
  return Add(fp, smiles, m.name());
}

int
GFPList::AddMolecules(const std::vector<Molecule*>& molecules, bool store_metadata) {
  if (!SetStoreMetadata(store_metadata)) {
    return 0;
  }
  Reserve(molecules.size());

  for (Molecule* molecule : molecules) {
    if (molecule == nullptr) {
      std::cerr << "GFPList::AddMolecules:null molecule pointer\n";
      return 0;
    }
    if (!Add(*molecule, store_metadata)) {
      return 0;
    }
  }

  return 1;
}

int
GFPList::ReadFile(const char* fname, int size_hint) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    std::cerr << "GFPList::ReadFile:cannot open '" << fname << "'\n";
    return 0;
  }

  if (size_hint <= 0) {
    size_hint = CountTdts(input, fname);
    if (size_hint < 0) {
      return 0;
    }
  }

  IW_TDT tdt;
  if (!tdt.next(input)) {
    std::cerr << "GFPList::ReadFile:no TDT records in '" << fname << "'\n";
    return 0;
  }

  if (!BuildContextIfNeeded(tdt)) {
    return 0;
  }

  Reserve(size_hint);

  if (!Add(tdt)) {
    return 0;
  }

  while (tdt.next(input)) {
    if (!Add(tdt)) {
      return 0;
    }
  }

  return 1;
}

float
GFPList::Distance(int i, int j) const {
  float similarity = 0.0f;
  for (const ActiveComponent& active : _context->active_components()) {
    switch (active.kind) {
      case ComponentKind::kMolecularProperties:
        similarity +=
            active.weight * _molecular_properties[i].Similarity(_molecular_properties[j]);
        break;
      case ComponentKind::kFixedBinary: {
        const int nfixed = _context->nfixed_binary();
        similarity += active.weight * _fixed_binary[i * nfixed + active.index].Tanimoto(
                                          _fixed_binary[j * nfixed + active.index]);
        break;
      }
      case ComponentKind::kSparse: {
        const int nsparse = _context->nsparse();
        similarity += active.weight * _sparse[i * nsparse + active.index].tanimoto(
                                          _sparse[j * nsparse + active.index]);
        break;
      }
      case ComponentKind::kFixedCounted: {
        const int nfc = _context->nfixed_counted();
        similarity += active.weight * _fixed_counted[i * nfc + active.index].tanimoto(
                                          _fixed_counted[j * nfc + active.index]);
        break;
      }
    }
  }

  if (similarity > 1.0f && similarity < 1.0001f) {
    similarity = 1.0f;
  }

  return 1.0f - similarity;
}

float
GFPList::Distance(const GFPFingerprint& fp, int j) const {
  if (fp.context_hash() != _context->context_hash()) {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float similarity = 0.0f;
  for (const ActiveComponent& active : _context->active_components()) {
    switch (active.kind) {
      case ComponentKind::kMolecularProperties:
        similarity += active.weight *
                      fp.molecular_properties().Similarity(_molecular_properties[j]);
        break;
      case ComponentKind::kFixedBinary: {
        const int nfixed = _context->nfixed_binary();
        similarity +=
            active.weight * fp.fixed_binary(active.index)
                                .Tanimoto(_fixed_binary[j * nfixed + active.index]);
        break;
      }
      case ComponentKind::kSparse: {
        const int nsparse = _context->nsparse();
        similarity +=
            active.weight *
            fp.sparse(active.index).tanimoto(_sparse[j * nsparse + active.index]);
        break;
      }
      case ComponentKind::kFixedCounted: {
        const int nfc = _context->nfixed_counted();
        similarity +=
            active.weight * fp.fixed_counted(active.index)
                                .tanimoto(_fixed_counted[j * nfc + active.index]);
        break;
      }
    }
  }

  if (similarity > 1.0f && similarity < 1.0001f) {
    similarity = 1.0f;
  }

  return 1.0f - similarity;
}

std::vector<NearestNeighbour>
GFPList::NearestNeighbours(int query, int k) const {
  if (k <= 0 || query < 0 || query >= size()) {
    return {};
  }

  NearestNeighbourMaxHeap heap;
  for (int i = 0; i < size(); ++i) {
    if (i == query) {
      continue;
    }

    MaybeAddNearestNeighbour(heap, k, NearestNeighbour{i, Distance(query, i)});
  }

  return SortedNearestNeighbours(heap);
}

std::vector<NearestNeighbour>
GFPList::NearestNeighbours(const GFPFingerprint& query, int k) const {
  if (k <= 0 || query.context_hash() != _context->context_hash()) {
    return {};
  }

  NearestNeighbourMaxHeap heap;
  for (int i = 0; i < size(); ++i) {
    MaybeAddNearestNeighbour(heap, k, NearestNeighbour{i, Distance(query, i)});
  }

  return SortedNearestNeighbours(heap);
}

std::vector<NearestNeighbour>
GFPList::NearestNeighboursWithinDistance(int query, float max_distance) const {
  std::vector<NearestNeighbour> result;
  if (query < 0 || query >= size() || max_distance < 0.0f) {
    return result;
  }

  for (int i = 0; i < size(); ++i) {
    if (i == query) {
      continue;
    }

    const float distance = Distance(query, i);
    if (distance <= max_distance) {
      result.push_back(NearestNeighbour{i, distance});
    }
  }

  std::sort(result.begin(), result.end(),
            [](const NearestNeighbour& lhs, const NearestNeighbour& rhs) {
              if (lhs.distance != rhs.distance) {
                return lhs.distance < rhs.distance;
              }
              return lhs.index < rhs.index;
            });

  return result;
}

std::vector<NearestNeighbour>
GFPList::NearestNeighboursWithinDistance(const GFPFingerprint& query,
                                         float max_distance) const {
  std::vector<NearestNeighbour> result;
  if (query.context_hash() != _context->context_hash() || max_distance < 0.0f) {
    return result;
  }

  for (int i = 0; i < size(); ++i) {
    const float distance = Distance(query, i);
    if (distance <= max_distance) {
      result.push_back(NearestNeighbour{i, distance});
    }
  }

  std::sort(result.begin(), result.end(),
            [](const NearestNeighbour& lhs, const NearestNeighbour& rhs) {
              if (lhs.distance != rhs.distance) {
                return lhs.distance < rhs.distance;
              }
              return lhs.index < rhs.index;
            });

  return result;
}

}  // namespace gfp_context
