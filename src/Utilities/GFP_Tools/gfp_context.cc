#include "Utilities/GFP_Tools/gfp_context.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"

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

bool
TagStartsWith(const IWString& tag, const char* prefix) {
  const int n = static_cast<int>(strlen(prefix));
  if (tag.length() < n) {
    return false;
  }

  return strncmp(tag.rawchars(), prefix, n) == 0;
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

MolecularProperties::~MolecularProperties() {
  delete [] _property;
}

MolecularProperties::MolecularProperties(const MolecularProperties& rhs) {
  *this = rhs;
}

MolecularProperties&
MolecularProperties::operator=(const MolecularProperties& rhs) {
  if (this == &rhs) {
    return *this;
  }

  delete [] _property;
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

  delete [] _property;

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
    delete [] _property;
    _nproperties = nproperties;
    _property = new int[_nproperties];
  }

  const unsigned char* bits = fp.bits();
  for (int i = 0; i < _nproperties; ++i) {
    _property[i] = bits[i];
  }

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

float
FixedBitVector::Tanimoto(const FixedBitVector& rhs) const {
  return FixedBinaryTanimoto(*this, rhs);
}

void
GFPFingerprint::FreeArrays() {
  delete [] _fixed_binary;
  delete [] _sparse;
  delete [] _fixed_counted;

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
          std::cerr << "GFPFingerprint::Build:cannot parse sparse '" << component.tag << "'\n";
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

int
GFPContext::BuildFromTdt(const IW_TDT& tdt) {
  _components.clear();
  _active.clear();
  _nfixed_binary = 0;
  _nsparse = 0;
  _nfixed_counted = 0;
  _has_molecular_properties = false;
  _context_hash = 0;

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
    } else if (TagStartsWith(tag, "FP")) {
      _components.push_back(Component{ComponentKind::kFixedBinary, tag, _nfixed_binary, 1.0f});
      ++_nfixed_binary;
    } else if (TagStartsWith(tag, "NC")) {
      _components.push_back(Component{ComponentKind::kSparse, tag, _nsparse, 1.0f});
      ++_nsparse;
    } else if (TagStartsWith(tag, "FC")) {
      _components.push_back(Component{ComponentKind::kFixedCounted, tag, _nfixed_counted, 1.0f});
      ++_nfixed_counted;
    }
  }

  if (_components.empty()) {
    std::cerr << "GFPContext::BuildFromTdt:no fingerprint components\n";
    return 0;
  }

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
      _active.push_back(ActiveComponent{component.kind, component.index, component.weight});
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

float
GFPContext::Distance(const GFPFingerprint& lhs, const GFPFingerprint& rhs) const {
  if (lhs.context_hash() != _context_hash || rhs.context_hash() != _context_hash) {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float similarity = 0.0f;
  for (const ActiveComponent& active : _active) {
    switch (active.kind) {
      case ComponentKind::kMolecularProperties:
        similarity += active.weight * lhs.molecular_properties().Similarity(rhs.molecular_properties());
        break;
      case ComponentKind::kFixedBinary:
        similarity += active.weight * lhs.fixed_binary(active.index).Tanimoto(rhs.fixed_binary(active.index));
        break;
      case ComponentKind::kSparse:
        similarity += active.weight * lhs.sparse(active.index).tanimoto(rhs.sparse(active.index));
        break;
      case ComponentKind::kFixedCounted:
        similarity += active.weight * lhs.fixed_counted(active.index).tanimoto(rhs.fixed_counted(active.index));
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

int
GFPList::Reserve(int size_hint) {
  if (size_hint <= 0) {
    return 1;
  }

  _molecular_properties.reserve(size_hint);
  _fixed_binary.reserve(size_hint * _context->nfixed_binary());
  _sparse.reserve(size_hint * _context->nsparse());
  _fixed_counted.reserve(size_hint * _context->nfixed_counted());
  _smiles.reserve(size_hint);
  _ids.reserve(size_hint);

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
          std::cerr << "GFPList::Add:cannot parse fixed counted '" << component.tag << "'\n";
          return 0;
        }
        break;
      }
    }
  }

  _smiles.push_back(smiles);
  _ids.push_back(id);

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
        similarity += active.weight * _molecular_properties[i].Similarity(_molecular_properties[j]);
        break;
      case ComponentKind::kFixedBinary: {
        const int nfixed = _context->nfixed_binary();
        similarity += active.weight *
            _fixed_binary[i * nfixed + active.index].Tanimoto(_fixed_binary[j * nfixed + active.index]);
        break;
      }
      case ComponentKind::kSparse: {
        const int nsparse = _context->nsparse();
        similarity += active.weight *
            _sparse[i * nsparse + active.index].tanimoto(_sparse[j * nsparse + active.index]);
        break;
      }
      case ComponentKind::kFixedCounted: {
        const int nfc = _context->nfixed_counted();
        similarity += active.weight *
            _fixed_counted[i * nfc + active.index].tanimoto(_fixed_counted[j * nfc + active.index]);
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
  std::vector<NearestNeighbour> result;
  if (k <= 0 || query < 0 || query >= size()) {
    return result;
  }

  result.reserve(size() > 1 ? std::min(k, size() - 1) : 0);
  for (int i = 0; i < size(); ++i) {
    if (i == query) {
      continue;
    }

    result.push_back(NearestNeighbour{i, Distance(query, i)});
  }

  auto compare = [](const NearestNeighbour& lhs, const NearestNeighbour& rhs) {
    if (lhs.distance != rhs.distance) {
      return lhs.distance < rhs.distance;
    }
    return lhs.index < rhs.index;
  };

  if (static_cast<int>(result.size()) > k) {
    std::nth_element(result.begin(), result.begin() + k, result.end(), compare);
    result.resize(k);
  }

  std::sort(result.begin(), result.end(), compare);

  return result;
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

}  // namespace gfp_context
