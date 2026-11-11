#include <algorithm>

#include "Foundational/iwmisc/proto_support.h"

#include "Utilities/GFP_Tools/bit_subset.h"

namespace bit_subset {

using std::cerr;

BitSubset::BitSubset() {
  _flatten_sparse_counted = 0;

  _nfixed = 0;
  _fixed_mask = nullptr;
  _nsparse = 0;
  _sparse_subset = nullptr;

}

BitSubset::~BitSubset() {
  if (_fixed_mask != nullptr) {
    delete [] _fixed_mask;
    _nfixed = 0;
  }

  if (_sparse_subset != nullptr) {
    delete [] _sparse_subset;
    _nsparse = 0;
  }
}

void
FormMask(std::unordered_set<uint32_t>& subset,
         const int nbits,
         IWDYFP& mask) {
  mask.allocate_space_for_bits(nbits);
  mask.clear();
  for (const auto iter : subset) {
    mask.set_bit(iter);
  }
}

int
BitSubset::InitialiseGfpKnown(const IW_General_Fingerprint& gfp) {
  if (number_fixed_size_counted_fingerprints() > 0) {
    cerr << "BitSubset::InitialiseGfpKnown:fixed size counted fingerprints are not supported\n";
    return 0;
  }

  _nfixed = number_fingerprints();
  _nsparse = number_sparse_fingerprints();
  if (_nfixed == 0 && _nsparse == 0) {
    cerr << "BitSubset::InitialiseGfpKnown:no fingerprints\n";
    return 0;
  }

  _fixed_mask = new IWDYFP[_nfixed];
  _sparse_subset = new Subset*[_nsparse];

  for (int i = 0; i < _nfixed; ++i) {
    const IWString& tag = fixed_fingerprint_tag(i);
    const auto iter =  _tag_to_bit_subset.find(tag);
    if (iter == _tag_to_bit_subset.end()) {
      cerr << "BitSubset::InitialiseGfpKnown:no fixed '" << tag << "' tag, nfixed " << _nfixed << '\n';
      return 0;
    }
    const IWDYFP& fp = gfp[i];
    FormMask(iter->second, fp.nbits(), _fixed_mask[i]);
  }

  for (int i = 0; i < _nsparse; ++i) {
    const IWString& tag = sparse_fingerprint_tag(i);
    auto iter = _tag_to_bit_subset.find(tag);
    if (iter == _tag_to_bit_subset.end()) {
      cerr << "BitSubset::InitialiseGfpKnown:no sparse '" << tag << "' tag, nsparse " << _nsparse << '\n';
      return 0;
    }
    _sparse_subset[i] = &(iter->second);
  }
  
  return 1;
}

int
BitSubset::Build(IWString& fname) {
  std::optional<GfpBitSubset::GfpBitSubset> proto =
                   iwmisc::ReadBinaryProto<GfpBitSubset::GfpBitSubset>(fname);
  if (! proto) {
    cerr << "BitSubset::Build:cannot read '" << fname << "'\n";
    return 0;
  }

  return Build(*proto);
}

int
BitSubset::Build(const GfpBitSubset::GfpBitSubset& proto) {
  for (const auto& [key, value] : proto.xref()) {
    IWString tag(key);
    Subset subset;

    for (const auto iter : value.bits()) {
      subset.insert(iter);
    }
    _tag_to_bit_subset.emplace(tag, std::move(subset));
  }

  _flatten_sparse_counted = proto.params().flatten_sparse_fingerprints();

  return _tag_to_bit_subset.size();
}

BitXref::BitXref() {
  _flatten_sparse_counted = 0;

  _nfixed = 0;
  _fixed = nullptr;
  _nsparse = 0;
  _sparse = nullptr;

  _nfixed_counted = 0;
  _fixed_counted = nullptr;

  _has_properties = false;
  _properties = nullptr;

  _highest_feature_number = 0;
}

BitXref::~BitXref() {
  if (_fixed != nullptr) {
    delete [] _fixed;
    _nfixed = 0;
  }

  if (_sparse != nullptr) {
    delete [] _sparse;
    _nsparse = 0;
  }

  if (_fixed_counted != nullptr) {
    delete [] _fixed_counted;
    _nfixed_counted = 0;
  }
}

int
BitXref::Build(IWString& fname) {
  std::optional<GfpBitSubset::GfpBitToFeature> proto = iwmisc::ReadBinaryProto<GfpBitSubset::GfpBitToFeature>(fname);
  if (! proto) {
    cerr << "BitXref::Build:cannot read '" << fname << "'\n";
    return 0;
  }

  return Build(*proto);
}

int
BitXref::Build(const GfpBitSubset::GfpBitToFeature& proto) {
  for (const auto& [key, value] : proto.xref()) {
    IWString tag(key);
    BitToFeature bit_to_feature;

    for (auto [key, value] : value.bit_to_feature()) {
      bit_to_feature.emplace(static_cast<gfp_bit_type>(key), static_cast<uint32_t>(value));
    } 
    _tag_to_bit_to_feature.emplace(tag, std::move(bit_to_feature));
  }

  _flatten_sparse_counted = proto.params().flatten_sparse_fingerprints();
  return 1;
}

struct FeatureCount {
  uint32_t feature;
  uint32_t count;
};

// Return the number of features set in `gfp`.
int
NumberFeatures(const IW_General_Fingerprint& gfp) {
  int nfixed = number_fingerprints();
  int nsparse = number_sparse_fingerprints();
  int nfeatures = 0;
  for (int i = 0; i < nfixed; ++i) {
    nfeatures += gfp[i].nset();
  }
  for (int i = 0; i < nsparse; ++i) {
    nfeatures += gfp.sparse_fingerprint(i).nset();
  }

  // An upper bound, positions with a count of zero are not features.
  const int nfixed_counted = number_fixed_size_counted_fingerprints();
  for (int i = 0; i < nfixed_counted; ++i) {
    nfeatures += gfp.fixed_size_counted_fingerprint(i).size();
  }

  nfeatures += gfp.molecular_properties_integer().nproperties();

  return nfeatures;
}

int
BitSubset::MakeSubset(IW_General_Fingerprint& gfp) {
  if (_nfixed == 0 && _nsparse == 0 && ! InitialiseGfpKnown(gfp)) {
    cerr << "BitSubset::MakeSubset:cannot initialise\n";
    return -1;  // Special return value indicating fatal error.
  }

  int rc = 0;
  for (int i = 0; i < _nfixed; ++i) {
    IWDYFP& fp = gfp[i];
    const int initial_nset = fp.nset();
    int notused = 0;
    fp.iwand(_fixed_mask[i], notused);
    rc += initial_nset - fp.nset();
  }

  for (int i = 0; i < _nsparse; ++i) {
    Sparse_Fingerprint & sparse_fingerprint = gfp.sparse_fingerprint(i);

    rc += sparse_fingerprint.ReduceToSubset(*_sparse_subset[i]);
  }

  return rc;
}


int
BitXref::WriteSvmlFeatures(const IW_General_Fingerprint& gfp,
                        IWString_and_File_Descriptor& output) {
  if (_nfixed == 0 && _nsparse == 0 && _nfixed_counted == 0 && ! _has_properties &&
      ! InitialiseGfpKnown(gfp)) {
    cerr << "BitXref::WriteSvmlFeatures:cannot initialise\n";
    return -1;  // Special return value indicating fatal error.
  }

//cerr << "BitXref::WriteSvmlFeatures writing " << _nfixed << " fixed and " << _nsparse << " sparse fingerprints\n";

  // The maximum number of features.
  const int nfeatures = NumberFeatures(gfp);

  int ndx = 0;
  std::unique_ptr<FeatureCount[]> feature_counts(new FeatureCount[nfeatures]);

  if (_has_properties) {
    const Molecular_Properties_Integer& mpi = gfp.molecular_properties_integer();
    const uint8_t* values = mpi.rawdata();
    for (int i = 0; i < mpi.nproperties(); ++i) {
      if (values[i] == 0) {
        continue;
      }
      const auto iter = _properties->find(i);
      if (iter == _properties->cend()) {
        continue;
      }
      feature_counts[ndx].feature = iter->second;
      feature_counts[ndx].count = values[i];
      ndx++;
    }
  }

  for (int i = 0; i < _nfixed; ++i) {
    const BitToFeature* b2f = _fixed[i];
    const IWDYFP & fp = gfp[i];
    int j = 0;
    int bit;
    while ((bit = fp.next_on_bit(j)) >= 0) {
      const auto iter = b2f->find(bit);
      if (iter == b2f->cend()) {
        continue;
      }
      feature_counts[ndx].feature = iter->second;
      feature_counts[ndx].count = 1;
      ndx++;
    }
  }

  for (int i = 0; i < _nsparse; ++i) {
    const BitToFeature* b2f = _sparse[i];
    const Sparse_Fingerprint& fp = gfp.sparse_fingerprint(i);
    int j = 0;
    uint32_t b = 0;
    int c = 0;
    while (fp.next_bit_set(j, b, c)) {
      const auto iter = b2f->find(b);
      if (iter == b2f->end()) {
        continue;
      }
      feature_counts[ndx].feature = iter->second;
      feature_counts[ndx].count = c;
      ndx++;
    }
  }

  for (int i = 0; i < _nfixed_counted; ++i) {
    const BitToFeature* b2f = _fixed_counted[i];
    const auto& fp = gfp.fixed_size_counted_fingerprint(i);
    for (int j = 0; j < fp.size(); ++j) {
      const uint32_t c = fp.count(j);
      if (c == 0) {
        continue;
      }
      const auto iter = b2f->find(j);
      if (iter == b2f->cend()) {
        continue;
      }
      feature_counts[ndx].feature = iter->second;
      feature_counts[ndx].count = c;
      ndx++;
    }
  }

  std::sort(feature_counts.get(), feature_counts.get() + ndx,
            [](const FeatureCount& fc1, const FeatureCount& fc2) {
    return fc1.feature < fc2.feature;
  });

  for (int i = 0; i < ndx; ++i) {
    output << ' ' << feature_counts[i].feature << ':' << feature_counts[i].count;
  }

  return ndx;  // The number of features written.
}

int
BitXref::InitialiseGfpKnown(const IW_General_Fingerprint& gfp) {
  _nfixed = number_fingerprints();
  _nsparse = number_sparse_fingerprints();
  _nfixed_counted = number_fixed_size_counted_fingerprints();
  _has_properties = gfp.molecular_properties_integer().active();
  _fixed = new BitToFeature*[_nfixed];
  _sparse = new BitToFeature*[_nsparse];
  _fixed_counted = new BitToFeature*[_nfixed_counted];
  for (int i = 0; i < _nfixed; ++i) {
    const IWString tag = fixed_fingerprint_tag(i);
    const auto iter = _tag_to_bit_to_feature.find(tag);
    if (iter == _tag_to_bit_to_feature.cend()) {
      cerr << "BitXref::InitialiseGfpKnown:no fixed '" << tag << "' fingerprint in input\n";
      return 0;
    }
    _fixed[i] = &(iter->second);
  }

  for (int i = 0; i < _nsparse; ++i) {
    const IWString tag = sparse_fingerprint_tag(i);
    const auto iter = _tag_to_bit_to_feature.find(tag);
    if (iter == _tag_to_bit_to_feature.cend()) {
      cerr << "BitXref::InitialiseGfpKnown:no sparse '" << tag << "' fingerprint in input\n";
      for (auto& [k, v] : _tag_to_bit_to_feature) {
        cerr << "key '" << k << "'\n";
      }
      return 0;
    }
    _sparse[i] = &(iter->second);
  }

  for (int i = 0; i < _nfixed_counted; ++i) {
    const IWString tag = fixed_size_counted_fingerprint_tag(i);
    const auto iter = _tag_to_bit_to_feature.find(tag);
    if (iter == _tag_to_bit_to_feature.cend()) {
      cerr << "BitXref::InitialiseGfpKnown:no fixed counted '" << tag << "' fingerprint in input\n";
      return 0;
    }
    _fixed_counted[i] = &(iter->second);
  }

  if (_has_properties) {
    const IWString tag = property_tag();
    const auto iter = _tag_to_bit_to_feature.find(tag);
    if (iter == _tag_to_bit_to_feature.cend()) {
      cerr << "BitXref::InitialiseGfpKnown:no molecular properties '" << tag << "' in input\n";
      return 0;
    }
    _properties = &(iter->second);
  }

  // Every fingerprint in the input is known. The cross reference must not
  // contain any others, or features that were used to build the model will
  // be silently absent from the output.
  const int nfingerprints = _nfixed + _nsparse + _nfixed_counted + (_has_properties ? 1 : 0);
  if (_tag_to_bit_to_feature.size() != static_cast<size_t>(nfingerprints)) {
    cerr << "BitXref::InitialiseGfpKnown:cross reference contains "
         << _tag_to_bit_to_feature.size() << " fingerprints, but the input has "
         << nfingerprints << '\n';
    for (const auto& [tag, _] : _tag_to_bit_to_feature) {
      bool in_input = false;
      for (int i = 0; i < _nfixed; ++i) {
        in_input = in_input || (tag == fixed_fingerprint_tag(i));
      }
      for (int i = 0; i < _nsparse; ++i) {
        in_input = in_input || (tag == sparse_fingerprint_tag(i));
      }
      for (int i = 0; i < _nfixed_counted; ++i) {
        in_input = in_input || (tag == fixed_size_counted_fingerprint_tag(i));
      }
      in_input = in_input || (_has_properties && tag == property_tag());
      if (!in_input) {
        cerr << "BitXref::InitialiseGfpKnown:fingerprint '" << tag
             << "' is in the cross reference, but not in the input\n";
      }
    }
    return 0;
  }

  return nfingerprints;
}

uint32_t
BitXref::HighestFeatureNumber() const {
  uint32_t result = 0;
  for (const auto& [tag, bit_to_feature] : _tag_to_bit_to_feature) {
    for (const auto& [_,  feature] : bit_to_feature) {
      if (feature > result) {
        result = feature;
      }
    }
  }

  return result;
}

template <typename T>
int
BitXref::PopulateFeatureVector(const IW_General_Fingerprint& gfp,
                               std::vector<T>& features) const {
  std::fill(features.begin(), features.end(), 0.0);
  int rc = 0;

  if (_has_properties) {
    const Molecular_Properties_Integer& mpi = gfp.molecular_properties_integer();
    const uint8_t* values = mpi.rawdata();
    for (int i = 0; i < mpi.nproperties(); ++i) {
      if (values[i] == 0) {
        continue;
      }
      const auto iter = _properties->find(i);
      if (iter == _properties->end()) {
        continue;
      }
      assert(iter->second < features.size());
      features[iter->second] = static_cast<float>(values[i]);
      rc++;
    }
  }

  for (int i = 0; i < _nfixed; ++i) {
    IWDYFP& fp = gfp[i];
    int j = 0;
    int bit;
    while ((bit = fp.next_on_bit(j)) >= 0) {
      const auto iter = _fixed[i]->find(bit);
      if (iter == _fixed[i]->end()) {
        continue;
      }
      assert(iter->second < features.size());
      features[iter->second] = 1.0f;
      rc++;
    }
  }

  for (int i = 0; i < _nsparse; ++i) {
    const Sparse_Fingerprint & sfp = gfp.sparse_fingerprint(i);
    int j = 0;
    unsigned int bit;
    int count;
    while (sfp.next_bit_set(j, bit, count)) {
      const auto iter = _sparse[i]->find(bit);
      if (iter == _sparse[i]->end()) {
        continue;
      }
      assert(iter->second < features.size());
      features[iter->second] = static_cast<float>(count);
      rc++;
    }
  }

  for (int i = 0; i < _nfixed_counted; ++i) {
    const auto& fp = gfp.fixed_size_counted_fingerprint(i);
    for (int j = 0; j < fp.size(); ++j) {
      const int count = fp.count(j);
      if (count == 0) {
        continue;
      }
      const auto iter = _fixed_counted[i]->find(j);
      if (iter == _fixed_counted[i]->end()) {
        continue;
      }
      assert(iter->second < features.size());
      features[iter->second] = static_cast<float>(count);
      rc++;
    }
  }

  return rc;
}
template int BitXref::PopulateFeatureVector(const IW_General_Fingerprint& gfp, std::vector<float>&) const;

int
BitXref::WriteDsv(const IW_General_Fingerprint& gfp, char output_separator,
                 IWString_and_File_Descriptor& output) {
  if (_highest_feature_number <= 0) {
    _highest_feature_number = HighestFeatureNumber();
    if (_highest_feature_number <= 0) {
      cerr << "BitXref::WriteDsv:highest feature number not discerned\n";
      return -1;  // Special return value indicating fatal error.
    }
  }

  std::vector<int> features(_highest_feature_number + 1);
  int rc = PopulateFeatureVector(gfp, features);

  output.reserve(_highest_feature_number * 2);
  output << features[0];
  for (uint32_t i = 1; i <= _highest_feature_number; ++i) {
    output << output_separator << features[i];
    output.write_if_buffer_holds_more_than(1024);  // Help pipe handling.
  }

  return rc;
}

int
BitSubset::DebugPrint(std::ostream& output) const {
  output << "BitXref::DebugPrint: _nfixed " << _nfixed << " _nsparse " << _nsparse << "\n";

  return 1;
}

}  // namespace bit_subset
