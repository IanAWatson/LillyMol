#include "Utilities/GFP_Tools/train_test_split_activity.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>

#include "absl/container/flat_hash_map.h"

#include "Foundational/iwmisc/iw_tabular_data.h"
#include "Foundational/iwstring/absl_hash.h"

namespace train_test_split {

using std::cerr;

namespace {

char
SeparatorForFileName(const IWString& fname) {
  if (fname.ends_with(".csv")) {
    return ',';
  }

  return ' ';
}

int
ParseBucketisation(const const_IWSubstring& value, ActivityBucketisation& bucketisation) {
  if (value == "quantile" || value == "quantiles" || value == "q") {
    bucketisation = ActivityBucketisation::kQuantile;
    return 1;
  }
  if (value == "width" || value == "equal_width" || value == "range" || value == "w") {
    bucketisation = ActivityBucketisation::kEqualWidth;
    return 1;
  }

  cerr << "Unrecognised activity bucketisation '" << value << "'\n";
  return 0;
}

int
ApplyDirective(const const_IWSubstring& token, ActivityProfileOptions& options) {
  if (token.starts_with("buckets=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(8);
    if (! value.numeric_value(options.buckets) || options.buckets == 0 ||
        options.buckets > std::numeric_limits<uint8_t>::max()) {
      cerr << "Invalid activity bucket count '" << token << "'\n";
      return 0;
    }
    return 1;
  }

  if (token.starts_with("bucket=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(7);
    return ParseBucketisation(value, options.bucketisation);
  }

  if (token.starts_with("bucketisation=")) {
    const_IWSubstring value(token);
    value.remove_leading_chars(14);
    return ParseBucketisation(value, options.bucketisation);
  }

  return ParseBucketisation(token, options.bucketisation);
}

uint64_t
AbsDiff(uint64_t v1, uint64_t v2) {
  return v1 > v2 ? v1 - v2 : v2 - v1;
}

}  // namespace

int
ActivityProfileSet::AddDirective(const const_IWSubstring& directive) {
  resizable_array_p<const_IWSubstring> tokens;
  if (directive.split(tokens, ':') == 0) {
    cerr << "ActivityProfileSet::AddDirective:empty directive\n";
    return 0;
  }

  ActivityProfileOptions options = _defaults;
  const_IWSubstring first(*tokens[0]);

  if (first.contains('=')) {
    for (int i = 0; i < tokens.number_elements(); ++i) {
      if (! ApplyDirective(*tokens[i], _defaults)) {
        return 0;
      }
    }
    return 1;
  }

  options.fname = first;
  for (int i = 1; i < tokens.number_elements(); ++i) {
    if (! ApplyDirective(*tokens[i], options)) {
      return 0;
    }
  }

  if (options.fname.empty()) {
    cerr << "ActivityProfileSet::AddDirective:no file name in '" << directive << "'\n";
    return 0;
  }

  _pending.push_back(std::move(options));
  return 1;
}

int
ActivityProfileSet::Build(const std::vector<std::string>& id,
                          const std::vector<uint32_t>& weight) {
  _profile.reserve(_pending.size());
  for (const ActivityProfileOptions& options : _pending) {
    ActivityProfile profile;
    if (! profile.Build(options, id, weight)) {
      return 0;
    }
    _profile.push_back(std::move(profile));
  }

  return 1;
}

int
ActivityProfile::Build(const ActivityProfileOptions& options,
                       const std::vector<std::string>& id,
                       const std::vector<uint32_t>& weight) {
  if (id.size() != weight.size()) {
    cerr << "ActivityProfile::Build:size mismatch " << id.size() << " vs " << weight.size() << '\n';
    return 0;
  }

  _options = options;

  IW_Tabular_Data<double> reader;
  reader.set_has_header(1);
  reader.set_first_column_is_identifier(1);

  const char sep = SeparatorForFileName(_options.fname);
  if (! reader.build(_options.fname.null_terminated_chars(), sep)) {
    cerr << "ActivityProfile::Build:cannot read activity data from '" << _options.fname << "'\n";
    return 0;
  }
  if (reader.ncols() == 0) {
    cerr << "ActivityProfile::Build:no activity columns in '" << _options.fname << "'\n";
    return 0;
  }

  absl::flat_hash_map<IWString, uint32_t> id_to_ndx;
  id_to_ndx.reserve(id.size());
  for (uint32_t i = 0; i < id.size(); ++i) {
    id_to_ndx[IWString(id[i])] = i;
  }

  std::vector<double> activity(id.size(), 0.0);
  std::vector<uint8_t> seen(id.size(), 0);
  uint32_t matched = 0;

  for (int i = 0; i < reader.nrows(); ++i) {
    const IWString& activity_id = reader.ids()[i];
    const auto iter = id_to_ndx.find(activity_id);
    if (iter == id_to_ndx.end()) {
      continue;
    }

    const uint32_t ndx = iter->second;
    if (seen[ndx]) {
      cerr << "ActivityProfile::Build:duplicate activity for '" << activity_id << "'\n";
      return 0;
    }
    seen[ndx] = 1;
    activity[ndx] = reader.value(i, 0);
    ++matched;
  }

  if (matched != id.size()) {
    cerr << "ActivityProfile::Build:have " << id.size() << " items but only read " <<
            matched << " activity values from '" << _options.fname << "'\n";
    return 0;
  }

  _scaled.resize(id.size());
  const auto [min_iter, max_iter] = std::minmax_element(activity.begin(), activity.end());
  const double minval = *min_iter;
  const double maxval = *max_iter;
  if (minval == maxval) {
    std::fill(_scaled.begin(), _scaled.end(), 0);
  } else {
    const double range = maxval - minval;
    for (uint32_t i = 0; i < activity.size(); ++i) {
      int scaled = static_cast<int>(std::lround(100.0 * (activity[i] - minval) / range));
      scaled = std::clamp(scaled, 0, 100);
      _scaled[i] = static_cast<uint8_t>(scaled);
    }
  }

  return AssignBuckets(activity);
}

int
ActivityProfile::AssignBuckets(const std::vector<double>& activity) {
  _bucket.resize(activity.size());
  _total_bucket_count.assign(_options.buckets, 0);
  _train_bucket_count.assign(_options.buckets, 0);

  if (_options.bucketisation == ActivityBucketisation::kEqualWidth) {
    return AssignEqualWidthBuckets();
  }

  return AssignQuantileBuckets(activity);
}

int
ActivityProfile::AssignEqualWidthBuckets() {
  for (uint32_t i = 0; i < _scaled.size(); ++i) {
    _bucket[i] = std::min<uint32_t>(_options.buckets - 1,
                                    (static_cast<uint32_t>(_scaled[i]) * _options.buckets) / 101);
  }

  return 1;
}

int
ActivityProfile::AssignQuantileBuckets(const std::vector<double>& activity) {
  std::vector<std::pair<double, uint32_t>> sorted;
  sorted.reserve(activity.size());
  for (uint32_t i = 0; i < activity.size(); ++i) {
    sorted.emplace_back(activity[i], i);
  }

  std::sort(sorted.begin(), sorted.end());

  uint32_t begin = 0;
  while (begin < sorted.size()) {
    uint32_t end = begin + 1;
    while (end < sorted.size() && sorted[end].first == sorted[begin].first) {
      ++end;
    }

    const uint32_t bucket = std::min<uint32_t>(_options.buckets - 1,
                                (begin * _options.buckets) / sorted.size());
    for (uint32_t i = begin; i < end; ++i) {
      _bucket[sorted[i].second] = bucket;
    }
    begin = end;
  }

  return 1;
}

uint64_t
ActivityProfile::Penalty(const std::vector<uint64_t>& train_bucket_count,
                         uint64_t train_weight) const {
  uint64_t penalty = 0;
  for (uint32_t i = 0; i < _options.buckets; ++i) {
    penalty += AbsDiff(train_bucket_count[i] * _total_weight,
                       _total_bucket_count[i] * train_weight);
  }

  return penalty;
}

int
ActivityProfile::InitialiseSplit(const std::vector<int>& in_train,
                                 const std::vector<uint32_t>& weight) {
  if (in_train.size() != _bucket.size() || weight.size() != _bucket.size()) {
    cerr << "ActivityProfile::InitialiseSplit:size mismatch\n";
    return 0;
  }

  std::fill(_total_bucket_count.begin(), _total_bucket_count.end(), 0);
  std::fill(_train_bucket_count.begin(), _train_bucket_count.end(), 0);
  _total_weight = 0;
  _train_weight = 0;

  for (uint32_t i = 0; i < _bucket.size(); ++i) {
    const uint32_t b = _bucket[i];
    _total_bucket_count[b] += weight[i];
    _total_weight += weight[i];
    if (in_train[i]) {
      _train_bucket_count[b] += weight[i];
      _train_weight += weight[i];
    }
  }

  _penalty = Penalty(_train_bucket_count, _train_weight);
  return 1;
}

uint64_t
ActivityProfile::RecomputePenalty(const std::vector<int>& in_train,
                                  const std::vector<uint32_t>& weight) {
  InitialiseSplit(in_train, weight);
  return _penalty;
}

int64_t
ActivityProfile::DeltaForSwap(uint32_t out_of_train, uint32_t into_train,
                              const std::vector<uint32_t>& weight) const {
  std::vector<uint64_t> train_bucket_count(_train_bucket_count);
  uint64_t train_weight = _train_weight;

  train_bucket_count[_bucket[out_of_train]] -= weight[out_of_train];
  train_bucket_count[_bucket[into_train]] += weight[into_train];
  train_weight = train_weight - weight[out_of_train] + weight[into_train];

  const uint64_t new_penalty = Penalty(train_bucket_count, train_weight);
  return static_cast<int64_t>(new_penalty) - static_cast<int64_t>(_penalty);
}

void
ActivityProfile::ApplySwap(uint32_t out_of_train, uint32_t into_train,
                           const std::vector<uint32_t>& weight) {
  _train_bucket_count[_bucket[out_of_train]] -= weight[out_of_train];
  _train_bucket_count[_bucket[into_train]] += weight[into_train];
  _train_weight = _train_weight - weight[out_of_train] + weight[into_train];
  _penalty = Penalty(_train_bucket_count, _train_weight);
}

int
ActivityProfileSet::InitialiseSplit(const std::vector<int>& in_train,
                                    const std::vector<uint32_t>& weight) {
  for (ActivityProfile& profile : _profile) {
    if (! profile.InitialiseSplit(in_train, weight)) {
      return 0;
    }
  }

  return 1;
}

uint64_t
ActivityProfileSet::RecomputePenalty(const std::vector<int>& in_train,
                                     const std::vector<uint32_t>& weight) {
  uint64_t result = 0;
  for (ActivityProfile& profile : _profile) {
    result += profile.RecomputePenalty(in_train, weight);
  }

  return result;
}

int64_t
ActivityProfileSet::DeltaForSwap(uint32_t out_of_train, uint32_t into_train,
                                 const std::vector<uint32_t>& weight) const {
  int64_t result = 0;
  for (const ActivityProfile& profile : _profile) {
    result += profile.DeltaForSwap(out_of_train, into_train, weight);
  }

  return result;
}

void
ActivityProfileSet::ApplySwap(uint32_t out_of_train, uint32_t into_train,
                              const std::vector<uint32_t>& weight) {
  for (ActivityProfile& profile : _profile) {
    profile.ApplySwap(out_of_train, into_train, weight);
  }
}

uint64_t
ActivityProfileSet::current_penalty() const {
  uint64_t result = 0;
  for (const ActivityProfile& profile : _profile) {
    result += profile.current_penalty();
  }

  return result;
}

}  // namespace train_test_split
