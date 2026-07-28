#ifndef UTILITIES_GFP_TOOLS_TRAIN_TEST_SPLIT_ACTIVITY_H_
#define UTILITIES_GFP_TOOLS_TRAIN_TEST_SPLIT_ACTIVITY_H_

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "Foundational/iwstring/iwstring.h"

namespace train_test_split {

enum class ActivityBucketisation {
  kQuantile,
  kEqualWidth,
};

struct ActivityProfileOptions {
  IWString fname;
  uint32_t buckets = 10;
  ActivityBucketisation bucketisation = ActivityBucketisation::kQuantile;
};

class ActivityProfile {
  private:
    ActivityProfileOptions _options;
    std::vector<uint8_t> _scaled;
    std::vector<uint32_t> _bucket;
    std::vector<uint64_t> _total_bucket_count;
    std::vector<uint64_t> _train_bucket_count;
    uint64_t _total_weight = 0;
    uint64_t _train_weight = 0;
    uint64_t _penalty = 0;

    int AssignBuckets(const std::vector<double>& activity);
    int AssignQuantileBuckets(const std::vector<double>& activity);
    int AssignEqualWidthBuckets();

    uint64_t Penalty(const std::vector<uint64_t>& train_bucket_count,
                     uint64_t train_weight) const;

  public:
    int Build(const ActivityProfileOptions& options,
              const std::vector<std::string>& id,
              const std::vector<uint32_t>& weight);

    int InitialiseSplit(const std::vector<int>& in_train,
                        const std::vector<uint32_t>& weight);
    uint64_t RecomputePenalty(const std::vector<int>& in_train,
                              const std::vector<uint32_t>& weight);
    int64_t DeltaForSwap(uint32_t out_of_train, uint32_t into_train,
                         const std::vector<uint32_t>& weight) const;
    void ApplySwap(uint32_t out_of_train, uint32_t into_train,
                   const std::vector<uint32_t>& weight);

    uint32_t bucket(uint32_t i) const {
      return _bucket[i];
    }
    uint8_t scaled(uint32_t i) const {
      return _scaled[i];
    }
    uint64_t current_penalty() const {
      return _penalty;
    }
    uint32_t nbuckets() const {
      return _options.buckets;
    }
    const IWString& fname() const {
      return _options.fname;
    }
};

class ActivityProfileSet {
  private:
    std::vector<ActivityProfileOptions> _pending;
    std::vector<ActivityProfile> _profile;
    ActivityProfileOptions _defaults;

  public:
    int AddDirective(const const_IWSubstring& directive);
    int Build(const std::vector<std::string>& id,
              const std::vector<uint32_t>& weight);

    int InitialiseSplit(const std::vector<int>& in_train,
                        const std::vector<uint32_t>& weight);
    uint64_t RecomputePenalty(const std::vector<int>& in_train,
                              const std::vector<uint32_t>& weight);
    int64_t DeltaForSwap(uint32_t out_of_train, uint32_t into_train,
                         const std::vector<uint32_t>& weight) const;
    void ApplySwap(uint32_t out_of_train, uint32_t into_train,
                   const std::vector<uint32_t>& weight);

    bool active() const {
      return !_profile.empty();
    }
    bool has_pending() const {
      return !_pending.empty();
    }
    int number_profiles() const {
      return _profile.size();
    }
    uint64_t current_penalty() const;
};

}  // namespace train_test_split

#endif  // UTILITIES_GFP_TOOLS_TRAIN_TEST_SPLIT_ACTIVITY_H_
