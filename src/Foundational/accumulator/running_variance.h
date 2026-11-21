#ifndef FOUNDATIONAL_ACCUMULATOR_RUNNING_VARIANCE_H_
#define FOUNDATIONAL_ACCUMULATOR_RUNNING_VARIANCE_H_

#include <algorithm>
#include <cstdint>

namespace iwaccumulator {

// Central moments avoid subtracting two nearly equal sums of squares. Counts
// belong to the enclosing accumulator; this state supplements its public sums.
class RunningVariance {
 private:
  double _mean = 0.0;
  double _m2 = 0.0;

 public:
  void Add(double value, uint64_t previous_count) {
    if (previous_count == 0) {
      _mean = value;
      _m2 = 0.0;
      return;
    }
    const double delta = value - _mean;
    _mean += delta / static_cast<double>(previous_count + 1);
    _m2 += delta * (value - _mean);
  }

  void AddCopies(double value, uint64_t copies, uint64_t previous_count) {
    if (previous_count == 0) {
      _mean = value;
      _m2 = 0.0;
      return;
    }
    const double total = static_cast<double>(previous_count) + copies;
    const double delta = value - _mean;
    _mean += delta * (static_cast<double>(copies) / total);
    _m2 += delta * delta * (static_cast<double>(previous_count) / total) * copies;
  }

  void Merge(const RunningVariance& rhs, uint64_t count, uint64_t rhs_count) {
    if (rhs_count == 0) return;
    if (count == 0) {
      *this = rhs;
      return;
    }
    const double total = static_cast<double>(count) + rhs_count;
    const double delta = rhs._mean - _mean;
    _mean += delta * (static_cast<double>(rhs_count) / total);
    _m2 += rhs._m2 + delta * delta * (static_cast<double>(count) / total) * rhs_count;
  }

  RunningVariance Without(const RunningVariance& rhs, uint64_t count,
                          uint64_t rhs_count) const {
    if (rhs_count == 0) return *this;
    RunningVariance result;
    const uint64_t remaining = count - rhs_count;
    if (remaining == 0) return result;
    const double delta = rhs._mean - _mean;
    result._mean = _mean - delta * (static_cast<double>(rhs_count) / remaining);
    result._m2 = std::max(_m2 - rhs._m2 - delta * delta *
                         (static_cast<double>(count) / remaining) * rhs_count, 0.0);
    return result;
  }

  double Variance(uint64_t count) const {
    return count > 1 ? std::max(_m2, 0.0) / static_cast<double>(count - 1) : 0.0;
  }
};

// Empty base optimisation keeps integer accumulators free of extra state.
template <typename Summer>
class VarianceState {};

}  // namespace iwaccumulator
#endif  // FOUNDATIONAL_ACCUMULATOR_RUNNING_VARIANCE_H_
