#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "gtest/gtest.h"

#include "accumulator.h"
#include "accumulator_v2.h"

namespace {

template <typename A, typename T>
struct Case {
  using Acc = A;
  using Value = T;
};

template <typename T>
class AccumulatorTest : public testing::Test {};

using AccumulatorTypes = testing::Types<
    Case<Accumulator<float>, float>, Case<Accumulator<double>, double>,
    Case<Accumulator_V2<unsigned int, float>, float>,
    Case<Accumulator_V2<unsigned int, double>, double>>;
TYPED_TEST_SUITE(AccumulatorTest, AccumulatorTypes);

TYPED_TEST(AccumulatorTest, IdenticalValuesHaveExactlyZeroVariance) {
  using T = typename TypeParam::Value;
  for (T value : {T(0.1), T(7.1), T(1.0e20)}) {
    typename TypeParam::Acc single, copies, array;
    for (int i = 0; i < 10; ++i) single.extra(value);
    copies.extra(value, 10);
    std::vector<T> values(10, value);
    array.extra(values.data(), values.size());
    for (const auto* acc : {&single, &copies, &array}) {
      EXPECT_EQ(acc->n(), 10u);
      EXPECT_DOUBLE_EQ(acc->variance(), 0.0);
      EXPECT_DOUBLE_EQ(acc->average(), static_cast<double>(value));
      EXPECT_EQ(acc->minval(), value);
      EXPECT_EQ(acc->maxval(), value);
    }
  }
}

TYPED_TEST(AccumulatorTest, InvalidDataDoesNotBecomeZeroVariance) {
  using T = typename TypeParam::Value;
  typename TypeParam::Acc acc;
  acc.extra(std::numeric_limits<T>::quiet_NaN());
  acc.extra(std::numeric_limits<T>::quiet_NaN());
  EXPECT_TRUE(std::isnan(acc.variance()));
}

TYPED_TEST(AccumulatorTest, SmallSpreadAroundLargeOffset) {
  using T = typename TypeParam::Value;
  const T first = T(10000.1);
  const T second = T(10000.2);
  const double difference = static_cast<double>(second) - first;
  typename TypeParam::Acc acc;
  acc.extra(first);
  acc.extra(second);
  EXPECT_NEAR(acc.variance(), difference * difference / 2.0, 1.0e-10);
}

TYPED_TEST(AccumulatorTest, ArrayOfOneAndEmptyInputs) {
  using T = typename TypeParam::Value;
  typename TypeParam::Acc acc;
  EXPECT_EQ(acc.extra(static_cast<const T*>(nullptr), 0), 0u);
  EXPECT_EQ(acc.extra(T(-10), 0), 0u);
  const std::array<T, 1> value{T(7.1)};
  EXPECT_EQ(acc.extra(value.data(), 1), 1u);
  EXPECT_EQ(acc.minval(), value[0]);
  EXPECT_EQ(acc.maxval(), value[0]);
  EXPECT_EQ(acc.extra(T(-100), 0), 1u);
  EXPECT_EQ(acc.minval(), value[0]);
  EXPECT_EQ(acc.maxval(), value[0]);
}

TYPED_TEST(AccumulatorTest, CopiesMergeAssignmentAndSelfMerge) {
  using T = typename TypeParam::Value;
  typename TypeParam::Acc a, b, sequential;
  a.extra(T(10000), 3);
  b.extra(T(10004), 2);
  for (T value : {T(10000), T(10000), T(10000), T(10004), T(10004)}) {
    sequential.extra(value);
  }
  typename TypeParam::Acc merged;
  merged.extra(a);
  merged += b;
  EXPECT_EQ(merged.n(), sequential.n());
  EXPECT_DOUBLE_EQ(merged.average(), sequential.average());
  EXPECT_NEAR(merged.variance(), 4.8, 1.0e-10);
  EXPECT_NEAR(merged.variance(), sequential.variance(), 1.0e-10);
  typename TypeParam::Acc copy(merged), assigned;
  assigned = merged;
  for (auto* acc : {&copy, &assigned}) {
    acc->extra(T(10002));
    EXPECT_EQ(acc->n(), 6u);
    EXPECT_NEAR(acc->variance(), 3.8666666666666667, 1.0e-10);
    acc->reset();
    EXPECT_EQ(acc->n(), 0u);
    acc->extra(T(2));
    acc->extra(T(4));
    EXPECT_DOUBLE_EQ(acc->variance(), 2.0);
  }
  a.extra(a);
  EXPECT_EQ(a.n(), 6u);
  EXPECT_DOUBLE_EQ(a.variance(), 0.0);
}

TYPED_TEST(AccumulatorTest, SubtractionAndAliasingPreserveCentralMoments) {
  using T = typename TypeParam::Value;
  typename TypeParam::Acc all, removed, remaining;
  for (T value : {T(1), T(3), T(7), T(9)}) all.extra(value);
  removed.extra(T(7));
  removed.extra(T(9));
  ASSERT_TRUE(all.subtract_data(removed, remaining));
  EXPECT_EQ(remaining.n(), 2u);
  EXPECT_DOUBLE_EQ(remaining.average(), 2.0);
  EXPECT_DOUBLE_EQ(remaining.variance(), 2.0);
  auto source_alias = all;
  ASSERT_TRUE(source_alias.subtract_data(removed, source_alias));
  EXPECT_DOUBLE_EQ(source_alias.variance(), 2.0);
  auto removed_alias = removed;
  ASSERT_TRUE(all.subtract_data(removed_alias, removed_alias));
  EXPECT_DOUBLE_EQ(removed_alias.variance(), 2.0);
  ASSERT_TRUE(all.subtract_data(all, remaining));
  EXPECT_EQ(remaining.n(), 0u);
  remaining.extra(T(3));
  remaining.extra(T(5));
  EXPECT_DOUBLE_EQ(remaining.variance(), 2.0);
}

TEST(AccumulatorPrecision, DoublePrecisionLargeOffsetVariance) {
  Accumulator<double> acc;
  acc.extra(1.0e12 + 1.0);
  acc.extra(1.0e12 + 2.0);
  acc.extra(1.0e12 + 3.0);
  EXPECT_DOUBLE_EQ(acc.variance(), 1.0);
}

TEST(AccumulatorPrecision, FloatSquaresAreComputedInDoublePrecision) {
  const float value = 1.0e20f;
  const double square = static_cast<double>(value) * value;
  for (int mode : {0, 1, 2}) {
    Accumulator<float> acc;
    if (mode == 0) {
      acc.extra(value);
      acc.extra(value);
    } else if (mode == 1) {
      acc.extra(value, 2);
    } else {
      const std::array<float, 2> values{value, value};
      acc.extra(values.data(), values.size());
    }
    accumulator::AccumulatorData proto;
    ASSERT_TRUE(acc.BuildProto(proto));
    EXPECT_TRUE(std::isfinite(proto.x2sum()));
    EXPECT_DOUBLE_EQ(proto.x2sum(), square * 2.0);
    EXPECT_DOUBLE_EQ(acc.variance(), 0.0);
  }
}

TEST(AccumulatorPrecision, V2IntegerValuesAreSquaredAfterPromotion) {
  Accumulator_V2<unsigned int, int> acc;
  acc.extra(50000);
  acc.extra(50001);
  EXPECT_DOUBLE_EQ(acc.sum_of_squares(), 50000.0 * 50000.0 + 50001.0 * 50001.0);
  EXPECT_DOUBLE_EQ(acc.variance(), 0.5);
}

TEST(AccumulatorInteger, ArrayAndZeroCopiesKeepExtremaCorrect) {
  Accumulator_Int<int> acc;
  acc.extra(100, 0);
  const std::array<int, 1> first{7};
  acc.extra(first.data(), 1);
  EXPECT_EQ(acc.minval(), 7);
  EXPECT_EQ(acc.maxval(), 7);
  acc.extra(9);
  EXPECT_EQ(acc.sum(), 16);
  EXPECT_EQ(acc.sum_of_squares(), 130);
  EXPECT_DOUBLE_EQ(acc.variance(), 2.0);
}

}  // namespace
