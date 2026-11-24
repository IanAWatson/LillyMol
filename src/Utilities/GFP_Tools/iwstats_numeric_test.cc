#include <array>
#include <cmath>

#include "gtest/gtest.h"

#include "Enrichment.h"
#include "Metric.h"
#include "bsquared.h"

namespace {
template <typename T>
class NumericInterface : public testing::Test {};
using NumericTypes = testing::Types<float, double>;
TYPED_TEST_SUITE(NumericInterface, NumericTypes);

TYPED_TEST(NumericInterface, PerfectAndReversedBSquared) {
  const std::array<TypeParam, 4> good{4, 3, 2, 1};
  const std::array<TypeParam, 4> bad{1, 2, 3, 4};
  std::array<TypeParam, 4> work;
  double result;
  ASSERT_TRUE(compute_b_squared(good.data(), work.data(), good.size(), result));
  EXPECT_NEAR(result, 1.0, 1e-12);
  ASSERT_TRUE(compute_b_squared(bad.data(), work.data(), bad.size(), result));
  EXPECT_NEAR(result, -1.0, 1e-12);
}

TYPED_TEST(NumericInterface, CumulativeBSquaredOverload) {
  std::array<TypeParam, 4> values{4, 3, 2, 1};
  std::array<TypeParam, 4> work;
  std::array<double, 4> result;
  ASSERT_TRUE(compute_b_squared(values.data(), work.data(), values.size(), result.data()));
  for (double value : result) EXPECT_NEAR(value, 1.0, 1e-12);
}

TYPED_TEST(NumericInterface, EnrichmentAcceptsBothElementTypes) {
  std::array<TypeParam, 4> activities{4, 1, 3, 0};
  Enrichment enrichment;
  ASSERT_TRUE(enrichment.buildFromSortedArray(activities.data(), activities.size(),
                                             static_cast<TypeParam>(2)));
  EXPECT_EQ(enrichment.getNbrActives(), 2u);
  EXPECT_EQ(enrichment.getNbrCompounds(), 4u);
  EXPECT_NEAR(calculateROC()(enrichment), 0.75, 1e-12);
}

TEST(DoubleNumericInterface, EnrichmentRetainsSmallDifferences) {
  const double activities[]{1.00000001, 1.00000002, 1.00000003};
  Enrichment enrichment;
  ASSERT_TRUE(enrichment.buildFromSortedArray(activities, 3, 1.000000015));
  EXPECT_EQ(enrichment.getNbrActives(), 2u);
}
}  // namespace
