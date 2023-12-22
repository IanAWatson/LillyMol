// Tests for highest ring number

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "highest_ring_number.h"

namespace lillymol {
namespace {
struct SmilesResult {
  IWString smiles;
  // the method being tested returns std::optional.
  int valid_result_expected;
  int expected;
};

class HighestRingTester : public testing::TestWithParam<SmilesResult> {
};

TEST_P(HighestRingTester, HighestRingTester) {
  const auto& params = GetParam();
  if (params.valid_result_expected) {
    EXPECT_THAT(HighestRingNumber(params.smiles), testing::Optional(params.expected));
  } else {
    EXPECT_EQ(HighestRingNumber(params.smiles), std::nullopt) << params.smiles <<
    " expected " << params.expected;
  }
}

INSTANTIATE_TEST_SUITE_P(HighestRingTester, HighestRingTester, testing::Values(
  SmilesResult{"", 1, 0},
  SmilesResult{"C", 1, 0},
  SmilesResult{"CC", 1, 0},
  SmilesResult{"C1CC1", 1, 1},
  SmilesResult{"C9CC9", 1, 9},
  SmilesResult{"C%10CC%10", 0, 0}
));


}  // namespace

}  // namespace lillymol
