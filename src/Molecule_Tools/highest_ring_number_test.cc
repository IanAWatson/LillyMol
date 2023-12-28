// Tests for highest ring number

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "highest_ring_number.h"

namespace lillymol {

namespace {

using testing::IsEmpty;

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

struct SmilesUnbalanced {
  IWString smiles;
  int expected_return_code;
  resizable_array<int> expected_rings;
};

class UnbalancedRingNumbersTester : public testing::TestWithParam<SmilesUnbalanced> {
  protected:
    resizable_array<int> ring_numbers;
};

TEST_P(UnbalancedRingNumbersTester, UnbalancedRingNumbersTester) {
  const auto& params = GetParam();
  int rc = UnbalancedRingNumbers(params.smiles, ring_numbers);
  EXPECT_EQ(rc, params.expected_return_code);

  if (rc == 0) {
    // If the result is a failure, there may be partial results accumulated, so do not check this.
    // EXPECT_THAT(ring_numbers, IsEmpty()) << params.smiles;
    return;
  }

  EXPECT_EQ(ring_numbers, params.expected_rings) << params.smiles;
}
INSTANTIATE_TEST_SUITE_P(UnbalancedRingNumbersTester, UnbalancedRingNumbersTester, testing::Values(
  SmilesUnbalanced{"", 0, {}},
  // Unbalanced ring numbers are preceded by % signs
  SmilesUnbalanced{"c1", 0, {}},
  // Failure if there are not 2 digits after the %
  SmilesUnbalanced{"c1%", 0, {}},
  SmilesUnbalanced{"c1% ", 0, {}},
  SmilesUnbalanced{"c1%q", 0, {}},
  SmilesUnbalanced{"c1%3", 0, {}},
  SmilesUnbalanced{"c1%3 ", 0, {}},
  SmilesUnbalanced{"c1%33qq%", 0, {}},
  SmilesUnbalanced{"c1%10%11", 2, {10, 11}},
  SmilesUnbalanced{"c1%11%10", 2, {11, 10}}
));


}  // namespace

}  // namespace lillymol
