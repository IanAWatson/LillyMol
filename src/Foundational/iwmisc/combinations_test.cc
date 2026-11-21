// Tester for the combinations class.

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "combinations.h"

namespace {

using combinations::Combinations;

using testing::ElementsAreArray;

TEST(Combinations, Test1NoVariability) {
  std::vector<uint32_t> count({1});
  Combinations<uint32_t> combinations(count);
  std::vector<uint32_t> state;
  state.resize(1, 0);
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0}));
  EXPECT_FALSE(combinations.Next(state));
}

TEST(Combinations, Test1WithVariability) {
  std::vector<uint32_t> count({2});
  Combinations<uint32_t> combinations(count);
  std::vector<uint32_t> state;
  state.resize(1, 0);
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0}));
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({1}));
  EXPECT_FALSE(combinations.Next(state));
}

TEST(Combinations, Test2NoVariability) {
  std::vector<uint32_t> count({1, 1});
  Combinations<uint32_t> combinations(count);
  std::vector<uint32_t> state;
  state.resize(2, 0);
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0, 0}));
  EXPECT_FALSE(combinations.Next(state));
}

TEST(Combinations, Test2VaryFirst) {
  std::vector<uint32_t> count({2, 1});
  Combinations<uint32_t> combinations(count);
  std::vector<uint32_t> state;
  state.resize(2, 0);
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0, 0}));
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({1, 0}));
  EXPECT_FALSE(combinations.Next(state));
}

TEST(Combinations, Test2VarySecond) {
  std::vector<uint32_t> count({1, 2});
  Combinations<uint32_t> combinations(count);
  std::vector<uint32_t> state;
  state.resize(2, 0);
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0, 0}));
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0, 1}));
  EXPECT_FALSE(combinations.Next(state));
}

TEST(Combinations, Test2VaryBoth) {
  std::vector<uint32_t> count({2, 2});
  Combinations<uint32_t> combinations(count);
  std::vector<uint32_t> state;
  state.resize(2, 0);
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0, 0}));
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0, 1}));
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({1, 0}));
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({1, 1}));
  EXPECT_FALSE(combinations.Next(state));
}

}  // namespace

TEST(Combinations, EmptyChoiceSetHasNoCombinations) {
  Combinations<uint32_t> combinations(std::vector<uint32_t>{2, 0, 3});
  std::vector<uint32_t> state;
  EXPECT_FALSE(combinations.Next(state));
  EXPECT_FALSE(combinations.Next(state));
  combinations.Reset();
  EXPECT_FALSE(combinations.Next(state));
}

TEST(Combinations, SizesStateAndRemainsExhaustedUntilReset) {
  Combinations<uint32_t> combinations(std::vector<uint32_t>{2});
  std::vector<uint32_t> state;
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0}));
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({1}));
  EXPECT_FALSE(combinations.Next(state));
  EXPECT_FALSE(combinations.Next(state));
  combinations.Reset();
  ASSERT_TRUE(combinations.Next(state));
  EXPECT_THAT(state, ElementsAreArray({0}));
}
