// Test code for the Matcher object.

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "matcher.h"

namespace {

TEST(TestMatcher, TestMatchAny) {
  iwmatcher::Matcher<int> matcher;
  ASSERT_TRUE(matcher.ok());
  EXPECT_TRUE(matcher.matches(0));
  EXPECT_TRUE(matcher.matches(1));
}

TEST(TestMatcher, TestMatchSingleValue) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(5);
  ASSERT_TRUE(matcher.ok());
  EXPECT_FALSE(matcher.matches(0));
  EXPECT_TRUE(matcher.matches(5));
  EXPECT_FALSE(matcher.matches(1));
  EXPECT_TRUE(matcher.matches(5));
}

TEST(TestMatcher, TestMatchMin) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(2);
  ASSERT_TRUE(matcher.ok());
  EXPECT_FALSE(matcher.matches(0));
  EXPECT_FALSE(matcher.matches(1));
  EXPECT_TRUE(matcher.matches(2));
  EXPECT_TRUE(matcher.matches(3));
}

TEST(TestMatcher, TestMatchMax) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(2);
  ASSERT_TRUE(matcher.ok());
  EXPECT_TRUE(matcher.matches(0));
  EXPECT_TRUE(matcher.matches(1));
  EXPECT_TRUE(matcher.matches(2));
  EXPECT_FALSE(matcher.matches(3));
}

TEST(TestMatcher, TestMatchBadBuildMinMax) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(1);
  matcher.set_min(2);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild1MaxViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(1);
  matcher.add(2);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild2MaxViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(1);
  matcher.add(0);
  EXPECT_TRUE(matcher.ok());
  matcher.add(2);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild1MinViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(2);
  matcher.add(1);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestMatchBadBuild2MinViolates) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(1);
  matcher.add(2);
  EXPECT_TRUE(matcher.ok());
  matcher.add(0);
  EXPECT_FALSE(matcher.ok());
}

TEST(TestMatcher, TestAddIfNotAlreadPresent) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  EXPECT_EQ(matcher.add_if_not_already_present(3), 0);
  matcher.add(4);
  EXPECT_EQ(matcher.add_if_not_already_present(4), 0);
  EXPECT_TRUE(matcher.ok());
}

TEST(TestMatcher, TestFromSeveralValues) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.add(4);
  matcher.add(5);

  EXPECT_FALSE(matcher.matches(1));
  EXPECT_FALSE(matcher.matches(2));
  EXPECT_TRUE(matcher.matches(3));
  EXPECT_TRUE(matcher.matches(4));
  EXPECT_TRUE(matcher.matches(5));
  EXPECT_FALSE(matcher.matches(6));
  EXPECT_FALSE(matcher.matches(7));
}

TEST(TestMatcher, TestMatchesMinOnly) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.set_min(1);
  ASSERT_TRUE(matcher.ok());
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatcher, TestMatchesMaxOnly) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.set_max(5);
  ASSERT_TRUE(matcher.ok());
  EXPECT_FALSE(matcher.matches(6));
  EXPECT_TRUE(matcher.matches(3));
}

TEST(TestMatchesMaxOnly, TestAdjustToAccommodateIndividualValue) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  EXPECT_FALSE(matcher.matches(2));
  EXPECT_TRUE(matcher.matches(3));
  matcher.adjust_to_accommodate(2);
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatchesMaxOnly, TestAdjustToAccommodateMin) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_min(3);
  EXPECT_FALSE(matcher.matches(2));
  matcher.adjust_to_accommodate(2);
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatchesMaxOnly, TestAdjustToAccommodateMax) {
  iwmatcher::Matcher<int> matcher;
  matcher.set_max(3);
  EXPECT_FALSE(matcher.matches(4));
  matcher.adjust_to_accommodate(2);
  EXPECT_TRUE(matcher.matches(2));
}

TEST(TestMatches, TestAnyValueFunctor) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.add(4);
  EXPECT_FALSE(matcher.AnyValue([] (const int v) {
    return v == 2;
  }));
  EXPECT_TRUE(matcher.AnyValue([] (const int v) {
    return v == 3;
  }));
}

TEST(TestMatcher, TestUPdateValues) {
  iwmatcher::Matcher<int> matcher;
  matcher.add(3);
  matcher.add(4);
  EXPECT_FALSE(matcher.matches(5));
  matcher.UpdateValues([](int v) {
    return v + 1;
  });
  EXPECT_FALSE(matcher.matches(3));
  EXPECT_TRUE(matcher.matches(4));
  EXPECT_TRUE(matcher.matches(5));
  EXPECT_FALSE(matcher.matches(6));
}

}  // namespace

TEST(MatcherOwnership, CopiesAreIndependent) {
  iwmatcher::Matcher<int> original;
  original.add(1);
  original.add(2);
  iwmatcher::Matcher<int> copy(original);
  EXPECT_TRUE(copy == original);
  original.reset();
  EXPECT_TRUE(copy.matches(2));
  EXPECT_FALSE(copy.matches(3));
  copy = copy;
  EXPECT_TRUE(copy.matches(2));
}

TEST(MatcherOwnership, AssignmentReplacesOldValues) {
  iwmatcher::Matcher<int> dest;
  dest.add(1);
  dest.add(2);
  iwmatcher::Matcher<int> source;
  source.add(9);
  dest = source;
  EXPECT_TRUE(dest == source);
  EXPECT_TRUE(dest.matches(9));
  EXPECT_FALSE(dest.matches(2));
  source.reset();
  dest = source;
  EXPECT_TRUE(dest.match_any());
  EXPECT_EQ(dest.number_elements(), 0);
}

TEST(MatcherOwnership, CopiesBoundsAndEmptyMatchers) {
  iwmatcher::Matcher<int> original;
  iwmatcher::Matcher<int> empty(original);
  EXPECT_TRUE(empty.match_any());
  original.set_min(2);
  original.set_max(5);
  iwmatcher::Matcher<int> bounded(original);
  EXPECT_TRUE(bounded == original);
  EXPECT_FALSE(bounded.matches(1));
  EXPECT_TRUE(bounded.matches(3));
  EXPECT_FALSE(bounded.matches(6));
}

TEST(MatcherAccessors, MinimumExplicitValue) {
  iwmatcher::Matcher<int> matcher;
  int value = 99;
  EXPECT_FALSE(matcher.min_explicit_value(value));
  matcher.add(7);
  matcher.add(0);
  ASSERT_TRUE(matcher.min_explicit_value(value));
  EXPECT_EQ(value, 0);
}

TEST(MatcherOwnership, AssignmentGrowsAndShrinksValueStorage) {
  iwmatcher::Matcher<int> dest;
  dest.add(1);
  dest.add(2);
  iwmatcher::Matcher<int> source;
  source.add(4);
  source.add(5);
  source.add(6);
  dest = source;
  EXPECT_TRUE(dest == source);
  source.resize(0);
  source.add(8);
  source.add(9);
  dest = source;
  EXPECT_TRUE(dest == source);
  EXPECT_FALSE(dest.matches(6));
  dest = dest;
  EXPECT_TRUE(dest.matches(9));
}
