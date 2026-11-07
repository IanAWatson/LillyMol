#include "cmdline_v2.h"

#include <cmath>
#include <initializer_list>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace {

class Arguments {
 private:
  std::vector<std::string> _strings;
  std::vector<char*> _argv;

 public:
  Arguments(std::initializer_list<const char*> args) : _strings(args.begin(), args.end()) {
    for (std::string& arg : _strings) {
      _argv.push_back(arg.data());
    }
    _argv.push_back(nullptr);
  }

  int argc() const { return static_cast<int>(_strings.size()); }
  char** argv() { return _argv.data(); }
};

template <typename T>
class FiniteFloatingPointValues : public testing::Test {};

using FloatingPointTypes = testing::Types<float, double>;
TYPED_TEST_SUITE(FiniteFloatingPointValues, FloatingPointTypes);

TYPED_TEST(FiniteFloatingPointValues, RejectsNonFiniteAndOverflowWithoutChangingDestination) {
  for (const char* text : {"nan", "NaN", "NAN", "nan(payload)", "-nan",
                           "inf", "+inf", "-inf", "INF", "infinity",
                           "-infinity", "1e9999", "-1e9999"}) {
    SCOPED_TRACE(text);
    Arguments args{"tool", "-number", text};
    Command_Line_v2 cl(args.argc(), args.argv(), "-number=s");
    ASSERT_TRUE(cl.good());
    TypeParam value = 17;
    EXPECT_FALSE(cl.value("number", value));
    EXPECT_EQ(value, static_cast<TypeParam>(17));
    IWString raw;
    ASSERT_TRUE(cl.value("number", raw));
    EXPECT_EQ(raw, text);
  }
}

TYPED_TEST(FiniteFloatingPointValues, AcceptsFiniteValuesIncludingSignedZero) {
  Arguments args{"tool", "-n", "0", "-n", "-0", "-n", "-1.25", "-n", "2.5e2",
                  "-n", "1e-30", "-n", "1e30"};
  Command_Line_v2 cl(args.argc(), args.argv(), "-n=f");
  ASSERT_TRUE(cl.good());
  const TypeParam expected[] = {0, -0.0, -1.25, 250, static_cast<TypeParam>(1e-30),
                                static_cast<TypeParam>(1e30)};
  for (int i = 0; i < 6; ++i) {
    TypeParam value = 17;
    ASSERT_TRUE(cl.value('n', value, i));
    EXPECT_EQ(value, expected[i]);
    if (i == 1) {
      EXPECT_TRUE(std::signbit(value));
    }
  }
}

TEST(CommandLineV2, NumericQualifiersRejectNonFiniteValuesAtConstruction) {
  for (const char* spec : {"-n=f", "-n=float", "-n=fraction"}) {
    for (const char* text : {"nan", "NaN", "inf", "-inf", "infinity", "1e9999"}) {
      SCOPED_TRACE(spec);
      SCOPED_TRACE(text);
      Arguments args{"tool", "-n", text};
      Command_Line_v2 cl(args.argc(), args.argv(), spec);
      EXPECT_FALSE(cl.good());
    }
  }
}

TEST(CommandLineV2, FractionBoundariesStillWork) {
  Arguments args{"tool", "-n", "0", "-n", "0.5", "-n", "1"};
  Command_Line_v2 cl(args.argc(), args.argv(), "-n=fraction");
  ASSERT_TRUE(cl.good());
  double value = -1;
  ASSERT_TRUE(cl.value('n', value));
  EXPECT_DOUBLE_EQ(value, 0);
  ASSERT_TRUE(cl.value('n', value, 1));
  EXPECT_DOUBLE_EQ(value, 0.5);
  ASSERT_TRUE(cl.value('n', value, 2));
  EXPECT_DOUBLE_EQ(value, 1);
}

TEST(CommandLineV2, FloatingPointRangeDependsOnDestinationType) {
  Arguments args{"tool", "-n", "1e100"};
  Command_Line_v2 cl(args.argc(), args.argv(), "-n=f");
  ASSERT_TRUE(cl.good());
  double d = 0;
  ASSERT_TRUE(cl.value('n', d));
  EXPECT_DOUBLE_EQ(d, 1e100);
  float f = 17;
  EXPECT_FALSE(cl.value('n', f));
  EXPECT_FLOAT_EQ(f, 17);
}

}  // namespace
