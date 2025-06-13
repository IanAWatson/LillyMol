// Tester for the string classes
#include <stdlib.h>

#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_set>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "iwstring.h"

namespace {

TEST(TestIWString, TestAsString) {
  IWString s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}
TEST(TestConstIWSubstring, TestAsString) {
  const_IWSubstring s("hello");
  const std::string as_string = s.AsString();
  EXPECT_EQ(as_string, "hello");
}
TEST(TestAppendOperator, TestStd) {
  IWString s("hello");
  std::string_view space(" ");
  std::string world("world");
  s << space << world;
  EXPECT_EQ(s, "hello world");
}
TEST(TestEqualsOperator, TestStd) {
  IWString hello = "hello";
  std::string world;
  EXPECT_FALSE(iwstring::Equals(hello, world));
}

TEST(TestIWString, TestExpandEnvironmentVariablesNothing) {
  const IWString hello("hello world");
  EXPECT_EQ(hello.ExpandEnvironmentVariables(), hello);
}
TEST(TestIWString, TestExpandEnvironmentVariablesTooShort) {
  const IWString hello("${}");
  EXPECT_EQ(hello.ExpandEnvironmentVariables(), hello);
}
TEST(TestIWString, TestExpandEnvironmentVariablesNotSet) {
  const IWString hello("${NOTSET}");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}
TEST(TestIWString, TestExpandEnvironmentVariablesEmpty) {
  const IWString hello("xx${}foobar");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}
TEST(TestIWString, TestExpandEnvironmentVariablesNotClosed) {
  const IWString hello("xx${xyfoobar");
  EXPECT_FALSE(hello.ExpandEnvironmentVariables());
}

struct EnvData {
  // Pairs of SHELL_VAR=>value
  std::unordered_map<std::string, std::string> to_set;
  // something like "hello ${world}"
  IWString input;
  // something like "hello world"
  IWString expected;
};

class TestExpandEnv: public testing::TestWithParam<EnvData> {
  protected:
    // Set of shell variables to be unset
    std::unordered_set<std::string> _to_clear;

    void TearDown();
};

void
TestExpandEnv::TearDown() {
  for (const auto& vname : _to_clear) {
    ::unsetenv(vname.c_str());
  }
}

TEST_P(TestExpandEnv, TestExpandEnv) {
  const auto params = GetParam();
  for (const auto& [vname, value] : params.to_set) {
    _to_clear.emplace(vname);
    ::setenv(vname.c_str(), value.c_str(), 1);
  }

  std::optional<IWString> expanded = params.input.ExpandEnvironmentVariables();
  EXPECT_EQ(expanded, params.expected) << " expected " << params.expected <<
        " got '" << expanded->AsString();
}
INSTANTIATE_TEST_SUITE_P(TestExpandEnv, TestExpandEnv, testing::Values(
  EnvData{{}, "hello", "hello"},
  EnvData{{{"world", "world"}}, "hello ${world}", "hello world"},
  EnvData{{{"hello", "hello"}}, "${hello} world", "hello world"},
  EnvData{{{"hello", "welcome"}}, "${hello} world", "welcome world"},
  EnvData{{{"a", "abcdefghi"}}, "${a} world", "abcdefghi world"},
  EnvData{{{"a", "abcdefghi"}}, "xxx ${a} world", "xxx abcdefghi world"},
  EnvData{{{"abcdefghi", "a"}}, "${abcdefghi} world", "a world"},
  EnvData{{{"abcdefghi", "a"}}, "xxx ${abcdefghi} y", "xxx a y"},
  EnvData{{{"mm93", "marc marquez"}, 
           {"vr46", "valentino rossi"}},
           "motogp ${mm93} and ${vr46} greats",
           "motogp marc marquez and valentino rossi greats"},
  EnvData{{{"mm93", "marc marquez"}}, "hello $mm93 motogp", "hello $mm93 motogp"},
  EnvData{{{"mm93", "marquez marquez"}}, "hello $mm93}", "hello $mm93}"}
));

TEST(TestIWString, TestEnsureEndsWithEmpty) {
  IWString s;
  EXPECT_TRUE(s.EnsureEndsWith('a'));
  EXPECT_EQ(s, 'a');
}

template <typename T>
class EnsureEndsWithTest : public testing::Test {
 protected:
  void TestEmpty();
  void TestAlreadyEndsWith();
  void TestMustBeAdded();
};

template <typename T>
void
EnsureEndsWithTest<T>::TestEmpty() {
  IWString s;
  T extra("a");
  EXPECT_EQ(s.EnsureEndsWith(extra), 1);
  EXPECT_EQ(s, "a");
}

template <typename T>
void
EnsureEndsWithTest<T>::TestAlreadyEndsWith() {
  IWString s('a');
  T extra("a");
  EXPECT_EQ(s.EnsureEndsWith(extra), 0);
  EXPECT_EQ(s, 'a');
}

template <typename T>
void
EnsureEndsWithTest<T>::TestMustBeAdded() {
  IWString s('a');
  T extra("b");
  EXPECT_EQ(s.EnsureEndsWith(extra), 1);
  EXPECT_EQ(s, "ab");
}

using MyTypes = ::testing::Types<const char*, const IWString&, const const_IWSubstring&>;
TYPED_TEST_SUITE_P(EnsureEndsWithTest);

TYPED_TEST_P(EnsureEndsWithTest, StartsEmpty) {
  // Inside a test, refer to TypeParam to get the type parameter.
  // TypeParam n = 0;
  // Maybe something could be done to combine the const char* type into the template?
  // std::is_pointer...
  // if (std::is_integral<TypeParam>::value) {
  // }

  // You will need to use `this` explicitly to refer to fixture members.
  this->TestEmpty();
}

TYPED_TEST_P(EnsureEndsWithTest, AlreadyEndsWith) {
  // TypeParam n = 0;

  this->TestAlreadyEndsWith();
}

TYPED_TEST_P(EnsureEndsWithTest, MustBeAdded) { 
  // TypeParam n = 0;

  this->TestMustBeAdded();
}

REGISTER_TYPED_TEST_SUITE_P(EnsureEndsWithTest,
                            StartsEmpty, AlreadyEndsWith, MustBeAdded);

INSTANTIATE_TYPED_TEST_SUITE_P(My, EnsureEndsWithTest, MyTypes);

// typed test suites are too complicated...

TEST(TestIWString, TestStartsWithChar) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with('c'));
  EXPECT_FALSE(foo.starts_with('b'));
  EXPECT_FALSE(foo.starts_with("cc"));
}

TEST(TestConst_IWSubstring, TestStartsWithChar) {
  const char* s = "c";
  const_IWSubstring foo(s);
  EXPECT_TRUE(foo.starts_with('c'));
  EXPECT_FALSE(foo.starts_with('b'));
  EXPECT_FALSE(foo.starts_with("cc"));
}

TEST(TestIWString, TestStartsWithCharStar) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with("c"));
  EXPECT_TRUE(foo.starts_with("c", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_FALSE(foo.starts_with("cx"));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with("a"));
  EXPECT_TRUE(foo.starts_with("ab"));
  EXPECT_TRUE(foo.starts_with("abc"));
  EXPECT_TRUE(foo.starts_with("abc", 1));
  EXPECT_TRUE(foo.starts_with("abc", 2));
  EXPECT_TRUE(foo.starts_with("abc", 3));
  EXPECT_FALSE(foo.starts_with("abcd", 4));
  EXPECT_FALSE(foo.starts_with("x"));
}

TEST(TestConst_IWSubstring, TestStartsWithCharStar) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with("c"));
  EXPECT_TRUE(foo.starts_with("c", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_TRUE(foo.starts_with("cx", 1));
  EXPECT_FALSE(foo.starts_with("x"));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with("a"));
  EXPECT_TRUE(foo.starts_with("ab"));
  EXPECT_TRUE(foo.starts_with("abc"));
  EXPECT_FALSE(foo.starts_with("abcd"));
  EXPECT_FALSE(foo.starts_with("x"));
}

TEST(TestIWString, TestStartsWithIWString) {
  IWString foo('c');
  EXPECT_TRUE(foo.starts_with(foo));
  IWString bar("cx");
  EXPECT_TRUE(bar.starts_with(bar));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with(foo));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "a";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "ab";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "abc";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
  bar = "abcd";
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
}

TEST(TestConst_IWSubstring, TestStartsWithIWString) {
  const_IWSubstring foo("c");
  EXPECT_TRUE(foo.starts_with(foo));
  const_IWSubstring bar("cx");
  EXPECT_TRUE(bar.starts_with(bar));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));

  foo = "abc";
  EXPECT_TRUE(foo.starts_with(foo));
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "a";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "ab";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_FALSE(bar.starts_with(foo));
  bar = "abc";
  EXPECT_TRUE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
  bar = "abcd";
  EXPECT_FALSE(foo.starts_with(bar));
  EXPECT_TRUE(bar.starts_with(foo));
}

TEST(TextNextWord, TestSpaces) {
  const_IWSubstring hello_world("hello  world foo");
  int i = 0;
  const_IWSubstring token;
  ASSERT_TRUE(hello_world.nextword(token, i));
  EXPECT_EQ(token, "hello");
  ASSERT_TRUE(hello_world.nextword(token, i));
  EXPECT_EQ(token, "world");
  ASSERT_TRUE(hello_world.nextword(token, i));
  EXPECT_EQ(token, "foo");
}

TEST(TestStrcmp, TestAll) {
  IWString foo("hello");
  IWString bar("hello");
  
  EXPECT_EQ(foo.strcmp(bar), 0) << "hello.strcmp(hello)";

  EXPECT_EQ(bar.strncmp(foo, 4), 0) << "hello.strncmp(hello, 4";

  bar = "HELLO";

  EXPECT_EQ(foo.strcasecmp(bar), 0) << "HELLO.strcasecmp(hello)";

  EXPECT_EQ(bar.strcasecmp(foo), 0) << "HELLO.strcasecmp(hello)";
}

TEST(TestStrncat, TestAll) {
  IWString foo("hello");

  foo.strncat(" world", 6);

  EXPECT_EQ(foo, "hello world") << "hello.strncat( world, 6)";
}

TEST(TestTestSplitIntoDirectiveAndValue, TestAll) {
  IWString foo("mm=93");

  const_IWSubstring directive;
  int value;
  ASSERT_TRUE(foo.split_into_directive_and_value(directive, '=', value));

  EXPECT_EQ(directive, "mm");
  EXPECT_EQ(value, 93);
}

TEST(TestAppendNumber, TestAll) {
  IWString foo;
  foo.append_number(0);
  EXPECT_EQ(foo, "0") << "append number 0";

  foo.append_number(-8);

  EXPECT_EQ(foo, "0-8") << "append number -8";

  foo.append_number(123456789);

  EXPECT_EQ(foo, "0-8123456789") << "append number 123456789";

  foo = "";

  foo.append_number(-9);
  EXPECT_EQ(foo, "-9") << "append -9";

  foo.append_number(1000);

  EXPECT_EQ(foo, "-91000") << "append 1000";
}


TEST(TestIndex, TestAll) {
  IWString foo("hello");

  auto i = foo.rindex(' ');
  EXPECT_EQ(i, -1) << foo << " rindex space";

  foo = " hello";

  i = foo.rindex(' ');
  EXPECT_EQ(i, 0) << foo << " rindex";

  const_IWSubstring xx(foo);
  i = xx.rindex(' ');

  EXPECT_EQ(i, 0) << xx << " rindex";

  foo = "hello";
  i = foo.rindex('l');
  EXPECT_EQ(i, 3) << foo << " rindex l";

  i = foo.index('l');

  EXPECT_EQ(i, 2) << foo << " index l";

  i = foo.index('h');
  EXPECT_EQ(i, 0) << foo << " index h";
}


TEST(TestStartsWith, TestAll) {
  IWString foo;
  foo = "abcdef";

  EXPECT_TRUE(foo.starts_with ("abc")) << foo << " starts with abc";

  IWString bar("abcd");
  EXPECT_TRUE(foo.starts_with(bar)) << bar << " starts with abcd";

  const_IWSubstring bb = substr(bar, 0, 2);
  EXPECT_TRUE(foo.starts_with(bb));
}

TEST(TestCharacterEquality, TestAll) {
  IWString foo;

  foo = "f";

  EXPECT_EQ(foo, 'f');
  EXPECT_EQ('f', foo);

  EXPECT_TRUE(foo == 'f');
  EXPECT_TRUE('f' == foo);
  EXPECT_TRUE(foo == foo);

  EXPECT_FALSE(foo != 'f');
  EXPECT_FALSE('f' != foo);
  
  const_IWSubstring bar = "b";
  EXPECT_TRUE(bar == 'b') << bar << " == b";
  EXPECT_TRUE('b' == bar) << " b == " << bar;

  EXPECT_FALSE(bar != 'b');
  EXPECT_FALSE('b' != bar);
}

TEST(TestChop, TestAll) {
  IWString foo ("a;sldkfjas;dlfkj");

  foo.chop();
  EXPECT_EQ(foo, "a;sldkfjas;dlfk");

  foo.chop(2);
  EXPECT_EQ(foo, "a;sldkfjas;dl");
}

}  // namespace
