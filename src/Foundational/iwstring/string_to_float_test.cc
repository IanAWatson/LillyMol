// Tests for functions in string_to_float

#include <cmath>

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"

#include "iwstring.h"

namespace {

TEST(TestDragonBox, TestAll) {
  IWString foo;
  float x = 3.14;
  foo.append_number_dragonbox(x);
  EXPECT_EQ(foo, "3.14e0") << "got '" << foo << "'\n";

  x = std::sqrt(x);

  foo.resize_keep_storage(0);
  foo.append_number_dragonbox(x);
  EXPECT_EQ(foo, "3.14e0") << "got '" << foo << "'\n";
}

}  // namespace
