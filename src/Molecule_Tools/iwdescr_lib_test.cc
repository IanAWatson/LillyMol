#include <cmath>
#include <cstdlib>
#include <string>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Tools/iwdescr_lib.h"

namespace {

TEST(IWDescr, InitialiseAllAndProcess) {
  if (getenv("LILLYMOL_HOME") == nullptr) {
    GTEST_SKIP() << "LILLYMOL_HOME is not defined";
  }

  IWDescr iwdescr;
  ASSERT_TRUE(iwdescr.InitialiseAll());
  ASSERT_GT(iwdescr.number_descriptors(), 0);

  std::unordered_set<std::string> names;
  for (int i = 0; i < iwdescr.number_descriptors(); ++i) {
    const std::string name = iwdescr.descriptor_name(i).AsString();
    EXPECT_FALSE(name.empty()) << i;
    EXPECT_TRUE(names.insert(name).second) << name;
  }

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCO"));
  std::vector<float> values(iwdescr.number_descriptors());
  ASSERT_TRUE(iwdescr.Process(m, values.data()));

  int finite_values = 0;
  for (float value : values) {
    finite_values += std::isfinite(value);
  }
  EXPECT_GT(finite_values, 0);
}

}  // namespace
