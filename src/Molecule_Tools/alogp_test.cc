// Tests for alogp

#include <filesystem>
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/alogp.h"

namespace {

class TestALogP : public testing::Test {
  protected:
    IWString _smiles;

    Molecule _m;

  // protected functions.
    void SetUp() override;

    Charge_Assigner _charge_assigner;

    alogp::ALogP _alogp;

  public:
};

void
TestALogP::SetUp() {
  const char* test_srcdir = getenv("TEST_SRCDIR");
  // Some diagnostic stuff to see the relationship between TEST_SRCDIR
  // and where our files actually show up.
#ifdef DEBUG_FILE_PATHS
  std::string qq(test_srcdir);
  qq += "/../alogp_test.runfiles/charge_assigner/";
  for (auto const& dir_entry : std::filesystem::directory_iterator{qq}) {
    std::cerr << dir_entry << '\n';
  }
#endif

  IWString queries_file(test_srcdir);
  queries_file << "/../alogp_test.runfiles/charge_assigner/queries";

  IWString cmd;
  cmd << "F:" << queries_file;

  if (! _charge_assigner.build(cmd)) {
    std::cerr << "TestALogP::SetUp:cannot initialise charge assigner '" << cmd << "'\n";
  } else {
    std::cerr << "Charge assigner initialised '" << cmd << "'\n";
    _charge_assigner.set_apply_charges_to_molecule(1);
  }
}

TEST_F(TestALogP, TestMethane) {
  _smiles = "C";
  ASSERT_TRUE(_m.build_from_smiles(_smiles));
  _charge_assigner.process(_m);
  std::optional<double> a = _alogp.LogP(_m);
  ASSERT_NE(a, std::nullopt);
  EXPECT_NEAR(*a, 0.636, 0.001);
}

struct SmilesExpected {
  IWString smiles;
  float alogp;
};

class TestAlogpP: public testing::TestWithParam<SmilesExpected> {
  protected:
    Molecule _m;
    alogp::ALogP _alogp;
};

TEST_P(TestAlogpP, TestMolecules) {
  const auto& params = GetParam();

  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  std::optional<double> a = _alogp.LogP(_m);
  ASSERT_NE(a, std::nullopt);
  EXPECT_NEAR(*a, params.alogp, 0.001) << params.smiles << " expected " << params.alogp;
}
INSTANTIATE_TEST_SUITE_P(TestAlogpP, TestAlogpP, testing::Values(
  SmilesExpected{"CC", 1.026},
  SmilesExpected{"CCC", 1.416},
  SmilesExpected{"CCCC", 1.806},
  SmilesExpected{"C1CC1", 1.170},
  SmilesExpected{"CC(C)C", 1.662},
  SmilesExpected{"CC(C)(C)C", 2.052},
  SmilesExpected{"c1ccccc1", 1.687},
  SmilesExpected{"CO", -0.392}
));

}  // namespace
