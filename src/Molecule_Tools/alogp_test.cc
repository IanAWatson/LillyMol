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

    // We can compare type atom types.
    IWString _expected_smiles;

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

  _alogp.set_label_with_atom_type(1);
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
  IWString labelled_smiles;
};

class TestAlogpP: public testing::TestWithParam<SmilesExpected> {
  protected:
    Molecule _mol;
    alogp::ALogP _alogp;
};

TEST_P(TestAlogpP, TestMolecules) {
  const auto& params = GetParam();

  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  _alogp.set_label_with_atom_type(1);

  std::optional<double> a = _alogp.LogP(_mol);
  ASSERT_NE(a, std::nullopt);
  EXPECT_NEAR(*a, params.alogp, 0.001) << params.smiles << " expect " <<
        static_cast<float>(params.alogp) <<
        " comp " << *a << 
        ' ' << _mol.aromatic_smiles();
  if (! params.labelled_smiles.empty()) {
    EXPECT_EQ(_mol.aromatic_smiles(), params.labelled_smiles) << params.labelled_smiles <<
        " got " << _mol.aromatic_smiles();
  }
}
INSTANTIATE_TEST_SUITE_P(TestAlogpP, TestAlogpP, testing::Values(
  SmilesExpected{"CC", 1.026, ""},
  SmilesExpected{"CCC", 1.416, ""},
  SmilesExpected{"CCCC", 1.806, ""},
  SmilesExpected{"C1CC1", 1.170, ""},
  SmilesExpected{"CC(C)C", 1.662, ""},
  SmilesExpected{"CC(C)(C)C", 2.052, ""},
  SmilesExpected{"c1ccccc1", 1.687, ""},
  SmilesExpected{"CO", -0.392, "[3CH3][50OH]"},
  SmilesExpected{"c1ccccc1O", 1.392, "[18cH]1[18cH][18cH][18cH][18cH][23c]1[50OH]"},
  SmilesExpected{"c1ccccc1S", 1.975, "[18cH]1[18cH][18cH][18cH][18cH][23c]1[50OH]"},
  SmilesExpected{"c1ccccc1F", 1.826, "[18cH]1[18cH][18cH][18cH][18cH][14c]1[62F]"},
  SmilesExpected{"c1ccccc1Cl", 2.340, "[18cH]1[18cH][18cH][18cH][18cH][15c]1[63Cl]"},
  SmilesExpected{"c1ccccc1Br", 2.449, "[18cH]1[18cH][18cH][18cH][18cH][16c]1[64Br]"},
  SmilesExpected{"c1ccccc1I", 2.291, "[18cH]1[18cH][18cH][18cH][18cH][17c]1[65I]"},
  SmilesExpected{"CC(F)(F)F", 1.569, "[1CH3][4C]([62F])([62F])[62F]"},
  SmilesExpected{"COC", 0.263, "[3CH3][51O][3CH3]"},
  SmilesExpected{"C12=C(C=CC=C1)C=CN2", 2.168, "[19c]12[19c]([18cH][18cH][18cH][18cH]1)[18cH][18cH][44nH]2"},
  SmilesExpected{"c1ccccc1c1ccccc1", 3.354, "[18cH]1[18cH][18cH][18cH][18cH][24c]1[68SH]"},
  SmilesExpected{"OC(C)C", 0.387, "[50OH][4CH]([1CH3])[1CH3]"}
));

}  // namespace
