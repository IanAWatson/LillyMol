// Tests for alogp

#include <filesystem>
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/alogp.h"

namespace {

#ifdef OLD_VERSION_NOT_USED__
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
#endif

struct SmilesExpected {
  IWString smiles;
  float alogp;
  IWString labelled_smiles;
};

class TestAlogpP: public testing::TestWithParam<SmilesExpected> {
  protected:
    Molecule _mol;
    alogp::ALogP _alogp;

    Charge_Assigner _charge_assigner;

  protected:
    void SetUp();
};

void
TestAlogpP::SetUp() {
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

  _alogp.set_use_alcohol_for_acid(1);
  _alogp.set_label_with_atom_type(1);
  _alogp.set_rdkit_charged_nitrogen(1);
}

TEST_P(TestAlogpP, TestMolecules) {
  const auto& params = GetParam();

  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));

  _charge_assigner.process(_mol);

  std::optional<double> a = _alogp.LogP(_mol);
  ASSERT_NE(a, std::nullopt);
  EXPECT_NEAR(*a, params.alogp, 0.001) << params.smiles <<
        " expect " << static_cast<float>(params.alogp) <<
        "\ncomp " << *a << 
        ' ' << _mol.aromatic_smiles();
  if (! params.labelled_smiles.empty()) {
    EXPECT_EQ(_mol.aromatic_smiles(), params.labelled_smiles) << params.labelled_smiles <<
        " got\n" << _mol.aromatic_smiles() << ' ' << params.alogp << ' ' <<
        _mol.name();
  }
}
INSTANTIATE_TEST_SUITE_P(TestAlogpP, TestAlogpP, testing::Values(
  SmilesExpected{"CC ethane", 1.026, ""},
  SmilesExpected{"CCC propane", 1.416, ""},
  SmilesExpected{"CCCC butane", 1.806, ""},
  SmilesExpected{"C1CC1 cyclopropane", 1.170, ""},
  SmilesExpected{"CC(C)C isobutane", 1.662, ""},
  SmilesExpected{"CC(C)(C)C neopentane", 2.052, ""},
  SmilesExpected{"c1ccccc1 benzene", 1.687, ""},
  SmilesExpected{"CO methanol", -0.392, "[3CH3][50OH]"},
  SmilesExpected{"c1ccccc1O phenol", 1.392, "[18cH]1[18cH][18cH][18cH][18cH][23c]1[50OH]"},
  SmilesExpected{"c1ccccc1S benzenethiol", 1.975, "[18cH]1[18cH][18cH][18cH][18cH][24c]1[68SH]"},
  SmilesExpected{"c1ccccc1F fluoro-benzene", 1.826, "[18cH]1[18cH][18cH][18cH][18cH][14c]1[62F]"},
  SmilesExpected{"c1ccccc1Cl chloro-benzene", 2.340, "[18cH]1[18cH][18cH][18cH][18cH][15c]1[63Cl]"},
  SmilesExpected{"c1ccccc1Br bromo-benzene", 2.449, "[18cH]1[18cH][18cH][18cH][18cH][16c]1[64Br]"},
  SmilesExpected{"c1ccccc1I iodo-benzene", 2.291, "[18cH]1[18cH][18cH][18cH][18cH][17c]1[65I]"},
  SmilesExpected{"CC(F)(F)F 1,1,1-trofluoroethane", 1.569, "[1CH3][4C]([62F])([62F])[62F]"},
  SmilesExpected{"COC dimethyl ether", 0.263, "[3CH3][51O][3CH3]"},
  SmilesExpected{"C12=C(C=CC=C1)C=CN2 indole", 2.168, "[19c]12[19c]([18cH][18cH][18cH][18cH]1)[18cH][18cH][44nH]2"},
  SmilesExpected{"c1ccccc1c1ccccc1 biphenyl", 3.354, "[18cH]1[18cH][18cH][18cH][18cH][20c]1[20c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"OC(C)C isopropyl alcohol", 0.387, "[50OH][4CH]([1CH3])[1CH3]"},
  SmilesExpected{"o1cccc1 furan", 1.280, "[49o]1[18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"s1cccc1 thiophene", 1.748, "[70s]1[18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC(=O)O acetic acid", 0.091, "[1CH3][5C](=[57O])[50O-]"},
  SmilesExpected{"CC=O acetaldehyde", 0.205, "[1CH3][5CH]=[57O]"},
  SmilesExpected{"CC#N acetonitrile", 0.530, "[1CH3][7C]#[42N]"},
  SmilesExpected{"C=C methene", 0.802, "[6CH2]=[6CH2]"},
  SmilesExpected{"CC=N ethanimine", 0.656, "[1CH3][5CH]=[40NH]"},
  SmilesExpected{"CCc1ccccc1 ethylbenzene", 2.249, "[1CH3][10CH2][21c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC(C)c1ccccc1 isopropylbenezene", 2.810, "[1CH3][11CH]([1CH3])[21c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC(C)(C)c1ccccc1 t-butyl benzene", 2.984, "[1CH3][12C]([1CH3])([1CH3])[21c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"C=Cc1ccccc1 styrene", 2.330, "[6CH2]=[26CH][21c]1[18cH][18cH][18cH][18cH][18cH]1"},

  // We only get concordance with RDKit if we assume 2 Hydrogens on the N+.
  SmilesExpected{"CN methylamine", -0.425, "[3CH3][34NH3+]"},
  SmilesExpected{"CNC dimethyl methylamine", -0.164, "[3CH3][35NH2+][3CH3]"},
  SmilesExpected{"n1ccccc1 pyridine", 1.082, "[44n]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"[nH]1cccc1 pyrole", 1.015, "[44nH]1[18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"Nc1ccccc1 aniline", 1.269, "[36NH2][22c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CNc1ccccc1 N-methylaniline", 1.728, "[3CH3][37NH][22c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CC=NC N-methylethanimine", 0.707, "[1CH3][5CH]=[39N][3CH3]"},
  SmilesExpected{"CN(C)C trimethylamine", 0.178, "[3CH3][40NH+]([3CH3])[3CH3]"},
  SmilesExpected{"CN(C)c1ccccc1 N,N-dimethylaniline", 1.753, "[3CH3][41N]([3CH3])[22c]1[18cH][18cH][18cH][18cH][18cH]1"},
  SmilesExpected{"CNC(=O)C N-METHYLACETAMIDE", -0.248, "[3CH3][35NH][5C](=[57O])[1CH3]"},
  SmilesExpected{"CNC(=O)NC 1,3-DIMETHYLUREA", -0.455, "[3CH3][35NH][5C](=[59O])[35NH][3CH3]"},
  SmilesExpected{"O=C1NN=CN1 CHEMBL1865594", -0.902, "[56O]=[25c]1[44nH][44n][18cH][44nH]1"},
  SmilesExpected{"C1(=S)C=CSS1 CHEMBL368700", 2.539, "[44n]1[18cH][18cH][18cH][18cH][18cH]1"}
));

}  // namespace
