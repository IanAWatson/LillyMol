// Tests for xlogp3

#include <filesystem>
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/xlogp3.h"

namespace {

struct SmilesExpected {
  IWString smiles;
  float xlogp;
  IWString labelled_smiles;
};

class TestXlogpP: public testing::TestWithParam<SmilesExpected> {
  protected:
    Molecule _mol;

  protected:
};

TEST_P(TestXlogpP, TestMolecules) {
  const auto& params = GetParam();

  ASSERT_TRUE(_mol.build_from_smiles(params.smiles)) << "bad smiles " << params.smiles;

  const int matoms = _mol.natoms();

  std::unique_ptr<int[]> status = std::make_unique<int[]>(matoms);
  std::optional<double> a = xlogp3::XLogP(_mol, status.get());
  ASSERT_NE(a, std::nullopt);
  EXPECT_NEAR(*a, params.xlogp, 0.001) << params.smiles <<
        " expect " << static_cast<float>(params.xlogp) <<
        "\ncomp " << *a << 
        ' ' << _mol.aromatic_smiles();
  if (! params.labelled_smiles.empty()) {
    EXPECT_EQ(_mol.aromatic_smiles(), params.labelled_smiles) << params.labelled_smiles <<
        " got\n" << _mol.aromatic_smiles() << ' ' << params.xlogp << ' ' <<
        _mol.name();
  }
}
INSTANTIATE_TEST_SUITE_P(TestXlogpP, TestXlogpP, testing::Values(
  // xlogp 3.2.2 1.30
  SmilesExpected{"CC ethane", 1.5791, ""},
  // xlogp 3.2.2 2.13
  SmilesExpected{"c1ccccc1 benzene", 1.8941, ""}
));

}  // namespace
