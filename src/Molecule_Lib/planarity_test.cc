#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/planarity.h"
#include "Molecule_Lib/smiles.h"

namespace iwplanarity {

std::ostream&
operator<<(std::ostream& os, const PlanarityStatus status) {
  switch (status) {
    case PlanarityStatus::kPlanar:
      return os << "planar";
    case PlanarityStatus::kNonPlanar:
      return os << "nonplanar";
    case PlanarityStatus::kError:
      return os << "error";
  }

  return os << "unknown";
}

std::ostream&
operator<<(std::ostream& os, const PlanarityResult& result) {
  return os << result.status
            << " obstruction_bonds " << result.obstruction_bonds.size();
}

}  // namespace iwplanarity

namespace {

struct PlanarityCase {
  IWString smiles;
  iwplanarity::PlanarityStatus expected_status;
  bool expect_obstruction_bonds;
};

class TestPlanarity : public testing::TestWithParam<PlanarityCase> {
 protected:
  Molecule BuildMolecule(const IWString& smiles) {
    Molecule m;
    EXPECT_TRUE(m.build_from_smiles(smiles)) << smiles;
    return m;
  }
};

TEST_P(TestPlanarity, TestStatus) {
  const PlanarityCase& test_case = GetParam();

  Molecule m = BuildMolecule(test_case.smiles);
  ASSERT_GT(m.natoms(), 0) << test_case.smiles;

  const iwplanarity::PlanarityResult result = iwplanarity::Planarity(m);

  EXPECT_EQ(result.status, test_case.expected_status) << test_case.smiles << " status "
            << test_case.expected_status << " cmp " << result.status;

  if (test_case.expect_obstruction_bonds) {
    EXPECT_FALSE(result.obstruction_bonds.empty()) << test_case.smiles;
  } else {
    EXPECT_TRUE(result.obstruction_bonds.empty()) << test_case.smiles;
  }
}

INSTANTIATE_TEST_SUITE_P(
    TestPlanarity,
    TestPlanarity,
    testing::Values(
        PlanarityCase{
            "C1CCCCC1 cyclohexane",
            iwplanarity::PlanarityStatus::kPlanar,
            false,
        },
        PlanarityCase{
            "c1ccccc1 benzene",
            iwplanarity::PlanarityStatus::kPlanar,
            false,
        },
        PlanarityCase{
            "C12C3C4C1C5C2C3C45 cubane",  // not crossed, draw one square inside the other and join bonds.
            iwplanarity::PlanarityStatus::kPlanar,
            false,
        },
        PlanarityCase{
          "C123C45C16C24C356 K5",
          iwplanarity::PlanarityStatus::kNonPlanar,
          true,
        },
        PlanarityCase{
          "C12C3C4C1C3C24 K3,3",
          iwplanarity::PlanarityStatus::kNonPlanar,
          true,
        }
    ));

}  // namespace
