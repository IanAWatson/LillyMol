// Tests for MFormula

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/molecule.h"

#include "mformula.h"

namespace {

struct SmilesFormula {
  IWString smiles;
  int carbon;
  int nitrogen;
  int oxygen;
  int fluorine;
  int phosphorus;
  int sulphur;
  int chlorine;
  int bromine;
  int iodine;

  int natoms;

//SmilesFormula();
};

#ifdef HAS_CONSTRUCTOR
SmilesFormula::SmilesFormula() {
  carbon = 0;
  nitrogen = 0;
  oxygen = 0;
  fluorine = 0;
  phosphorus = 0;
  sulphur = 0;
  chlorine = 0;
  bromine = 0;
  iodine = 0;
  natoms = 0;
}
#endif

class TestMFormula: public testing::TestWithParam<SmilesFormula> {
  protected:
    mformula::MFormula _formula;
};

TEST_P(TestMFormula, Tests) {
  const auto params = GetParam();
  EXPECT_EQ(_formula.Build(params.smiles), params.natoms);
  EXPECT_EQ(_formula.Carbon(), params.carbon);
  EXPECT_EQ(_formula.Nitrogen(), params.nitrogen) << params.smiles;
  EXPECT_EQ(_formula.Oxygen(), params.oxygen);
  EXPECT_EQ(_formula.Fluorine(), params.fluorine);
  EXPECT_EQ(_formula.Phosphorus(), params.phosphorus);
  EXPECT_EQ(_formula.Sulphur(), params.sulphur) << params.smiles;
  EXPECT_EQ(_formula.Chlorine(), params.chlorine) << params.smiles;
  EXPECT_EQ(_formula.Bromine(), params.bromine);
  EXPECT_EQ(_formula.Iodine(), params.iodine);

  EXPECT_TRUE(_formula.initialised());
}
INSTANTIATE_TEST_SUITE_P(TestMFormula, TestMFormula, testing::Values(
  SmilesFormula{"C", 1, 0, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"c", 1, 0, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"N", 0, 1, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"n", 0, 1, 0, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"O", 0, 0, 1, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"o", 0, 0, 1, 0, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"F", 0, 0, 0, 1, 0, 0, 0, 0, 0,   1},
  SmilesFormula{"P", 0, 0, 0, 0, 1, 0, 0, 0, 0,   1},
  SmilesFormula{"S", 0, 0, 0, 0, 0, 1, 0, 0, 0,   1},
  SmilesFormula{"s", 0, 0, 0, 0, 0, 1, 0, 0, 0,   1},
  SmilesFormula{"Cl", 0, 0, 0, 0, 0, 0, 1, 0, 0,   1},
  SmilesFormula{"Br", 0, 0, 0, 0, 0, 0, 0, 1, 0,   1},
  SmilesFormula{"I", 0, 0, 0, 0, 0, 0, 0, 0, 1,   1},

  SmilesFormula{"CC", 2, 0, 0, 0, 0, 0, 0, 0, 0,   2},
  SmilesFormula{"CCC", 3, 0, 0, 0, 0, 0, 0, 0, 0,   3},
  SmilesFormula{"CN", 1, 1, 0, 0, 0, 0, 0, 0, 0,   2},
  SmilesFormula{"CCl", 1, 0, 0, 0, 0, 0, 1, 0, 0,   2},
  SmilesFormula{"ClC", 1, 0, 0, 0, 0, 0, 1, 0, 0,   2},
  SmilesFormula{"[Cl][C]", 1, 0, 0, 0, 0, 0, 1, 0, 0,   2},
  SmilesFormula{"[4NH2+][12C]", 1, 1, 0, 0, 0, 0, 0, 0, 0,   2}
));

// The positions in the fixed counted fingerprint, which must be kept in sync with
// the constants in mformula.cc.
constexpr int kHydrogenOnHeteroatom = 13;
constexpr int kHydrogenAromatic = 14;
constexpr int kHydrogen = 15;
constexpr int kRingAtom = 16;

std::vector<int>
FixedCounted(const IWString& smiles, bool include_hydrogens) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles(smiles));

  mformula::MFormula formula;
  formula.set_include_hydrogens(include_hydrogens);
  EXPECT_TRUE(formula.Build(m));

  std::vector<int> result(mformula::kMFOther + 1);
  EXPECT_TRUE(formula.ToFixedCountedFingerprint(result.data(), result.size()));
  return result;
}

TEST(TestMFormulaHydrogens, HydrogenCountsAreOptional) {
  // Phenol has an OH, 5 aromatic CH atoms and no other carbon with a hydrogen.
  const std::vector<int> with = FixedCounted("Oc1ccccc1", true);
  EXPECT_EQ(with[kHydrogenOnHeteroatom], 1);
  EXPECT_EQ(with[kHydrogenAromatic], 5);
  EXPECT_EQ(with[kHydrogen], 0);
  EXPECT_EQ(with[kRingAtom], 6);

  const std::vector<int> without = FixedCounted("Oc1ccccc1", false);
  EXPECT_EQ(without[kHydrogenOnHeteroatom], 0);
  EXPECT_EQ(without[kHydrogenAromatic], 0);
  EXPECT_EQ(without[kHydrogen], 0);

  // Nothing else changes.
  for (int i = 0; i < 13; ++i) {
    EXPECT_EQ(with[i], without[i]) << i;
  }
  EXPECT_EQ(with[kRingAtom], without[kRingAtom]);
  EXPECT_EQ(with[mformula::kMFOther], without[mformula::kMFOther]);
}

// The point of leaving out the hydrogen counts: a molecule that contains the
// query, but with substituted atoms, has at least every count of the query.
TEST(TestMFormulaHydrogens, SubstitutedSupersetHasAtLeastTheQueryCounts) {
  const std::vector<int> query = FixedCounted("Oc1ccccc1", false);
  const std::vector<int> candidate = FixedCounted("Cc1ccccc1O", false);
  for (size_t i = 0; i < query.size(); ++i) {
    EXPECT_GE(candidate[i], query[i]) << i;
  }

  // With the hydrogen counts the aromatic CH count of the candidate is lower.
  const std::vector<int> query_h = FixedCounted("Oc1ccccc1", true);
  const std::vector<int> candidate_h = FixedCounted("Cc1ccccc1O", true);
  EXPECT_LT(candidate_h[kHydrogenAromatic], query_h[kHydrogenAromatic]);
}

}  // namespace
