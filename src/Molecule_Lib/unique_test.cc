// Tester for unique smiles generation.

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "aromatic.h"
#include "molecule.h"
#include "smiles.h"

namespace {

// A smiles and an expected unique smiles.
struct SmilesAndExpected {
  IWString smiles;
  IWString expected;
};

class TestUniqueSmiles: public testing::TestWithParam<SmilesAndExpected> {
  protected:
    void SetUp();
};

void
TestUniqueSmiles::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestUniqueSmiles, Lots) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.smiles));
//std::cerr << "usmi " << m.unique_smiles() << std::endl;
  EXPECT_EQ(m.unique_smiles(), params.expected);
}

INSTANTIATE_TEST_SUITE_P(RandomlySelectedMolecules, TestUniqueSmiles, testing::Values(
   SmilesAndExpected{"C", "C"},
   SmilesAndExpected{"CC", "CC"},
   SmilesAndExpected{"OC", "OC"},
   SmilesAndExpected{"CO", "OC"},
   SmilesAndExpected{"O=C(O)C12CC3(CC(C1)CC(C2)C3)", "O=C(O)C12CC3CC(CC(C3)C1)C2"},
   SmilesAndExpected{"S1C=CC=C1", "s1cccc1"}
));

class TestUniqueKekuleSmiles: public testing::TestWithParam<SmilesAndExpected> {
  protected:
    Molecule _m;

    void SetUp();
};

void
TestUniqueKekuleSmiles::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestUniqueKekuleSmiles, Kekule) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.UniqueKekuleSmiles(), params.expected);
}

INSTANTIATE_TEST_SUITE_P(AromaticMolecules, TestUniqueKekuleSmiles, testing::Values(
  SmilesAndExpected{"C", "C"},
  SmilesAndExpected{"CC", "CC"},
  SmilesAndExpected{"C1CC1", "C1CC1"},
  SmilesAndExpected{"c1ccccc1", "c1=cc=cc=c1"},
  SmilesAndExpected{"N1=CC=CC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=NC=CC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CN=CC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CC=NC=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CC=CN=C1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=CC=CC=N1", "[n]1=cc=cc=c1"},
  SmilesAndExpected{"C1=C2C=CC=CC2=NC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C12=CC=CN=C2C=CC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C1=C2C=CC=CC2=CN=C1", "[n]1=cc2=c(c=c1)c=cc=c2"},
  SmilesAndExpected{"C1=CC=CC2=CC=NC=C12", "[n]1=cc2=c(c=c1)c=cc=c2"}
));

class TestUniqueKekuleSmilesRandom: public testing::TestWithParam<SmilesAndExpected> {
  protected:
    Molecule _m;

  void SetUp();
};

void
TestUniqueKekuleSmilesRandom::SetUp() {
  reset_aromatic_file_scope_variables();
}

TEST_P(TestUniqueKekuleSmilesRandom, Kekule) {
  const auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(_m.UniqueKekuleSmiles(), params.expected);
  constexpr int kNtest = 10;
  for (int i = 0; i < kNtest; ++i) {
    const IWString& smi = _m.random_smiles();
    Molecule mcopy;
    ASSERT_TRUE(mcopy.build_from_smiles(smi));
    ASSERT_EQ(mcopy.unique_smiles(), _m.unique_smiles());
    EXPECT_EQ(mcopy.UniqueKekuleSmiles(), params.expected);
  }
}

INSTANTIATE_TEST_SUITE_P(KekuleForms, TestUniqueKekuleSmilesRandom, testing::Values(
  SmilesAndExpected{"C1=C2C=CC=CC2=NC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C12=CC=CN=C2C=CC=C1", "[n]1=cc=cc2=c1c=cc=c2"},
  SmilesAndExpected{"C1=C2C=CC=CC2=CN=C1", "[n]1=cc2=c(c=c1)c=cc=c2"},
  SmilesAndExpected{"C1=CC=CC2=CC=NC=C12", "[n]1=cc2=c(c=c1)c=cc=c2"},
  SmilesAndExpected{"N1C=CC2=C1C=CC=C2", "[nH]1c=cc2=cc=cc=c12"},  // indole
  SmilesAndExpected{"N1=CN=CC=C1", "c1=c[n]=c[n]=c1"},   // pyrimidine
  SmilesAndExpected{"N1=CNC2=C1C=NC=N2", "[n]1=cc2=c([n]=c1)[nH]c=[n]2"},
  SmilesAndExpected{"C1=CC=CC2=C1CC3=C(O2)C=CC=C3", "O1c2=c(Cc3=cc=cc=c13)c=cc=c2"},  // xanthene
  SmilesAndExpected{"C1=CC=CC2=C1NC3=C2C=CC=C3", "[nH]1c2=c(c=cc=c2)c2=c1c=cc=c2"}  // carbazole
));

// Bracket-specified and inferred H counts must give the same canonical text
// when the specified count equals the normal valence-derived count. Ranking
// already considers those atoms equivalent; emission must agree with ranking.
class TestRedundantBracketHydrogens : public TestUniqueSmiles {};

TEST_P(TestRedundantBracketHydrogens, CanonicalAndRoundTrip) {
  const auto params = GetParam();
  Molecule reference;
  ASSERT_TRUE(reference.build_from_smiles(params.expected));
  const IWString expected = reference.unique_smiles();
  Molecule original;
  ASSERT_TRUE(original.build_from_smiles(params.smiles));
  const IWString formula = original.molecular_formula();
  std::vector<int> known(original.natoms());
  for (atom_number_t a = 0; a < original.natoms(); ++a) {
    known[a] = original.implicit_hydrogens_known(a);
  }
  EXPECT_EQ(original.unique_smiles(), expected);
  // Canonical output must not change the stored hydrogen specifications.
  for (atom_number_t a = 0; a < original.natoms(); ++a) {
    EXPECT_EQ(original.implicit_hydrogens_known(a), known[a]);
  }
  for (int i = 0; i < 100; ++i) {
    Molecule molecule;
    ASSERT_TRUE(molecule.build_from_smiles(original.random_smiles()));
    EXPECT_EQ(molecule.unique_smiles(), expected);
    EXPECT_EQ(molecule.molecular_formula(), formula);
    Molecule roundtrip;
    ASSERT_TRUE(roundtrip.build_from_smiles(molecule.unique_smiles()));
    EXPECT_EQ(roundtrip.unique_smiles(), expected);
    EXPECT_EQ(roundtrip.molecular_formula(), formula);
  }
}

INSTANTIATE_TEST_SUITE_P(BracketHydrogens, TestRedundantBracketHydrogens, testing::Values(
    SmilesAndExpected{"[NH]1C(=Nc2ccccc2)NCC1", "N1C(=Nc2ccccc2)NCC1"},
    SmilesAndExpected{"C[P+](O)([OH])OC", "C[P+](O)(O)OC"},
    SmilesAndExpected{"[CH4]", "C"},
    SmilesAndExpected{"[NH3]", "N"},
    SmilesAndExpected{"[OH2]", "O"},
    SmilesAndExpected{"[SH2]", "S"},
    SmilesAndExpected{"C[OH]", "CO"},
    SmilesAndExpected{"C[NH2]", "CN"},
    SmilesAndExpected{"C[CH2]C", "CCC"},
    SmilesAndExpected{"[cH]1ccccc1", "c1ccccc1"},
    SmilesAndExpected{"[nH]1cccc1", "[nH]1cccc1"},
    SmilesAndExpected{"C[13CH2]C", "C[13CH2]C"},
    SmilesAndExpected{"[NH4+]", "[NH4+]"},
    SmilesAndExpected{"F[C@H](Cl)C", "F[C@H](Cl)C"},
    SmilesAndExpected{"C[CH2:17]C", "C[CH2:17]C"},
    SmilesAndExpected{"O[2H]", "O[2H]"},
    SmilesAndExpected{"[CH2]", "[CH2]"},
    SmilesAndExpected{"[C]C", "[C]C"},
    SmilesAndExpected{"C[N]C", "C[N]C"}
));

TEST(UniqueSmilesHydrogens, CanonicalSubsetOmitsRedundantBrackets) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C[CH2]C"));
  Smiles_Information subset;
  const int include[] = {1, 1, 1};
  EXPECT_EQ(m.unique_smiles(subset, include), "CCC");
  EXPECT_TRUE(m.implicit_hydrogens_known(1));
}

TEST(UniqueSmilesHydrogens, OrdinarySmilesRetainsSpecifiedHydrogens) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C[CH2]C"));
  EXPECT_TRUE(m.smiles().contains("[CH2]"));
  m.unique_smiles();
  const IWString random = m.random_smiles();
  EXPECT_TRUE(random.contains("[CH2]"));
  EXPECT_TRUE(m.implicit_hydrogens_known(1));
}

TEST(UniqueSmilesHydrogens, RadicalsRemainDistinctFromNormalValence) {
  for (const auto& [radical, saturated] : std::vector<std::pair<const char*, const char*>>{
      {"[CH2]", "C"}, {"C[N]C", "CNC"}, {"[C]C", "CC"}}) {
    Molecule r, s;
    ASSERT_TRUE(r.build_from_smiles(radical));
    ASSERT_TRUE(s.build_from_smiles(saturated));
    EXPECT_NE(r.unique_smiles(), s.unique_smiles());
  }
}

}  // namespace
