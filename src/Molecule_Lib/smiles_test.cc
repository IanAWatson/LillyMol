// Tests for smiles reading and generation

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "molecule.h"
#include "smiles.h"

namespace {

// For tests that need an input smiles and an expected output.
struct SmilesAndExpected {
  IWString input_smiles;
  IWString expected_smiles;
};

// For tests that need an input, a parameter to set, and an expected output.
struct SmilesAndParam {
  IWString input_smiles;
  int int_param;
  IWString expected_smiles;
};

class TestSmilesCharges : public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestSmilesCharges, MultiplePlusCharges) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  set_write_formal_charge_as_consecutive_signs(params.int_param);
  EXPECT_EQ(m.smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestMultipleCharges, TestSmilesCharges, testing::Values(
   SmilesAndParam{"C", 0, "C"},
   SmilesAndParam{"C", 1, "C"},
   SmilesAndParam{"[Fe+++]", 1, "[Fe+++]"},
   SmilesAndParam{"[Fe+++]", 0, "[Fe+3]"},
   SmilesAndParam{"[P---]", 1, "[P---]"},
   SmilesAndParam{"[P---]", 0, "[P-3]"},
   SmilesAndParam{"[Zn--]", 1, "[Zn--]"},
   SmilesAndParam{"[Zn--]", 0, "[Zn-2]"}
));

class TestAromaticSmiles : public testing::TestWithParam<SmilesAndExpected> {
  public:
    void SetUp() {
      set_include_aromaticity_in_smiles(1);
    }
    void TearDown() {
      set_include_aromaticity_in_smiles(0);
    }
};

TEST_P(TestAromaticSmiles, AromaticSmiles) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  EXPECT_EQ(m.smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestAromaticSmiles, TestAromaticSmiles, testing::Values(
  SmilesAndExpected{"C", "C"},
  SmilesAndExpected{"C1=CC=CC=C1", "c1ccccc1"}
));

class TestHOnAromaticNP: public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestHOnAromaticNP, TestHOnAromaticNP) {
  const auto params = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  set_include_implicit_hydrogens_on_aromatic_n_and_p(params.int_param);
  EXPECT_EQ(m.unique_smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestHOnAromaticNP, TestHOnAromaticNP, testing::Values(
  SmilesAndParam{"C", 1, "C"},
  SmilesAndParam{"C1=CC=CN1", 1, "[nH]1cccc1"},
  SmilesAndParam{"C1=CC=CN1", 0, "n1cccc1"}
));

class TestAddHToIsotopes : public testing::TestWithParam<SmilesAndParam> {
};

TEST_P(TestAddHToIsotopes, TestAddHToIsotopes) {
  const auto params = GetParam();
  set_unset_implicit_hydrogens_known_if_possible(params.int_param);
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(params.input_smiles));
  EXPECT_EQ(m.unique_smiles(), params.expected_smiles);
}

INSTANTIATE_TEST_SUITE_P(TestAddHToIsotopes, TestAddHToIsotopes, testing::Values(
  SmilesAndParam{"C", 1, "C"},
  SmilesAndParam{"[2C]", 1, "[2C]"},
  SmilesAndParam{"[2C]C", 1, "[2C]C"},
  SmilesAndParam{"C[2C]C", 1, "C[2C]C"}
//SmilesAndParam{"[9C]1=CC=CC=C1", 1, "[9cH]1ccccc1"}  not sure why this does not work, investigate.
));


// Ordinary coordinate SMILES remain supported after retiring boxed coordinates.
class TestSmilesCoordinates : public testing::Test {
 protected:
  void SetUp() override {
    lillymol::set_include_coordinates_with_smiles(1);
  }
  void TearDown() override {
    lillymol::set_include_coordinates_with_smiles(0);
  }
};

TEST_F(TestSmilesCoordinates, RoundTripCoordinates) {
  Molecule molecule;
  ASSERT_TRUE(molecule.build_from_smiles(
      "C{{1.25,-2.5,3}}[NH3+]{{-4,5.5,-6.25}}"));
  ASSERT_EQ(molecule.natoms(), 2);
  EXPECT_FLOAT_EQ(molecule.x(0), 1.25f);
  EXPECT_FLOAT_EQ(molecule.y(0), -2.5f);
  EXPECT_FLOAT_EQ(molecule.z(0), 3.0f);
  EXPECT_FLOAT_EQ(molecule.x(1), -4.0f);
  EXPECT_FLOAT_EQ(molecule.y(1), 5.5f);
  EXPECT_FLOAT_EQ(molecule.z(1), -6.25f);

  Molecule roundtrip;
  ASSERT_TRUE(roundtrip.build_from_smiles(molecule.smiles()));
  ASSERT_EQ(roundtrip.natoms(), molecule.natoms());
  for (int i = 0; i < molecule.natoms(); ++i) {
    EXPECT_FLOAT_EQ(roundtrip.x(i), molecule.x(i));
    EXPECT_FLOAT_EQ(roundtrip.y(i), molecule.y(i));
    EXPECT_FLOAT_EQ(roundtrip.z(i), molecule.z(i));
  }
}

struct RemoveHTest {
  IWString starting_smiles;
  IWString expected_smiles;
};

// Tests to see if we can convert an explicit Hydrogen to an
// implicit Hydrogen when a chiral centre is involved.
class TestRemoveH : public testing::TestWithParam<RemoveHTest> {
  protected:
    Molecule _mol;
};

TEST_P(TestRemoveH, Test) {
  const auto params = GetParam();
  ASSERT_TRUE(_mol.build_from_smiles(params.starting_smiles));
  _mol.remove_all(1);
  EXPECT_EQ(_mol.unique_smiles(), params.expected_smiles) << "Mismatch " << _mol.unique_smiles() << " vs " << params.expected_smiles << '\n';
}

INSTANTIATE_TEST_SUITE_P(TestRemoveH, TestRemoveH, testing::Values(
  RemoveHTest{"C", "C"},
  RemoveHTest{"[H][C@](F)(N)C", "F[C@@H](N)C"},
  RemoveHTest{"F[C@]([H])(C)N", "F[C@@H](N)C"},
  RemoveHTest{"N[C@](C)([H])F", "F[C@@H](N)C"},
  RemoveHTest{"C[C@](N)(F)[H]", "F[C@@H](N)C"},

  RemoveHTest{"[H][C@@](F)(N)C", "F[C@H](N)C"},
  RemoveHTest{"F[C@@]([H])(C)N", "F[C@H](N)C"},
  RemoveHTest{"N[C@@](C)([H])F", "F[C@H](N)C"},
  RemoveHTest{"C[C@@](N)(F)[H]", "F[C@H](N)C"}

));

struct SmilesResult {
  IWString smiles;
  int expected_ok;
  IWString echo_smiles;
  std::vector<atomic_number_t> atomic_number;
};

class TestAnyLength : public testing::TestWithParam<SmilesResult> {
  protected:
    Molecule _mol;

  void SetUp() {
    set_atomic_symbols_can_have_arbitrary_length(1);
    set_include_aromaticity_in_smiles(0);
  }
  void TearDown() {
    set_atomic_symbols_can_have_arbitrary_length(0);
  }
};
TEST_P(TestAnyLength, Test) {
  const auto params = GetParam();
  if (_mol.build_from_smiles(params.smiles)) {
    EXPECT_EQ(params.expected_ok, 1) << params.smiles << " not ok";
  } else {
    EXPECT_EQ(params.expected_ok, 0) << params.smiles << " should not be ok";
  }

  if (! params.expected_ok) {
    return;
  }

  EXPECT_EQ(params.echo_smiles, _mol.smiles()) << params.smiles << ' ' << _mol.smiles();

  for (int i = 0; i < _mol.natoms(); ++i) {
    EXPECT_EQ(_mol.atomic_number(i), params.atomic_number[i]) << params.smiles;
  }

}
INSTANTIATE_TEST_SUITE_P(TestAnyLength, TestAnyLength, testing::Values(
  SmilesResult{"C", 1, "C", {6}},
  SmilesResult{"1", 0, "", {}},
  SmilesResult{"[]", 0, "", {}},
  SmilesResult{"[1]", 0, "", {}},
  SmilesResult{"[C]", 1, "[C]", {6}},
  SmilesResult{"[CH4]", 1, "C", {6}},
  SmilesResult{"[NH3]", 1, "N", {7}},
  SmilesResult{"[OH2]", 1, "O", {8}},
  SmilesResult{"[F]", 1, "[F]", {9}},
  SmilesResult{"[1F]", 1, "[1F]", {9}},
  SmilesResult{"[K]", 1, "[K]", {19}},
  SmilesResult{"[Fr]", 1, "[Fr]", {87}},
  SmilesResult{"[Ca]", 1, "[Ca]", {20}},
  SmilesResult{"[Na]", 1, "[Na]", {11}},

  SmilesResult{"[1]", 0, "", {}},
  SmilesResult{"[12]", 0, "", {}},
  SmilesResult{"[+]", 0, "", {}},
  SmilesResult{"[-]", 0, "", {}},

  SmilesResult{"[O-]", 1, "[O-]", {8}},
  SmilesResult{"[N+]", 1, "[N+]", {7}},
  SmilesResult{"[NH+]", 1, "[NH+]", {7}},

  SmilesResult{"[Foo]", 1, "[Foo]", {kInvalidAtomicNumber}},
  SmilesResult{"[1Foo]", 1, "[1Foo]", {kInvalidAtomicNumber}},
  SmilesResult{"[1Foo+]", 1, "[1Foo+]", {kInvalidAtomicNumber}},
  SmilesResult{"[1FooH+]", 1, "[1FooH+]", {kInvalidAtomicNumber}},
  SmilesResult{"c2nnn[nH]2", 1, "C1=NN=NN1", {6, 7, 7, 7, 7}}
));


}  // namespace
