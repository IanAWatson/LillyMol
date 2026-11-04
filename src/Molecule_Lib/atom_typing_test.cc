#include <cstdint>
#include <vector>

#include "gtest/gtest.h"

#include "aromatic.h"
#include "atom_typing.h"

namespace {

struct SmartsForAtomTypeTestCase {
  const char* smiles;
  atom_number_t atom;
  uint32_t atom_type;
  const char* expected;
};

class SmartsForAtomTypeTest
    : public testing::TestWithParam<SmartsForAtomTypeTestCase> {
 protected:
  void SetUp() override {
    set_global_aromaticity_type(Daylight);
  }
};

TEST_P(SmartsForAtomTypeTest, GeneratesExpectedSmarts) {
  const SmartsForAtomTypeTestCase& test = GetParam();
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(test.smiles));
  ASSERT_GE(test.atom, 0);
  ASSERT_LT(test.atom, m.natoms());

  const IWString smarts =
      lillymol::SmartsForAtomType(m, test.atom, test.atom_type);

  EXPECT_EQ(smarts.AsString(), test.expected);
}

INSTANTIATE_TEST_SUITE_P(
    AtomProperties, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_Z, "[#6]"},
        SmartsForAtomTypeTestCase{"N", 0, IWATTYPE_USP_Z, "[#7]"},
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_A, "[A]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0, IWATTYPE_USP_A, "[a]"},
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_Z | IWATTYPE_USP_A,
                                  "[#6]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A, "[c]"},
        SmartsForAtomTypeTestCase{"Br", 0, IWATTYPE_USP_Y, "[#17]"},
        SmartsForAtomTypeTestCase{"I", 0, IWATTYPE_USP_Y, "[#17]"},
        SmartsForAtomTypeTestCase{"F", 0, IWATTYPE_USP_Y, "[#9]"},
        SmartsForAtomTypeTestCase{"C", 0, IWATTYPE_USP_N, "[*]"}));

INSTANTIATE_TEST_SUITE_P(
    Connections, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"CC(C)(C)C", 1, IWATTYPE_USP_C, "[D4]"},
        SmartsForAtomTypeTestCase{"CC(C)(C)C", 1,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A |
                                      IWATTYPE_USP_C,
                                  "[#6;D4]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A |
                                      IWATTYPE_USP_C,
                                  "[c;D2]"},
        SmartsForAtomTypeTestCase{"Cc1ccccc1", 1,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_A |
                                      IWATTYPE_USP_C,
                                  "[c;D3]"}));

INSTANTIATE_TEST_SUITE_P(
    RingsAndUnsaturation, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"C1CC1", 0, IWATTYPE_USP_R, "[x2]"},
        SmartsForAtomTypeTestCase{"C1CC1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_R, "[#6;x2]"},
        SmartsForAtomTypeTestCase{"C1CC1", 0, IWATTYPE_USP_S, "[r3]"},
        SmartsForAtomTypeTestCase{"C1CC1", 0, IWATTYPE_USP_L, "[r3]"},
        SmartsForAtomTypeTestCase{"C1CC2(C1)CCCC2", 2, IWATTYPE_USP_S, "[r4]"},
        SmartsForAtomTypeTestCase{"C1CC2(C1)CCCC2", 2, IWATTYPE_USP_L, "[r5]"},
        SmartsForAtomTypeTestCase{"CC", 0, IWATTYPE_USP_U, "[G0]"},
        SmartsForAtomTypeTestCase{"C=C", 0, IWATTYPE_USP_U, "[G1]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0, IWATTYPE_USP_U, "[G1]"},
        SmartsForAtomTypeTestCase{"C=C", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_B, "[#6;G1]"},
        SmartsForAtomTypeTestCase{"c1ccccc1", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_B, "[#6]"}));

INSTANTIATE_TEST_SUITE_P(
    ChargeAndCompositeTypes, SmartsForAtomTypeTest,
    testing::Values(
        SmartsForAtomTypeTestCase{"[NH4+]", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_O, "[#7;+]"},
        SmartsForAtomTypeTestCase{"[O-]C=O", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_O, "[#8;-]"},
        SmartsForAtomTypeTestCase{"O", 0,
                                  IWATTYPE_USP_Z | IWATTYPE_USP_O, "[#8]"},
        SmartsForAtomTypeTestCase{
            "C1CC1", 0,
            IWATTYPE_USP_Z | IWATTYPE_USP_A | IWATTYPE_USP_C |
                IWATTYPE_USP_R | IWATTYPE_USP_U | IWATTYPE_USP_S |
                IWATTYPE_USP_L | IWATTYPE_USP_O,
            "[#6;D2;x2;G0;r3;r3]"}));

TEST(CagedRingAtomType, TagIncludesModifier) {
  Atom_Typing_Specification spec;
  ASSERT_TRUE(spec.build("UST:ZW"));
  IWString tag;
  ASSERT_TRUE(spec.append_to_tag(tag));
  EXPECT_TRUE(tag.contains('W'));
}

TEST(CagedRingAtomType, ModifierAddsToExistingTypes) {
  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles("NC1CC2CCC1C2"));
  Atom_Typing_Specification base, combined;
  ASSERT_TRUE(base.build("UST:Z"));
  ASSERT_TRUE(combined.build("UST:ZW"));
  std::vector<uint32_t> before(mol.natoms()), after(mol.natoms());
  ASSERT_TRUE(base.assign_atom_types(mol, before.data()));
  ASSERT_TRUE(combined.assign_atom_types(mol, after.data()));
  EXPECT_EQ(after[0], before[0]);
  for (int i = 1; i < mol.natoms(); ++i) {
    EXPECT_EQ(after[i], before[i] + 24737);
  }
}

struct CageCase {
  const char* _smiles;
  int _marked;
};

class CagedRingAtomType : public testing::TestWithParam<CageCase> {};

TEST_P(CagedRingAtomType, RandomAtomOrderPreservesTypes) {
  set_global_aromaticity_type(Daylight);
  const auto& test = GetParam();
  Molecule original;
  ASSERT_TRUE(original.build_from_smiles(test._smiles));
  for (int i = 0; i < original.natoms(); ++i) {
    original.set_isotope(i, i + 1);
  }
  Atom_Typing_Specification spec;
  ASSERT_TRUE(spec.build("UST:W"));
  std::vector<uint32_t> expected(original.natoms());
  ASSERT_TRUE(spec.assign_atom_types(original, expected.data()));
  int marked = 0;
  for (uint32_t value : expected) {
    EXPECT_TRUE(value == 1 || value == 24738);
    marked += value == 24738;
  }
  EXPECT_EQ(marked, test._marked) << test._smiles;

  for (int trial = 0; trial < 200; ++trial) {
    const IWString smiles = original.random_smiles();
    SCOPED_TRACE(smiles.AsString());
    Molecule reordered;
    ASSERT_TRUE(reordered.build_from_smiles(smiles));
    ASSERT_EQ(reordered.natoms(), original.natoms());
    std::vector<uint32_t> actual(reordered.natoms());
    ASSERT_TRUE(spec.assign_atom_types(reordered, actual.data()));
    for (int i = 0; i < reordered.natoms(); ++i) {
      const int index = reordered.isotope(i) - 1;
      ASSERT_GE(index, 0);
      ASSERT_LT(index, original.natoms());
      EXPECT_EQ(actual[i], expected[index]) << "original atom " << index;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    RingTopologies, CagedRingAtomType,
    testing::Values(
        CageCase{"CCC", 0},
        CageCase{"C1CCCCC1", 0},
        CageCase{"C1CCC2CCCCC2C1", 0},
        CageCase{"C1CCC2(CC1)CCCC2", 0},
        CageCase{"c1ccc2ccccc2c1", 0},
        CageCase{"C1CC2CCC1C2", 7},
        CageCase{"C1CCC2CCC1C2", 8},
        CageCase{"C1C2CC3CC1CC(C2)C3", 10},
        CageCase{"C12C3C4C1C5C2C3C45", 0},
        CageCase{"CC1CC2CCC1C2", 7},
        CageCase{"C13C(CCCC3)C2CCC1C2", 7},
        CageCase{"C1CC2CCC1C23CCCCC3", 7},
        CageCase{"C1CC2CCC1C2C3CCCCC3", 7}));

}  // namespace
