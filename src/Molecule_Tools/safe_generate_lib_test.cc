// Tests for safe_generate_library

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Tools/safe_generate_lib.h"

namespace {

using safe_generate::BondsOk;
using safe_generate::SafeFragment;

TEST(TestBondsOk, Test1) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("N#[1C][1CH]1[1CH]([1NH][1CH3])[1CH]([1S][1CH2][1CH]2NC3C(N2)CCCC3)[1CH]([1CH]2CCCCC2)C=C1"));
  EXPECT_TRUE(BondsOk(m));
}

struct SmilesResult {
  IWString smiles;
  int result;
};

class TestBondsOk: public testing::TestWithParam<SmilesResult> {
  protected:
    Molecule _m;
};
TEST_P(TestBondsOk, Tests) {
  auto params = GetParam();
  ASSERT_TRUE(_m.build_from_smiles(params.smiles));
  EXPECT_EQ(BondsOk(_m), params.result) << params.smiles;
}
INSTANTIATE_TEST_SUITE_P(TestBondsOk, TestBondsOk, testing::Values(
  SmilesResult{"[1Cl]-[1CH3]", 0},
  SmilesResult{"[1Br]-[1CH3]", 0},
  SmilesResult{"[1I]-[1CH3]", 0},
  SmilesResult{"[1F]-[1CH3]", 1},
  SmilesResult{"[1F]-[1c]1ccccc1", 1},
  SmilesResult{"[1Cl]-[1c]1ccccc1", 1},
  SmilesResult{"[1Br]-[1c]1ccccc1", 1},
  SmilesResult{"[1I]-[1c]1ccccc1", 1},
  SmilesResult{"C[1NH][1F]", 0},
  SmilesResult{"C[1NH][1Cl]", 0},
  SmilesResult{"C[1NH][1Br]", 0},
  SmilesResult{"C[1NH][1I]", 0},
  SmilesResult{"C[1O][1F]", 0},
  SmilesResult{"C[1O][1Cl]", 0},
  SmilesResult{"C[1O][1Br]", 0},
  SmilesResult{"C[1O][1I]", 0},
  SmilesResult{"C[1O]-[1O]C", 0},
  SmilesResult{"C[1NH]-[1NH]C", 0},
  SmilesResult{"C[1S]-[1S]C", 0},
  SmilesResult{"C[1S]-[1CH2]C", 1},
  SmilesResult{"C[1S](=O)-[1NH]C", 1},
  SmilesResult{"C[1S][1NH]C", 0},
  SmilesResult{"C[1S]-[1P](=O)OC", 0},
  SmilesResult{"C[1O][1O]C", 0},
  // In a ring is not checked
  SmilesResult{"C1[1O][1O]C1", 1},
  SmilesResult{"C[1NH][1NH]C", 0},
  SmilesResult{"C1[1NH][1NH]C1", 1},
  SmilesResult{"N#[1C][1P]1(=O)O[1P](=O)([1OH])OC[1CH]([1OH])[1CH]([1CH2][1OH])O1", 1},
  SmilesResult{"[1Cl][1CH]1[1N]([1C]2=NO[1C]([1CH]3NCCCC3)([1CH3])O2)CCC1", 0}
));



TEST(TestBondsOk, Test2) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("N#CC1C(NC)C(SCC2NC3C(N2)CCCC3)C(C2CCCCC2)C=C1"));
  EXPECT_TRUE(BondsOk(m));
}

TEST(TestSameRingOpenings, Test1) {
  SafeFragment f1, f2;
  ASSERT_TRUE(f1.Build("CC[1N]%12CCC[1C]%13O"));
  EXPECT_EQ(f1.natoms(), 8);
  EXPECT_EQ(f1.nrings(), 0);
  ASSERT_TRUE(f2.Build("[1O]%17C[1N]%18O"));
  EXPECT_EQ(f2.natoms(), 4);
  IWString new_smiles;
  EXPECT_TRUE(f1.SameNumbers(f2, new_smiles));
  EXPECT_EQ(new_smiles, "[1O]%12C[1N]%13O");
}



TEST(TestAttachmentIsotope, SafeInputIgnoresOtherLabels) {
  SafeFragment fragment;
  ASSERT_TRUE(fragment.Build("[13CH3][7CH2]%10CC[7CH2]%11[15NH2]", 7));
  EXPECT_EQ(fragment.ncon(), 2);
  EXPECT_EQ(fragment.distance(), 3);
  EXPECT_EQ(fragment.mol().isotope(0), 13u);
  EXPECT_EQ(fragment.mol().isotope(5), 15u);

  SafeFragment legacy;
  ASSERT_TRUE(legacy.Build("[13CH3][7CH2]%10CC[7CH2]%11[15NH2]"));
  EXPECT_EQ(legacy.distance(), 1);
}

TEST(TestAttachmentIsotope, LibraryOnlyAddsConnectionsAtSelectedIsotope) {
  dicer_data::DicerFragment proto;
  proto.set_smi("[13CH3][7CH2]CC[7CH2][15NH2]");
  proto.set_par("labelled fragment");
  SafeFragment fragment;
  ASSERT_TRUE(fragment.Build(proto, 7));
  EXPECT_EQ(fragment.ncon(), 2);
  EXPECT_EQ(fragment.distance(), 3);
  EXPECT_EQ(fragment.smiles(), "[13CH3][7CH2]%11CC[7CH2]%12[15NH2]");

  SafeFragment legacy;
  ASSERT_TRUE(legacy.Build(proto));
  EXPECT_EQ(legacy.ncon(), 4);
  EXPECT_EQ(legacy.distance(), 1);
}

TEST(TestAttachmentIsotope, MatchWholeNumber) {
  dicer_data::DicerFragment proto;
  proto.set_smi("[13CH3][1CH2][15NH2]");
  SafeFragment fragment;
  ASSERT_TRUE(fragment.Build(proto, 1));
  EXPECT_EQ(fragment.ncon(), 1);
  EXPECT_EQ(fragment.distance(), 0);
  EXPECT_EQ(fragment.smiles(), "[13CH3][1CH2]%11[15NH2]");
}

TEST(TestAttachmentIsotope, RejectMissingMarker) {
  SafeFragment safe;
  EXPECT_FALSE(safe.Build("[13CH3][1CH2]%10C", 7));
  dicer_data::DicerFragment proto;
  proto.set_smi("[13CH3][1CH2]C");
  SafeFragment library;
  EXPECT_FALSE(library.Build(proto, 7));
}

TEST(TestAttachmentIsotope, BondChecksOnlySelectedLabels) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("[13NH2][15NH2]"));
  EXPECT_FALSE(BondsOk(m));
  EXPECT_TRUE(BondsOk(m, 7));
  ASSERT_TRUE(m.build_from_smiles("[7NH2][15NH2]"));
  EXPECT_TRUE(BondsOk(m, 7));
  ASSERT_TRUE(m.build_from_smiles("[7NH2][7NH2]"));
  EXPECT_FALSE(BondsOk(m, 7));
}

TEST(TestAttachmentIsotope, RemoveOnlyAttachmentLabels) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("[13CH3][7CH2][2CH2][15NH2]"));
  safe_generate::RemoveAttachmentIsotopes(m, 7);
  EXPECT_EQ(m.isotope(0), 13u);
  EXPECT_EQ(m.isotope(1), 0u);
  EXPECT_EQ(m.isotope(2), 2u);
  EXPECT_EQ(m.isotope(3), 15u);
  safe_generate::RemoveAttachmentIsotopes(m, 0);
  EXPECT_EQ(m.number_isotopic_atoms(), 0);
}

TEST(TestFragmentMetadata, LibraryRingCountAndSize) {
  dicer_data::DicerFragment proto;
  proto.set_smi("[7CH]1CC1");
  SafeFragment fragment;
  ASSERT_TRUE(fragment.Build(proto, 7));
  EXPECT_EQ(fragment.nrings(), 1);
  EXPECT_EQ(fragment.largest_ring_size(), 3);
}

TEST(TestSafeRingNumbers, RejectsInternalTwoDigitRings) {
  SafeFragment source;
  EXPECT_FALSE(source.Build("[7CH]%10%12CCCCC%10", 7));
  dicer_data::DicerFragment proto;
  proto.set_smi("[7CH]%11CCCCC%11");
  SafeFragment replacement;
  EXPECT_FALSE(replacement.Build(proto, 7));
}

TEST(TestSafeRingNumbers, RejectsMalformedLabels) {
  for (const char* smiles : {"[1C]%", "[1C]%1", "[1C]%1x"}) {
    SafeFragment fragment;
    EXPECT_FALSE(fragment.Build(smiles, 1)) << smiles;
  }
}

TEST(TestSafeRingNumbers, RejectsConnectionCountMismatch) {
  SafeFragment source, replacement;
  ASSERT_TRUE(source.Build("[7CH2]%10C[7CH2]%12", 7));
  ASSERT_TRUE(replacement.Build("[7CH2]%11C", 7));
  IWString smiles;
  EXPECT_FALSE(source.SameNumbers(replacement, smiles));
}

}  // namespace
