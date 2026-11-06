// Tests for fragment_molecule

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Molecule_Tools/fragment_molecule.h"

namespace {


TEST(TestFragmentMolecule, TestNoBreakCC) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CC"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 0);
  EXPECT_EQ(bonds.size(), 0u);
}

TEST(TestFragmentMolecule, TestBreakCC) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CCC(C)C"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 0);
  fragmenter.set_break_carbon_carbon_bonds(1);
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 4);
  EXPECT_EQ(bonds.size(), 4u);
}

TEST(TestFragmentMolecule, TestCAmide) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CNC(=O)C"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 2);
  EXPECT_EQ(bonds.size(), 2u);
}

TEST(TestFragmentMolecule, TestSAmide1) {
  Molecule m;
  EXPECT_TRUE(m.build_from_smiles("CNS(=O)C"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 2);
  EXPECT_EQ(bonds.size(), 2u);
}

TEST(TestFragmentMolecule, IsotopeProtectsOnlyBondsWithBothEndpointsLabelled) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C[7NH][7CH2]O"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  fragmenter.set_preserve_isotope(7);
  resizable_array<int> bonds;
  ASSERT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 2);
  EXPECT_THAT(bonds, testing::UnorderedElementsAre(m.which_bond(0, 1), m.which_bond(2, 3)));
  EXPECT_EQ(fragmenter.bonds_protected_by_isotope(), 1u);

  bonds.resize_keep_storage(0);
  Molecule different;
  ASSERT_TRUE(different.build_from_smiles("C[7NH][13CH2]O"));
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(different, bonds), 3);
  EXPECT_EQ(fragmenter.bonds_protected_by_isotope(), 0u);
}

TEST(TestFragmentMolecule, IsotopeProtectionIsDisabledByDefault) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C[7NH][7CH2]O"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 3);
  EXPECT_EQ(fragmenter.bonds_protected_by_isotope(), 0u);
}

TEST(TestFragmentMolecule, IsotopeCanProtectEveryCandidate) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("[7CH3][7OH]"));
  fragment_molecule::MoleculeFragmenter fragmenter;
  fragmenter.set_preserve_isotope(7);
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 0);
  EXPECT_EQ(fragmenter.bonds_protected_by_isotope(), 1u);
}

TEST(TestFragmentMolecule, ExclusionQueriesSupplementDefaults) {
  const char* argv[] = {"fragmenter", "-Q", "SMARTS:COC"};
  Command_Line cl(3, const_cast<char**>(argv), "Q:");
  fragment_molecule::MoleculeFragmenter fragmenter;
  ASSERT_TRUE(fragmenter.Initialise(cl));
  fragmenter.set_add_queries_to_default_rules(1);
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("COCN"));
  resizable_array<int> bonds;
  ASSERT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 1);
  EXPECT_EQ(bonds[0], m.which_bond(2, 3));
}

TEST(TestFragmentMolecule, ExclusionOverridesCf3Rule) {
  const char* argv[] = {"fragmenter", "-Q", "SMARTS:FC(F)(F)C"};
  Command_Line cl(3, const_cast<char**>(argv), "Q:");
  fragment_molecule::MoleculeFragmenter fragmenter;
  ASSERT_TRUE(fragmenter.Initialise(cl));
  fragmenter.set_add_queries_to_default_rules(1);
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("FC(F)(F)CO"));
  resizable_array<int> bonds;
  ASSERT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 1);
  EXPECT_EQ(bonds[0], m.which_bond(4, 5));
}

TEST(TestFragmentMolecule, IsotopeExcludesExplicitQueryCuts) {
  const char* argv[] = {"fragmenter", "-q", "SMARTS:CC"};
  Command_Line cl(3, const_cast<char**>(argv), "q:");
  fragment_molecule::MoleculeFragmenter fragmenter;
  ASSERT_TRUE(fragmenter.Initialise(cl));
  fragmenter.set_preserve_isotope(7);
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("[7CH3][7CH3]"));
  resizable_array<int> bonds;
  EXPECT_EQ(fragmenter.IdentifyBreakableBonds(m, bonds), 0);
  EXPECT_EQ(fragmenter.bonds_protected_by_isotope(), 1u);
}

TEST(TestFragmentMolecule, NoQueryOptionsIsValid) {
  const char* argv[] = {"fragmenter"};
  Command_Line cl(1, const_cast<char**>(argv), "q:Q:");
  fragment_molecule::MoleculeFragmenter fragmenter;
  EXPECT_TRUE(fragmenter.Initialise(cl));
}

}  // namespace
