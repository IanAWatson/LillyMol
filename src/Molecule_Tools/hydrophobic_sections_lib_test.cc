#include "Molecule_Tools/hydrophobic_sections_lib.h"

#include <vector>

#include "Molecule_Lib/molecule.h"
#include "gtest/gtest.h"

namespace {

std::vector<int>
CarbonMask(const Molecule& m) {
  std::vector<int> result(m.natoms());
  for (int atom = 0; atom < m.natoms(); ++atom) {
    result[atom] = m.atomic_number(atom) == 6;
  }
  return result;
}

TEST(HydrophobicSections, SizeThresholdAndConsecutiveLabels) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("COCCOCCC"));
  const auto mask = CarbonMask(m);
  hydrophobic::Sections sections;
  ASSERT_TRUE(sections.Build(m, mask.data(), 2));
  EXPECT_EQ(sections.sizes(), (std::vector<int>{2, 3}));
  EXPECT_EQ(sections.membership(), (std::vector<int>{0, 0, 1, 1, 0, 2, 2, 2}));
  EXPECT_EQ(sections.atoms_in_sections(), 5);
  EXPECT_EQ(sections.largest_section(), 3);
  EXPECT_EQ(mask, CarbonMask(m));
}

TEST(HydrophobicSections, RingsAndDisconnectedComponents) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C1CCCCC1.O.CC"));
  const auto mask = CarbonMask(m);
  hydrophobic::Sections sections;
  ASSERT_TRUE(sections.Build(m, mask.data()));
  EXPECT_EQ(sections.sizes(), (std::vector<int>{6, 2}));
  EXPECT_EQ(sections.number_sections(), 2);
  EXPECT_EQ(sections.atoms_in_sections(), 8);
  ASSERT_TRUE(sections.Build(m, sections.membership().data(), 3));
  EXPECT_EQ(sections.sizes(), (std::vector<int>{6}));
  EXPECT_EQ(sections.atoms_in_sections(), 6);
}

TEST(HydrophobicSections, EmptyAndAllExcludedResetResults) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("CCC"));
  const auto mask = CarbonMask(m);
  hydrophobic::Sections sections;
  ASSERT_TRUE(sections.Build(m, mask.data()));
  ASSERT_TRUE(sections.Build(m, mask.data(), 4));
  EXPECT_EQ(sections.number_sections(), 0);
  EXPECT_EQ(sections.atoms_in_sections(), 0);
  EXPECT_EQ(sections.largest_section(), 0);
  EXPECT_EQ(sections.membership(), (std::vector<int>{0, 0, 0}));
  Molecule empty;
  ASSERT_TRUE(sections.Build(empty, nullptr));
  EXPECT_TRUE(sections.membership().empty());
}

TEST(HydrophobicSections, InvalidArguments) {
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles("C"));
  hydrophobic::Sections sections;
  const int mask[] = {1};
  EXPECT_FALSE(sections.Build(m, nullptr));
  EXPECT_FALSE(sections.Build(m, mask, 0));
}

}  // namespace
