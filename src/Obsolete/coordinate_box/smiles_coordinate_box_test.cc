// Archived SMILES integration tests; not part of the active build.
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "coordinate_box.h"
#include "Molecule_Lib/smiles.h"

namespace {
using testing::ElementsAre;
using testing::FloatNear;

class TestBoxedCoordinates : public testing::Test {
    void SetUp() {
      set_append_coordinate_box_after_each_atom(1);
    }
    void TearDown() {
      set_append_coordinate_box_after_each_atom(0);
    }
};

TEST_F(TestBoxedCoordinates, TestBoxedCoordinates0) {
  const_IWSubstring smiles("C{{B0:0}}");
  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));
  EXPECT_EQ(m.x(0), 0.0f);
  EXPECT_EQ(m.y(0), 0.0f);
  EXPECT_EQ(m.z(0), 0.0f);
}

TEST_F(TestBoxedCoordinates, TestBoxedCoordinates1) {
  Coordinates coords(1.0f, 2.0f, 0.5f);
  coordinate_box::ConcentricBox box;
  const coordinate_box::LayerPosition layer_position = box.Position(coords);

  IWString smiles("CC");
  smiles << "{{B" << layer_position << "}}";

  constexpr float kAbsDiff = 0.001;

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));
  EXPECT_THAT(m.atom(1).ToResizableArray(),
              ElementsAre(FloatNear(coords.x(), kAbsDiff),
              FloatNear(coords.y(), kAbsDiff),
              FloatNear(coords.z(), kAbsDiff)));
}

TEST_F(TestBoxedCoordinates, TestBoxedCoordinates2) {
  constexpr int matoms = 10;
  std::vector<Coordinates> coords(matoms);
  for (int i = 0; i < matoms; ++i) {
    coords[i].setxyz(i, 2 * i, i + 8);
  }

  coordinate_box::ConcentricBox box;

  IWString smiles;
  for (int i = 0; i < matoms; ++i) {
    smiles << 'C';
    const coordinate_box::LayerPosition layer_position = box.Position(coords[i]);
    smiles << "{{B" << layer_position << "}}";
  }

  Molecule m;
  ASSERT_TRUE(m.build_from_smiles(smiles));

  constexpr float kAbsDiff = 0.001;

  ASSERT_EQ(m.natoms(), matoms);

  for (int i = 0; i < matoms; ++i) {
    EXPECT_THAT(m.atom(i).ToResizableArray(),
                ElementsAre(FloatNear(coords[i].x(), kAbsDiff),
                            FloatNear(coords[i].y(), kAbsDiff),
                            FloatNear(coords[i].z(), kAbsDiff)));
  }
}

}  // namespace
