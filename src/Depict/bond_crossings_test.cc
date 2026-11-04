#include "Depict/bond_crossings.h"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "Molecule_Lib/molecule.h"
#include "gtest/gtest.h"

namespace {
using bond_crossings::Analyse;
using bond_crossings::Options;
using bond_crossings::Result;
using bond_crossings::Status;
using XY = std::array<double, 2>;
using Edge = std::array<int, 2>;

void Build(Molecule& mol, const std::vector<XY>& points, const std::vector<Edge>& edges) {
  for (const auto& p : points) {
    const int atom = mol.add(new Atom(6));
    mol.setxyz(atom, p[0], p[1], 0);
  }
  for (const auto& edge : edges) {
    ASSERT_TRUE(mol.add_bond(edge[0], edge[1], SINGLE_BOND));
  }
}

const std::vector<Edge> kSquare = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
const std::vector<XY> kBowtie = {{1, 1}, {3, 3}, {1, 3}, {3, 1}};

TEST(BondCrossings, CleanRingAndDistances) {
  Molecule mol;
  Build(mol, {{1, 1}, {3, 1}, {3, 3}, {1, 3}}, kSquare);
  Result result;
  ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
  EXPECT_EQ(result.crossings(), 0);
  EXPECT_DOUBLE_EQ(result.score(), 0);
  ASSERT_EQ(result.systems().size(), 1u);
  EXPECT_EQ(result.ring_atoms(), 4);
  EXPECT_DOUBLE_EQ(result.min_distance(), 2);
  EXPECT_DOUBLE_EQ(result.min_distance_per_atom(), 0.5);
  EXPECT_DOUBLE_EQ(result.max_distance(), std::sqrt(8.0));
  EXPECT_DOUBLE_EQ(result.max_distance_per_atom(), std::sqrt(8.0) / 4);
  EXPECT_EQ(result.topological_diameter(), 2);
  EXPECT_DOUBLE_EQ(result.extent_per_topology(), std::sqrt(2.0));
}

TEST(BondCrossings, OneCrossingCountsOnePairAndTwoBonds) {
  Molecule mol;
  Build(mol, kBowtie, kSquare);
  Result result;
  ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
  EXPECT_EQ(result.crossings(), 1);
  EXPECT_EQ(result.ring_ring(), 1);
  EXPECT_EQ(result.ring_chain(), 0);
  EXPECT_EQ(result.crossed_bonds(), 2);
  EXPECT_DOUBLE_EQ(result.score(), 1);
  EXPECT_DOUBLE_EQ(result.score_per_atom(), 0.25);
  EXPECT_DOUBLE_EQ(result.worst_density(), 0.25);
}

TEST(BondCrossings, RingChainCongestionIsCubicAndChainChainIsIgnored) {
  const std::vector<XY> positions = {{0, 0}, {4, 0}, {4, 4}, {0, 4}, {2, 6}, {1, 6}, {3, 6}};
  for (int count = 1; count <= 3; ++count) {
    Molecule mol;
    std::vector<XY> points(positions.begin(), positions.begin() + 4 + count);
    auto edges = kSquare;
    edges.push_back({0, 4});
    if (count >= 2) edges.push_back({1, 5});
    if (count >= 3) edges.push_back({0, 6});
    Build(mol, points, edges);
    Result result;
    ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
    EXPECT_EQ(result.ring_ring(), 0);
    EXPECT_EQ(result.ring_chain(), count);
    EXPECT_EQ(result.crossed_bonds(), count + 1);
    EXPECT_DOUBLE_EQ(result.score(), count * count * count);
    EXPECT_DOUBLE_EQ(result.worst_density(), count * count * count / 4.0);
  }
}

TEST(BondCrossings, AcyclicCrossingsNeverScore) {
  Molecule mol;
  Build(mol, kBowtie, {{0, 1}, {1, 2}, {2, 3}});
  Options options;
  options.set_diagnostics(true);
  Result result;
  ASSERT_EQ(Analyse(mol, options, result), Status::kOk);
  EXPECT_TRUE(result.systems().empty());
  EXPECT_EQ(result.crossings(), 0);
  EXPECT_EQ(result.overlaps(), 0);
  EXPECT_EQ(result.near_misses(), 0);
  EXPECT_DOUBLE_EQ(result.score(), 0);
  EXPECT_DOUBLE_EQ(result.score_per_system(), 0);
}

TEST(BondCrossings, CrossingsInSeparateSystemsAreLessSevere) {
  Molecule mol;
  auto points = kBowtie;
  for (const auto& p : kBowtie) points.push_back({p[0] + 10, p[1]});
  auto edges = kSquare;
  for (const auto& e : kSquare) edges.push_back({e[0] + 4, e[1] + 4});
  edges.push_back({3, 4});
  Build(mol, points, edges);
  Result result;
  ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
  ASSERT_EQ(result.systems().size(), 2u);
  EXPECT_EQ(result.crossings(), 2);
  EXPECT_EQ(result.systems_crossed(), 2);
  EXPECT_EQ(result.max_crossings(), 1);
  EXPECT_EQ(result.inter_system(), 0);
  EXPECT_DOUBLE_EQ(result.score(), 2);
}

TEST(BondCrossings, SpiroSystemsStaySeparateWithoutLosingBondsAtSharedAtom) {
  Molecule mol;
  Build(mol, {{1, 1}, {3, 3}, {1, 3}, {3, 1}, {-1, -1}, {1, -1}, {-1, 1}},
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {4, 5}, {5, 6}, {6, 0}});
  Result result;
  ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
  ASSERT_EQ(result.systems().size(), 2u);
  EXPECT_EQ(result.ring_atoms(), 7);
  for (const auto& system : result.systems()) {
    EXPECT_EQ(system.atoms(), 4);
    EXPECT_EQ(system.ring_ring(), 1);
  }
  EXPECT_DOUBLE_EQ(result.score(), 2);
  EXPECT_EQ(result.inter_system(), 0);
}

TEST(BondCrossings, FusedAndBridgedRingsBelongToOneSystem) {
  for (const char* smiles : {"C1CCC2CCCCC2C1", "C1CC2CCC1C2", "C12C3C4C1C5C2C3C45"}) {
    Molecule mol;
    ASSERT_TRUE(mol.build_from_smiles(smiles));
    for (int i = 0; i < mol.natoms(); ++i) {
      mol.setxyz(i, 2 + std::cos(i * 0.7), 2 + std::sin(i * 0.7), 0);
    }
    Result result;
    ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
    ASSERT_EQ(result.systems().size(), 1u);
    EXPECT_GE(result.systems()[0].rings(), 2);
    EXPECT_EQ(result.ring_atoms(), mol.natoms());
    EXPECT_EQ(result.ring_chain(), 0);
  }
}

TEST(BondCrossings, InterSystemCrossingsAreCountedOnceWithLinearPenalty) {
  Molecule mol;
  Build(mol, {{1, 1}, {5, 1}, {3, 5}, {1, 4}, {5, 4}, {3, 0}},
        {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}});
  Result result;
  ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
  EXPECT_EQ(result.inter_system(), 6);
  EXPECT_EQ(result.ring_ring(), 0);
  EXPECT_EQ(result.ring_chain(), 0);
  EXPECT_DOUBLE_EQ(result.score(), 6);
}

TEST(BondCrossings, RigidTransformsAndScalePreserveCounts) {
  Molecule mol;
  Build(mol, kBowtie, kSquare);
  Result original, transformed;
  ASSERT_EQ(Analyse(mol, Options(), original), Status::kOk);
  for (int i = 0; i < mol.natoms(); ++i) {
    const double x = mol.x(i), y = mol.y(i);
    mol.setxyz(i, 10 + 2.5 * (x * std::cos(0.37) - y * std::sin(0.37)),
                  -7 + 2.5 * (x * std::sin(0.37) + y * std::cos(0.37)), 0);
  }
  ASSERT_EQ(Analyse(mol, Options(), transformed), Status::kOk);
  EXPECT_EQ(transformed.crossings(), original.crossings());
  EXPECT_DOUBLE_EQ(transformed.score(), original.score());
  EXPECT_NEAR(transformed.min_distance(), original.min_distance() * 2.5, 1e-5);
}

TEST(BondCrossings, DiagnosticsDoNotChangePrimaryScore) {
  Options options;
  options.set_diagnostics(true);
  for (double y : {-0.01, 0.0}) {
    Molecule mol;
    auto points = std::vector<XY>{{0, 0}, {4, 0}, {4, 4}, {0, 4}, {-1, y}, {5, y}};
    auto edges = kSquare;
    edges.push_back({0, 4});
    edges.push_back({4, 5});
    Build(mol, points, edges);
    Result basic, diagnostic;
    ASSERT_EQ(Analyse(mol, Options(), basic), Status::kOk);
    ASSERT_EQ(Analyse(mol, options, diagnostic), Status::kOk);
    EXPECT_EQ(diagnostic.crossings(), 0);
    EXPECT_DOUBLE_EQ(diagnostic.score(), basic.score());
    if (y == 0) EXPECT_GT(diagnostic.overlaps(), 0);
    else EXPECT_GT(diagnostic.near_misses(), 0);
    EXPECT_GT(diagnostic.min_nonbonded(), 0);
  }
}

TEST(BondCrossings, AtomsOnRingBondsAndDegeneraciesAreOptional) {
  Options options;
  options.set_diagnostics(true);
  for (const XY& point : std::vector<XY>{{2, 0}, {0, 0}}) {
    Molecule mol;
    Build(mol, {{0, 0}, {4, 0}, {4, 4}, {0, 4}, point},
          {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}});
    Result result;
    ASSERT_EQ(Analyse(mol, options, result), Status::kOk);
    EXPECT_EQ(result.crossings(), 0);
    if (point[0] == 2) EXPECT_EQ(result.atoms_on_ring_bonds(), 1);
    else {
      EXPECT_EQ(result.degenerate_bonds(), 1);
      EXPECT_EQ(result.coincident_atoms(), 1);
    }
  }
}

TEST(BondCrossings, InvalidInputsAndResultReuse) {
  Molecule mol;
  Build(mol, kBowtie, kSquare);
  Result result;
  ASSERT_EQ(Analyse(mol, Options(), result), Status::kOk);
  EXPECT_GT(result.score(), 0);
  mol.setz(0, 1);
  EXPECT_EQ(Analyse(mol, Options(), result), Status::kNot2D);
  EXPECT_EQ(result.crossings(), 0);
  EXPECT_EQ(result.atoms(), 0);
  mol.setxyz(0, std::numeric_limits<float>::quiet_NaN(), 1, 0);
  EXPECT_EQ(Analyse(mol, Options(), result), Status::kInvalidCoordinates);
  Molecule disconnected;
  Build(disconnected, {{1, 1}, {2, 2}}, {});
  EXPECT_EQ(Analyse(disconnected, Options(), result), Status::kDisconnected);
  Molecule empty;
  EXPECT_EQ(Analyse(empty, Options(), result), Status::kEmpty);
  Molecule line;
  Build(line, {{0, 0}, {1, 0}}, {{0, 1}});
  EXPECT_EQ(Analyse(line, Options(), result), Status::kNot2D);
  Options invalid;
  invalid.set_near_miss_fraction(-1);
  EXPECT_EQ(Analyse(line, invalid, result), Status::kInvalidOptions);
}

}  // namespace
