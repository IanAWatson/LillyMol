// Tests for functions that happen to reside in ss_ring_base.
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

#include "molecule.h"
#include "substructure.h"
#include "target.h"

namespace {

using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;

class TestSubstituent : public testing::Test
{
  protected:

  protected:
    std::string _string_proto;
    SubstructureSearch::SubstructureQuery _proto;
    IWString _smiles;
    Molecule _mol;
    Substructure_Query _query;
    Substructure_Results _sresults;
};

TEST_F(TestSubstituent, TestAtomCountMatches) {
  _string_proto = R"pb(
query {
  ring_specifier {
    aromatic: false
    base {
      substituent {
        natoms: 1;
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_EQ(_sresults.embedding(0)->item(0), 3);
}

TEST_F(TestSubstituent, TestAtomCountNoMatch) {
  _string_proto = R"pb(
query {
  ring_specifier {
    aromatic: false
    base {
      substituent {
        natoms: 2;
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, TestAtomCountMultipleMatches) {
  _string_proto = R"pb(
query {
  ring_specifier {
    aromatic: false
    base {
      substituent {
        natoms: 1;
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CC1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 2u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(0));
  EXPECT_THAT(*_sresults.embedding(1), UnorderedElementsAre(4));
}

TEST_F(TestSubstituent, TestSSMatches) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_specifier {
    aromatic: false
    base {
      substituent {
        required_smarts: "[CD1]";
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CCC1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({0, 1, 5}));
}

TEST_F(TestSubstituent, TestSSMatchesNumericQualifier) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_specifier {
    aromatic: false
    base {
      substituent {
        required_smarts: "2[CD1]";
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CC(C)C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({0, 1, 2}));
}

TEST_F(TestSubstituent, TestSSDisqualifying) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_specifier {
    aromatic: false
    base {
      substituent {
        disqualifying_smarts: "O";
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "CC(C)C1CC1C";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({0, 1, 2, 6}));
}

TEST_F(TestSubstituent, TestNrings0) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 0
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CC)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({8, 9}));
}

TEST_F(TestSubstituent, TestNrings0NoMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 0
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}


TEST_F(TestSubstituent, TestNrings1Match) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({8, 9, 10, 11, 12}));
}

TEST_F(TestSubstituent, TestNrings1NoMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 2
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, TestNrings1NoMatchTooManyRings) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c12ccccc1cc(CCC1CC1C1CC1)cc2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, TestNrings1NoSubstituents) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
      }
    }
  }
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "S1N=NC2=C1SC=C2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, TestNringsNoMatchNoRing) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 2
    rings_in_system: 2
    base {
      substituent {
        nrings: 1
      }
    }
  }
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "N1=C2C=NNC2=C(C)N1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, TestNringsSeparateRings) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 1
    rings_in_system: 1
    base {
      substituent {
        nrings: 2
        required_smarts: "[/IWfsid1].[/IWfsid2]"
        set_global_id: 2
      }
    }
  }
  smarts: "[/IWgid2]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c1cncc1C(C1CC1)(C1CC1)";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({5, 6, 7, 8, 9, 10, 11}));
}
TEST_F(TestSubstituent, TestNringsFusedRingsShouldNotMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_system_specifier {
    aromatic_ring_count: 1
    rings_in_system: 1
    base {
      substituent {
        nrings: 2
        required_smarts: "[/IWfsid1R].[/IWfsid2R]"
        set_global_id: 2
      }
    }
  }
  smarts: "[/IWgid2]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "c1cncc1CC12CC1C2";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, InterRingRegion1) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    natoms: 1
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3}));
}

TEST_F(TestSubstituent, InterRingRegion1NoMatch) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    natoms: 1
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1CCC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, InterRingRegionOkLength2) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    length: 2
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3}));
}

TEST_F(TestSubstituent, InterRingRegionOkLength3) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    length: 3
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3, 7}));
}

TEST_F(TestSubstituent, InterRingRegionOkRequiredLength3) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    required_length: [2, 3, 3]
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3, 7}));
}

TEST_F(TestSubstituent, InterRingRegionRequiredLengthNotMatched) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    required_length: [2, 3, 4]
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, InterRingRegionRingConnections) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)CC1CC1CCCC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({11, 12, 13}));
}

TEST_F(TestSubstituent, InterRingRegionRingHitsNeededOk) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 2
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({10, 14}));
}

TEST_F(TestSubstituent, InterRingRegionRingHitsNeededTooSmall) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 1
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, InterRingRegionRingHitsNeededTooLarge) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 3
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 0u);
}

TEST_F(TestSubstituent, InterRingRegionRingRequiredSmarts) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    ring_connections: 2
    hits_needed: 2
    required_smarts: "[R]!@[AD2]!@[R]"
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({10, 14}));
}

TEST_F(TestSubstituent, InterRingRegionRingDisqualifyingSmarts) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    disqualifying_smarts: "[R]!@[AD2]!@[R]"
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();
  _smiles = "C1CC1C(C1CC1)C1CC1CC1CC1CC1CC1";
  ASSERT_TRUE(_mol.build_from_smiles(_smiles));
  EXPECT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAreArray({3}));
}

struct RegionSelectionCase {
  const char* smiles;
  const char* criteria;
  std::vector<int> expected;
};

class TestInterRingSelection : public testing::TestWithParam<RegionSelectionCase> {};

TEST_P(TestInterRingSelection, SelectsExpectedAtoms) {
  const auto& params = GetParam();
  SubstructureSearch::SubstructureQuery proto;
  const std::string text = std::string(
      "query { compress_embeddings: true inter_ring_region { set_global_id: 1 ") +
      params.criteria + " } smarts: \"[/IWgid1]\" }";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(text, &proto));
  Substructure_Query query;
  ASSERT_TRUE(query.ConstructFromProto(proto));
  Molecule mol;
  ASSERT_TRUE(mol.build_from_smiles(params.smiles));
  Substructure_Results results;
  const int nhits = query.substructure_search(&mol, results);
  if (params.expected.empty()) {
    EXPECT_EQ(nhits, 0);
    return;
  }
  ASSERT_EQ(nhits, 1);
  EXPECT_THAT(*results.embedding(0), UnorderedElementsAreArray(params.expected));

  // Preserve original atom identities across random SMILES orderings.
  for (int i = 0; i < mol.natoms(); ++i) {
    mol.set_isotope(i, i + 1);
  }
  for (int iteration = 0; iteration < 5; ++iteration) {
    Molecule reordered;
    ASSERT_TRUE(reordered.build_from_smiles(mol.random_smiles()));
    Substructure_Results reordered_results;
    ASSERT_EQ(query.substructure_search(&reordered, reordered_results), 1u);
    std::vector<int> original_atoms;
    for (atom_number_t atom : *reordered_results.embedding(0)) {
      original_atoms.push_back(reordered.isotope(atom) - 1);
    }
    EXPECT_THAT(original_atoms, UnorderedElementsAreArray(params.expected));
  }
}

INSTANTIATE_TEST_SUITE_P(InterRingSelection, TestInterRingSelection, testing::Values(
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1", "", {3, 4, 8, 9, 10, 11}},
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1",
                        "selection: MAX_ATOMS", {8, 9, 10, 11}},
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1",
                        "selection: LONGEST_LENGTH", {8, 9, 10, 11}},
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1",
                        "selection: MAX_ATOMS max_natoms: 3", {3, 4}},
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1",
                        "selection: MAX_ATOMS min_natoms: 5", {}},
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1",
                        "selection: MAX_ATOMS hits_needed: 2", {8, 9, 10, 11}},
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1",
                        "selection: MAX_ATOMS hits_needed: 1", {}},
    RegionSelectionCase{"C1CC1CCC1CC1CCC1CC1",
                        "selection: MAX_ATOMS", {3, 4, 8, 9}},
    RegionSelectionCase{"C1CC1CCC1CC1CCC1CC1",
                        "selection: LONGEST_LENGTH", {3, 4, 8, 9}},
    RegionSelectionCase{"C1CC1NCC1CC1CCCCC1CC1",
                        "selection: MAX_ATOMS required_smarts: \"[N]\"", {3, 4}},
    RegionSelectionCase{"C1CC1CCC1CC1CCCCC1CC1",
                        "selection: LONGEST_LENGTH max_length: 3", {3, 4}},
    RegionSelectionCase{"C1CC1CCC1CC1", "selection: MAX_ATOMS", {3, 4}},
    // A pendant carbon must not count towards region size or receive a global id.
    RegionSelectionCase{"C1CC1C(C)CC1CC1",
                        "selection: MAX_ATOMS natoms: 2 hits_needed: 1", {3, 5}},
    RegionSelectionCase{"C1CC1CCC1CC1", "selection: ALL hits_needed: 1", {3, 4}}
));

TEST_F(TestSubstituent, InterRingSelectionIsPerSpecification) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    max_natoms: 2
    selection: MAX_ATOMS
    set_global_id: 1
  }
  inter_ring_region {
    min_natoms: 3
    selection: LONGEST_LENGTH
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles("C1CC1CCC1CC1CCCCC1CC1"));
  ASSERT_EQ(_query.substructure_search(&_mol, _sresults), 1u);
  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(3, 4, 8, 9, 10, 11));
}

TEST_F(TestSubstituent, InterRingSelectionSerializes) {
  SubstructureSearch::InterRingAtoms proto;
  proto.set_selection(SubstructureSearch::InterRingAtoms::LONGEST_LENGTH);
  proto.set_min_natoms(3);
  proto.set_set_global_id(7);
  InterRingAtoms specification;
  ASSERT_TRUE(specification.ConstructFromProto(proto));
  SubstructureSearch::InterRingAtoms output;
  // BuildProto currently reports incomplete support for other legacy fields.
  specification.BuildProto(output);
  EXPECT_EQ(output.selection(), SubstructureSearch::InterRingAtoms::LONGEST_LENGTH);
  EXPECT_EQ(output.min_natoms(), 3u);
  EXPECT_EQ(output.set_global_id(), 7u);
}

TEST_F(TestSubstituent, InterRingSelectionCriteriaDiffer) {
  // The branched region has 7 atoms and length 6; the linear region has
  // 6 atoms and length 7. Ranking should choose different components.
  ASSERT_TRUE(_mol.build_from_smiles(
      "C(CCC1CC1)(CCC1CC1)CCC1CC1.C1CC1CCCCCCC1CC1"));
  for (const auto selection : {SubstructureSearch::InterRingAtoms::MAX_ATOMS,
                               SubstructureSearch::InterRingAtoms::LONGEST_LENGTH}) {
    SubstructureSearch::SubstructureQuery proto;
    auto* query = proto.add_query();
    query->set_smarts("[/IWgid1]");
    query->set_compress_embeddings(true);
    auto* region = query->add_inter_ring_region();
    region->set_set_global_id(1);
    region->set_selection(selection);
    Substructure_Query q;
    ASSERT_TRUE(q.ConstructFromProto(proto));
    Substructure_Results results;
    ASSERT_EQ(q.substructure_search(&_mol, results), 1u);
    const bool max_atoms = selection == SubstructureSearch::InterRingAtoms::MAX_ATOMS;
    EXPECT_EQ(results.embedding(0)->number_elements(), max_atoms ? 7 : 6);
    for (atom_number_t atom : *results.embedding(0)) {
      EXPECT_EQ(atom < 16, max_atoms);
    }
  }
}

struct ProtoSmilesExpected {
  std::string proto;
  IWString smiles;
  uint32_t expected;
};

std::ostream&
operator<<(std::ostream& output, const ProtoSmilesExpected& pse) {
  output << pse.proto << ' ' << pse.smiles << " expect " << pse.expected;

  return output;
}

class TestSubstructureP: public testing::TestWithParam<ProtoSmilesExpected> {
  protected:
    SubstructureSearch::SubstructureQuery _proto;
    Substructure_Query _query;
    Molecule _mol;
    Substructure_Results _sresults;
};

TEST_P(TestSubstructureP, Tests) {
  const auto params = GetParam();
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(params.proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles(params.smiles));
  // cerr << "Testing " << params.smiles << " expecting " << params.expected << '\n';
  EXPECT_EQ(_query.substructure_search(_mol, _sresults), params.expected) << params;
}
INSTANTIATE_TEST_SUITE_P(TestSubstructureP, TestSubstructureP, testing::Values(
  ProtoSmilesExpected{R"pb(
query {
  ring_specifier {
    base {
      substituent {
        natoms: 1;
        heteroatom_count: 0
        set_global_id: 1
        hits_needed: 2
      }
    }
  }
  smarts: "[/IWgid1C]"
}
)pb", "C1CC1(C)C", 2},

  ProtoSmilesExpected{R"pb(
query {
  substituent {
    max_natoms: 3
    heteroatom_count: 0
    no_other_substituents_allowed: true
  }
  smarts: "[CD3T2](=O)-[ND2]"
}
)pb", "NCCC(=O)NCCC", 0}

));


}  // namespace
