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

TEST_F(TestSubstituent, SubsetMaskDoesNotAlterOuterTargetOrMolecule) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  ring_specifier {
    base {
      substituent {
        required_smarts: "2[CD1]"
        disqualifying_smarts: "O"
        set_global_id: 1
      }
    }
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles("[13CH3]C(C)C1CC1C[C@H](O)F"));
  const IWString original = _mol.unique_smiles();
  ASSERT_EQ(_mol.chiral_centres(), 1);
  Molecule_to_Match target(&_mol);
  Substructure_Query oxygen_query;
  ASSERT_TRUE(oxygen_query.create_from_smarts("O-C-F"));
  for (int i = 0; i < 10; ++i) {
    Substructure_Results results;
    ASSERT_EQ(_query.substructure_search(target, results), 1u);
    EXPECT_THAT(*results.embedding(0), UnorderedElementsAre(0, 1, 2));
    // The atoms masked in the inner search remain visible to the outer target.
    Substructure_Results oxygen_results;
    EXPECT_EQ(oxygen_query.substructure_search(target, oxygen_results), 1u);
    EXPECT_EQ(_mol.unique_smiles(), original);
    EXPECT_EQ(_mol.chiral_centres(), 1);
    EXPECT_EQ(_mol.isotope(0), 13u);
  }
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

TEST_F(TestSubstituent, DisqualifyingSmartsRejectsOnlyTheMatchingSubstituent) {
  struct TestCase {
    const char* smiles;
    std::vector<int> expected;
  };
  const std::vector<TestCase> cases = {
      {"OCC1CC1", {}},          // The only substituent contains O.
      {"NCC1CC1", {}},          // The first disqualifying query matches.
      {"CCC1CC1", {0, 1}},      // Neither disqualifying query matches.
      {"OCC1CC1CC", {5, 6}},    // Reject one substituent, retain the other.
      {"C1OC1CC", {3, 4}},      // O in the ring must not reject its ethyl group.
  };
  for (const char* specification : {"ring_specifier", "ring_system_specifier"}) {
    for (bool required : {false, true}) {
      std::string text = "query { compress_embeddings: true ";
      text += specification;
      text += R"pb( { base { substituent {
        disqualifying_smarts: "N"
        disqualifying_smarts: "O"
        set_global_id: 1
      )pb";
      if (required) {
        text += "required_smarts: \"C\" ";
      }
      text += R"pb( } } } smarts: "[/IWgid1]" })pb";
      for (const TestCase& test : cases) {
        SCOPED_TRACE(test.smiles);
        SCOPED_TRACE(specification);
        SCOPED_TRACE(required);
        SubstructureSearch::SubstructureQuery proto;
        ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(text, &proto));
        Substructure_Query query;
        ASSERT_TRUE(query.ConstructFromProto(proto));
        Molecule molecule;
        ASSERT_TRUE(molecule.build_from_smiles(test.smiles));
        ASSERT_TRUE(molecule.valence_ok());
        Substructure_Results results;
        EXPECT_EQ(query.substructure_search(molecule, results),
                  test.expected.empty() ? 0u : 1u);
        if (!test.expected.empty()) {
          ASSERT_EQ(results.number_embeddings(), 1);
          EXPECT_THAT(*results.embedding(0), UnorderedElementsAreArray(test.expected));
        }
      }
    }
  }
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

TEST_F(TestSubstituent, InterRingSubsetMaskDoesNotAlterOuterTargetOrMolecule) {
  _string_proto = R"pb(
query {
  compress_embeddings: true
  inter_ring_region {
    # Count the two three-connected ring attachment carbons in this region.
    required_smarts: "2[CD3]"
    disqualifying_smarts: "O"
    set_global_id: 1
  }
  smarts: "[/IWgid1]"
}
)pb";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_mol.build_from_smiles("C1CC1[13CH2]C2CC2OC3CC3[C@H](O)F"));
  ASSERT_TRUE(_mol.valence_ok());
  const IWString original = _mol.unique_smiles();
  ASSERT_EQ(_mol.chiral_centres(), 1);
  Molecule_to_Match target(&_mol);
  Substructure_Query oxygen_query;
  ASSERT_TRUE(oxygen_query.create_from_smarts("O"));
  for (int i = 0; i < 10; ++i) {
    Substructure_Results results;
    ASSERT_EQ(_query.substructure_search(target, results), 1u);
    // Retain the carbon linker, rejecting the separate oxygen linker.
    EXPECT_THAT(*results.embedding(0), UnorderedElementsAre(3));
    Substructure_Results oxygen_results;
    EXPECT_EQ(oxygen_query.substructure_search(target, oxygen_results), 2u);
    EXPECT_EQ(_mol.unique_smiles(), original);
    EXPECT_EQ(_mol.chiral_centres(), 1);
    EXPECT_EQ(_mol.isotope(3), 13u);
  }
}

TEST_F(TestSubstituent, InterRingPositiveElementsRestrictNumericCountsToRegion) {
  // The whole molecule has nine ring carbons and two carbon-13 atoms.
  // The accepted region has two ring attachment carbons and one carbon-13.
  ASSERT_TRUE(_mol.build_from_smiles("C1CC1[13CH2]C2CC2OC3CC3[13C@H](O)F"));
  const IWString original = _mol.unique_smiles();
  Molecule_to_Match target(&_mol);
  for (const char* required : {"2[#6R]", "3[#6v4]", "1[13C]", "2[CD>2]", "3[C,N]"}) {
    SCOPED_TRACE(required);
    std::string text = R"pb(query {
      compress_embeddings: true
      inter_ring_region {
        disqualifying_smarts: "O"
        set_global_id: 1
        required_smarts: ")pb";
    text += required;
    text += R"pb(" } smarts: "[/IWgid1]" })pb";
    SubstructureSearch::SubstructureQuery proto;
    ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(text, &proto));
    Substructure_Query query;
    ASSERT_TRUE(query.ConstructFromProto(proto));
    Substructure_Results results;
    ASSERT_EQ(query.substructure_search(target, results), 1u);
    EXPECT_THAT(*results.embedding(0), UnorderedElementsAre(3));
    EXPECT_EQ(_mol.unique_smiles(), original);
  }
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
