#include <iostream>
#include <string>

#include "aromatic.h"
#include "substructure.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using std::cerr;
using std::endl;

using testing::UnorderedElementsAre;
using testing::Contains;

class TestSubstructure : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    std::string _string_proto;

    IWString _smiles;

    SubstructureSearch::SubstructureQuery _proto;

    Substructure_Query _query;

    Substructure_Results _sresults;

    Molecule _m;

    static constexpr int _ntest = 10;

  protected:
    void _WriteQuery(const char * fname) {
      IWString tmp(fname);
      _query.write_msi(tmp);
    }

    const Set_of_Atoms _FirstAtomEachEmbedding()
    {
      Set_of_Atoms to_be_returned;
      for (const auto * e : _sresults.embeddings()) {
        to_be_returned.add(e->item(0));
      }

      return to_be_returned;
    }

    bool _DoPerumationsTests(const int expected);
};

void
TestSubstructure::SetUp()
{
  set_global_aromaticity_type(Daylight);
}

bool TestSubstructure::_DoPerumationsTests(const int expected)
{
  const IWString initial_smiles = _m.smiles();

  for (int i = 0; i < _ntest; ++i) {
    const IWString smiles = _m.random_smiles();
    if (! _m.build_from_smiles(smiles)) {
      cerr << "_DoPerumationsTests:invalid smiles '" << smiles << "'\n";
      return false;
    }

    const int hits = _query.substructure_search(_m, _sresults);

    if (expected == hits)
      continue;

    cerr << "_DoPerumationsTests:permutation test failure, expected " << expected << " got " << hits << " hits\n";
    cerr << "Initial smiles " << initial_smiles << endl;
    cerr << "Random  smiles " << smiles << endl;

    return false;
  }

  return true;
}

TEST_F(TestSubstructure, TestOneEmbeddingPerStartAtom)
{
  _string_proto = R"(
    query {
      one_embedding_per_start_atom: true
      smarts: "C(F)(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F.FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestUniqueEmbeddingsOnly)
{
  _string_proto = R"(
    query {
      unique_embeddings_only: true
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMaxMatchesToFine)
{
  _string_proto = R"(
    query {
      max_matches_to_find: 3
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3u);

  EXPECT_TRUE(_DoPerumationsTests(3));
}

TEST_F(TestSubstructure, TestSubtractFromRc)
{
  _string_proto = R"(
    query {
      subtract_from_rc: 2
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 22u);

  EXPECT_TRUE(_DoPerumationsTests(22));
}

TEST_F(TestSubstructure, TestHitsNeeded)
{
  _string_proto = R"(
    query {
      min_hits_needed: 15
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 24u);

  EXPECT_TRUE(_DoPerumationsTests(24));
}

TEST_F(TestSubstructure, TestNoSymmetricHits)
{
  _string_proto = R"(
    query {
      perceive_symmetric_equivalents: false
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestEmbeddingsDoNotOverlap1)
{
  _string_proto = R"(
    query {
      embeddings_do_not_overlap: true
      smarts: "FC(F)(F)F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();


  _smiles = "FC(F)(F)F";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

// Hmmm, was going to use benzene as a target molecule, with the query cc, but then
// I realised this is undefined. There might be either 2 or 3 possible embeddings that
// do not violate the overlapping embeddings constraint. Doing it optimally would be
// hard, so I just do not run that test.
TEST_F(TestSubstructure, TestEmbeddingsDoNotOverlapAmbiguous)
{
  _string_proto = R"(
    query {
      embeddings_do_not_overlap: true
      smarts: "cc"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  // currently the software returns 3, but there is no guarantee of that.
//EXPECT_EQ(_query.substructure_search(_m, _sresults), 3);
}

TEST_F(TestSubstructure, TestEmbeddingsDoNotOverlapEasy)
{
  _string_proto = R"(
    query {
      embeddings_do_not_overlap: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();


  _smiles = "CC(=O)NC(=O)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAllHitsInSameFragment)
{
  _string_proto = R"(
    query {
      all_hits_in_same_fragment: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();


  _smiles = "CC(=O)NC(=O)C.CC(=O)NC(=O)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestOnlyMatchLargestFragment)
{
  _string_proto = R"(
    query {
      only_match_largest_fragment: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(=O)NC(=O)C.CC(=O)NC(=O)CC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  Set_of_Atoms matched;
  for (const auto* e : _sresults.embeddings()) {
    matched += *e;
  }

  EXPECT_THAT(matched, UnorderedElementsAre(8, 9, 10, 10, 11, 12));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestOnlyMatchLargestFragmentNoMatchesInFirst)
{
  _string_proto = R"(
    query {
      only_match_largest_fragment: true
      smarts: "O=C-N"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(=O)NC(=O)CC.CCCCCCCCCCCCCCCCCCCCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  _DoPerumationsTests(0);
}

TEST_F(TestSubstructure, TestImplicitRingCondition1)
{
  _string_proto = R"(
    query {
      implicit_ring_condition: 1
      smarts: "C-C-C-C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 8u);

  EXPECT_TRUE(_DoPerumationsTests(8));
}

TEST_F(TestSubstructure, TestImplicitRingCondition0)
{
  _string_proto = R"(
    query {
      implicit_ring_condition: 0
      smarts: "C-C-C-C"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "C1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestDistanceBetweenHitsMax)
{
  _string_proto = R"(
    query {
      max_distance_between_hits: 3
      smarts: "F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1c(F)cc(F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(0, 3));

//EXPECT_TRUE(_DoPerumationsTests(2)); // permutation tests not done because of ordering effects
}

// Here again results are ambiguous and not tested.
// The algorithm works by removing hits that are too close
// to a previous match. but since the order in which matches
// are discovered is undefined, there is no real way of 
// telling how the winnowing will work.
TEST_F(TestSubstructure, TestDistanceBetweenHitsMin)
{
  _string_proto = R"(
    query {
      min_distance_between_hits: 4
      smarts: "F"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "Fc1c(F)cc(F)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, Contains(6));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

// This might be flaky. I think this should work with distance_between_hits:5
// but it did not. Should be investigated.
TEST_F(TestSubstructure, TestDistanceBetweenHitsMatchedAtom)
{
  _string_proto = R"(
    query {
      distance_between_hits_ncheck: 2
      max_distance_between_hits: 5
      smarts: "OC(=O)c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)cc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(0, 7));

//EXPECT_TRUE(_DoPerumationsTests(2));  permutation tests not done
}

TEST_F(TestSubstructure, TestAttachedHeteroatomCount)
{
  _string_proto = R"(
    query {
      attached_heteroatom_count: 0
      smarts: "OC(=O)c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  const Set_of_Atoms matches = _FirstAtomEachEmbedding();

  EXPECT_THAT(matches, UnorderedElementsAre(0, 7));

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestRingAtomsMatched)
{
  _string_proto = R"(
    query {
      ring_atoms_matched: 1
      smarts: "OC(=O)c"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 3u);

  EXPECT_TRUE(_DoPerumationsTests(3));
}

TEST_F(TestSubstructure, TestHeteroatomsMatched)
{
  _string_proto = R"(
    query {
      heteroatoms_matched: 2
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4u);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 3
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 9
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4u);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestNatomsNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        natoms: 9
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNatomsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        natoms: 19
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4u);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestNringsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 2
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNringsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 1
      }
      smarts: "OC(=*)[#6]"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "OC(=O)c1c(C(=O)O)cc(C(=O)O)nc1.CC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 4u);

  EXPECT_TRUE(_DoPerumationsTests(4));
}

TEST_F(TestSubstructure, TestNcon)
{
  _string_proto = R"(
    query {
      ncon: 2
      smarts: "FC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "FCC.FC(C)C.FC(C)(C)C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  const Set_of_Atoms matched = _FirstAtomEachEmbedding();

  EXPECT_EQ(matched[0], 3);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNcon2)
{
  _string_proto = R"(
    query {
      ncon: 2
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC.CNC.CNCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_THAT(*_sresults.embedding(0), UnorderedElementsAre(6, 7));

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestHeteroatomsInMoleculeMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        heteroatoms_in_molecule: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestNringsNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNrings0OK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 0
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCN";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestNrings1)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        nrings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestFusedRingsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        fused_rings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestFusedRingsZeroOk)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        fused_rings: 0
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestFusedRings2)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        fused_rings: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CNC1CC2C1CC2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}

TEST_F(TestSubstructure, TestStronglyFusedRingsNo)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        min_strongly_fused_rings: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1CCC12C3C4C1C5C2C3C45";  // Cubane has no rings that share more than 1 bond.

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestStronglyFusedRingsOK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        min_strongly_fused_rings: 2
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "N1CCC1C1C3CC2CC(CC1C2)C3";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 2u);

  EXPECT_TRUE(_DoPerumationsTests(2));
}


TEST_F(TestSubstructure, TestIsolatedRingsNotIsolated)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        isolated_rings: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "c12cccc1c(CN)ccn2";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingsNoRing)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        isolated_rings: 1
      }
      smarts: "NC"
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingsOK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        isolated_rings: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestIsolatedRingObjectsNoRings)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        ring_systems: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingObjectsNotMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        ring_systems: 2
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCc1cccc2ccncc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestIsolatedRingObjectsMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        ring_systems: 1
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCc1cccc2ccncc12";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAromaticAtomsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        aromatic_atoms: 0
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAromaticAtoms0NoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        aromatic_atoms: 0
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestAromaticAtomsOK)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        aromatic_atoms: 6
      }
      smarts: "CN";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCc1ccccc1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestHeteroatomsMatches)
{
  _string_proto = R"(
    query {
      heteroatoms: [9, 17]
      attached_heteroatom_count: 1
      heteroatoms_matched: 1
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestHeteroatomsNoMatches)
{
  _string_proto = R"(
    query {
      heteroatoms: [9, 17]
      attached_heteroatom_count: 1
      heteroatoms_matched: 1
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementHitsNeededMatches)
{
  _string_proto = R"(
    query {
      element_hits_needed {
        atomic_number: 7
        hits_needed: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestElementHitsNeededNoMatches)
{
  _string_proto = R"(
    query {
      element_hits_needed {
        atomic_number: 7
        hits_needed: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementsNeededNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 7
          hits_needed: 1
        }
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementsNeededNumberNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 7
          hits_needed: 2
        }
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestElementsNeededNumberMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        elements_needed {
          atomic_number: 7
          hits_needed: 1
        }
        elements_needed {
          atomic_number: 6
          hits_needed: 3
        }
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNonAromaticRingsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        non_aromatic_rings: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNonAromaticRingsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        non_aromatic_rings: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNonAromaticRingsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        non_aromatic_rings: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNumberIsotopicAtomsNone)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_isotopic_atoms: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNumberIsotopicAtomsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_isotopic_atoms: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1[3CH2]C1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNumberIsotopicAtomsMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_isotopic_atoms: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1[3CH2]C1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNumberFragmentsNoMatch)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_fragments: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestNumberFragmentsMatch1)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_fragments: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNumberFragmentsMatch2)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        number_fragments: 2
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCC1CC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAtomsInSpinachAllSpinach)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        atoms_in_spinach: 7
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCC.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestAtomsInSpinachMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        atoms_in_spinach: 3
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestInterRingAtomsNoneMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        inter_ring_atoms: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestInterRingAtomsZero)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        inter_ring_atoms: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C1CC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestInterRingAtomsOne)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        inter_ring_atoms: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CC1CC1.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsNoMatch)
{
  _string_proto = R"(
    query {
      unmatched_atoms: 1
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C.[Br-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestUnmatchedAtomsMatch)
{
  _string_proto = R"(
    query {
      unmatched_atoms: 5
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C.[Br-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargeZeroMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargeZeroNoMatches)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: 0
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargePositive)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: 1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "[NH3+]CC1CCC1C";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNetFormalChargeNegative)
{
  _string_proto = R"(
    query {
      required_molecular_properties {
        net_formal_charge: -1
      }
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1C.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedFailOne)
{
  _string_proto = R"(
    query {
      min_fraction_atoms_matched: 0.34
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedMatchOne)
{
  _string_proto = R"(
    query {
      min_fraction_atoms_matched: 0.330
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedMatchTwo)
{
  _string_proto = R"(
    query {
      min_fraction_atoms_matched: 0.660
      smarts: "N-CC.CCO";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestMinFractionAtomsMatchedMaxViolated)
{
  _string_proto = R"(
    query {
      max_fraction_atoms_matched: 0.330
      smarts: "N-CC";
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestRespectInitialAlllNumbering)
{
  _string_proto = R"(
    query {
      respect_initial_atom_numbering: true
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 6
        }
        single_bond: 5
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCC1CCC1CO.[Cl-]";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  const Set_of_Atoms * e = _sresults.embedding(0);

  EXPECT_THAT(*e, testing::ElementsAre(-1, -1, -1, 1, -1, 0));

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNoMatchedAtomsBetweenMatches)
{
  _string_proto = R"(
    query {
      unique_embeddings_only: true
      no_matched_atoms_between {
        a1: 1
        a2: 4
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 8
        }
        double_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 6
        }
        single_bond: 3
      }
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 8
        }
        double_bond: 4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "COC(=O)CCCCC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestNoMatchedAtomsBetweenNoMatch)
{
  _string_proto = R"(
    query {
      unique_embeddings_only: true
      no_matched_atoms_between {
        a1: 1
        a2: 4
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 8
        }
        double_bond: 1
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 8
        }
      }
      query_atom {
        id: 4
        atom_properties {
          atomic_number: 6
        }
        single_bond: 3
      }
      query_atom {
        id: 5
        atom_properties {
          atomic_number: 8
        }
        double_bond: 4
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "CC(=O)OCCCCC(=O)O";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, TestLinkAtomsMatches)
{
  _string_proto = R"(
    query {
      link_atoms {
        a1: 1
        a2: 2
        distance: 3
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCCF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);

  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, TestLinkAtomsNoMatches)
{
  _string_proto = R"(
    query {
      link_atoms {
        a1: 1
        a2: 2
        min_distance: 4
      }
      query_atom {
        id: 0
        atom_properties {
          atomic_number: 7
        }
      }
      query_atom {
        id: 1
        atom_properties {
          atomic_number: 6
        }
        single_bond: 0
      }
      query_atom {
        id: 2
        atom_properties {
          atomic_number: 6
        }
      }
      query_atom {
        id: 3
        atom_properties {
          atomic_number: 9
        }
        single_bond: 2
      }
    }
  )";

  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));

  EXPECT_TRUE(_query.ConstructFromProto(_proto)) << "Cannot parse proto " << _proto.ShortDebugString();

  _smiles = "NCCCCCF";

  ASSERT_TRUE(_m.build_from_smiles(_smiles));

  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  EXPECT_TRUE(_DoPerumationsTests(0));
}

TEST_F(TestSubstructure, SmilesSingleComponent) {
  ASSERT_TRUE(_query.create_from_smiles("C"));

  ASSERT_TRUE(_m.build_from_smiles("CC"));
  EXPECT_GT(_query.substructure_search(_m, _sresults), 0u);
}

// Each '.' separated component of the smiles is a separate root atom, and
// all of them must match. Previously only the first component was used, so a
// molecule with no oxygen matched.
TEST_F(TestSubstructure, SmilesMultipleComponentsAllMustMatch) {
  ASSERT_TRUE(_query.create_from_smiles("C.O"));

  ASSERT_TRUE(_m.build_from_smiles("CC"));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 0u);

  ASSERT_TRUE(_m.build_from_smiles("CCO"));
  EXPECT_GT(_query.substructure_search(_m, _sresults), 0u);
}

TEST_F(TestSubstructure, NeutralMatchedAtomsRuntimeConstraint) {
  struct Example {
    const char* smarts;
    const char* smiles;
    int unrestricted;
    int neutral;
  };
  for (const Example& example : {
           Example{"[#7]", "N", 1, 1},
           Example{"[#7]", "[NH4+]", 1, 0},
           Example{"[#8]", "[O-]C", 1, 0},
           Example{"[#7]~[#6]~[#8]", "[NH3+]C[O-]", 1, 0},
           Example{"[#7]~[#6]~[#8]", "NCO", 1, 1},
           Example{"CC", "CC[NH3+]", 2, 2},
           Example{"CC", "CC.[Cl-]", 2, 2},
           Example{"[N+]", "[NH4+]", 1, 0},
           Example{"[N,O]", "[NH4+].[OH-]", 2, 0},
           Example{"[#7].[#8]", "[NH4+].O", 1, 0},
           Example{"[$(C[N+])]", "C[NH3+]", 1, 1},
           Example{"[#7]", "[NH4+].N", 2, 1}}) {
    SCOPED_TRACE(example.smarts);
    SCOPED_TRACE(example.smiles);
    // Test absent, explicit false and true, including count-only searches.
    for (int flag : {-1, 0, 1}) {
      for (bool save : {false, true}) {
        SubstructureSearch::SubstructureQuery proto;
        auto* component = proto.add_query();
        component->set_smarts(example.smarts);
        component->set_save_matched_atoms(save);
        if (flag >= 0) {
          component->set_require_neutral_matched_atoms(flag);
        }
        Substructure_Query query;
        ASSERT_TRUE(query.ConstructFromProto(proto));
        Molecule original;
        ASSERT_TRUE(original.build_from_smiles(example.smiles));
        const uint32_t expected = flag == 1 ? example.neutral : example.unrestricted;
        for (int i = 0; i < _ntest; ++i) {
          Molecule molecule;
          ASSERT_TRUE(molecule.build_from_smiles(original.random_smiles()));
          Substructure_Results results;
          EXPECT_EQ(query.substructure_search(molecule, results), expected);
        }
      }
    }
  }
}

TEST_F(TestSubstructure, NeutralMatchedAtomsChecksBeforeMatchLimit) {
  _string_proto = R"(query {
    smarts: "[#7]"
    require_neutral_matched_atoms: true
    max_matches_to_find: 1
  })";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_m.build_from_smiles("[NH4+].N"));
  ASSERT_EQ(_query.substructure_search(_m, _sresults), 1u);
  ASSERT_EQ(_sresults.embeddings().size(), 1u);
  EXPECT_EQ(_m.formal_charge(_sresults.embeddings()[0]->front()), 0);
}

TEST_F(TestSubstructure, NeutralMatchedAtomsDoesNotConstrainEnvironment) {
  _string_proto = R"(query {
    smarts: "C"
    require_neutral_matched_atoms: true
    environment {
      smarts: "[N+]"
      attachment { attachment_point: 0 btype: SS_SINGLE_BOND }
    }
  })";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_m.build_from_smiles("C[NH3+]"));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);
  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, NeutralMatchedAtomsIgnoresExcludedQueryAtoms) {
  _string_proto = R"(query {
    require_neutral_matched_atoms: true
    query_atom { id: 0 atom_smarts: "[C]" }
    query_atom {
      id: 1 atom_smarts: "[N+]" single_bond: 0 include_in_embedding: false
    }
  })";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  ASSERT_TRUE(_m.build_from_smiles("C[NH3+]"));
  EXPECT_EQ(_query.substructure_search(_m, _sresults), 1u);
  EXPECT_TRUE(_DoPerumationsTests(1));
}

TEST_F(TestSubstructure, NeutralMatchedAtomsProtoRoundTripAndComponents) {
  _string_proto = R"(
    query { smarts: "[#7]" require_neutral_matched_atoms: true }
    query { smarts: "[#8]" require_neutral_matched_atoms: false }
  )";
  ASSERT_TRUE(google::protobuf::TextFormat::ParseFromString(_string_proto, &_proto));
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  auto written = _query.BuildProto();
  ASSERT_EQ(written.query_size(), 2);
  EXPECT_TRUE(written.query(0).require_neutral_matched_atoms());
  EXPECT_FALSE(written.query(1).require_neutral_matched_atoms());
  Substructure_Query restored;
  ASSERT_TRUE(restored.ConstructFromProto(written));
  ASSERT_TRUE(_m.build_from_smiles("[NH4+]"));
  EXPECT_EQ(restored.substructure_search(_m, _sresults), 0u);
  ASSERT_TRUE(_m.build_from_smiles("[OH-]"));
  EXPECT_EQ(restored.substructure_search(_m, _sresults), 1u);
}

TEST_F(TestSubstructure, ProtoRoundTripPreservesSignedGlobalIdMatch) {
  auto* atom = _proto.add_query()->add_query_atom();
  atom->set_id(0);
  atom->set_atom_smarts("[C]");
  atom->set_global_id_match(7);
  ASSERT_TRUE(_query.ConstructFromProto(_proto));
  const auto written = _query.BuildProto();
  ASSERT_EQ(written.query_size(), 1);
  ASSERT_EQ(written.query(0).query_atom_size(), 1);
  EXPECT_EQ(written.query(0).query_atom(0).global_id_match(), 7);
}

}  // namespace
