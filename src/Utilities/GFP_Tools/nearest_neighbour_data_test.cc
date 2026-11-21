#include "nearest_neighbour_data.h"

#include <unistd.h>

#include <cstdio>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "Foundational/data_source/tfdatarecord.h"

#include "gtest/gtest.h"

namespace nearest_neighbour_data {
namespace {
class NearestNeighbourDataTest : public testing::Test {
 protected:
  std::string _fname =
      testing::TempDir() + "nearest_neighbour_data." + std::to_string(getpid());
  NearestNeighbourData _data;

  void
  TearDown() override {
    std::remove(_fname.c_str());
  }

  nnbr::NearNeighbours
  Record(const std::string& id, const std::vector<float>& distances) {
    nnbr::NearNeighbours result;
    result.set_name(id);
    result.set_smiles("CC");
    for (size_t i = 0; i < distances.size(); ++i) {
      auto* nbr = result.add_nbr();
      nbr->set_id("nbr" + std::to_string(i));
      nbr->set_dist(distances[i]);
      nbr->set_smi("CCC");
    }
    return result;
  }

  void
  Write(const std::vector<nnbr::NearNeighbours>& protos) {
    iw_tf_data_record::TFDataWriter writer;
    ASSERT_TRUE(writer.Open(_fname.c_str()));
    for (const auto& proto : protos) {
      ASSERT_TRUE(writer.WriteSerializedProto(proto));
    }
    ASSERT_TRUE(writer.Close());
  }
  nnbr::NearNeighboursIndices
  IndexedRecord(const std::string& id,
                const std::vector<std::pair<uint32_t, float>>& neighbours) {
    nnbr::NearNeighboursIndices result;
    result.set_name(id);
    result.set_smiles("CC");
    for (const auto& [index, distance] : neighbours) {
      auto* nbr = result.add_nbr();
      nbr->set_id(index);
      nbr->set_dist(distance);
    }
    return result;
  }

  void
  WriteIndexed(const std::vector<nnbr::NearNeighboursIndices>& protos) {
    iw_tf_data_record::TFDataWriter writer;
    ASSERT_TRUE(writer.Open(_fname.c_str()));
    for (const auto& proto : protos) {
      ASSERT_TRUE(writer.WriteSerializedProto(proto));
    }
    ASSERT_TRUE(writer.Close());
  }

};

TEST_F(NearestNeighbourDataTest, RetainsConfiguredCountAndLooksUpByName) {
  auto proto = Record("a", {0.1f, 0.2f, 0.4f});
  proto.set_cluster(2);
  Write({proto, Record("b", {})});
  ASSERT_TRUE(_data.SetNeighboursToRetain(2));
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  EXPECT_EQ(_data.size(), 2u);
  const auto* result = _data.Lookup("a");
  ASSERT_NE(result, nullptr);
  EXPECT_EQ(result->nbr_size(), 2);
  EXPECT_EQ(result->nbr(1).id(), "nbr1");
  EXPECT_EQ(result->nbr(1).smi(), "CCC");
  EXPECT_EQ(result->cluster(), 2u);
  ASSERT_TRUE(_data.NearestDistance("a"));
  EXPECT_FLOAT_EQ(*_data.NearestDistance("a"), 0.1f);
  ASSERT_NE(_data.Lookup("b"), nullptr);
  EXPECT_FALSE(_data.NearestDistance("b"));
  EXPECT_EQ(_data.Lookup("missing"), nullptr);
  EXPECT_FALSE(_data.NearestDistance("missing"));
  EXPECT_FALSE(_data.SetNeighboursToRetain(3));
  EXPECT_FALSE(_data.SetNeighboursToRetain(2));
  EXPECT_EQ(_data.neighbours_to_retain(), 2);
}

TEST_F(NearestNeighbourDataTest, SkipsSelfButPreservesNonSelfZeroDistance) {
  auto proto = Record("a", {0.0f, 0.0f, 0.2f});
  proto.mutable_nbr(0)->set_id("a");
  Write({proto});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  ASSERT_TRUE(_data.NearestDistance("a"));
  EXPECT_FLOAT_EQ(*_data.NearestDistance("a"), 0.0f);
  EXPECT_EQ(_data.Lookup("a")->nbr(0).id(), "nbr1");
  EXPECT_EQ(_data.Lookup("a")->nbr_size(), 1);
}

TEST_F(NearestNeighbourDataTest, SelfOnlyBecomesEmpty) {
  auto proto = Record("a", {0.0f});
  proto.mutable_nbr(0)->set_id("a");
  Write({proto});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  ASSERT_NE(_data.Lookup("a"), nullptr);
  EXPECT_FALSE(_data.NearestDistance("a"));
}

TEST_F(NearestNeighbourDataTest, CompactCopyDoesNotRetainDiscardedStorage) {
  auto proto = Record("a", std::vector<float>(1000, 0.1f));
  Write({proto});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  ASSERT_NE(_data.Lookup("a"), nullptr);
  EXPECT_LT(_data.Lookup("a")->SpaceUsedLong(), proto.SpaceUsedLong() / 10);
}

TEST_F(NearestNeighbourDataTest, InvalidRetentionCountsAndEmptyFile) {
  EXPECT_FALSE(_data.SetNeighboursToRetain(0));
  EXPECT_FALSE(_data.SetNeighboursToRetain(-1));
  Write({});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  EXPECT_EQ(_data.size(), 0u);
  EXPECT_FALSE(_data.SetNeighboursToRetain(1));
}

TEST_F(NearestNeighbourDataTest, FailedReadPreservesExistingData) {
  Write({Record("a", {0.1f})});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  Write({Record("b", {0.2f}), Record("b", {0.3f})});
  EXPECT_FALSE(_data.Read(_fname.c_str()));
  EXPECT_EQ(_data.size(), 1u);
  EXPECT_NE(_data.Lookup("a"), nullptr);
  EXPECT_EQ(_data.Lookup("b"), nullptr);
  Write({Record("b", {0.2f})});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  EXPECT_EQ(_data.Lookup("a"), nullptr);
  EXPECT_NE(_data.Lookup("b"), nullptr);
}

TEST_F(NearestNeighbourDataTest, RejectsMissingIdsAndInvalidDistances) {
  Write({Record("", {0.1f})});
  EXPECT_FALSE(_data.Read(_fname.c_str()));
  for (float distance : {-0.1f, std::numeric_limits<float>::infinity(),
                         std::numeric_limits<float>::quiet_NaN()}) {
    Write({Record("a", {distance})});
    EXPECT_FALSE(_data.Read(_fname.c_str()));
    EXPECT_EQ(_data.size(), 0u);
  }
  EXPECT_TRUE(_data.SetNeighboursToRetain(2));
}

TEST_F(NearestNeighbourDataTest, RejectsMalformedProtoAndMissingFile) {
  {
    iw_tf_data_record::TFDataWriter writer;
    ASSERT_TRUE(writer.Open(_fname.c_str()));
    ASSERT_TRUE(writer.Write("\xff", 1));
    ASSERT_TRUE(writer.Close());
  }
  EXPECT_FALSE(_data.Read(_fname.c_str()));
  EXPECT_FALSE(_data.Read((_fname + ".missing").c_str()));
}
TEST_F(NearestNeighbourDataTest, IndexedResolvesForwardAndBackwardReferences) {
  auto c = IndexedRecord("c", {});
  c.set_smiles("CCC");
  WriteIndexed({IndexedRecord("a", {{0, 0.0f}, {2, 0.1f}, {1, 0.2f}, {3, 0.4f}}),
                IndexedRecord("b", {{0, 0.0f}}), c,
                IndexedRecord("d", {{3, 0.0f}})});
  ASSERT_TRUE(_data.SetNeighboursToRetain(2));
  ASSERT_TRUE(_data.ReadIndexed(_fname.c_str()));
  ASSERT_EQ(_data.size(), 4u);
  const auto* a = _data.Lookup("a");
  ASSERT_NE(a, nullptr);
  ASSERT_EQ(a->nbr_size(), 2);
  EXPECT_EQ(a->nbr(0).id(), "c");
  EXPECT_EQ(a->nbr(0).smi(), "CCC");
  EXPECT_EQ(a->nbr(1).id(), "b");
  EXPECT_FLOAT_EQ(*_data.NearestDistance("a"), 0.1f);
  ASSERT_TRUE(_data.NearestDistance("b"));
  EXPECT_FLOAT_EQ(*_data.NearestDistance("b"), 0.0f);
  EXPECT_EQ(_data.Lookup("b")->nbr(0).id(), "a");
  EXPECT_FALSE(_data.NearestDistance("c"));
  ASSERT_NE(_data.Lookup("d"), nullptr);
  EXPECT_FALSE(_data.NearestDistance("d"));
  EXPECT_EQ(_data.Lookup("missing"), nullptr);
  EXPECT_FALSE(_data.SetNeighboursToRetain(3));
}

TEST_F(NearestNeighbourDataTest, IndexedFailurePreservesPreviouslyLoadedNamedData) {
  Write({Record("original", {0.2f})});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  WriteIndexed({IndexedRecord("a", {{1, 0.1f}})});
  EXPECT_FALSE(_data.ReadIndexed(_fname.c_str()));
  EXPECT_EQ(_data.size(), 1u);
  EXPECT_EQ(_data.Lookup("a"), nullptr);
  ASSERT_TRUE(_data.NearestDistance("original"));
  EXPECT_FLOAT_EQ(*_data.NearestDistance("original"), 0.2f);
  WriteIndexed({IndexedRecord("a", {})});
  ASSERT_TRUE(_data.ReadIndexed(_fname.c_str()));
  EXPECT_EQ(_data.Lookup("original"), nullptr);
  EXPECT_NE(_data.Lookup("a"), nullptr);
  Write({Record("named", {0.3f})});
  ASSERT_TRUE(_data.Read(_fname.c_str()));
  EXPECT_NE(_data.Lookup("named"), nullptr);
}

TEST_F(NearestNeighbourDataTest, IndexedRejectsInvalidIdentifiersAndDistances) {
  WriteIndexed({IndexedRecord("", {})});
  EXPECT_FALSE(_data.ReadIndexed(_fname.c_str()));
  WriteIndexed({IndexedRecord("a", {}), IndexedRecord("a", {})});
  EXPECT_FALSE(_data.ReadIndexed(_fname.c_str()));
  for (float distance : {-0.1f, std::numeric_limits<float>::infinity(),
                         std::numeric_limits<float>::quiet_NaN()}) {
    WriteIndexed({IndexedRecord("a", {{1, distance}}), IndexedRecord("b", {})});
    EXPECT_FALSE(_data.ReadIndexed(_fname.c_str()));
    EXPECT_EQ(_data.size(), 0u);
  }
  EXPECT_TRUE(_data.SetNeighboursToRetain(2));
}

TEST_F(NearestNeighbourDataTest, IndexedCompactCopyDiscardsUnretainedNeighbours) {
  auto a = IndexedRecord("a", std::vector<std::pair<uint32_t, float>>(1000, {1, 0.1f}));
  WriteIndexed({a, IndexedRecord("b", {})});
  ASSERT_TRUE(_data.ReadIndexed(_fname.c_str()));
  ASSERT_NE(_data.Lookup("a"), nullptr);
  EXPECT_EQ(_data.Lookup("a")->nbr_size(), 1);
  EXPECT_LT(_data.Lookup("a")->SpaceUsedLong(), a.SpaceUsedLong() / 10);
}

TEST_F(NearestNeighbourDataTest, IndexedEmptyFileLocksRetentionCount) {
  WriteIndexed({});
  ASSERT_TRUE(_data.ReadIndexed(_fname.c_str()));
  EXPECT_EQ(_data.size(), 0u);
  EXPECT_FALSE(_data.SetNeighboursToRetain(2));
}

TEST_F(NearestNeighbourDataTest, IndexedRejectsMalformedProtoAndMissingFile) {
  {
    iw_tf_data_record::TFDataWriter writer;
    ASSERT_TRUE(writer.Open(_fname.c_str()));
    ASSERT_TRUE(writer.Write("\xff", 1));
    ASSERT_TRUE(writer.Close());
  }
  EXPECT_FALSE(_data.ReadIndexed(_fname.c_str()));
  EXPECT_FALSE(_data.ReadIndexed((_fname + ".missing").c_str()));
}

}  // namespace
}  // namespace nearest_neighbour_data
