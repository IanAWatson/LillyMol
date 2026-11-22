#include "tsubstructure_fp.h"

#include <sstream>
#include <string>
#include <vector>

#include "Foundational/iwbits/dy_fingerprint.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "gtest/gtest.h"

namespace {
class TsubFingerprintTest : public testing::Test {
 protected:
  Molecule _m;
  TSubstructure_FP _generator;

  void
  SetUp() override {
    ASSERT_TRUE(_m.build_from_smiles("C example"));
    _generator.set_default_fingerprint_nbits(8);
  }

  std::string
  Generate(const char* tag, const std::vector<int>& hits) {
    IWString tag_with_delimiter(tag);
    tag_with_delimiter << '<';
    _generator.set_fingerprint_tag(tag_with_delimiter);
    std::ostringstream output;
    EXPECT_TRUE(_generator.do_fingerprint_output(_m, hits.size(), hits.data(),
                                                 static_cast<std::ostream&>(output)));
    const std::string text = output.str();
    const size_t first = text.find(std::string(tag) + '<');
    EXPECT_NE(first, std::string::npos);
    if (first == std::string::npos) {
      return {};
    }
    return text.substr(first, text.find('\n', first) - first + 1);
  }

  IWString
  Payload(const std::string& line) {
    const size_t open = line.find('<');
    const size_t close = line.find('>');
    if (open == std::string::npos || close == std::string::npos || close <= open) {
      ADD_FAILURE() << line;
      return IWString();
    }
    return IWString(line.data() + open + 1, close - open - 1);
  }

  std::vector<unsigned char>
  Counts(const std::string& line) {
    IWString encoded = Payload(line);
    encoded.truncate_at_first(';');
    std::vector<unsigned char> bytes(encoded.length(), 0);
    unsigned int nbytes = 0;
    EXPECT_TRUE(du_ascii2bin(encoded.data(), encoded.length(), bytes.data(), nbytes));
    bytes.resize(nbytes);
    return bytes;
  }
};

TEST_F(TsubFingerprintTest, BinaryExtraReplicatesAreIndependent) {
  _generator.set_extra_bit_total_hits(3);
  IWString payload = Payload(Generate("FPTEST", {2}));
  payload.truncate_at_first(';');
  IW_Bits_Base bits;
  ASSERT_TRUE(bits.construct_from_daylight_ascii_representation(payload));
  EXPECT_EQ(bits.nset(), 4);
  EXPECT_TRUE(bits.is_set(0));
  for (int bit : {8, 9, 10}) {
    EXPECT_TRUE(bits.is_set(bit));
  }
}

TEST_F(TsubFingerprintTest, BinaryNineReplicatesAndOneExtra) {
  _generator.set_bit_replicates(9);
  _generator.set_extra_bit_total_hits(1);
  IWString payload = Payload(Generate("FPTEST", {2}));
  payload.truncate_at_first(';');
  IW_Bits_Base bits;
  ASSERT_TRUE(bits.construct_from_daylight_ascii_representation(payload));
  EXPECT_EQ(bits.nbits(), 80);
  EXPECT_EQ(bits.nset(), 10);
  EXPECT_TRUE(bits.is_set(72));
}

TEST_F(TsubFingerprintTest, BinaryAggregateIndicatesAnyMatch) {
  _generator.set_each_query_writes_bit(0);
  _generator.set_bit_replicates(3);
  for (const std::vector<int>& hits : {std::vector<int>{2, 0}, {0, 0}}) {
    IWString payload = Payload(Generate("FPTEST", hits));
    payload.truncate_at_first(';');
    IW_Bits_Base bits;
    ASSERT_TRUE(bits.construct_from_daylight_ascii_representation(payload));
    EXPECT_EQ(bits.nset(), hits[0] ? 3 : 0);
  }
}

TEST_F(TsubFingerprintTest, SparseExtraFeaturesHaveDistinctNumbers) {
  _generator.set_bit_replicates(2);
  _generator.set_extra_bit_total_hits(3);
  Sparse_Fingerprint_Creator expected;
  for (unsigned int bit : {0u, 1u, 4u, 5u, 6u}) {
    expected.hit_bit(bit, 2);
  }
  std::ostringstream output;
  expected.write_fingerprint(IWString("NCTEST<"), static_cast<std::ostream&>(output));
  EXPECT_EQ(Generate("NCTEST", {2, 0}), output.str());
}

TEST_F(TsubFingerprintTest, CountedExtraFeatureUsesAByteIndex) {
  _generator.set_extra_bit_total_hits(1);
  const std::string line = Generate("FCTEST", {2});
  const auto bytes = Counts(line);
  ASSERT_EQ(bytes.size(), 8u);
  EXPECT_EQ(bytes[0], 2);
  EXPECT_EQ(bytes[1], 2);
  for (size_t i = 2; i < bytes.size(); ++i) {
    EXPECT_EQ(bytes[i], 0);
  }
  EXPECT_NE(line.find(";2;2;2;2>"), std::string::npos);
}

TEST_F(TsubFingerprintTest, CountedPaddingCrossesWordBoundary) {
  _generator.set_extra_bit_total_hits(1);
  const auto bytes = Counts(Generate("FCTEST", std::vector<int>(8, 1)));
  ASSERT_EQ(bytes.size(), 16u);
  for (int i = 0; i < 8; ++i) {
    EXPECT_EQ(bytes[i], 1);
  }
  EXPECT_EQ(bytes[8], 8);
  for (size_t i = 9; i < bytes.size(); ++i) {
    EXPECT_EQ(bytes[i], 0);
  }
}

TEST_F(TsubFingerprintTest, CountedReplicatesSaturateAt255) {
  _generator.set_bit_replicates(2);
  _generator.set_extra_bit_total_hits(3);
  const std::string line = Generate("FCTEST", {2, 300});
  const auto bytes = Counts(line);
  ASSERT_EQ(bytes.size(), 8u);
  const std::vector<unsigned char> expected{2, 2, 255, 255, 255, 255, 255, 0};
  EXPECT_EQ(bytes, expected);
  EXPECT_NE(line.find(";7;7;7;7>"), std::string::npos);
}

TEST_F(TsubFingerprintTest, CountedAggregateCountsQueries) {
  _generator.set_each_query_writes_bit(0);
  _generator.set_bit_replicates(2);
  const auto bytes = Counts(Generate("FCTEST", {2, 300, 0}));
  ASSERT_EQ(bytes.size(), 8u);
  EXPECT_EQ(bytes[0], 2);
  EXPECT_EQ(bytes[1], 2);
  for (size_t i = 2; i < bytes.size(); ++i) {
    EXPECT_EQ(bytes[i], 0);
  }
}

TEST_F(TsubFingerprintTest, CountedNoMatchesRemainZero) {
  _generator.set_extra_bit_total_hits(2);
  const std::string line = Generate("FCTEST", {0, 0});
  const auto bytes = Counts(line);
  for (unsigned char count : bytes) {
    EXPECT_EQ(count, 0);
  }
  EXPECT_NE(line.find(";4;0;4;0>"), std::string::npos);
}
}  // namespace
