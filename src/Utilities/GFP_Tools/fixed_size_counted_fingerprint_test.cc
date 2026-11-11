#include "Utilities/GFP_Tools/fixed_size_counted_fingerprint.h"

#include "Foundational/iwstring/iwstring.h"

#include <sstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

Fixed_Size_Counted_Fingerprint_uchar
Fingerprint(std::initializer_list<unsigned char> values) {
  Fixed_Size_Counted_Fingerprint_uchar result;
  EXPECT_TRUE(result.construct_from_array(values.begin(), values.size()));
  return result;
}

TEST(TestFixedSizeCountedFingerprint, CopyConstructor) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar copy(fp);

  EXPECT_EQ(copy.size(), fp.size());
  EXPECT_EQ(copy.nset(), fp.nset());
  EXPECT_FLOAT_EQ(copy.tanimoto(fp), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, CopyAssignment) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar copy = Fingerprint({1, 0});

  copy = fp;

  EXPECT_EQ(copy.size(), fp.size());
  EXPECT_EQ(copy.nset(), fp.nset());
  EXPECT_FLOAT_EQ(copy.tanimoto(fp), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, MoveConstructor) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar moved(std::move(fp));

  Fixed_Size_Counted_Fingerprint_uchar expected = Fingerprint({1, 1, 0, 1});
  EXPECT_EQ(moved.size(), expected.size());
  EXPECT_EQ(moved.nset(), expected.nset());
  EXPECT_FLOAT_EQ(moved.tanimoto(expected), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, MoveAssignment) {
  Fixed_Size_Counted_Fingerprint_uchar fp = Fingerprint({1, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar moved = Fingerprint({1, 0});

  moved = std::move(fp);

  Fixed_Size_Counted_Fingerprint_uchar expected = Fingerprint({1, 1, 0, 1});
  EXPECT_EQ(moved.size(), expected.size());
  EXPECT_EQ(moved.nset(), expected.nset());
  EXPECT_FLOAT_EQ(moved.tanimoto(expected), 1.0f);
}

TEST(TestFixedSizeCountedFingerprint, VectorReallocation) {
  std::vector<Fixed_Size_Counted_Fingerprint_uchar> fingerprints;
  fingerprints.reserve(1);

  fingerprints.push_back(Fingerprint({1, 0, 1, 0}));
  fingerprints.push_back(Fingerprint({0, 1, 0, 1}));
  fingerprints.push_back(Fingerprint({1, 0, 0, 1}));

  Fixed_Size_Counted_Fingerprint_uchar expected0 = Fingerprint({1, 0, 1, 0});
  Fixed_Size_Counted_Fingerprint_uchar expected1 = Fingerprint({0, 1, 0, 1});
  Fixed_Size_Counted_Fingerprint_uchar expected2 = Fingerprint({1, 0, 0, 1});

  EXPECT_FLOAT_EQ(fingerprints[0].tanimoto(expected0), 1.0f);
  EXPECT_FLOAT_EQ(fingerprints[1].tanimoto(expected1), 1.0f);
  EXPECT_FLOAT_EQ(fingerprints[2].tanimoto(expected2), 1.0f);
}

// Counts larger than one, and a size that is not a multiple of 8. The
// fingerprint is built as the formula fingerprint generator does, so that nset is
// the sum of the counts, which is also what reading a record computes.
TEST(TestFixedSizeCountedFingerprint, WriteAndReadPreservesCounts) {
  const std::vector<int> counts = {0, 5, 255, 37, 124, 0, 1, 2, 3, 99, 0, 0, 7, 200, 0, 1, 129, 64};
  Fixed_Size_Counted_Fingerprint_uchar fp;
  ASSERT_TRUE(fp.construct_from_array_of_ints(counts.data(), counts.size()));
  ASSERT_EQ(fp.nset(), 927);

  std::ostringstream output;
  const const_IWSubstring tag("FCT<");
  ASSERT_TRUE(fp.write_daylight_ascii_representation(output, tag, 1));

  std::string record = output.str();
  ASSERT_EQ(record.back(), '\n');
  record.pop_back();

  Fixed_Size_Counted_Fingerprint_uchar read;
  ASSERT_TRUE(read.construct_from_tdt_record(const_IWSubstring(record.c_str(), record.size())));

  EXPECT_EQ(read.size(), fp.size());
  EXPECT_EQ(read.nset(), fp.nset());
  EXPECT_EQ(read.bits_in_common(fp), fp.nset());  // equal counts in every position
  EXPECT_FLOAT_EQ(read.tanimoto(fp), 1.0f);
}

}  // namespace
