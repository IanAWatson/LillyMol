#include "Utilities/GFP_Tools/gfp_context.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Utilities/GFP_Tools/fixed_size_counted_fingerprint.h"
#include "Utilities/GFP_Tools/sparsefp.h"

namespace {

using gfp_context::GFPList;
using testing::ElementsAre;
using testing::FloatNear;

IWString
FingerprintRecord(const IWString& tag, const std::vector<int>& bits) {
  fixed_bit_vector::FixedBitVector fp(128);
  for (int bit : bits) {
    fp.set_bit(bit);
  }

  IWString result;
  result << tag << fp.DaylightAsciiRepresentationIncludingNsetInfo() << ">\n";
  return result;
}

IWString
SparseRecord() {
  Sparse_Fingerprint_Creator sfc;
  sfc.hit_bit(100, 2);
  sfc.hit_bit(100000, 3);

  IWString result;
  sfc.daylight_ascii_form_with_counts_encoded("NCA<", result);
  result << '\n';
  return result;
}

IWString
MolecularPropertiesRecord() {
  return FingerprintRecord("MPR<", {0, 3, 5, 11, 17, 31, 47, 63});
}

IWString
FixedCountedRecord() {
  const unsigned char counts[] = {1, 0, 1, 0, 1, 0, 1, 0,
                                  1, 0, 1, 0, 1, 0, 1, 0};
  Fixed_Size_Counted_Fingerprint_uchar fp;
  fp.construct_from_array(counts, 16);

  std::ostringstream output;
  fp.write_daylight_ascii_representation(output, "FCA<", 1);

  return IWString(output.str().c_str());
}

std::string
TempFileName() {
  const char* tmpdir = std::getenv("TEST_TMPDIR");
  if (tmpdir == nullptr) {
    tmpdir = "/tmp";
  }

  std::string result(tmpdir);
  result.append("/gfp_context_test.gfp");
  return result;
}

std::string
CarbonTempFileName() {
  const char* tmpdir = std::getenv("TEST_TMPDIR");
  if (tmpdir == nullptr) {
    tmpdir = "/tmp";
  }

  std::string result(tmpdir);
  result.append("/gfp_context_carbon_test.gfp");
  return result;
}

std::string
WriteCarbonGfpFile() {
  const std::string fname = CarbonTempFileName();

  std::ofstream output(fname);
  output << "$SMI<C>\n";
  output << "PCN<Methane>\n";
  output << "FCTS<.E..........2;1;1;1;1>\n";
  output << "|\n";
  output << "$SMI<CC>\n";
  output << "PCN<Ethane>\n";
  output << "FCTS<.U..........2;1;1;1;1>\n";
  output << "|\n";

  return fname;
}

std::string
WriteTestGfpFile() {
  const std::string fname = TempFileName();

  IWString contents;
  contents << "$SMI<C>\n";
  contents << "PCN<methane>\n";
  contents << FingerprintRecord("FPA<", {0, 1, 2});
  contents << FingerprintRecord("FPB<", {10, 11});
  contents << SparseRecord();
  contents << MolecularPropertiesRecord();
  contents << FixedCountedRecord();
  contents << "|\n";
  contents << "$SMI<CC>\n";
  contents << "PCN<ethane>\n";
  contents << FingerprintRecord("FPA<", {0, 1});
  contents << FingerprintRecord("FPB<", {10, 11});
  contents << SparseRecord();
  contents << MolecularPropertiesRecord();
  contents << FixedCountedRecord();
  contents << "|\n";
  contents << "$SMI<O>\n";
  contents << "PCN<water>\n";
  contents << FingerprintRecord("FPA<", {20, 21});
  contents << FingerprintRecord("FPB<", {30, 31});
  contents << SparseRecord();
  contents << MolecularPropertiesRecord();
  contents << FixedCountedRecord();
  contents << "|\n";

  std::ofstream output(fname);
  output << contents;

  return fname;
}

TEST(TestGFPContext, ReadFileDiscoversContextAndComputesDistances) {
  const std::string fname = WriteTestGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  EXPECT_EQ(gfp.size(), 3);
  EXPECT_EQ(gfp.smiles(0), "C");
  EXPECT_EQ(gfp.id(1), "ethane");
  EXPECT_THAT(gfp.context().Tags(), ElementsAre("FPA<", "FPB<", "NCA<", "MPR<", "FCA<"));

  EXPECT_THAT(gfp.Distance(0, 0), FloatNear(0.0f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 15.0f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(1, 0), FloatNear(gfp.Distance(0, 1), 1.0e-6f));

  const uint64_t hash = gfp.context().context_hash();
  ASSERT_TRUE(gfp.mutable_context().SetWeight("FPB<", 0.0f));
  EXPECT_EQ(gfp.context().context_hash(), hash);
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 12.0f, 1.0e-6f));

  std::vector<IWString> only_fpa;
  only_fpa.emplace_back("FPA<");
  ASSERT_TRUE(gfp.mutable_context().UseOnly(only_fpa));
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 3.0f, 1.0e-6f));

  gfp.mutable_context().UseAll();
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(1.0f / 12.0f, 1.0e-6f));
}

TEST(TestGFPContext, NearestNeighbours) {
  const std::string fname = WriteTestGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  const auto nbrs = gfp.NearestNeighbours(0, 2);
  ASSERT_EQ(nbrs.size(), 2u);
  EXPECT_EQ(nbrs[0].index, 1);
  EXPECT_THAT(nbrs[0].distance, FloatNear(1.0f / 15.0f, 1.0e-6f));
  EXPECT_EQ(nbrs[1].index, 2);
  EXPECT_THAT(nbrs[1].distance, FloatNear(0.4f, 1.0e-6f));
}

TEST(TestGFPContext, FixedCountedCarbonDistance) {
  const std::string fname = WriteCarbonGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  EXPECT_EQ(gfp.size(), 2);
  EXPECT_THAT(gfp.context().Tags(), ElementsAre("FCTS<"));
  EXPECT_THAT(gfp.Distance(0, 0), FloatNear(0.0f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(0, 1), FloatNear(0.5f, 1.0e-6f));
  EXPECT_THAT(gfp.Distance(1, 0), FloatNear(0.5f, 1.0e-6f));
}

TEST(TestGFPContext, NearestNeighboursWithinDistance) {
  const std::string fname = WriteTestGfpFile();

  GFPList gfp;
  ASSERT_TRUE(gfp.ReadFile(fname.c_str()));

  const auto close = gfp.NearestNeighboursWithinDistance(0, 0.1f);
  ASSERT_EQ(close.size(), 1u);
  EXPECT_EQ(close[0].index, 1);
  EXPECT_THAT(close[0].distance, FloatNear(1.0f / 15.0f, 1.0e-6f));

  const auto none = gfp.NearestNeighboursWithinDistance(0, 0.01f);
  EXPECT_TRUE(none.empty());

  const auto all = gfp.NearestNeighboursWithinDistance(0, 1.0f);
  ASSERT_EQ(all.size(), 2u);
  EXPECT_EQ(all[0].index, 1);
  EXPECT_EQ(all[1].index, 2);
}

}  // namespace
