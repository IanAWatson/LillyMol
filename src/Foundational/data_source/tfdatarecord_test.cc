#include "Foundational/data_source/tfdatarecord.h"

#include <unistd.h>

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace {
class TFDataRecordTest : public testing::Test {
 protected:
  std::string _fname;

  void
  SetUp() override {
    const char* tempdir = std::getenv("TEST_TMPDIR");
    _fname = std::string(tempdir ? tempdir : "/tmp") + "/tfdatarecordXXXXXX";
    const int fd = mkstemp(_fname.data());
    ASSERT_GE(fd, 0);
    close(fd);
  }

  void
  TearDown() override {
    unlink(_fname.c_str());
  }

  void
  WriteBytes(const std::string& bytes) {
    std::ofstream output(_fname, std::ios::binary);
    output.write(bytes.data(), bytes.size());
  }

  std::string
  EncodedRecord() {
    iw_tf_data_record::TFDataWriter writer;
    EXPECT_TRUE(writer.Open(_fname.c_str()));
    EXPECT_GT(writer.Write("abcdefgh", 8), 0);
    writer.Close();
    std::ifstream input(_fname, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(input), {});
  }
};

TEST_F(TFDataRecordTest, PackedRecordsIncludingEmpty) {
  iw_tf_data_record::TFDataWriter writer;
  ASSERT_TRUE(writer.Open(_fname.c_str()));
  ASSERT_GT(writer.Write("", 0), 0);
  ASSERT_GT(writer.Write("a", 1), 0);
  ASSERT_GT(writer.Write("unaligned", 9), 0);
  writer.Close();
  iw_tf_data_record::TFDataReader reader(_fname.c_str());
  for (const std::string expected : {"", "a", "unaligned"}) {
    auto record = reader.Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->AsString(), expected);
  }
  EXPECT_FALSE(reader.Next().has_value());
  EXPECT_TRUE(reader.good());
  EXPECT_TRUE(reader.eof());
}

TEST_F(TFDataRecordTest, EmptyFileIsCleanEof) {
  iw_tf_data_record::TFDataReader reader(_fname.c_str());
  EXPECT_FALSE(reader.Next().has_value());
  EXPECT_TRUE(reader.good());
  EXPECT_TRUE(reader.eof());
}

TEST_F(TFDataRecordTest, RoundTripAcrossBufferBoundaries) {
  std::vector<std::string> records;
  iw_tf_data_record::TFDataWriter writer;
  ASSERT_TRUE(writer.Open(_fname.c_str()));
  for (int i = 0; i < 400; ++i) {
    // Exercise split headers, unaligned CRCs, empty records and buffer growth.
    const size_t length = i == 200 ? 12000 : (i * 17) % 67;
    records.emplace_back(length, 'A' + i % 26);
    ASSERT_GT(writer.Write(records.back().data(), records.back().size()), 0);
  }
  writer.Close();
  iw_tf_data_record::TFDataReader reader(_fname.c_str());
  for (const std::string& expected : records) {
    auto record = reader.Next();
    ASSERT_TRUE(record.has_value());
    EXPECT_EQ(record->AsString(), expected);
  }
  EXPECT_EQ(reader.items_read(), 400);
  EXPECT_FALSE(reader.Next().has_value());
  EXPECT_TRUE(reader.good());
  EXPECT_TRUE(reader.eof());
}

TEST_F(TFDataRecordTest, EveryTruncatedPrefixFails) {
  const std::string complete = EncodedRecord();
  for (size_t length = 1; length < complete.size(); ++length) {
    SCOPED_TRACE(length);
    WriteBytes(complete.substr(0, length));
    iw_tf_data_record::TFDataReader reader(_fname.c_str());
    EXPECT_FALSE(reader.Next().has_value());
    EXPECT_FALSE(reader.good());
  }
}

TEST_F(TFDataRecordTest, InvalidChecksumsFail) {
  const std::string complete = EncodedRecord();
  for (size_t offset : {size_t{8}, complete.size() - 1}) {
    std::string corrupt = complete;
    corrupt[offset] ^= 1;
    WriteBytes(corrupt);
    iw_tf_data_record::TFDataReader reader(_fname.c_str());
    EXPECT_FALSE(reader.Next().has_value());
    EXPECT_FALSE(reader.good());
  }
}
}  // namespace
