#include <filesystem>
#include <fstream>
#include <string>

#include "gtest/gtest.h"

#define NUMERIC_DATA_FILE_IMPLEMENTATION
#define ACTIVITY_DATA_IMPLEMENATION_H
#include "activity_data_from_file.h"
#include "iw_tabular_data.h"
#include "numeric_data_from_file.h"

namespace {
class TabularReaders : public testing::Test {
 protected:
  std::string _path;

  void
  SetUp() override {
    _path = (std::filesystem::path(testing::TempDir()) / "iwmisc_reader.txt").string();
  }

  void
  TearDown() override {
    std::filesystem::remove(_path);
    std::filesystem::remove(_path + ".out");
  }

  void
  Write(const char* contents) {
    std::ofstream output(_path);
    output << contents;
  }

  const_IWSubstring
  Filename() const {
    return const_IWSubstring(_path.c_str());
  }
};

TEST_F(TabularReaders, NumericHeaderAndDoublePrecision) {
  Write("# comment\n\nID Value\na 1.23456789012345\n");
  Numeric_Data_From_File<double> data;
  ASSERT_EQ(data.read_data(Filename()), 1);
  EXPECT_EQ(data.header(), "Value");
  double value;
  EXPECT_FALSE(data.get_data(const_IWSubstring("ID"), value));
  ASSERT_TRUE(data.get_data(const_IWSubstring("a"), value));
  EXPECT_DOUBLE_EQ(value, 1.23456789012345);
}

TEST_F(TabularReaders, NumericRejectsInvalidRecordsAfterHeader) {
  Write("ID Value\na invalid\n");
  Numeric_Data_From_File<float> data;
  EXPECT_FALSE(data.read_data(Filename()));
  EXPECT_TRUE(data.empty());
}

TEST_F(TabularReaders, IntegerReaderDoesNotRoundThroughFloat) {
  Write("a 16777217\n");
  Numeric_Data_From_File<int> data;
  ASSERT_TRUE(data.read_data(Filename()));
  int value;
  ASSERT_TRUE(data.get_data(const_IWSubstring("a"), value));
  EXPECT_EQ(value, 16777217);
}

TEST_F(TabularReaders, ActivityHeaderAndPrecision) {
  Write("ID Value\na 1.23456789012345\n");
  Activity_Data_From_File<double> data;
  ASSERT_EQ(data.read_activity_data(Filename()), 1);
  double value;
  EXPECT_FALSE(data.get_activity(const_IWSubstring("ID"), value));
  ASSERT_TRUE(data.get_activity(const_IWSubstring("a"), value));
  EXPECT_DOUBLE_EQ(value, 1.23456789012345);
}

TEST_F(TabularReaders, ActivityAutoSeparatorAndSelectedColumn) {
  _path += ".csv";
  Write("ID,Other,Value\na,99,1.25\n");
  const char* argv[] = {"test", "-A", "sep=auto", "-A", "col=3", "-A", _path.c_str()};
  Command_Line cl(7, const_cast<char**>(argv), "A:");
  Activity_Data_From_File<double> data;
  ASSERT_TRUE(data.construct_from_command_line(cl, 'A', 0));
  double value;
  ASSERT_TRUE(data.get_activity(const_IWSubstring("a"), value));
  EXPECT_DOUBLE_EQ(value, 1.25);
}

TEST_F(TabularReaders, TableHeaderIdentifiersMissingAndRebuild) {
  Write("Name A B\na 1 N/A\nb 3 4\n");
  IW_Tabular_Data<double> table;
  table.set_has_header(1);
  table.set_first_column_is_identifier(1);
  ASSERT_TRUE(table.build(_path.c_str(), ' '));
  EXPECT_EQ(table.nrows(), 2);
  EXPECT_EQ(table.ncols(), 2);
  EXPECT_EQ(table.ids()[1], "b");
  EXPECT_DOUBLE_EQ(table.value(0, 1), 0);
  EXPECT_DOUBLE_EQ(table.value(1, 1), 4);
  ASSERT_TRUE(table.build(_path.c_str(), ' '));
  EXPECT_DOUBLE_EQ(table.value(1, 0), 3);
}

TEST_F(TabularReaders, TableRejectsRaggedRows) {
  IW_Tabular_Data<int> table;
  Write("1 2\n3\n");
  EXPECT_FALSE(table.build(_path.c_str(), ' '));
  Write("1 2\n3 4 5\n");
  EXPECT_FALSE(table.build(_path.c_str(), ' '));
}

TEST_F(TabularReaders, TableCsvHeader) {
  Write("A,B\n1,2\n3,4\n");
  IW_Tabular_Data<int> table;
  table.set_has_header(1);
  ASSERT_TRUE(table.build(_path.c_str(), ','));
  EXPECT_EQ(table.value(1, 1), 4);
}

TEST_F(TabularReaders, RemoveColumnCompactsAllRows) {
  IW_Tabular_Data<int> table;
  ASSERT_TRUE(table.resize(2, 3));
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 3; ++c) {
      table.set(r, c, r * 10 + c);
    }
  }
  ASSERT_TRUE(table.remove_column(1));
  EXPECT_EQ(table.ncols(), 2);
  EXPECT_EQ(table.value(0, 0), 0);
  EXPECT_EQ(table.value(0, 1), 2);
  EXPECT_EQ(table.value(1, 0), 10);
  EXPECT_EQ(table.value(1, 1), 12);
  EXPECT_FALSE(table.remove_column(2));
}
}  // namespace

TEST_F(TabularReaders, TableWritesHeadersIdentifiersAndMissingValuesAfterRemoval) {
  Write("Name A B C\na 1 2 N/A\nb 3 4 5\n");
  IW_Tabular_Data<double> table;
  table.set_has_header(1);
  table.set_first_column_is_identifier(1);
  ASSERT_TRUE(table.build(_path.c_str(), ' '));
  ASSERT_TRUE(table.remove_column(1));
  {
    IWString_and_File_Descriptor output;
    ASSERT_TRUE(output.open((_path + ".out").c_str()));
    ASSERT_TRUE(table.do_write(output, ' '));
  }
  std::ifstream input(_path + ".out");
  std::string contents((std::istreambuf_iterator<char>(input)), {});
  EXPECT_EQ(contents, "Name A C\na 1 N/A\nb 3 5\n");
}

TEST_F(TabularReaders, NumericSubsequentFileMayHaveItsOwnHeader) {
  Write("ID Value\na 1\n");
  Numeric_Data_From_File<float> data;
  ASSERT_EQ(data.read_data(Filename()), 1);
  Write("ID Value\nb 2\n");
  EXPECT_EQ(data.read_data(Filename()), 2);
}
