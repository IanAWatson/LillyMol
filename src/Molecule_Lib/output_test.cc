#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "gtest/gtest.h"

#include "Foundational/cmdline/cmdline.h"
#include "molecule.h"
#include "output.h"

namespace {
namespace fs = std::filesystem;

// Redirect stdout for a scope, including tests that deliberately make it fail.
class StdoutRedirect {
 private:
  std::streambuf* _previous;
  std::ios::iostate _state;
 public:
  explicit StdoutRedirect(std::streambuf* buffer)
      : _previous(std::cout.rdbuf()), _state(std::cout.rdstate()) {
    std::cout.rdbuf(buffer);
    std::cout.clear();
  }
  ~StdoutRedirect() {
    std::cout.rdbuf(_previous);
    std::cout.clear(_state);
  }
};

class RejectWritesBuffer : public std::streambuf {
 public:
  std::streamsize xsputn(const char*, std::streamsize) override { return 0; }
  int_type overflow(int_type) override { return traits_type::eof(); }
};

class FailingBuffer : public std::stringbuf {
 public:
  int sync() override { return -1; }
};

class MoleculeOutputTest : public testing::Test {
 protected:
  fs::path _directory;
  Molecule _molecule;
  void SetUp() override {
    _directory = fs::path(testing::TempDir()) / "molecule_output_tests";
    fs::create_directories(_directory);
    ASSERT_TRUE(_molecule.build_from_smiles("CC example"));
  }
  void TearDown() override {
    lillymol_csv::set_csv_separator(',');
    fs::remove_all(_directory);
  }
  IWString Name(const char* basename) const {
    return IWString((_directory / basename).string());
  }
  std::string Read(const char* basename) const {
    std::ifstream input(_directory / basename);
    return std::string(std::istreambuf_iterator<char>(input), {});
  }
};

TEST_F(MoleculeOutputTest, MultipleFormats) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SDF));
  ASSERT_TRUE(output.new_stem(Name("both")));
  ASSERT_TRUE(output.write(_molecule));
  EXPECT_EQ(output.molecules_written(), 1);
  EXPECT_TRUE(output.do_close());
  EXPECT_NE(Read("both.smi").find("example"), std::string::npos);
  EXPECT_NE(Read("both.sdf").find("$$$$"), std::string::npos);
}

TEST_F(MoleculeOutputTest, PartialOpenFailureLeavesNoStreams) {
  fs::create_directory(_directory / "blocked.csv");
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_CSV));
  EXPECT_FALSE(output.new_stem(Name("blocked")));
  EXPECT_EQ(output.number_elements(), 0);
  EXPECT_TRUE(output.ok());
  EXPECT_FALSE(output.write(_molecule));
  EXPECT_TRUE(output.new_stem(Name("recovered")));
  EXPECT_TRUE(output.write(_molecule));
  EXPECT_TRUE(output.do_close());
  EXPECT_EQ(Read("recovered.csv"), "SMILES,ID\nCC,example\n");
}

TEST_F(MoleculeOutputTest, FirstFormatOpenFailure) {
  fs::create_directory(_directory / "blocked.smi");
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_CSV));
  EXPECT_FALSE(output.new_stem(Name("blocked")));
  EXPECT_FALSE(fs::exists(_directory / "blocked.csv"));
  EXPECT_EQ(output.number_elements(), 0);
}

TEST_F(MoleculeOutputTest, StdoutThenFileThenStdout) {
  std::ostringstream captured;
  StdoutRedirect redirect(captured.rdbuf());
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem("-"));
  ASSERT_TRUE(output.write(_molecule));
  EXPECT_EQ(&output.stream_for_type(FILE_TYPE_SMI), &std::cout);
  ASSERT_TRUE(output.new_stem(Name("file")));
  ASSERT_TRUE(output.write(_molecule));
  ASSERT_TRUE(output.new_stem("-"));
  ASSERT_TRUE(output.write(_molecule));
  EXPECT_EQ(output.molecules_written(), 3);
  ASSERT_TRUE(output.do_close());
  EXPECT_EQ(captured.str(), "CC example\nCC example\n");
  EXPECT_EQ(Read("file.smi"), "CC example\n");
  EXPECT_FALSE(output.active());
  EXPECT_FALSE(output.good());
  EXPECT_EQ(output.first_output_type(), FILE_TYPE_INVALID);
  EXPECT_FALSE(output.write(_molecule));
}

TEST_F(MoleculeOutputTest, StdoutFlushAndCloseReportFailure) {
  FailingBuffer buffer;
  StdoutRedirect redirect(&buffer);
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem("-"));
  ASSERT_TRUE(output.write(_molecule));
  EXPECT_FALSE(output.do_flush());
  EXPECT_FALSE(output.do_close());
  EXPECT_FALSE(output.active());
  EXPECT_TRUE(output.do_close());
}

TEST_F(MoleculeOutputTest, RepeatedFailedStdoutWritesAreRejected) {
  RejectWritesBuffer buffer;
  StdoutRedirect redirect(&buffer);
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem("-"));
  EXPECT_FALSE(output.write(_molecule));
  EXPECT_FALSE(output.write(_molecule));
  EXPECT_EQ(output.molecules_written(), 0);
}

TEST_F(MoleculeOutputTest, FileFlushAndCloseReportFailure) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem(Name("file")));
  output.stream_for_type(FILE_TYPE_SMI).setstate(std::ios::badbit);
  EXPECT_FALSE(output.do_flush());
  EXPECT_FALSE(output.do_close());
  EXPECT_FALSE(output.active());
  EXPECT_TRUE(output.do_close());
}

TEST_F(MoleculeOutputTest, NewStemReportsPreviousWriteFailure) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem(Name("first")));
  output.stream_for_type(FILE_TYPE_SMI).setstate(std::ios::badbit);
  EXPECT_FALSE(output.new_stem(Name("second")));
  EXPECT_FALSE(fs::exists(_directory / "second.smi"));
  EXPECT_TRUE(output.ok());
}

TEST_F(MoleculeOutputTest, CsvHeaderUsesConfiguredSeparator) {
  lillymol_csv::set_csv_separator('\t');
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_CSV));
  ASSERT_TRUE(output.new_stem(Name("table")));
  ASSERT_TRUE(output.write(_molecule));
  ASSERT_TRUE(output.do_close());
  EXPECT_EQ(Read("table.csv"), "SMILES\tID\nCC\texample\n");
}

TEST_F(MoleculeOutputTest, CsvHeaderFailureOnStdout) {
  RejectWritesBuffer buffer;
  StdoutRedirect redirect(&buffer);
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_CSV));
  EXPECT_FALSE(output.new_stem("-"));
  EXPECT_TRUE(output.ok());
}

TEST_F(MoleculeOutputTest, SplitFilesAndCsvHeaders) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_CSV));
  ASSERT_TRUE(output.set_molecules_per_file(2));
  ASSERT_TRUE(output.new_stem(Name("part")));
  EXPECT_TRUE(output.good());
  for (int i = 0; i < 3; ++i) {
    ASSERT_TRUE(output.write(_molecule));
  }
  EXPECT_EQ(output.molecules_written(), 3);
  ASSERT_TRUE(output.do_close());
  EXPECT_EQ(Read("part0.csv"), "SMILES,ID\nCC,example\nCC,example\n");
  EXPECT_EQ(Read("part1.csv"), "SMILES,ID\nCC,example\n");
}

TEST_F(MoleculeOutputTest, SplitFilesKeepSuffix) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.set_molecules_per_file(1));
  ASSERT_TRUE(output.new_stem(Name("part."), 1));
  ASSERT_TRUE(output.write(_molecule));
  ASSERT_TRUE(output.write(_molecule));
  ASSERT_TRUE(output.do_close());
  EXPECT_EQ(Read("part.0.smi"), "CC example\n");
  EXPECT_EQ(Read("part.1.smi"), "CC example\n");
}

TEST_F(MoleculeOutputTest, RejectInvalidSplitConfiguration) {
  Molecule_Output_Object output;
  EXPECT_FALSE(output.set_molecules_per_file(-1));
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem(Name("file")));
  ASSERT_TRUE(output.write(_molecule));
  EXPECT_FALSE(output.set_molecules_per_file(2));
  EXPECT_EQ(output.molecules_per_file(), 0);
  EXPECT_TRUE(output.write(_molecule));
}

TEST_F(MoleculeOutputTest, EnableSplittingAfterOpenBeforeWrite) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem(Name("part")));
  ASSERT_TRUE(output.set_molecules_per_file(1));
  ASSERT_TRUE(output.write(_molecule));
  ASSERT_TRUE(output.do_close());
  EXPECT_EQ(Read("part0.smi"), "CC example\n");
}

TEST_F(MoleculeOutputTest, DisableSplittingBeforeWrite) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.set_molecules_per_file(2));
  ASSERT_TRUE(output.new_stem(Name("part")));
  ASSERT_TRUE(output.set_molecules_per_file(0));
  ASSERT_TRUE(output.write(_molecule));
  ASSERT_TRUE(output.do_close());
  EXPECT_EQ(Read("part.smi"), "CC example\n");
}

TEST_F(MoleculeOutputTest, CannotSplitStdout) {
  std::ostringstream captured;
  StdoutRedirect redirect(captured.rdbuf());
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem("-"));
  EXPECT_FALSE(output.set_molecules_per_file(1));
  ASSERT_TRUE(output.do_close());
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.set_molecules_per_file(1));
  EXPECT_FALSE(output.new_stem("-"));
}

TEST_F(MoleculeOutputTest, OverwriteChecksResolveFileAliases) {
  IWString input = Name("input.smi");
  std::ofstream(input.AsString()) << "CC example\n";
  fs::create_hard_link(_directory / "input.smi", _directory / "hard.smi");
  fs::create_symlink(_directory / "input.smi", _directory / "symbolic.smi");
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  EXPECT_TRUE(output.would_use_name(Name("./input"), input.null_terminated_chars()));
  EXPECT_TRUE(output.would_use_name(Name("hard"), input.null_terminated_chars()));
  EXPECT_TRUE(output.would_use_name(Name("symbolic"), input.null_terminated_chars()));
  EXPECT_FALSE(output.would_use_name(Name("other"), input.null_terminated_chars()));
  EXPECT_FALSE(output.would_use_name("-", input.null_terminated_chars()));
}

TEST_F(MoleculeOutputTest, SplitOverwriteCheckUsesGeneratedFilename) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.set_molecules_per_file(2));
  EXPECT_TRUE(output.would_use_name(Name("part"), Name("part0.smi").null_terminated_chars()));
}

TEST_F(MoleculeOutputTest, OverwriteComparisonFailureIsReportedAsConflict) {
  fs::create_symlink("loop", _directory / "loop");
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  EXPECT_TRUE(output.would_use_name(Name("loop/file"), "input.smi"));
}

TEST_F(MoleculeOutputTest, MissingStreamDoesNotFallBackToStdout) {
  std::ostringstream captured;
  StdoutRedirect redirect(captured.rdbuf());
  Molecule_Output_Object output;
  output.stream_for_type(FILE_TYPE_SMI) << "unwanted";
  EXPECT_TRUE(captured.str().empty());
}

TEST_F(MoleculeOutputTest, BufferedWriteFailureIsDetectedAtClose) {
  if (!fs::exists("/dev/full")) {
    GTEST_SKIP() << "Requires /dev/full";
  }
  fs::create_symlink("/dev/full", _directory / "full.smi");
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.new_stem(Name("full")));
  output.write(_molecule);  // Failure can be delayed until the buffer is flushed.
  EXPECT_FALSE(output.do_close());
  EXPECT_FALSE(output.active());
}

TEST_F(MoleculeOutputTest, NamedSplitFiles) {
  Molecule_Output_Object output;
  ASSERT_TRUE(output.add_output_type(FILE_TYPE_SMI));
  ASSERT_TRUE(output.set_molecules_per_file(1));
  output.set_name_token_for_file_name(0);
  ASSERT_TRUE(output.new_stem(Name("unused")));
  _molecule.set_name(Name("named_first"));
  ASSERT_TRUE(output.write(_molecule));
  _molecule.set_name(Name("named_second"));
  ASSERT_TRUE(output.write(_molecule));
  ASSERT_TRUE(output.do_close());
  EXPECT_EQ(output.molecules_written(), 2);
  EXPECT_NE(Read("named_first.smi").find("CC "), std::string::npos);
  EXPECT_NE(Read("named_second.smi").find("CC "), std::string::npos);
}

TEST_F(MoleculeOutputTest, TypedStreamConstructorSetsType) {
  ofstream_and_type stream(FILE_TYPE_SMI);
  ASSERT_TRUE(stream.valid());
  IWString fname = Name("typed.smi");
  ASSERT_TRUE(stream.open(fname));
  ASSERT_TRUE(stream.write_molecule(&_molecule));
  stream.close();
  EXPECT_EQ(Read("typed.smi"), "CC example\n");
}
}  // namespace
