#include <filesystem>
#include <fstream>
#include <system_error>

#include "gtest/gtest.h"

#include "misc.h"

namespace {
namespace fs = std::filesystem;

class SameFileTest : public testing::Test {
 protected:
  fs::path _directory;
  void SetUp() override {
    _directory = fs::path(testing::TempDir()) / "same_file_tests";
    fs::create_directories(_directory);
    std::ofstream(_directory / "first") << "same content\n";
    std::ofstream(_directory / "second") << "same content\n";
  }
  void TearDown() override { fs::remove_all(_directory); }
  void Check(const fs::path& a, const fs::path& b, bool expected) {
    // Successful comparisons must clear an error supplied by the caller.
    std::error_code error = std::make_error_code(std::errc::invalid_argument);
    EXPECT_EQ(iwmisc::SameFile(IWString(a.string()), IWString(b.string()), error), expected);
    EXPECT_FALSE(error) << error.message();
    EXPECT_EQ(iwmisc::SameFile(IWString(a.string()), IWString(b.string())), expected);
  }
};

TEST_F(SameFileTest, ExistingFilesUseIdentityNotContents) {
  Check(_directory / "first", _directory / "first", true);
  Check(_directory / "first", _directory / "second", false);
}

TEST_F(SameFileTest, RelativeAndAbsolutePaths) {
  const auto relative = fs::relative(_directory / "first", fs::current_path());
  Check(relative, fs::absolute(_directory / "first"), true);
}

TEST_F(SameFileTest, LexicalAliasesOfNewPaths) {
  Check(_directory / "new", _directory / "." / "new", true);
  Check(_directory / "subdir" / ".." / "new", _directory / "new", true);
  Check(_directory / "uncreated" / "one", _directory / "uncreated" / "two", false);
}

TEST_F(SameFileTest, MissingBasenameAgainstAbsolutePath) {
  const fs::path relative("iwmisc_same_file_nonexistent_basename");
  ASSERT_FALSE(fs::exists(relative));
  Check(relative, fs::current_path() / relative, true);
}

TEST_F(SameFileTest, ExistingAndNewFilesDiffer) {
  Check(_directory / "first", _directory / "new", false);
  Check(_directory / "new", _directory / "first", false);
}

TEST_F(SameFileTest, HardLink) {
  fs::create_hard_link(_directory / "first", _directory / "hardlink");
  Check(_directory / "first", _directory / "hardlink", true);
}

TEST_F(SameFileTest, SymbolicLink) {
  fs::create_symlink(_directory / "first", _directory / "symlink");
  Check(_directory / "first", _directory / "symlink", true);
}

TEST_F(SameFileTest, ExistingSymlinkedParentForNewFiles) {
  fs::create_directory(_directory / "parent");
  fs::create_directory_symlink(_directory / "parent", _directory / "alias");
  Check(_directory / "parent" / "new", _directory / "alias" / "new", true);
}

TEST_F(SameFileTest, EmptyFilenameIsAnError) {
  std::error_code error;
  EXPECT_EQ(std::nullopt, iwmisc::SameFile(IWString(), IWString("file"), error));
  EXPECT_EQ(error, std::errc::invalid_argument);
  EXPECT_EQ(iwmisc::SameFile(IWString("file"), IWString()), std::nullopt);
}

TEST_F(SameFileTest, SymlinkLoopIsAnError) {
  fs::create_symlink("loop2", _directory / "loop1");
  fs::create_symlink("loop1", _directory / "loop2");
  std::error_code error;
  EXPECT_EQ(std::nullopt, iwmisc::SameFile(IWString((_directory / "loop1").string()),
                              IWString((_directory / "first").string()), error));
  EXPECT_TRUE(error);
  EXPECT_EQ(iwmisc::SameFile(IWString((_directory / "loop1").string()),
                            IWString((_directory / "first").string())), std::nullopt);
}

TEST_F(SameFileTest, DashIsAnOrdinaryFilename) {
  Check(fs::path("-"), fs::current_path() / "-", true);
}
}  // namespace
