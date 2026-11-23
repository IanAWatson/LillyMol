#include <filesystem>
#include <optional>
#include <string>

#include "misc.h"

namespace iwmisc {

namespace fs = std::filesystem;

std::optional<bool>
SameFile(const IWString& lhs, const IWString& rhs, std::error_code& error) {
  error.clear();
  if (lhs.empty() || rhs.empty() || lhs.contains('\0') || rhs.contains('\0')) {
    error = std::make_error_code(std::errc::invalid_argument);
    return std::nullopt;
  }
  const fs::path a(std::string(lhs.data(), lhs.length()));
  const fs::path b(std::string(rhs.data(), rhs.length()));
  // equivalent() reports missing files differently across standard libraries;
  // libc++ can use "operation not supported" instead of ENOENT. Determine
  // existence with status() first, retaining errors other than missing paths.
  const fs::file_status status_a = fs::status(a, error);
  if (error && error != std::errc::no_such_file_or_directory) {
    return std::nullopt;
  }
  error.clear();
  const fs::file_status status_b = fs::status(b, error);
  if (error && error != std::errc::no_such_file_or_directory) {
    return std::nullopt;
  }
  error.clear();
  if (fs::exists(status_a) && fs::exists(status_b)) {
    const bool same = fs::equivalent(a, b, error);
    if (error) {
      return std::nullopt;
    }
    return same;
  }
  const auto absolute_a = fs::absolute(a, error);
  if (error) {
    return std::nullopt;
  }
  const auto absolute_b = fs::absolute(b, error);
  if (error) {
    return std::nullopt;
  }
  const auto resolved_a = fs::weakly_canonical(absolute_a, error);
  if (error) {
    return std::nullopt;
  }
  const auto resolved_b = fs::weakly_canonical(absolute_b, error);
  if (error) {
    return std::nullopt;
  }
  return resolved_a == resolved_b;
}

std::optional<bool>
SameFile(const IWString& lhs, const IWString& rhs) {
  std::error_code error;
  return SameFile(lhs, rhs, error);
}

std::optional<IWString>
FileOrPath(const IWString& outer_file_name,
           const IWString& inner_file_name) {
  const std::string possible_file(inner_file_name.data(), inner_file_name.length());

  if (fs::exists(fs::path(possible_file))) {
    // std::cerr << "File name " << possible_file << " alread exists, good\n";
    return inner_file_name;
  }

  // Try query_file_name in same directory as proto_file_name.

  // cerr << "outer_file_name '" << outer_file_name << "'\n";
  fs::path new_name = fs::path(std::string(outer_file_name.data(), outer_file_name.length())).replace_filename(possible_file);
  // cerr << "With directory '" << new_name.string() << "'\n";

  if (fs::exists(new_name)) {
    return new_name.string();
  }

  return std::nullopt;
}

}  // namespace iwmisc
