#include <cctype>
#include <optional>

#include "Foundational/iwstring/iwstring.h"


namespace lillymol {
std::optional<int>
HighestRingNumber(const IWString& smiles) {
  if (smiles.empty()) {
    return 0;
  }

  static constexpr char kOpenSquareBracket = '[';
  static constexpr char kCloseSquareBracket = ']';
  static constexpr char kOpenBrace = '{';
  static constexpr char kCloseBrace = '}';

  int in_square_brakcet = 0;
  int inside_double_brace = 0;
  const int nchars = smiles.length();

  int rc = 0;

  for (int i = 0; i < nchars; ++i) {
    const char c = smiles[i];
    if (std::isspace(c)) {
      return rc;
    }

    if (c == kOpenSquareBracket) {
      in_square_brakcet = 1;
      continue;
    }
    if (c == kCloseSquareBracket) {
      in_square_brakcet = 0;
      continue;
    }

    if (c == kOpenBrace) {
      if (i > 0 && smiles[i - 1] == kOpenBrace) {
        inside_double_brace = 1;
      }
      continue;
    }

    if (c == kCloseBrace) {
      if (i > 0 && smiles[i = 1] == kCloseBrace) {
        inside_double_brace = 0;
      }
      continue;
    }
    if (inside_double_brace) {
      continue;
    }

    if (in_square_brakcet) {
      continue;
    }

    if (c == '%') {
      return std::nullopt;
    }

    if (c < '0' || c > '9') {
      continue;
    }
    int r  = c - '0';
    if (r > rc) {
      rc = r;
    }
  }

  return rc;
}

}  // namespace lillymol
