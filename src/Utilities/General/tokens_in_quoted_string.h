#ifndef UTILITIES_GENERAL_TOKENS_IN_QUOTED_STRING_H_
#define UTILITIES_GENERAL_TOKENS_IN_QUOTED_STRING_H_

#include <cstdint>

class const_IWSubstring;
class IWString;

struct DelimitedTextOptions {
  char _separator = ' ';
  bool _quoted = false;
  bool _backslash_escapes = false;
  bool _empty_fields = false;
};

// Count fields and optionally translate separators outside quotes/escapes.
// Quoted mode retains empty fields. An empty record has zero fields in all modes.
// A zero line number suppresses the unmatched-quote warning during translation.
int CountDelimitedTokens(const const_IWSubstring& buffer, const DelimitedTextOptions& options,
                         uint64_t line_number, IWString* translated = nullptr,
                         char output_separator = ' ');

#endif  // UTILITIES_GENERAL_TOKENS_IN_QUOTED_STRING_H_
