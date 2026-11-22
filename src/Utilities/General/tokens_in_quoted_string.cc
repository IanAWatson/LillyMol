#include "Foundational/iwstring/iwstring.h"

#include "tokens_in_quoted_string.h"

int
CountDelimitedTokens(const const_IWSubstring& buffer, const DelimitedTextOptions& options,
                     uint64_t line_number, IWString* translated, char output_separator) {
  if (translated != nullptr) {
    *translated = buffer;
  }
  if (buffer.empty()) {
    return 0;
  }

  const bool keep_empty = options._empty_fields || options._quoted;
  int count = keep_empty ? 1 : 0;
  bool in_word = false;
  bool in_quote = false;
  bool escaped = false;
  for (int i = 0; i < buffer.length(); ++i) {
    const char c = buffer[i];
    bool separator = false;
    if (escaped) {
      escaped = false;
    } else if (options._backslash_escapes && c == '\\') {
      escaped = true;
    } else if (options._quoted && c == '"') {
      in_quote = !in_quote;
    } else {
      separator = !in_quote && c == options._separator;
    }

    if (separator) {
      if (keep_empty) {
        ++count;
      }
      in_word = false;
      if (translated != nullptr) {
        (*translated)[i] = output_separator;
      }
    } else if (!in_word) {
      if (!keep_empty) {
        ++count;
      }
      in_word = true;
    }
  }
  if (in_quote && line_number != 0) {
    std::cerr << "tcount: warning unmatched quote on line " << line_number
              << ": '" << buffer << "'\n";
  }
  return count;
}
