
#include <cstdint>
#include <optional>

#include "Foundational/iwstring/iwstring.h"

#include "iwmtypes.h"
#include "istream_and_type.h"

namespace lillymol {
// Given a file name, return the number of molecules in that file.
std::optional<uint64_t> CountMoleculesInFile(const IWString& fname, FileType fype);
}  // namespace lillymol
