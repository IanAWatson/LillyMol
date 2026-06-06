#include <iostream>

#include "util.h"

#include "istream_and_type.h"

namespace lillymol {

using std::cerr;

std::optional<uint64_t>
CountMoleculesInFile(const IWString& fname, FileType ftype) {
  if (ftype == FILE_TYPE_INVALID) {
    ftype = discern_file_type_from_name(fname);
    if (ftype == FILE_TYPE_INVALID) {
      cerr << "CountMoleculesInFile:cannot discern file type '" << fname << "'\n";
      return std::nullopt;
    }
  }

  data_source_and_type<Molecule> input(ftype, fname);
  if (! input.good()) {
    return std::nullopt;
  }

  return input.molecules_remaining();
}

}  // namespace lillymol
