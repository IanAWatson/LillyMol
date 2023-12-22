#ifndef MOLECULE_TOOLS_HIGHEST_RING_NUMBER_H
#define MOLECULE_TOOLS_HIGHEST_RING_NUMBER_H

#include <optional>

#include "Foundational/iwstring/iwstring.h"

namespace lillymol {

// Examine `smiles` and return the highest ring number used in the smiles.
// If there are no rings, 0 will be returned.
// If there are too many rings, the % sign is encountered, nullopt is returned.
// There is no smiles interpretation.
std::optional<int> HighestRingNumber(const IWString& smiles);

}  // namespace lillymol

#endif
