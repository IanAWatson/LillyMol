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

// Given a SAFE smiles, identify the unbalanced ring closures.
// by convention, the unbalanced rings that have been added are 
// all 2 digits and will be preceded by a % sign. This makes
// parsing easy.
int UnbalancedRingNumbers(const IWString& smiles, resizable_array<int>& ring_numbers);

// Given a smiles with isotopic labels, convert those to ring openings.
// [nnnC]C becomes [nnnC]%rC where `r` is derived from `ring_number`
// Note that `new_smiles` will be appended to.
int IsotopeToRingOpening(const IWString& smiles, int& ring_number, IWString& new_smiles);

}  // namespace lillymol

#endif
