#include "Molecule_Tools/hydrophobic_sections_lib.h"

#include <algorithm>

#include "Molecule_Lib/molecule.h"

namespace hydrophobic {

bool
Sections::Build(const Molecule& m, const int* eligible, int minimum_size) {
  _membership.resize(m.natoms());
  _sizes.clear();
  _atoms_in_sections = 0;
  _largest_section = 0;
  if (minimum_size < 1 || (m.natoms() > 0 && eligible == nullptr)) {
    std::fill(_membership.begin(), _membership.end(), 0);
    return false;
  }
  for (int atom = 0; atom < m.natoms(); ++atom) {
    // Read before writing so a prior membership vector can serve as the mask.
    _membership[atom] = eligible[atom] ? -1 : 0;
  }

  std::vector<atom_number_t> component;
  for (atom_number_t root = 0; root < m.natoms(); ++root) {
    if (_membership[root] != -1) {
      continue;
    }
    const int label = _sizes.size() + 1;
    component.clear();
    component.push_back(root);
    _membership[root] = label;
    for (size_t head = 0; head < component.size(); ++head) {
      const atom_number_t atom = component[head];
      for (const Bond* bond : m[atom]) {
        const atom_number_t other = bond->other(atom);
        if (_membership[other] == -1) {
          _membership[other] = label;
          component.push_back(other);
        }
      }
    }
    const int size = component.size();
    if (size < minimum_size) {
      for (atom_number_t atom : component) {
        _membership[atom] = 0;
      }
      continue;
    }
    _sizes.push_back(size);
    _atoms_in_sections += size;
    _largest_section = std::max(_largest_section, size);
  }
  return true;
}

}  // namespace hydrophobic
