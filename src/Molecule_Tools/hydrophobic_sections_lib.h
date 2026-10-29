#ifndef MOLECULE_TOOLS_HYDROPHOBIC_SECTIONS_LIB_H_
#define MOLECULE_TOOLS_HYDROPHOBIC_SECTIONS_LIB_H_

#include <vector>

class Molecule;

namespace hydrophobic {

// Connected components of an atom-selection mask. Atom selection (queries,
// simple rules, etc.) and molecular preprocessing are the caller's choices.
class Sections {
 private:
  std::vector<int> _membership;
  std::vector<int> _sizes;
  int _atoms_in_sections = 0;
  int _largest_section = 0;

 public:
  // eligible must describe all m.natoms() atoms; any nonzero value selects an
  // atom. Small components are discarded, and remaining labels are consecutive
  // positive integers. Zero means unselected or discarded. Does not modify m or
  // eligible. Empty molecules are valid; invalid arguments return false.
  bool Build(const Molecule& m, const int* eligible, int minimum_size = 1);

  const std::vector<int>& membership() const { return _membership; }
  const std::vector<int>& sizes() const { return _sizes; }
  int number_sections() const { return _sizes.size(); }
  int atoms_in_sections() const { return _atoms_in_sections; }
  int largest_section() const { return _largest_section; }
};

}  // namespace hydrophobic

#endif  // MOLECULE_TOOLS_HYDROPHOBIC_SECTIONS_LIB_H_
