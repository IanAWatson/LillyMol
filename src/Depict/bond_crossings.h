#ifndef DEPICT_BOND_CROSSINGS_H_
#define DEPICT_BOND_CROSSINGS_H_

#include <vector>

class Molecule;

namespace bond_crossings {

enum class Status { kOk, kEmpty, kDisconnected, kNot2D, kInvalidCoordinates, kInvalidOptions };

class Options {
 private:
  bool _diagnostics = false;
  double _near_miss_fraction = 0.1;
 public:
  bool diagnostics() const { return _diagnostics; }
  void set_diagnostics(bool value) { _diagnostics = value; }
  double near_miss_fraction() const { return _near_miss_fraction; }
  void set_near_miss_fraction(double value) { _near_miss_fraction = value; }
};

class Result;
Status Analyse(Molecule& mol, const Options& options, Result& result);

class RingSystemResult {
 private:
  int _atoms = 0;
  int _rings = 0;
  int _ring_ring = 0;
  int _ring_chain = 0;
  friend Status Analyse(Molecule&, const Options&, Result&);
 public:
  RingSystemResult(int atoms, int rings) : _atoms(atoms), _rings(rings) {}
  int atoms() const { return _atoms; }
  int rings() const { return _rings; }
  int ring_ring() const { return _ring_ring; }
  int ring_chain() const { return _ring_chain; }
  int crossings() const { return _ring_ring + _ring_chain; }
  double score() const {
    const double count = crossings();
    return count * count * count;
  }
};

// Analyses existing XY coordinates only. Ring systems sharing only a spiro
// atom remain separate. All bonds are centre-line segments, irrespective of
// bond order or wedge/hash style. Shared endpoints are never crossings.
class Result {
 private:
  int _atoms = 0;
  int _ring_atoms = 0;
  int _ring_ring = 0;
  int _ring_chain = 0;
  int _inter_system = 0;
  int _crossed_bonds = 0;
  int _systems_crossed = 0;
  int _max_crossings = 0;
  int _topological_diameter = 0;
  double _score = 0;
  double _worst_density = 0;
  double _min_distance = 0;
  double _max_distance = 0;
  int _degenerate_bonds = 0;
  int _overlaps = 0;
  int _near_misses = 0;
  int _atoms_on_ring_bonds = 0;
  int _coincident_atoms = 0;
  double _min_nonbonded = 0;
  std::vector<RingSystemResult> _systems;
  friend Status Analyse(Molecule&, const Options&, Result&);
 public:
  int atoms() const { return _atoms; }
  int ring_atoms() const { return _ring_atoms; }
  int ring_ring() const { return _ring_ring; }
  int ring_chain() const { return _ring_chain; }
  int inter_system() const { return _inter_system; }
  int crossings() const { return _ring_ring + _ring_chain + _inter_system; }
  int crossed_bonds() const { return _crossed_bonds; }
  int systems_crossed() const { return _systems_crossed; }
  int max_crossings() const { return _max_crossings; }
  int topological_diameter() const { return _topological_diameter; }
  double score() const { return _score; }
  double score_per_atom() const { return _atoms ? _score / _atoms : 0; }
  double score_per_ring_atom() const { return _ring_atoms ? _score / _ring_atoms : 0; }
  double score_per_system() const { return _systems.empty() ? 0 : _score / _systems.size(); }
  double worst_density() const { return _worst_density; }
  double min_distance() const { return _min_distance; }
  double max_distance() const { return _max_distance; }
  double min_distance_per_atom() const { return _atoms ? _min_distance / _atoms : 0; }
  double max_distance_per_atom() const { return _atoms ? _max_distance / _atoms : 0; }
  double extent_per_topology() const {
    return _topological_diameter ? _max_distance / _topological_diameter : 0;
  }
  const std::vector<RingSystemResult>& systems() const { return _systems; }
  int degenerate_bonds() const { return _degenerate_bonds; }
  int overlaps() const { return _overlaps; }
  int near_misses() const { return _near_misses; }
  int atoms_on_ring_bonds() const { return _atoms_on_ring_bonds; }
  int coincident_atoms() const { return _coincident_atoms; }
  // Zero if no pair of nonbonded atoms exists.
  double min_nonbonded() const { return _min_nonbonded; }
};

// Resets result on every call. Does not generate coordinates, delete atoms,
// or change molecular connectivity. Like other LillyMol perceptions, ring and
// topological-distance caches may be populated. Requires finite coordinates,
// one connected component, and highest_coordinate_dimensionality() == 2.
Status Analyse(Molecule& mol, const Options& options, Result& result);
const char* StatusMessage(Status status);

}  // namespace bond_crossings
#endif
