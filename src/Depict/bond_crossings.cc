#include "Depict/bond_crossings.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

namespace bond_crossings {
namespace {

struct Point {
  double _x;
  double _y;
};

struct Segment {
  int _a1;
  int _a2;
  Point _p;
  Point _q;
  double _length;
  int _system = -1;
};

double SquaredDistance(Point p, Point q) {
  const double dx = p._x - q._x;
  const double dy = p._y - q._y;
  return dx * dx + dy * dy;
}

double Cross(Point a, Point b, Point c) {
  return (b._x - a._x) * (c._y - a._y) - (b._y - a._y) * (c._x - a._x);
}

bool Opposite(double a, double b, double tolerance) {
  return (a > tolerance && b < -tolerance) || (b > tolerance && a < -tolerance);
}

bool ProperCrossing(const Segment& a, const Segment& b, double tolerance) {
  return Opposite(Cross(a._p, a._q, b._p), Cross(a._p, a._q, b._q), tolerance * a._length) &&
         Opposite(Cross(b._p, b._q, a._p), Cross(b._p, b._q, a._q), tolerance * b._length);
}

bool BoxesClose(const Segment& a, const Segment& b, double margin) {
  return std::max(a._p._x, a._q._x) + margin >= std::min(b._p._x, b._q._x) &&
         std::max(b._p._x, b._q._x) + margin >= std::min(a._p._x, a._q._x) &&
         std::max(a._p._y, a._q._y) + margin >= std::min(b._p._y, b._q._y) &&
         std::max(b._p._y, b._q._y) + margin >= std::min(a._p._y, a._q._y);
}

double PointSegmentDistanceSquared(Point p, const Segment& b) {
  if (b._length == 0) {
    return SquaredDistance(p, b._p);
  }
  const double dx = b._q._x - b._p._x;
  const double dy = b._q._y - b._p._y;
  const double t = std::clamp(((p._x - b._p._x) * dx + (p._y - b._p._y) * dy) /
                                (b._length * b._length), 0.0, 1.0);
  return SquaredDistance(p, {b._p._x + t * dx, b._p._y + t * dy});
}

bool Overlap(const Segment& a, const Segment& b, double tolerance) {
  if (std::abs(Cross(a._p, a._q, b._p)) > tolerance * a._length ||
      std::abs(Cross(a._p, a._q, b._q)) > tolerance * a._length) {
    return false;
  }
  const double ux = (a._q._x - a._p._x) / a._length;
  const double uy = (a._q._y - a._p._y) / a._length;
  const double first = (b._p._x - a._p._x) * ux + (b._p._y - a._p._y) * uy;
  const double second = (b._q._x - a._p._x) * ux + (b._q._y - a._p._y) * uy;
  return std::min(a._length, std::max(first, second)) -
             std::max(0.0, std::min(first, second)) > tolerance;
}

// Atom labels alone are ambiguous at spiro atoms. Assign each ring bond from
// the rings themselves, merging fused/bridged rings but not spiro-only joins.
std::vector<RingSystemResult> AssignRingSystems(Molecule& mol, std::vector<Segment>& bonds) {
  std::map<int, int> fused_system;
  std::vector<std::vector<int>> atoms;
  std::vector<int> rings;
  for (int i = 0; i < mol.nrings(); ++i) {
    const Ring& ring = *mol.ringi(i);
    int system;
    if (ring.is_fused()) {
      auto [iter, inserted] = fused_system.emplace(ring.fused_system_identifier(), atoms.size());
      system = iter->second;
      if (inserted) {
        atoms.emplace_back();
        rings.push_back(0);
      }
    } else {
      system = atoms.size();
      atoms.emplace_back();
      rings.push_back(0);
    }
    ++rings[system];
    for (int j = 0; j < ring.number_elements(); ++j) {
      atoms[system].push_back(ring[j]);
      const int next = ring[(j + 1) % ring.number_elements()];
      bonds[mol.which_bond(ring[j], next)]._system = system;
    }
  }
  std::vector<RingSystemResult> result;
  for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
    auto& members = atoms[i];
    std::sort(members.begin(), members.end());
    members.erase(std::unique(members.begin(), members.end()), members.end());
    result.emplace_back(members.size(), rings[i]);
  }
  return result;
}

}  // namespace

const char* StatusMessage(Status status) {
  switch (status) {
    case Status::kOk: return "ok";
    case Status::kEmpty: return "empty molecule";
    case Status::kDisconnected: return "disconnected molecule; use -l to retain the largest fragment";
    case Status::kNot2D: return "highest_coordinate_dimensionality() is not 2; use -z i to skip";
    case Status::kInvalidCoordinates: return "non-finite coordinates";
    case Status::kInvalidOptions: return "near-miss fraction must be finite and non-negative";
  }
  return "unknown error";
}

Status Analyse(Molecule& mol, const Options& options, Result& result) {
  result = Result();
  if (!std::isfinite(options.near_miss_fraction()) || options.near_miss_fraction() < 0) {
    return Status::kInvalidOptions;
  }
  const int n = mol.natoms();
  if (n == 0) {
    return Status::kEmpty;
  }
  if (mol.number_fragments() != 1) {
    return Status::kDisconnected;
  }
  std::vector<Point> points;
  points.reserve(n);
  for (const Atom* atom : mol) {
    // Crossings use only XY, but validate Z as well before checking input
    // dimensionality so malformed coordinates are reported as invalid input.
    if (!std::isfinite(atom->x()) || !std::isfinite(atom->y()) || !std::isfinite(atom->z())) {
      return Status::kInvalidCoordinates;
    }
    points.push_back({atom->x(), atom->y()});
  }
  if (mol.highest_coordinate_dimensionality() != 2) {
    return Status::kNot2D;
  }
  result._atoms = n;
  std::vector<Segment> bonds;
  std::vector<double> lengths;
  for (const Bond* bond : mol.bond_list()) {
    const int a1 = bond->a1(), a2 = bond->a2();
    const double length = std::sqrt(SquaredDistance(points[a1], points[a2]));
    bonds.push_back({a1, a2, points[a1], points[a2], length});
    if (length > 0) {
      lengths.push_back(length);
    }
  }
  double scale = 1;
  if (!lengths.empty()) {
    const auto middle = lengths.begin() + lengths.size() / 2;
    std::nth_element(lengths.begin(), middle, lengths.end());
    scale = *middle;
  }
  const double tolerance = scale * 1.0e-7;
  const double tolerance_squared = tolerance * tolerance;
  double min_squared = std::numeric_limits<double>::infinity();
  double max_squared = 0;
  double nonbonded_squared = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double distance = SquaredDistance(points[i], points[j]);
      min_squared = std::min(min_squared, distance);
      max_squared = std::max(max_squared, distance);
      if (options.diagnostics()) {
        result._coincident_atoms += distance <= tolerance_squared;
        if (!mol.are_bonded(i, j)) {
          nonbonded_squared = std::min(nonbonded_squared, distance);
        }
      }
    }
  }
  result._min_distance = n > 1 ? std::sqrt(min_squared) : 0;
  result._max_distance = std::sqrt(max_squared);
  result._min_nonbonded = std::isfinite(nonbonded_squared) ? std::sqrt(nonbonded_squared) : 0;
  result._topological_diameter = n > 1 ? mol.longest_path() : 0;
  std::vector<int> labels(n);
  mol.label_atoms_by_ring_system(labels.data());
  result._ring_atoms = std::count_if(labels.begin(), labels.end(), [](int label) { return label > 0; });
  result._systems = AssignRingSystems(mol, bonds);
  std::vector<int> ring_bonds, chain_bonds;
  for (int i = 0; i < static_cast<int>(bonds.size()); ++i) {
    (bonds[i]._system >= 0 ? ring_bonds : chain_bonds).push_back(i);
    if (options.diagnostics()) {
      result._degenerate_bonds += bonds[i]._length <= tolerance;
    }
  }
  const double near_distance = scale * options.near_miss_fraction();
  std::vector<bool> crossed(bonds.size(), false);
  auto examine = [&](int i, int j) {
    const Segment& a = bonds[i];
    const Segment& b = bonds[j];
    if (a._a1 == b._a1 || a._a1 == b._a2 || a._a2 == b._a1 || a._a2 == b._a2 ||
        a._length <= tolerance || b._length <= tolerance ||
        !BoxesClose(a, b, options.diagnostics() ? std::max(tolerance, near_distance) : tolerance)) {
      return;
    }
    if (ProperCrossing(a, b, tolerance)) {
      crossed[i] = crossed[j] = true;
      if (b._system < 0) {
        ++result._ring_chain;
        ++result._systems[a._system]._ring_chain;
      } else if (a._system == b._system) {
        ++result._ring_ring;
        ++result._systems[a._system]._ring_ring;
      } else {
        ++result._inter_system;
      }
    } else if (options.diagnostics()) {
      if (Overlap(a, b, tolerance)) {
        ++result._overlaps;
      } else {
        const double distance = std::min({PointSegmentDistanceSquared(a._p, b),
            PointSegmentDistanceSquared(a._q, b), PointSegmentDistanceSquared(b._p, a),
            PointSegmentDistanceSquared(b._q, a)});
        if (distance > tolerance_squared && distance < near_distance * near_distance) {
          ++result._near_misses;
        }
      }
    }
  };
  // These two loops count each eligible pair once and never compare chain-chain.
  for (int i = 0; i < static_cast<int>(ring_bonds.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(ring_bonds.size()); ++j) {
      examine(ring_bonds[i], ring_bonds[j]);
    }
    for (int chain : chain_bonds) {
      examine(ring_bonds[i], chain);
    }
  }
  result._crossed_bonds = std::count(crossed.begin(), crossed.end(), true);
  result._score = result._inter_system;
  for (const auto& system : result._systems) {
    result._systems_crossed += system.crossings() > 0;
    result._max_crossings = std::max(result._max_crossings, system.crossings());
    result._score += system.score();
    result._worst_density = std::max(result._worst_density, system.score() / system.atoms());
  }
  if (options.diagnostics()) {
    for (int index : ring_bonds) {
      const Segment& bond = bonds[index];
      if (bond._length <= tolerance) {
        continue;
      }
      for (int i = 0; i < n; ++i) {
        if (i == bond._a1 || i == bond._a2) {
          continue;
        }
        const double projection = ((points[i]._x - bond._p._x) * (bond._q._x - bond._p._x) +
            (points[i]._y - bond._p._y) * (bond._q._y - bond._p._y)) / bond._length;
        if (projection > tolerance && projection < bond._length - tolerance &&
            PointSegmentDistanceSquared(points[i], bond) <= tolerance_squared) {
          ++result._atoms_on_ring_bonds;
        }
      }
    }
  }
  return Status::kOk;
}

}  // namespace bond_crossings
