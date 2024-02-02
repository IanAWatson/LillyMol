#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/xlogp3.h"

namespace xlogp3 {

using std::cerr;

struct PerMoleculeData {
  public:
    Molecule& mol;
    std::unique_ptr<atomic_number_t[]> z;
    int* lipophilic;
    int* aromatic;
    int* conjugated;
    int* ncon;
    int* ring_bond_count;
    int* attached_heteroatom_count;
    int* single_bond_count;
    int* double_bond_count;
    int* triple_bond_count;
    int* connected_to_conjugated;
    int* atom_assigned;
    // Useful for atoms that have only one attachement
    atom_number_t* first_atom_attached;

  // Private functions
    int IdentifyLipophilic();
    int IsLiphphilic(int* visited);

    int DoublyBonded(atom_number_t zatom, atomic_number_t target) const;
    int IsAmide(atom_number_t zatom) const;
    atom_number_t Exocyclic(atom_number_t zatom) const;
    atom_number_t DoublyBonded(atom_number_t zatom) const;

    int IdentifyGroups(float& result);
    int IsCyano(atom_number_t zatom, float& result);
    int IsDiazo(atom_number_t zatom, float& result);
    int IsNitro(atom_number_t zatom, float& result);
    int IsNOxide(atom_number_t zatom, float& result);

    int AromaticCarbon(atom_number_t zatom, float& result);
    int SaturatedCarbonCh0(atom_number_t zatom, float& result);
    int SaturatedCarbonCh1(atom_number_t zatom, float& result);
    int SaturatedCarbonCh2(atom_number_t zatom, float& result);
    int SaturatedCarbonCh3(atom_number_t zatom, float& result);
    int UnSaturatedCarbonH2(atom_number_t zatom, float& result);
    int UnSaturatedCarbonH1(atom_number_t zatom, float& result);
    int UnSaturatedCarbonH0(atom_number_t zatom, float& result);
    int SaturatedCarbon(atom_number_t zatom, float& result);
    int UnSaturatedCarbon(atom_number_t zatom, float& result);
    int Carbon(atom_number_t zatom, float& result);
    int Oxygen(atom_number_t zatom, float& result);
    int SaturatedOxygen(atom_number_t zatom, float& result);
    int UnSaturatedOxygen(atom_number_t zatom, float& result);
    int AromaticNitrogen5(atom_number_t zatom, float& result);
    int AromaticNitrogen6(atom_number_t zatom, float& result);
    int AromaticNitrogen(atom_number_t zatom, float& result);
    int SaturatedNitrogen(atom_number_t zatom, float& result);
    int UnSaturatedNitrogen(atom_number_t zatom, float& result);
    int Nitrogen(atom_number_t zatom, float& result);
    int Fluorine(atom_number_t zatom, float& result);
    int Phosphorus(atom_number_t zatom, float& result);
    int Sulphur(atom_number_t zatom, float& result);
    int SaturatedSulphur(atom_number_t zatom, float& result);
    int UnSaturatedSulphur(atom_number_t zatom, float& result);
    int Chlorine(atom_number_t zatom, float& result);
    int Bromine(atom_number_t zatom, float& result);
    int Iodine(atom_number_t zatom, float& result);

  public:
    PerMoleculeData(Molecule& m);
    ~PerMoleculeData();
};

PerMoleculeData::PerMoleculeData(Molecule& m) : mol(m) {
  const int matoms = m.natoms();

  z = m.AtomicNumbers();

  m.compute_aromaticity_if_needed();

  lipophilic = new_int(matoms);
  aromatic = new_int(matoms);
  conjugated = new_int(matoms);
  ncon = new int[matoms];
  ring_bond_count = new int[matoms];
  attached_heteroatom_count = new_int(matoms);
  single_bond_count = new_int(matoms);
  double_bond_count = new_int(matoms);
  triple_bond_count = new_int(matoms);
  connected_to_conjugated = new_int(matoms);
  first_atom_attached = new atom_number_t[matoms];

  atom_assigned = new_int(matoms);

  for (int i = 0; i < matoms; ++i) {
    ring_bond_count[i] = m.ring_bond_count(i);
    ncon[i] = m.ncon(i);

    if (ncon[i] == 0) {
      first_atom_attached[i] = kInvalidAtomNumber;
    } else {
      first_atom_attached[i] = m.other(i, 0);
    }
  }

  for (const Bond* b : m.bond_list()) {
    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    if (z[a1] != 6) {
      ++attached_heteroatom_count[a2];
    }
    if (z[a2] != 6) {
      ++attached_heteroatom_count[a1];
    }

    if (b->is_aromatic()) {
      ++aromatic[a1];
      ++aromatic[a2];
      ++conjugated[a1];
      ++conjugated[a2];
    } else if (b->is_single_bond()) {
      ++single_bond_count[a1];
      ++single_bond_count[a2];
    } else if (b->is_double_bond()) {
      ++double_bond_count[a1];
      ++double_bond_count[a2];
      ++conjugated[a1];
      ++conjugated[a2];
    } else if (b->is_triple_bond()) {
      ++triple_bond_count[a1];
      ++triple_bond_count[a2];
      ++conjugated[a1];
      ++conjugated[a2];
    }
  }

  for (const Bond* b : m.bond_list()) {
    if (b->is_aromatic()) {
      continue;
    }
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();
    if (conjugated[a1]) {
      ++connected_to_conjugated[a2];
    }
    if (conjugated[a2]) {
      ++connected_to_conjugated[a1];
    }
  }

  IdentifyLipophilic();
}

PerMoleculeData::~PerMoleculeData() {
}

constexpr int kVisited = 1;
constexpr int kNextShell = 2;
constexpr int kMarkNextShell = 3;

// Return 0 if moving outward from any atom for which
// visited[i] == kNextShell
// we encounter a conjugated atom.
int
PerMoleculeData::IsLiphphilic(int* visited) {
  const int matoms = mol.natoms();

  for (int i = 0; i < matoms; ++i) {
    if (visited[i] != kNextShell) {
      continue;
    }

    for (const Bond* b : mol[i]) {
      atom_number_t o = b->other(i);
      if (conjugated[o]) {
        return 0;
      }
      visited[o] = kMarkNextShell;
    }
  }

  for (int i = 0; i < matoms; ++i) {
    if (visited[i] == 0) {
      continue;
    }
    if (visited[i] == kVisited) {
      continue;
    }

    if (visited[i] == kNextShell) {
      visited[i] = kVisited;
    } else if (visited[i] == kMarkNextShell) {
      visited[i] = kNextShell;
    }
  }

  return 0;
}

// The paper defines a lipophilic atom as something that
// is a carbon and has no heteroatom, or conjugated atom within
// 3 single bonds.
// Avoid computing the a the distance matrix, since we only need
// distance 3.
int
PerMoleculeData::IdentifyLipophilic() {
  const int matoms = mol.natoms();
  std::unique_ptr<int[]> visited = std::make_unique<int[]>(matoms);

  for (int i = 0; i < matoms; ++i) {
    if (z[i] != 6) {
      continue;
    }
    if (conjugated[i]) {
      continue;
    }

    std::fill_n(visited.get(), matoms, 0);
    visited[i] = kNextShell;

    // the value to be assigned to lipophilic[i];
    int to_set = 1;
    for (int j = 0; j < 3; ++j) {
      if (! IsLiphphilic(visited.get())) {
        to_set = 0;
        break;
      }
    }
    lipophilic[i] = to_set;
  }

  return 1;
}

int
PerMoleculeData::IsNitro(atom_number_t zatom, float& result) {
  if (ncon[zatom] != 3) {
    return 0;
  }
  if (double_bond_count[zatom] != 2) {
    return 0;
  }
  if (attached_heteroatom_count[zatom] < 2) {
    return 0;
  }

  atom_number_t o1 = kInvalidAtomNumber;
  atom_number_t o2 = kInvalidAtomNumber;
  for (const Bond* b : mol[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }

    atom_number_t o = b->other(zatom);
    if (z[o] != 8) {
      return 0;
    }
    if (o1 == kInvalidAtomNumber) {
      o1 = o;
    } else {
      o2 = o;
    }
  }

  if (o2 == kInvalidAtomNumber) {
    return 0;
  }

  atom_assigned[zatom] = kNitro;
  atom_assigned[o1] = kNitro;
  atom_assigned[o2] = kNitro;

  result += 1.2442;

  return 1;
}

int
PerMoleculeData::IsCyano(atom_number_t zatom, float& result) {
  if (ncon[zatom] != 1) {
    return 0;
  }
  if (triple_bond_count[zatom] != 1) {
    return 0;
  }

  atom_number_t carbon = first_atom_attached[zatom];
  if (z[carbon] != 6) {
    return 0;
  }

  atom_assigned[zatom] = kCyano;
  atom_assigned[carbon] = kCyano;

  result += 0.0337;

  return 1;
}

int
PerMoleculeData::IsDiazo(atom_number_t zatom, float& result) {
  if (ncon[zatom] != 2) {
    return 0;
  }
  if (ring_bond_count[zatom]) {
    return 0;
  }

  if (double_bond_count[zatom] != 2) {
    return 0;
  }

  atom_number_t n1 = kInvalidAtomNumber;
  atom_number_t n2 = kInvalidAtomNumber;

  for (const Bond* b : mol[zatom]) {
    atom_number_t n = b->other(zatom);
    if (z[n] != 7) {
      return 0;
    }

    if (n1 == kInvalidAtomNumber) {
      n1 = n;
    } else {
      n2 = n;
    }
  }

  if (n2 == kInvalidAtomNumber) {
    return 0;
  }

  atom_assigned[zatom] = kDiazo;
  atom_assigned[n1] = kDiazo;
  atom_assigned[n2] = kDiazo;

  result += 0.5339;

  return 1;
}

int
PerMoleculeData::IsNOxide(atom_number_t zatom, float& result) {
  if (conjugated[zatom]) {
    return 0;
  }

  for (const Bond* b : mol[zatom]) {
    atom_number_t o = b->other(zatom);

    if (z[o] != 8) {
      continue;
    }
    if (ncon[o] != 1) {
      continue;
    }

    atom_assigned[zatom] = kNOxide;
    atom_assigned[o] = kNOxide;

    result += -1.2147;

    return 1;
  }

  return 0;
}

// Place into the atom_assigned array any atoms that are covered by
// the groups defined in the paper.
int
PerMoleculeData::IdentifyGroups(float& result) {
  const int matoms = mol.natoms();
  int rc = 0;
  for (int i = 0; i < matoms; ++i) {
    // All functional groups start with nitrogen
    if (z[i] != 7) {
      continue;
    }
    if (IsNitro(i, result)) {
      ++rc;
    } else if (IsCyano(i, result)) {
      ++rc;
    } else if (IsDiazo(i, result)) {
      ++rc;
    } else if (IsNOxide(i, result)) {
      ++rc;
    }
  }

  return rc;
}

// `zatom` is aromatic
atom_number_t
PerMoleculeData::Exocyclic(atom_number_t zatom) const {
  assert(aromatic[zatom] == 2);

  for (const Bond* b : mol[zatom]) {
    if (b->is_aromatic()) {
      continue;
    }
    return b->other(zatom);
  }

  return kInvalidAtomNumber;
}

int
PerMoleculeData::AromaticCarbon(atom_number_t zatom, float& result) {
  const int hcount = mol.hcount(zatom);
  if (hcount && attached_heteroatom_count[zatom] == 0) {
    atom_assigned[zatom] = kCarh;  // 20
    result += 0.3157;
    return 1;
  }

  if (hcount) {
    atom_assigned[zatom] = kCarhx;  // 19
    result += -0.1039;
    return 1;

  }

  if (ncon[zatom] == 3 && aromatic[zatom] == 3) {
    atom_assigned[zatom] = kCarar;  // 21
    result += 0.1038;
    return 1;
  }

  if (ncon[zatom] == 3 && aromatic[zatom] == 2 && attached_heteroatom_count[zatom] == 0) {
    atom_assigned[zatom] = kCar;  // 25
    result += 0.1911;
    return 1;
  }

  atom_number_t o = Exocyclic(zatom);
  if (o == kInvalidAtomNumber) {
    cerr << "Huh, no non-aromatic connection " << mol.smarts_equivalent_for_atom(zatom) << '\n';
    return 0;
  }

  // c=O and c=N are classified as sp2 in the xlogp world.
  if (ncon[zatom] == 3 && ring_bond_count[zatom] == 2 &&
      aromatic[zatom] == 2 && double_bond_count[zatom] == 1 &&
      ncon[o] == 1 && z[o] != 6) {
    atom_assigned[zatom] = kC2dx;  // 36
    result += -0.6093;
    return 1;
  }

  if (ncon[zatom] == 3 && ring_bond_count[zatom] == 2) {
    if (z[o] != 6 && attached_heteroatom_count[zatom] == 1) {
      atom_assigned[zatom] = kCary;  // 23
      result += -0.0112;
      return 1;
    }
  }

  if (z[o] == 6) {
    atom_assigned[zatom] = kCarx;  // 24
    result += -0.1874;
    return 1;
  }

  atom_assigned[zatom] = kCarxy;  // 22
  result += -0.1003;
  return 1;
}

int
PerMoleculeData::UnSaturatedCarbonH2(atom_number_t zatom, float& result) {
  atom_assigned[zatom] = kC2hc; // 29
  result += 0.3214;
  return 1;
}

atom_number_t
PerMoleculeData::DoublyBonded(atom_number_t zatom) const {
  for (const Bond* b : mol[zatom]) {
    if (b->is_double_bond()) {
      return b->other(zatom);
    }
  }

  return kInvalidAtomNumber;
}

int
PerMoleculeData::UnSaturatedCarbonH1(atom_number_t zatom, float& result) {
  atom_number_t o = DoublyBonded(zatom);
  assert(o != kInvalidAtomNumber);

  if (z[o] == 6 && attached_heteroatom_count[zatom] == 1) {
    atom_assigned[zatom] = kC2hcx;  // 27
    result += -0.0967;
    return 1;
  }

  if (z[o] == 6 && attached_heteroatom_count[zatom] == 0 &&
      ring_bond_count[zatom]) {
    atom_assigned[zatom] = kC2hcring;  // 28
    result += 0.4004;
    return 1;
  }

  if (z[o] == 6 && attached_heteroatom_count[zatom] == 0) {
    atom_assigned[zatom] = kC2hc;  // 29
    result += 0.3214;
    return 1;
  }

  if (z[o] != 6) {
    atom_assigned[zatom] = kC2hx;  // 30
    result += -0.8756;
    return 1;
  }

  cerr << "Unrecognised UnSaturatedCarbonH1 " << mol.smarts_equivalent_for_atom(zatom) << '\n';
  return 0;
}

int
PerMoleculeData::UnSaturatedCarbonH0(atom_number_t zatom, float& result) {
  atom_number_t o = DoublyBonded(zatom);
  assert(o != kInvalidAtomNumber);

  if (z[o] == 6 && attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC2cx;  // 31
    result += -0.2069;
    return 1;
  }
  if (z[o] == 6 && attached_heteroatom_count[zatom] == 0 &&
      ring_bond_count[zatom]) {
    atom_assigned[zatom] = kC2cring;  // 32
    result += -0.2084;
    return 1;
  }
  if (z[o] == 6) {
    atom_assigned[zatom] = kC2c;  // 33
    result += 0.4840;
    return 1;
  }

  // attached atom is heteroatom
  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC2xx;  // 34
    result += -0.8076;
    return 1;
  }
  if (ring_bond_count[zatom]) {
    atom_assigned[zatom] = kC2xring;  // 35
    result += -0.5034;
    return 1;
  }

  atom_assigned[zatom] = kC2dx;  // 36
  result += -0.6093;
  return 1;
}

int
PerMoleculeData::UnSaturatedCarbon(atom_number_t zatom, float& result) {
  if (ncon[zatom] == 2 && double_bond_count[zatom] == 2) {
    atom_assigned[zatom] = kC1d;  // 37
    result += -0.5879;
    return 1;
  }

  if (ncon[zatom] == 2 && triple_bond_count[zatom] == 1) {
    atom_assigned[zatom] = kC1;  // 38
    result += 0.1945;
    return 1;
  }

  const int hcount = mol.hcount(zatom);
  if (hcount == 2) {
    return UnSaturatedCarbonH2(zatom, result);
  }

  if (hcount == 1) {
    return UnSaturatedCarbonH1(zatom, result);
  }

  return UnSaturatedCarbonH0(zatom, result);
}

int
PerMoleculeData::SaturatedCarbonCh3(atom_number_t zatom, float& result) {
  assert(ncon[zatom] == 1);
  if (connected_to_conjugated[zatom] && attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kCh3XPi;  // 2
    result += -0.0753;
    return 1;
  }

  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kCh3X;  // 3
    result += 0.0402;
    return 1;
  }

  if (connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kCh3Pi;  // 4
    result += 0.5018;
    return 1;
  }

  atom_assigned[zatom] = kC33h;  // 5
  result += 0.5240;
  return 1;
}

int
PerMoleculeData::SaturatedCarbonCh2(atom_number_t zatom, float& result) {
  if (connected_to_conjugated[zatom] && attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC32hxpi;  // 7
    result += -0.2441;
    return 1;
  }
  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC32hx;  // 8
    result += -0.0821;
    return 1;
  }
  if (connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kC32hpi;  // 9
    result += 0.2718;
    return 1;
  }

  atom_assigned[zatom] = kC32h;  // 10
  result += 0.3436;
  return 1;
}

int
PerMoleculeData::SaturatedCarbonCh1(atom_number_t zatom, float& result) {
  if (connected_to_conjugated[zatom] && attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC3hxpi;  // 11
    result += -0.3711;
    return 1;
  }
  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC3hx;  // 12
    result += -0.1426;
    return 1;
  }
  if (connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kC3hpi;  // 13
    result += 0.0841;
    return 1;
  }

  atom_assigned[zatom] = kC3h;  // 14
  result += 0.1485;
  return 1;
}

int
PerMoleculeData::SaturatedCarbonCh0(atom_number_t zatom, float& result) {
  if (connected_to_conjugated[zatom] && attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC3xpi;  // 15
    result += -0.5475;
    return 1;
  }
  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kC3x;  // 16
    result += -0.4447;
    return 1;
  }
  if (connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kC3pi;  // 17
    result += 0.0885;
    return 1;
  }

  atom_assigned[zatom] = kC3;  // 18
  result += 0.0596;
  return 1;
}

int
PerMoleculeData::SaturatedCarbon(atom_number_t zatom, float& result) {
  if (lipophilic[zatom] && ncon[zatom] == 1) {
    atom_assigned[zatom] = kLipophilicCH3; // 1
    result += 0.7896;
    return 1;
  }
  if (lipophilic[zatom] && ncon[zatom] == 2) {
    atom_assigned[zatom] = kLipophilicCH3; // 1
    result += 0.7896;
    return 1;
  }

  if (ncon[zatom] == 1) {
    return SaturatedCarbonCh3(zatom, result);
  }
  if (ncon[zatom] == 2) {
    return SaturatedCarbonCh2(zatom, result);
  }
  if (ncon[zatom] == 3) {
    return SaturatedCarbonCh1(zatom, result);
  }
  if (ncon[zatom] == 4) {
    return SaturatedCarbonCh0(zatom, result);
  }

  return 0;
}

int
PerMoleculeData::Carbon(atom_number_t zatom, float& result) {
  if (aromatic[zatom]) {
    return AromaticCarbon(zatom, result);
  }
  if (conjugated[zatom]) {
    return UnSaturatedCarbon(zatom, result);
  }

  return SaturatedCarbon(zatom, result);
}

int
PerMoleculeData::AromaticNitrogen6(atom_number_t zatom, float& result) {
  if (attached_heteroatom_count[zatom] == 2) {
    atom_assigned[zatom] = kNarx2;  // 48
    result += -0.2167;
    return 1;
  }

  if (attached_heteroatom_count[zatom] == 1) {
    atom_assigned[zatom] = kNarx;  // 49
    result += -0.2974;
    return 1;
  }

  atom_assigned[zatom] = kNar;
  result += 0.0888;
  return 1;
}

int
PerMoleculeData::AromaticNitrogen5(atom_number_t zatom, float& result) {
  int hcount = mol.hcount(zatom);

  if (hcount && attached_heteroatom_count[zatom] == 1) {
    atom_assigned[zatom] = kNarhx;  // 51
    result += 0.3675;
    return 1;
  }

  if (hcount && attached_heteroatom_count[zatom] == 0) {
    atom_assigned[zatom] = kNarh;  // 52
    result += 0.2364;
    return 1;
  }

  if (attached_heteroatom_count[zatom] == 2) {
    atom_assigned[zatom] = kNarx2_2;  // 53
    result += 1.1022;
    return 1;
  }

  if (attached_heteroatom_count[zatom] == 1) {
    atom_assigned[zatom] = kNarx_2;  // 54
    result += 0.4854;
    return 1;
  }

  if (attached_heteroatom_count[zatom] == 0) {
    atom_assigned[zatom] = kNar_2;  // 55
    result += 0.3181;
    return 1;
  }

  cerr << "Unrecognised aromatic 5 " << mol.smarts_equivalent_for_atom(zatom) << '\n';
  return 0;
}

int
PerMoleculeData::AromaticNitrogen(atom_number_t zatom, float& result) {
  const Ring* r = mol.ring_containing_atom(zatom);
  assert (r != nullptr);

  if (r->size() == 5) {
    return AromaticNitrogen5(zatom, result);
  }
  if (r->size() == 6) {
    return AromaticNitrogen6(zatom, result);
  }

  return 0;
}

int
PerMoleculeData::IsAmide(atom_number_t zatom) const {
  for (const Bond* b : mol[zatom]) {
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t cs = b->other(zatom);
    if (z[cs] == 6) {
    } else if (z[cs] == 16) {
    } else {
      continue;
    }
    if (ncon[cs] < 2) {
      continue;
    }

    if (! conjugated[cs]) {
      continue;
    }

    for (const Bond* b2 : mol[cs]) {
      if (! b2->is_double_bond()) {
        continue;
      }
      atom_number_t o = b2->other(cs);

      if (z[o] == 8) {
        return 1;
      }
      if (z[o] == 16) {
        return 1;
      }
    }
  }

  return 0;
}

atom_number_t
DoublyBondedTo(const Molecule& m, atom_number_t zatom) {
  for (const Bond* b : m[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }
    return b->other(zatom);
  }

  cerr << "Huh, no doubly bonded atom\n";
  return kInvalidAtomNumber;
}

int
PerMoleculeData::UnSaturatedNitrogen(atom_number_t zatom, float& result) {
  if (mol.hcount(zatom)) {
    atom_assigned[zatom] = kN2h; // 56
    result += 0.6927;
    return 1;
  }
  
  atom_number_t d = DoublyBondedTo(mol, zatom);

  if (ring_bond_count[zatom] && z[d] == 6) {
    atom_assigned[zatom] = kN2cring; // 57
    result += 0.7974;
    return 1;
  }

  if (z[d] == 6) {
    atom_assigned[zatom] = kN2c; // 58
    result += 0.9794;
    return 1;
  }

  atom_assigned[zatom] = kN2c; // 59
  result += 0.2698;
  return 1;
}

int
PerMoleculeData::SaturatedNitrogen(atom_number_t zatom, float& result) {
  if (IsAmide(zatom)) {
    if (ncon[zatom] == 1) {
      atom_assigned[zatom] = kNam2h;  // 39
      result += -0.6414;
    } else if (ncon[zatom] == 2) {
      atom_assigned[zatom] = kNamh;  // 42
      result += -0.3333;
    } else if (ncon[zatom] == 3) {
      atom_assigned[zatom] = kNam;  // 45
      result += -0.1551;
    }
    return 1;
  }

  const int hcount = mol.hcount(zatom);

  if (hcount == 2 && ncon[zatom] == 1 && connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kN32hpi;  // 40
    result += -0.3637;
    return 1;
  }

  if (hcount == 2 && ncon[zatom] == 1) {
    atom_assigned[zatom] = kN32h;   // 41
    result += -0.7445;
    return 1;
  }

  if (hcount == 1 && ncon[zatom] == 2 && connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kN3hpi;   // 43
    result += 0.2172;
    return 1;
  }

  if (hcount == 1 && ncon[zatom] == 2) {
    atom_assigned[zatom] = kN3h;   // 44
    result += -0.2610;
    return 1;
  }

  if (ncon[zatom] == 3 && connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kN3pi;   // 46
    result += 0.3776;
    return 1;
  }

  if (ncon[zatom] == 3) {
    atom_assigned[zatom] = kN3;   // 47
    result += 0.1779;
    return 1;
  }

  cerr << "Unclassified saturated Nitrogen atom " << mol.smarts_equivalent_for_atom(zatom) << '\n';
  return 0;
}

int
PerMoleculeData::Nitrogen(atom_number_t zatom, float& result) {
  if (aromatic[zatom]) {
    return AromaticNitrogen(zatom, result);
  }

  if (conjugated[zatom]) {
    return UnSaturatedNitrogen(zatom, result);
  }

    return SaturatedNitrogen(zatom, result);
}

int
PerMoleculeData::UnSaturatedOxygen(atom_number_t zatom, float& result) {
  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kO2c;  // 65
    result += 0.7148;
    return 1;
  }

  atom_assigned[zatom] = kO2x;  // 66
  result += -0.5411;
  return 1;
}

int
PerMoleculeData::SaturatedOxygen(atom_number_t zatom, float& result) {
  if (ncon[zatom] == 1) {
    if (connected_to_conjugated[zatom]) {
      atom_assigned[zatom] = kO3hpi;
      result += -0.0381;
    } else {
      atom_assigned[zatom] = kO3h;
      result += -0.4802;
    }

    return 1;
  }

  // Ethers

  if (connected_to_conjugated[zatom]) {
    atom_assigned[zatom] = kO3pi;
    result += 0.2701;
    return 1;
  }

  atom_assigned[zatom] = kO3;
  result += 0.0059;
  return 1;
}

int
PerMoleculeData::Oxygen(atom_number_t zatom, float& result) {
  if (aromatic[zatom]) {
    atom_assigned[zatom] = kOar;
    result += 0.5238;
    return 1;
  }

  if (conjugated[zatom]) {
    return UnSaturatedOxygen(zatom, result);
  }

  return SaturatedOxygen(zatom, result);
}

int
PerMoleculeData::SaturatedSulphur(atom_number_t zatom, float& result) {
  if (ncon[zatom] == 1 && single_bond_count[zatom] == 1) {
    atom_assigned[zatom] = kS3h;  // 67
    result += 0.4927;
    return 1;
  }

  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kS3x;  // 69
    result += 0.4125;
    return 1;
  }

  atom_assigned[zatom] = kS3;  // 70
  result += 0.8300;
  return 1;
}

// Return the number of doubly bonded `target` atoms attached to `zatom`.
int
PerMoleculeData::DoublyBonded(atom_number_t zatom, atomic_number_t target) const {
  int rc = 0;
  for (const Bond* b : mol[zatom]) {
    if (! b->is_double_bond()) {
      continue;
    }
    atom_number_t o = b->other(zatom);
    if (z[o] == target) {
      ++rc;
    }
  }

  return rc;
}

int
PerMoleculeData::UnSaturatedSulphur(atom_number_t zatom, float& result) {
  int doubly_bonded_oxygens = DoublyBonded(zatom, 8);

  if (doubly_bonded_oxygens >= 2) {
    atom_assigned[zatom] = kSo2;
    result += 0.5729;
    return 1;
  }

  if (doubly_bonded_oxygens == 1) {
    atom_assigned[zatom] = kSo;
    result += 0.0525;
    return 1;
  }

  if (attached_heteroatom_count[zatom]) {
    atom_assigned[zatom] = kS2x;
    result += 1.2218;
    return 1;
  }

  atom_assigned[zatom] = kS2c;
  result += 1.3544;
  return 1;
}

int
PerMoleculeData::Sulphur(atom_number_t zatom, float& result) {
  if (aromatic[zatom]) {
    atom_assigned[zatom] = kSar;  // 68
    result += 1.1715;
    return 1;
  }

  if (conjugated[zatom]) {
    return UnSaturatedSulphur(zatom, result);
  }

  return SaturatedSulphur(zatom, result);
}

int
PerMoleculeData::Phosphorus(atom_number_t zatom, float& result) {
  atom_assigned[zatom] = kP3;
  result += -0.6694;
  return 1;
}

int
PerMoleculeData::Fluorine(atom_number_t zatom, float& result) {
  atom_number_t o = first_atom_attached[zatom];
  if (conjugated[o]) {
    atom_assigned[zatom] = kFpi;
    result += 0.4401;
  } else {
    atom_assigned[zatom] = kF;
    result += 0.5360;
  }

  return 1;
}

int
PerMoleculeData::Chlorine(atom_number_t zatom, float& result) {
  atom_number_t o = first_atom_attached[zatom];
  if (conjugated[o]) {
    atom_assigned[zatom] = kClpi;
    result += 0.9610;
  } else {
    atom_assigned[zatom] = kCl;
    result += 0.8036;
  }

  return 1;
}

int
PerMoleculeData::Bromine(atom_number_t zatom, float& result) {
  atom_number_t o = first_atom_attached[zatom];
  if (conjugated[o]) {
    atom_assigned[zatom] = kBrpi;
    result += 1.0295;
  } else {
    atom_assigned[zatom] = kBr;
    result += 0.9664;
  }

  return 1;
}

int
PerMoleculeData::Iodine(atom_number_t zatom, float& result) {
  atom_number_t o = first_atom_attached[zatom];
  if (conjugated[o]) {
    atom_assigned[zatom] = kIpi;
    result += 0.7801;
  } else {
    atom_assigned[zatom] = kI;
    result += 0.9071;
  }

  return 1;
}

int
XLogPAtom(PerMoleculeData& pmd, atom_number_t zatom, float& result) {
  atomic_number_t z = pmd.z[zatom];
  if (z == 6) {
    return pmd.Carbon(zatom, result);
  }
  if (z == 7) {
    return pmd.Nitrogen(zatom, result);
  }
  if (z == 8) {
    return pmd.Oxygen(zatom, result);
  }
  if (z == 9) {
    return pmd.Fluorine(zatom, result);
  }
  if (z == 15) {
    return pmd.Phosphorus(zatom, result);
  }
  if (z == 16) {
    return pmd.Sulphur(zatom, result);
  }
  if (z == 17) {
    return pmd.Chlorine(zatom, result);
  }
  if (z == 35) {
    return pmd.Bromine(zatom, result);
  }
  if (z == 53) {
    return pmd.Iodine(zatom, result);
  }

  return 0;
}

std::optional<double>
XLogP(Molecule& m, int* status) {
  const int matoms = m.natoms();
  if (matoms == 0) {
    return std::nullopt;
  }

  std::fill_n(status, matoms, 0);

  PerMoleculeData pmd(m);

  float result = 0.0f;
  pmd.IdentifyGroups(result);

  for (int i = 0; i < matoms; ++i) {
    if (pmd.atom_assigned[i]) {
      continue;
    }

    if (! XLogPAtom(pmd, i, result)) {
      return std::nullopt;
    }
  }

  return result;
}

std::optional<double>
XLogP(Molecule& m) {
  std::unique_ptr<int[]> status = std::make_unique<int[]>(m.natoms());

  return XLogP(m, status.get());
}

}  // namespace xlogp3
