#include <algorithm>

#include "ring_substitution.h"

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/path.h"
#include "Molecule_Lib/path_around_ring.h"

namespace ring_substitution {

RingSubstitutionGenerator::RingSubstitutionGenerator() {
  _positional_information_only = 1;
  _simple_atom_types = 0;
  _max_path_length = 24;
  _place_single_feature_bits = 0;
}

/*
  Make sure the ring atom types have gaps in them because we
  incrememt by 1 if the ring atom is aromatic
*/

#define RSTYPE_DOUBLE_BOND_OUTSIDE_RING 1
#define RSTYPE_SPIRO 2
#define RSTYPE_RING_JOIN 3
#define RSTYPE_TWO_SUBSTITUENTS 4
#define RSTYPE_SUBSTITUTED 5
#define RSTYPE_ANOTHER_RING 7

#define RSTYPE_SATURATED_CARBON 9
#define RSTYPE_UNSATURATED_CARBON 11
#define RSTYPE_SATURATED_HETEROATOM 13
#define RSTYPE_UNSATURATED_HETEROATOM 15

#define RSTYPE_METHYL 17
#define RSTYPE_TERMINAL_N 19
#define RSTYPE_NITRO 21
#define RSTYPE_SATURATED_NITROGEN 23
#define RSTYPE_SP2_NITROGEN 25
#define RSTYPE_HYDROXY 27
#define RSTYPE_ETHER 29
#define RSTYPE_SULPH 31
#define RSTYPE_FLUORINE 33
#define RSTYPE_HALOGEN 35

int
RingSubstitutionGenerator::Generate(const IWString& mname,
                  const resizable_array<int>& abstract_path,
                  int* tmp, Sparse_Fingerprint_Creator& sfpc) const {
  int n = abstract_path.number_elements();

// #define DEBUG_RING_SUBSTITUTION
#ifdef DEBUG_RING_SUBSTITUTION
  cerr << mname << '\n';
  for (int i = 0; i < n; i++) {
    cerr << " i = " << i << " abstract_path " << abstract_path[i] << '\n';
  }
#endif

  assert(n > 2);

  copy_vector(tmp, abstract_path.rawdata(), n);
  copy_vector(tmp + n, abstract_path.rawdata(), n);

  int jstop = n;
  if (jstop > _max_path_length) {
    jstop = _max_path_length;
  }

  for (int i = 0; i < n; i++) {
    int ta = tmp[i];

    if (0 == ta) {
      continue;
    }

    for (int j = 1; j < jstop; j++) {
      int k = j + i;  // index of next atom

      int tb = tmp[k];
      if (0 == tb) {
        continue;
      }

      unsigned int b = 40 + 30 + j * 30 * 30 * _max_path_length +
                       (ta * tb);  // 30 is larger than all atom types assigned

#ifdef DEBUG_RING_SUBSTITUTION
      cerr << " i = " << i << " ar " << ta << " j = " << j << " k = " << k << " ar " << tb
           << " b = " << b << '\n';
#endif
      sfpc.hit_bit(b);
    }
  }

  return 1;
}

int
RingSubstitutionGenerator::Generate(Molecule& m, const int* atype,
                  const Set_of_Atoms& par, int* tmp,
                  Sparse_Fingerprint_Creator& sfpc) const {
  resizable_array<int> abstract_path;

  const int n = par.number_elements();

  int number_features = 0;
  int the_feature = -1;

  for (int i = 0; i < n; i++) {
    atom_number_t j = par[i];

    abstract_path.add(atype[j]);

    if (0 == atype[j]) {
      continue;
    }

    if (0 == number_features) {
      the_feature = atype[j];
    }

    number_features++;

    if (_place_single_feature_bits) {
      sfpc.hit_bit(atype[j]);
    }
  }

  // cerr << "Abstract ring contains " << number_features << " features\n";

  if (1 == number_features) {
    sfpc.hit_bit(1383 * n + the_feature);
    return 1;
  }

  assert(abstract_path.number_elements() == n);

  return Generate(m.name(), abstract_path, tmp, sfpc);
}

int
RingSubstitutionGenerator::Generate(Molecule& m, const int* atype,
                           Sparse_Fingerprint_Creator& sfpc) const {
  const int nr = m.nrings();

  const int matoms = m.natoms();

  int* ring_already_done = new_int(nr);
  std::unique_ptr<int[]> free_ring_already_done(ring_already_done);

  int* tmp = new int[matoms + matoms];
  std::unique_ptr<int[]> free_tmp(tmp);

  m.compute_aromaticity_if_needed();

  // Do all the non-fused rings first. Includes spiro fused

  int rings_processed = 0;

  for (int i = 0; i < nr; i++) {
    if (ring_already_done[i]) {
      continue;
    }

    const Ring* ri = m.ringi(i);

    //  cerr << "Ring " << i << " has " <<
    //  ri->largest_number_of_bonds_shared_with_another_ring() << " shared bonds\n";

    if (ri->largest_number_of_bonds_shared_with_another_ring() > 0) {
      continue;
    }

    Generate(m, atype, *ri, tmp, sfpc);

    ring_already_done[i] = 1;
    rings_processed++;
  }

  if (nr == rings_processed) {
    return 1;
  }

  // Now any fused rings

  int* in_ring_system = new_int(matoms);
  std::unique_ptr<int[]> free_in_ring_system(in_ring_system);

  for (int i = 0; i < nr; i++) {
    if (ring_already_done[i]) {
      continue;
    }

    std::fill_n(in_ring_system, matoms, 0);

    const Ring* ri = m.ringi(i);

    ri->set_vector(in_ring_system, 1);

    rings_processed++;

    for (int j = i + 1; j < nr; j++) {
      if (ring_already_done[j]) {
        continue;
      }

      const Ring* rj = m.ringi(j);

      if (rj->fused_system_identifier() != ri->fused_system_identifier()) {
        continue;
      }

      rj->set_vector(in_ring_system, 1);
      ring_already_done[j] = 1;

      rings_processed++;
    }

    Set_of_Atoms s;
    if (!path_around_edge_of_ring_system(m, in_ring_system, 1, s))  // strongly fused
    {
      sfpc.hit_bit(count_occurrences_of_item_in_array(
          1, matoms, in_ring_system));  // unprocessed ring, hit bit according to size
      continue;
    }

    Generate(m, atype, s, tmp, sfpc);

    if (nr == rings_processed) {
      break;
    }
  }

  return 1;
}

int
RingSubstitutionGenerator::determine_substitution_type(Molecule& m,
                        atom_number_t zatom, const Atom& a) const {
  int acon = a.ncon();
  // cerr << "determine_substitution_type for atom type " <<
  // m.smarts_equivalent_for_atom(zatom) << '\n';

  if (2 == acon) {  // unsubstituted
    return 0;
  }

  if (m.nrings(zatom) > 1) {
    return RSTYPE_RING_JOIN;
  }

  if (4 == acon) {
    return RSTYPE_TWO_SUBSTITUENTS;
  }

  for (int i = 0; i < acon; i++) {
    const Bond* b = a[i];

    if (b->nrings()) {
      continue;
    }

    if (b->is_double_bond()) {
      return RSTYPE_DOUBLE_BOND_OUTSIDE_RING;
    }

    if (_positional_information_only) {
      return RSTYPE_SUBSTITUTED;
    }

    atom_number_t o = b->other(zatom);

    if (m.is_ring_atom(o)) {
      return RSTYPE_ANOTHER_RING;
    }

    atomic_number_t zo = m.atomic_number(o);

    int ocon = m.ncon(o);

    if (_simple_atom_types) {
      int unsaturation = m.nbonds(o) - ocon;
      if (6 == zo) {
        if (0 == unsaturation) {
          return RSTYPE_SATURATED_CARBON;
        } else {
          return RSTYPE_UNSATURATED_CARBON;
        }
      } else if (1 == ocon && (9 == zo || 17 == zo || 35 == zo || 53 == zo)) {
        return RSTYPE_HALOGEN;
      } else if (0 == unsaturation) {
        return RSTYPE_SATURATED_HETEROATOM;
      } else {
        return RSTYPE_UNSATURATED_HETEROATOM;
      }

      continue;
    }

    //  Full atom typing

    if (6 == zo) {
      if (1 == ocon) {
        return RSTYPE_METHYL;
      }

      int nbonds = m.nbonds(o);
      if (ocon == nbonds) {
        return RSTYPE_SATURATED_CARBON;
      } else {
        return RSTYPE_UNSATURATED_CARBON;
      }
    } else if (7 == zo) {
      if (1 == ocon) {
        return RSTYPE_TERMINAL_N;
      }

      int nbonds = m.nbonds(o);

      if (3 == ocon && 5 == nbonds) {
        return RSTYPE_NITRO;
      } else if (ocon == nbonds) {
        return RSTYPE_SATURATED_NITROGEN;
      } else {
        return RSTYPE_SP2_NITROGEN;
      }
    } else if (8 == zo) {
      if (1 == ocon) {
        return RSTYPE_HYDROXY;
      } else {
        return RSTYPE_ETHER;
      }
    } else if (16 == zo) {
      if (1 == ocon) {
        return RSTYPE_HYDROXY;
      } else if (2 == ocon) {
        return RSTYPE_ETHER;
      } else {
        return RSTYPE_SULPH;  // some other kind of state, who knows
      }
    } else if (9 == zo) {
      return RSTYPE_FLUORINE;
    } else if (17 == zo || 35 == zo || 53 == zo) {
      return RSTYPE_HALOGEN;
    }

    return RSTYPE_SUBSTITUTED;  // of not classified above
  }

  assert(nullptr == "Could not classify non-ring bond???");

  return 0;
}

int
is_spiro_fused(Molecule& m, atom_number_t zatom, const Atom& a) {
  int acon = a.ncon();

  if (4 != acon) {
    return 0;
  }

  if (2 != m.nrings(zatom)) {  // too hard and rare otherwise
    return 0;
  }

  for (int i = 0; i < acon; i++) {
    const Bond* b = a[i];

    if (0 == b->nrings()) {
      return 0;
    }
  }

  return 1;
}

int
has_double_bond_outside_ring(Molecule& m, atom_number_t zatom, const Atom& a) {
  int acon = a.ncon();

  if (3 != acon) {
    return 0;
  }

  for (int i = 0; i < acon; i++) {
    const Bond* b = a[i];
    //  cerr << " from atom " << zatom << " have bond to " << b->other(zatom) << " ring "
    //  << b->nrings() << " single? " << b->is_single_bond() << '\n';

    if (b->nrings()) {
      continue;
    }

    if (b->is_double_bond()) {
      return 1;
    }
  }

  return 0;
}

int
RingSubstitutionGenerator::assign_atom_types(Molecule& m, int* atype) const {
  m.ring_membership();

  int matoms = m.natoms();

  for (int i = 0; i < matoms; i++) {
    int nr = m.nrings(i);

    if (0 == nr) {
      continue;
    }

    const Atom* a = m.atomi(i);

    if (has_double_bond_outside_ring(m, i, *a)) {
      atype[i] = RSTYPE_DOUBLE_BOND_OUTSIDE_RING;
      if (m.is_aromatic(i)) {
        atype[i]++;
      }
      continue;
    }

    if (is_spiro_fused(m, i, *a)) {
      atype[i] = RSTYPE_SPIRO;
      continue;
    }

    atype[i] = determine_substitution_type(m, i, *a);

    if (0 == atype[i])
      ;
    else if (m.is_aromatic(i)) {
      atype[i]++;
    }
  }

  return 1;
}


int
RingSubstitutionGenerator::Generate(Molecule& m,
                 Sparse_Fingerprint_Creator& sfc) const {

  int* atype = new_int(m.natoms());
  std::unique_ptr<int[]> free_atype(atype);

  assign_atom_types(m, atype);
  
  return Generate(m, atype, sfc);
}

}  // namespace ring_substitution
