#ifndef MOLECULE_TOOLS_RING_SUBSTITUTION_H
#define MOLECULE_TOOLS_RING_SUBSTITUTION_H

#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Molecule_Lib/molecule.h"

namespace ring_substitution {

class RingSubstitutionGenerator {
  private:
    int _positional_information_only;
    int _simple_atom_types;
    int _full_atom_types;
    int _max_path_length;

    // We can lower the distances between molecules by including presence
    // or absence of the features
    int _place_single_feature_bits;

  // private functions.
    int Generate(Molecule& m, const int* atype, Sparse_Fingerprint_Creator& sfpc) const;
    int Generate(const IWString& mname,
                  const resizable_array<int>& abstract_path,
                  int* tmp, Sparse_Fingerprint_Creator& sfpc)const;
    int determine_substitution_type(Molecule& m,
                        atom_number_t zatom, const Atom& a) const;
    int assign_atom_types(Molecule& m, int* atype) const;
    int Generate(Molecule& m, const int* atype,
                  const Set_of_Atoms& par, int* tmp,
                  Sparse_Fingerprint_Creator& sfpc)const ;

  public:
    RingSubstitutionGenerator();

    void set_positional_information_only(int s) {
      _positional_information_only = s;
    }
    void set_simple_atom_types(int s) {
      _simple_atom_types = s;
    }
    void set_full_atom_types(int s) {
      _full_atom_types = s;
    }
    void set_max_path_length(int s) {
      _max_path_length = s;
    }
    void set_place_single_feature_bits(int s) {
      _place_single_feature_bits = s;
    }

    int Generate(Molecule& m,
                 Sparse_Fingerprint_Creator& sfc) const;
};

}  // namespace ring_substitution

#endif // MOLECULE_TOOLS_RING_SUBSTITUTION_H
