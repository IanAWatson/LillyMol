#ifndef MOLECULE_TOOLS_SAFE_GENERATE_LIB_H_
#define MOLECULE_TOOLS_SAFE_GENERATE_LIB_H_

#include <algorithm>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/path.h"

#include "Molecule_Tools/mformula.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

namespace safe_generate {

// Zero preserves the legacy convention that every nonzero isotope is a marker.
inline bool
IsAttachmentIsotope(isotope_t isotope, isotope_t attachment_isotope) {
  return isotope != 0 && (attachment_isotope == 0 || isotope == attachment_isotope);
}

// Remove attachment markers, retaining all other isotopes in explicit-marker mode.
void RemoveAttachmentIsotopes(Molecule& m, isotope_t attachment_isotope);

class SafeFragment {
  private:
    // the smiles of the fragment
    // Individual dot separated tokens from something like
    // [1Cl]%10.[1C]%11.[1C]%121=NO[1C]%13%11O1.[1C]%131CCCCN1.C1CC[1C]%10[1C]%12C1
    IWString _smiles;

    // The indices within _smiles where the first digit after the % signs are.
    resizable_array<int> _first_digit;

    // A molecule build from _smiles by dropping the %nn characters
    Molecule _m;

    mformula::MFormula _mformula;

    // Number of atoms in the fragment.
    int _natoms;

    isotope_t _attachment_isotope = 0;

    // Number of rings in the fragment.
    int _nrings;

    // A list of ring sizes.
    resizable_array<int> _ring_sizes;

    // Number of connections (ring openings) in the fragment.
    int _ncon;

    // If two or more connections, the number of bonds between the closest two.
    int _distance;

    // In a molecule consisting of individual SafeFragment's each fragment
    // can be marked as ok to select or not - depending on substructure
    // matches in the parent molecule.
    int _ok_to_select;

  // Private functions
    int ProcessSquareBracket(const IWString& smi, int &i, int& next_ring,
                     IWString& destination);

  public:
    SafeFragment();

    int DebugPrint(std::ostream& output) const;

    const IWString& name() const {
      return _m.name();
    }
    void set_name(const IWString& s) {
      _m.set_name(s);
    }

    int natoms() const {
      return _natoms;
    }

    int ncon() const {
      return _ncon;
    }

    int nrings() const {
      return _nrings;
    }

    int largest_ring_size() const {
      int result = 0;
      for (int size : _ring_sizes) {
        result = std::max(result, size);
      }
      return result;
    }

    int distance() const {
      return _distance;
    }

    int Build(const const_IWSubstring& buffer, isotope_t attachment_isotope = 0);
    int Build(const dicer_data::DicerFragment& proto, isotope_t attachment_isotope = 0);

    int ok_to_select() const {
      return _ok_to_select;
    }
    void set_ok_to_select(int s) {
      _ok_to_select = s;
    }

    const IWString& smiles () const {
      return _smiles;
    }

    Molecule& mol() {
      return _m;
    }

    // Formula difference
    uint32_t FormulaDifference(const SafeFragment& rhs) const {
      return _mformula.Diff(rhs._mformula);
    }

    // Place into `new_smiles` the smiles of `f2` with the ring numbers
    // from `this`.
    int SameNumbers(const SafeFragment& f2, IWString& new_smiles) const;
};

// Examines all bonds in `m` and returns false if we find bonds that should
// be rejected.
int BondsOk(Molecule& m, isotope_t attachment_isotope = 0);

}  // namespace safe_generate

#endif  // MOLECULE_TOOLS_SAFE_GENERATE_LIB_H_
