#include <iostream> 
#include <random> 

#include "Molecule_Lib/molecule.h"

#include "safe_generate_lib.h"

namespace safe_generate {

constexpr char kPercent = '%';
constexpr char kOpenSquareBracket = '[';
constexpr char kCloseSquareBracket = ']';

using std::cerr;

int
IsSulfonamide(const Molecule& m, atom_number_t s, atom_number_t n) {
  assert(m.atomic_number(s) == lillymol::kSulphur);

  for (const Bond* b : m[s]) {
    atom_number_t o = b->other(s);
    if (o == n) {
      continue;
    }

    if (b->is_double_bond() && m.atomic_number(o) == 8) {
      return 1;
    }
  }

  return 0;
}

int
BondsOk(Molecule& m, isotope_t attachment_isotope) {
  m.compute_aromaticity_if_needed();

  for (const Bond* b : m.bond_list()) {
    // We don't form new rings, so ignore all ring bonds.
    if (b->nrings()) {
      continue;
    }
    if (! b->is_single_bond()) {
      continue;
    }

    atom_number_t a1 = b->a1();
    atom_number_t a2 = b->a2();

    // Bonds are only formed with isotopes
    if (!IsAttachmentIsotope(m.isotope(a1), attachment_isotope) ||
        !IsAttachmentIsotope(m.isotope(a2), attachment_isotope)) {
      continue;
    }

    atomic_number_t z1 = m.atomic_number(a1);
    atomic_number_t z2 = m.atomic_number(a2);
    // Carbon-Carbon always OK.
    if (z1 == 6 && z2 == 6) {
      continue;
    }

    // All NN or OO bonds fail
    if (z1 == 7 && z2 == 7) {
      return 0;
    }
    if (z1 == 8 && z2 == 8) {
      return 0;
    }
    // S-O bonds fail
    if (z1 == 8 && z2 == 16) {
      return 0;
    }
    if (z1 == 16 && z2 == 8) {
      return 0;
    }

    // One atom must be a carbon.
    if (z1 > z2) {
      std::swap(a1, a2);
      std::swap(z1, z2);
    }

    // All #6-F bonds are OK.
    if (z2 == 9) {
      if (z1 == 6) {
        continue;
      } else {
        return 0;
      }
    }

    // #6-#7 and #6-#8 are OK.
    if (z2 == 7 || z2 == 8) {
      continue;
    }

    // reject alkyl halides.
    if (z2 == 17 || z2 == 35 || z2 == 53) {
      if (m.is_aromatic(a1)) {  // aromatic attachments OK.
        continue;
      } else {
        return 0;
      }
    }

    if (z2 == lillymol::kSulphur) {
      if (z1 == 7) {
        if (IsSulfonamide(m, a2, a1)) {
          continue;
        } else {
          return 0;
        }
      }

      if (z1 == 6 || z1 == 8) {
        continue;
      } else {
        return 0;
      }
    }

    // P-O bonds are OK.
    if (z2 == 15) {
      if (z1 == 8 || z1 == 6) {
        continue;
      } else {
        return 0;
      }
    }

    return 0;
  }

  return 1;
}

// Keep the original all-isotope removal behavior when no marker is specified.
void
RemoveAttachmentIsotopes(Molecule& m, isotope_t attachment_isotope) {
  if (attachment_isotope == 0) {
    m.transform_to_non_isotopic_form();
    return;
  }
  for (int i = 0; i < m.natoms(); ++i) {
    if (m.isotope(i) == attachment_isotope) {
      m.set_isotope(i, 0);
    }
  }
}

SafeFragment::SafeFragment() {
  _natoms = 0;
  _nrings = 0;
  _ncon = 0;
  _distance = 0;
  _ok_to_select = 1;
}

int
SafeFragment::DebugPrint(std::ostream& output) const {
  output << "SafeFragment smiles " << _smiles << " natoms " << _natoms << " ncon " << _ncon << '\n';
  output << "first";;
  for (int f : _first_digit) {
    output << ' ' << f << " '" << _smiles[f] << "'";
  }
  output << '\n';

  return output.good();
}

int
SafeFragment::Build(const const_IWSubstring& smi, isotope_t attachment_isotope) {
  _attachment_isotope = attachment_isotope;
  const int nchars = smi.length();
  // cerr << "SafeFragment::Build from '" << smi << "'\n";
  // Must have at least an atomic symbol followed by % and a two digit ring number.
  if (nchars < 4) {
    cerr << "SafeFragment::Build:too short '" << smi << "'\n";
    return 0;
  }

  // All two-digit labels must be SAFE connections to other fragments.
  // Internal two-digit ring closures are outside this tool's input convention.
  const_IWSubstring smiles(smi);
  smiles.truncate_at_first(' ');
  bool seen[100] = {};
  IWString tmp;
  tmp.reserve(smiles.length());
  for (int i = 0; i < smiles.length(); ++i) {
    if (smiles[i] != kPercent) {
      tmp << smiles[i];
      continue;
    }
    if (i + 2 >= smiles.length() ||
        !std::isdigit(static_cast<unsigned char>(smiles[i + 1])) ||
        !std::isdigit(static_cast<unsigned char>(smiles[i + 2]))) {
      cerr << "SafeFragment::Build:invalid %nn ring label '" << smiles << "'\n";
      return 0;
    }
    const int label = (smiles[i + 1] - '0') * 10 + smiles[i + 2] - '0';
    if (seen[label]) {
      cerr << "SafeFragment::Build:internal two-digit ring closures are unsupported\n";
      return 0;
    }
    seen[label] = true;
    _first_digit << (i + 1);
    i += 2;
  }

  if (! _m.build_from_smiles(tmp)) {
    cerr << "SafeFragment::Build:invalid smiles '" << tmp << "' from '" << smi << "'\n";
    return 0;
  }

  _smiles = smi;
  _smiles.truncate_at_first(' ');
  _natoms = _m.natoms();
  _nrings = _m.nrings();
  _ncon = _first_digit.number_elements();

  for (const Ring* r : _m.sssr_rings()) {
    _ring_sizes.add_if_not_already_present(r->size());
  }

  int niso = 0;
  for (int i = 0; i < _natoms; ++i) {
    niso += IsAttachmentIsotope(_m.isotope(i), _attachment_isotope);
  }
  if (niso == 0) {
    cerr << "SafeFragment::Build:no attachment isotopes '" << _smiles << "' from '" << smi << "'\n";
    return 0;
  }

  if (niso > 1) {
    _distance = _natoms + 1;
    for (int i = 0; i < _natoms; ++i) {
      if (!IsAttachmentIsotope(_m.isotope(i), _attachment_isotope)) {
        continue;
      }

      for (int j = i + 1; j < _natoms; ++j) {
        if (!IsAttachmentIsotope(_m.isotope(j), _attachment_isotope)) {
          continue;
        }
        const int d = _m.bonds_between(i, j);
        if (d < _distance) {
          _distance = d;
        }
      }
    }
  }

  _mformula.Build(_m);

  return 1;
}

int
SafeFragment::ProcessSquareBracket(const IWString& smi, int &i, int& next_ring,
                     IWString& destination) {
  assert(smi[i] == kOpenSquareBracket);
  destination << kOpenSquareBracket;
  ++i;

  // Should not happen.
  if (i == smi.length()) {
    return 0;
  }

  // Read the full isotope number: e.g. marker 1 must not match isotope 13.
  const int isotope_start = i;
  int isotope_end = i;
  while (isotope_end < smi.length() &&
         std::isdigit(static_cast<unsigned char>(smi[isotope_end]))) {
    ++isotope_end;
  }
  isotope_t isotope = 0;
  if (isotope_end > isotope_start) {
    const_IWSubstring digits(smi.data() + isotope_start, isotope_end - isotope_start);
    if (!digits.numeric_value(isotope)) {
      return 0;
    }
  }
  const bool got_isotope = IsAttachmentIsotope(isotope, _attachment_isotope);

  for ( ;i < smi.length(); ++i) {
    destination << smi[i];
    if (smi[i] == kCloseSquareBracket) {
      break;
    }
  }

  if (got_isotope) {
    // Skip over any ring closures already present.
    for ( ; (i + 1) < smi.length() && std::isdigit(smi[i + 1]); ++i) {
      destination << smi[i + 1];
    }
    if (next_ring > 99) {
      cerr << "SafeFragment::ProcessSquareBracket:no available two-digit ring label\n";
      return 0;
    }
    destination << '%';
    _first_digit << destination.length();
    destination << next_ring;
    ++next_ring;
  }

  return 1;
}

int
SafeFragment::Build(const dicer_data::DicerFragment& proto, isotope_t attachment_isotope) {
  _attachment_isotope = attachment_isotope;
  IWString smi = proto.smi();
  if (smi.empty()) {
    cerr << "SafeFragment::Build:empty smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }
  if (smi.contains('%')) {
    cerr << "SafeFragment::Build:internal two-digit ring closures are unsupported\n";
    return 0;
  }
  if (! _m.build_from_smiles(smi)) {
    cerr << "SafeFragment::Build:invalid smiles " << proto.ShortDebugString() << '\n';
    return 0;
  }

  _m.set_name(proto.par());

  _mformula.Build(_m);

  _natoms = _m.natoms();
  _nrings = _m.nrings();
  for (const Ring* ring : _m.sssr_rings()) {
    _ring_sizes.add_if_not_already_present(ring->size());
  }
  // cerr << "SafeFragment::build:smiles " << smi << " natoms " << _natoms << '\n';

  _ncon = 0;
  for (int i = 0; i < _natoms; ++i) {
    _ncon += IsAttachmentIsotope(_m.isotope(i), _attachment_isotope);
  }
  if (_ncon == 0) {
    cerr << "SafeFragment::Build:no attachment isotopes '" << smi << "' from '" <<
            proto.ShortDebugString() << "'\n";
    return 0;
  }

  if (_ncon > 1) {
    _distance = _natoms + 1;
    for (int i = 0; i < _natoms; ++i) {
      if (!IsAttachmentIsotope(_m.isotope(i), _attachment_isotope)) {
        continue;
      }
      for (int j = i + 1; j < _natoms; ++j) {
        if (!IsAttachmentIsotope(_m.isotope(j), _attachment_isotope)) {
          continue;
        }
        const int d = _m.bonds_between(i, j);
        if (d < _distance) {
          _distance = d;
        }
      }
    }
  }

  // The smiles that is read does not have the %nn ring opening/closing motifs.
  // We insert those after the isotopic labels.

  IWString tmp;
  tmp.reserve(smi.length() + 9);
  int next_ring = 11;
  for (int i = 0; i < smi.length(); ++i) {
    const char c = smi[i];
    if (c == kOpenSquareBracket) {
      if (!ProcessSquareBracket(smi, i, next_ring, tmp)) {
        return 0;
      }
    } else {
      tmp << c;
    }
  }

  // cerr << "from " << smi << " transform to '" << tmp << " atoms " << _natoms << '\n';
  _smiles = tmp;

  return 1;
}

int
SafeFragment::SameNumbers(const SafeFragment& f2, IWString& new_smiles) const {
  if (_first_digit.size() != f2._first_digit.size()) {
    cerr << "SafeFragment::SameNumbers:attachment count mismatch\n";
    return 0;
  }
  new_smiles = f2._smiles;

  // Now copy our ring numbers to `new_smiles`.
  const int nf = _first_digit.number_elements();
  for (int i = 0; i < nf; ++i) {
    new_smiles[f2._first_digit[i]] = _smiles[_first_digit[i]];
    new_smiles[f2._first_digit[i] + 1] = _smiles[_first_digit[i] + 1];
  }

  return 1;
}

}  // namespace safe_generate
