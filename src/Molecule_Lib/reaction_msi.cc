// Reading and writing reactions in the msi (.rxn) file format.
//
// The member functions here are declared in iwreaction.h. They were moved out of iwreaction.cc,
// which now has the mechanics of performing reactions. The reaction proto form is in
// reaction_proto.cc.

#include <algorithm>
#include <iostream>
#include <memory>
#include <unordered_map>

using std::cerr;

#include "Foundational/iwmisc/misc.h"

#include "chiral_centre.h"
#include "istream_and_type.h"
#include "iwreaction.h"
#include "misc2.h"
#include "molecule.h"
#include "molecule_to_query.h"
#include "output.h"
#include "path.h"
#include "target.h"

static int warn_no_query_in_reaction_object = 1;
void
set_warn_no_query_in_reaction_object(int s) {
  warn_no_query_in_reaction_object = s;
}

/*
  An element transformation will look like

    (A C change_element "3 S")

  which means change atom 3 to a sulphur
*/

int
Reaction_Change_Element::construct_from_msi_attribute(const msi_attribute* att) {
  const IWString& rce = att->stringval();

  if (2 != rce.nwords()) {
    cerr << "Reaction_Change_Element::construct_from_msi_attribute: must have two string "
            "attributes\n";
    cerr << (*att) << '\n';
    cerr << rce << '\n';
    return 0;
  }

  const_IWSubstring token;
  rce.word(0, token);

  if (!token.numeric_value(_atom) || _atom < 0) {
    cerr << "Reaction_Change_Element::construct_from_msi_attribute: atom specifiers must "
            "be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  (void)rce.word(1, token);  // get the new element type

  _element = get_element_from_symbol_no_case_conversion(token);
  if (nullptr == _element) {
    cerr
        << "Reaction_Change_Element::construct_from_msi_attribute: unrecognised element '"
        << token << "'\n";
    cerr << (*att) << '\n';
    return 0;
  }

  return 1;
}

int
Reaction_Change_Element::write_msi(std::ostream& os, const const_IWSubstring& ind,
                                   const const_IWSubstring& attribute_name) const {
  os << ind << "  (A C " << attribute_name << " \"" << _atom << ' ' << _element->symbol()
     << "\")\n";

  return os.good();
}

int
Reaction_Formal_Charge::construct_from_msi_attribute(const msi_attribute* att) {
  if (2 != att->number_int_values()) {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: must have two int "
            "attributes\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _atom = att->int_multi_value(0);

  if (_atom < 0) {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: atom specifiers must "
            "be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _fc = att->int_multi_value(1);

  return 1;
}

int
Reaction_Change_Formal_Charge::construct_from_msi_attribute(const msi_attribute* att) {
  if (2 != att->number_int_values()) {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: must have two int "
            "attributes\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _atom = att->int_multi_value(0);

  if (_atom < 0) {
    cerr << "Reaction_Formal_Charge::construct_from_msi_attribute: atom specifiers must "
            "be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _delta = att->int_multi_value(1);

  return 1;
}

int
Reaction_Formal_Charge::write_msi(std::ostream& os, const const_IWSubstring& ind,
                                  const const_IWSubstring& attribute_name) const {
  os << ind << "  (A I " << attribute_name << " \"" << _atom << ' ' << _fc << "\")\n";

  return os.good();
}

int
Reaction_Change_Formal_Charge::write_msi(std::ostream& os, const const_IWSubstring& ind,
                                         const const_IWSubstring& attribute_name) const {
  os << ind << "  (A I " << attribute_name << " \"" << _atom << ' ' << _delta << "\")\n";

  return os.good();
}

int
Reaction_Place_Isotope::construct_from_msi_attribute(const msi_attribute* att) {
  if (2 != att->number_int_values()) {
    cerr << "Reaction_Place_Isotope::construct_from_msi_attribute: must have two int "
            "attributes\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _atom << att->int_multi_value(0);

  if (_atom[0] < 0) {
    cerr << "Reaction_Place_Isotope::construct_from_msi_attribute: atom specifiers must "
            "be non negative ints\n";
    cerr << (*att) << '\n';
    return 0;
  }

  _isotope = att->int_multi_value(1);

  return 1;
}

int
Reaction_Place_Isotope::write_msi(std::ostream& os, const const_IWSubstring& ind,
                                  const const_IWSubstring& attribute_name) const {
  os << ind << "  (A I " << attribute_name << " (" << _atom[0] << ' ' << _isotope
     << "))\n";

  return os.good();
}

static Bond*
parse_bond_attribute(const msi_attribute* att, bond_type_t bt = SINGLE_BOND) {
  if (2 != att->number_int_values()) {
    cerr << "A bond specification must have two int values\n";
    return 0;
  }

  int a1 = att->int_multi_value(0);
  int a2 = att->int_multi_value(1);

  if (a1 == a2) {
    cerr << "Identical atoms found at ends of bond " << a1 << " and " << a2 << '\n';
    cerr << (*att);
    return 0;
  }

  Bond* rc = new Bond(a1, a2, bt);

  assert(rc->ok());

  return rc;
}

static int
parse_atom_values(const msi_attribute* att, resizable_array<int>& zatoms) {
  if (1 == att->number_int_values()) {
    int tmp;
    if (!att->value(tmp) || tmp < 0) {
      cerr << "Atom values must be whole non-negative numbers\n";
      return 0;
    }

    zatoms.add(tmp);

    return 1;
  }

  int nv = att->number_int_values();
  for (int i = 0; i < nv; i++) {
    int j = att->int_multi_value(i);

    if (zatoms.contains(j)) {
      cerr << "Duplicate specification " << j << "\n";
      return 0;
    }

    zatoms.add(j);
  }

  return 1;
}

template <typename T>
int
process_attribute(const msi_attribute* att, resizable_array_p<T>& things) {
  T* t = new T;
  if (!t->construct_from_msi_attribute(att)) {
    cerr << "Scaffold_Reaction_Site::construct_from_msi_object: cannot parse " << (*att);
    delete t;
    return 0;
  }

  things.add(t);

  return 1;
}

template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Change_Element>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Formal_Charge>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Change_Formal_Charge>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Place_Isotope>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Increment_Isotope>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Invert_Isotope>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Stereo_Centre>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Dihedral_Angle>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Bond_Length>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_Bond_Angle>&);
template int process_attribute(const msi_attribute*,
                               resizable_array_p<Reaction_3D_Replace>&);

int
Reaction_Site::_parse_wedge_bonds(const msi_attribute* att) {
  Reaction_Wedge_Bond* w = new Reaction_Wedge_Bond;

  if (!w->construct_from_msi_attribute(att)) {
    cerr << "Reaction_Site::_parse_wedge_bonds: invalid attribute " << (*att) << '\n';
    delete w;
    return 0;
  }

  _wedge_bonds_to_place.add(w);

  return 1;
}

int
Reaction_Site::_parse_recompute_implicit_hydrogens(const msi_attribute* att) {
  int n = att->number_int_values();
  if (0 == n) {
    cerr << "Reaction_Site::_parse_recompute_implicit_hydrogens: no values\n";
    return 0;
  }

  if (_unfix_implicit_hydrogens.number_elements() + n <
      _unfix_implicit_hydrogens.elements_allocated()) {
    _unfix_implicit_hydrogens.resize(_unfix_implicit_hydrogens.elements_allocated() + n);
  }

  for (int i = 0; i < n; i++) {
    int j = att->int_multi_value(i);
    _unfix_implicit_hydrogens.add(j);
  }

  return 1;
}

int
Reaction_Site::_read_inactive_query(const msi_attribute* att) {
  const IWString& smarts = att->stringval();

  return _read_inactive_query(smarts);
}

int
Reaction_Site::_read_inactive_query(const IWString& smarts) {
  if (smarts.starts_with("Q:")) {  // single query file
    return _read_inactive_from_query_file(smarts);
  }

  const_IWSubstring tmp(smarts);

  if (smarts.starts_with("F:") || smarts.starts_with("Q:")) {
    return process_cmdline_token(' ', tmp, _inactive, 0);
  }

  Substructure_Query* q = new Substructure_Query;

  if (!q->create_from_smarts(smarts)) {
    cerr << "Reaction_Site::_read_inactive_query: invalid smarts '" << smarts << "'\n";
    delete q;
    return 0;
  }

  _inactive.add(q);

  return 1;
}

int
Reaction_Site::_read_inactive_from_query_file(const IWString& fname) {
  IWString tmp(fname);

  assert(tmp.starts_with("Q:"));

  tmp.remove_leading_chars(2);

  Substructure_Query* q = new Substructure_Query;

  if (!q->read(tmp)) {
    cerr << "Reaction_Site::_read_inactive_from_query_file: cannot read query file '"
         << fname << "'\n";
    delete q;

    return 0;
  }

  _inactive.add(q);

  return 1;
}

int
Reaction_Site::_read_query_file_from_molecule(IWString& fname) {
  FileType input_type = discern_file_type_from_name(fname);

  if (FILE_TYPE_INVALID == input_type) {
    input_type = FILE_TYPE_SMI;
  }

  data_source_and_type<Molecule> input(input_type, fname.null_terminated_chars());

  if (!input.good()) {
    cerr << "Reaction_Site::_read_query_file_from_molecule:cannot open '" << fname
         << "'\n";
    return 0;
  }

  Molecule* m = input.next_molecule();

  if (nullptr == m) {
    cerr << "Reaction_Site::_read_query_file_from_molecule:cannot read molecule from '"
         << fname << "'\n";
    return 0;
  }

  std::unique_ptr<Molecule> free_m(m);

  Molecule_to_Query_Specifications mqs;

  // Duplicate this until we eliminate the inheritance.
  if (!Substructure_Query::create_from_molecule(*m, mqs)) {
    cerr << "Reaction_Site::_read_query_file_from_molecule:cannot build query from "
         << m->smiles() << '\n';
    return 0;
  }

  Substructure_Query* q = new Substructure_Query();
  q->create_from_molecule(*m, mqs);
  _queries.add(q);

  return 1;
}

int
Reaction_Site::_read_query_file(const IWString& fname,
                                const int query_files_in_current_directory,
                                const IWString& reaction_directory) {
  if (fname.starts_with("M:")) {
    IWString tmp(fname);
    tmp.remove_leading_chars(2);
    return _read_query_file_from_molecule(tmp);
  }

  IWString query_fname;

  // cerr << "Parsing query file attribute, same " << query_files_in_current_directory <<
  // " rxn dir " << reaction_directory << " name '" << fname << "'\n";
  if (query_files_in_current_directory) {
    query_fname = fname;
  } else if (fname.starts_with('/')) {
    query_fname = fname;
  } else if (fname.contains("${")) {
    // shell variables will be expanded
    query_fname = fname;
  } else if (fname.starts_with("$")) {
    // Old style shell variable expansion, starts with $NAME
    IWString temp(fname);
    expand_environment_variables(temp.c_str(), query_fname);
  } else if (0 == reaction_directory.length()) {
    query_fname = fname;
  } else {
    query_fname = reaction_directory;
    query_fname << '/' << fname;
  }

  if (!Substructure_Query::read(query_fname)) {
    cerr << "Reaction_Site::_read_query_file_from_molecule: cannot process query '"
         << query_fname << "'\n";
    cerr << "In directory " << query_files_in_current_directory << " dir '"
         << reaction_directory << "'\n";
    return 0;
  }

  return 1;
}

/*
  Note that we'd really like things to be independent of the order of items
  in the input msi_object, but we really haven't achieved that. For example
  the NAME_OF_ONE_EMBEDDING_PER_START_ATOM_ATTRIBUTE attribute may fail if
  it appears AFTER the smiles. Fix someday if it becomes a problem
*/

int
Reaction_Site::construct_from_msi_object(const msi_object& msi,
                                         int query_files_in_current_directory,
                                         const IWString& reaction_directory) {
  int nmsi = msi.number_elements();

  int found_query = 0;  // must be specified only once

  _unique_id = msi.object_id();

  std::unordered_map<int, int> atom_map_to_unique_id;

  for (int i = 0; i < nmsi; i++) {
    const msi_object* obj = msi[i];

    if ("Query" == obj->name() || "query" == obj->name()) {
      if (!Substructure_Query::construct_from_msi_object(*obj)) {
        cerr << "Reaction_Site::construct_from_msi_object: query interpretation failed\n";
        cerr << (*obj);
        return 0;
      }

      found_query = 1;
    } else if (NAME_OF_NO_REACTION_OBJECT ==
               obj->name())  // will be parsed by the owning sidechain object
      ;
    else {
      cerr << "Reaction_Site::construct_from_msi_object: WHAT KIND OF OBJECT IS THIS '"
           << obj->name() << "'\n";
      return 0;
    }
  }

  if (_match_via_atom_map) {
    const int h = highest_initial_atom_number();

    for (int i = 0; i <= h; ++i) {
      const Substructure_Atom* a = query_atom_with_initial_atom_number(i);
      if (nullptr == a) {
        cerr << "Reaction_Site::construct_from_msi_object:no query atom with initial "
                "atom number "
             << i << '\n';
        return 0;
      }

      const int amap = a->atom_map_number();
      if (amap >= 0) {
        atom_map_to_unique_id[amap] = i;
      }
    }
  }

  int one_embedding = 0;  // NAME_OF_ONE_EMBEDDING_PER_START_ATOM_ATTRIBUTE present or not
  int unique_embeddings = 0;
  int ignore_symmetry = 0;
  int max_matches_allowed = 0;
  int embeddings_do_not_overlap = 0;

  int nat = msi.attribute_count();
  for (int i = 0; i < nat; i++) {
    const msi_attribute* att = msi.attribute(i);

    if (NAME_OF_QUERY_FILE_ATTRIBUTE == att->name()) {
      if (found_query) {
        cerr << "Reaction_Site::construct_from_msi_object: already found query\n";
        abort();
        return 0;
      }

      if (!_read_query_file(att->stringval(), query_files_in_current_directory,
                            reaction_directory)) {
        cerr << "Reaction_Site::construct_from_msi_object:cannot process query file '"
             << att->stringval() << "'\n";
        return 0;
      }

      found_query = 1;
    } else if (NAME_OF_ONE_EMBEDDING_PER_START_ATOM_ATTRIBUTE == att->name()) {
      att->value(one_embedding);
    } else if (NAME_OF_UNIQUE_EMBEDDINGS_ONLY_ATTRIBUTE == att->name()) {
      att->value(unique_embeddings);
    } else if (NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE == att->name()) {
      att->value(embeddings_do_not_overlap);
    } else if (NAME_OF_IGNORE_SYMMETRY_RELATED_MATCHES_ATTRIBUTE == att->name()) {
      att->value(ignore_symmetry);
    } else if (NAME_OF_REACTION_MAX_MATCHES_ATTRIBUTE == att->name()) {
      if (!att->value(max_matches_allowed) || max_matches_allowed < 1) {
        cerr << "IWReaction::_construct_from_msi_object: the '"
             << NAME_OF_REACTION_MAX_MATCHES_ATTRIBUTE
             << "' must be a whole positive number\n";
        return 0;
      }
    } else if (NAME_OF_QUERY_SMARTS_ATTRIBUTE == att->name()) {
      if (found_query) {
        cerr << "Reaction_Site::construct_from_msi_object: already found query\n";
        abort();
        return 0;
      }

      const_IWSubstring s;
      att->value(s);

      if (!Substructure_Query::create_from_smarts(s)) {
        cerr << "Reaction_Site::construct_from_msi_object: cannot parse smarts '" << s
             << "'\n";
        return 0;
      }

      _smarts << new IWString(s);

      found_query = 1;
    } else if (NAME_OF_INACTIVE_ATTRIBUTE == att->name()) {
      if (!_read_inactive_query(att)) {
        cerr << "Reaction_Site::construct_from_msi_object: invalid "
             << NAME_OF_INACTIVE_ATTRIBUTE << " attribite " << (*att) << '\n';
        return 0;
      }
    } else if (NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE == att->name()) {
      int tmp;
      if (!att->value(tmp)) {
        cerr << "Reaction_Site::construct_from_msi_object: invalid "
             << NAME_OF_EACH_COMPONENT_SEARCH_ATTRIBUTE << " attribute " << (*att)
             << '\n';
        return 0;
      }

      Substructure_Query::set_each_component_search(tmp);
    } else if (NAME_OF_BOND_TO_BREAK_ATTRIBUTE == att->name()) {
      Bond* b = parse_bond_attribute(att);
      if (nullptr == b) {
        return 0;
      }

      _bonds_to_be_broken.add(b);
    } else if (NAME_OF_SINGLE_BOND_ATTRIBUTE == att->name()) {
      Bond* b = parse_bond_attribute(att);
      if (nullptr == b) {
        return 0;
      }

      _bonds_to_be_made.add(b);
    } else if (NAME_OF_AROMATIC_BOND_ATTRIBUTE == att->name()) {
      Bond* b = parse_bond_attribute(att, AROMATIC_BOND);
      if (nullptr == b) {
        return 0;
      }

      _bonds_to_be_made.add(b);
    } else if (NAME_OF_DOUBLE_BOND_ATTRIBUTE == att->name()) {
      Bond* b = parse_bond_attribute(att, DOUBLE_BOND);
      if (nullptr == b) {
        return 0;
      }

      _bonds_to_be_made.add(b);
    } else if (NAME_OF_TRIPLE_BOND_ATTRIBUTE == att->name()) {
      Bond* b = parse_bond_attribute(att, TRIPLE_BOND);
      if (nullptr == b) {
        return 0;
      }

      _bonds_to_be_made.add(b);
    } else if (NAME_OF_BOND_TO_BREAK_ATTRIBUTE == att->name()) {
      Bond* b = parse_bond_attribute(att);
      if (nullptr == b) {
        return 0;
      }

      _bonds_to_be_broken.add(b);
    } else if (NAME_OF_REMOVE_ATOM_ATTRIBUTE == att->name()) {
      if (!parse_atom_values(att, _atoms_to_be_removed)) {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    } else if (NAME_OF_REMOVE_FRAGMENT_ATTRIBUTE == att->name()) {
      if (!parse_atom_values(att, _fragments_to_be_removed)) {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    } else if (NAME_OF_KEEP_FRAGMENT_ATTRIBUTE == att->name()) {
      if (!parse_atom_values(att, _fragments_to_be_kept)) {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    } else if (NAME_OF_RECOMPUTE_IMPLICIT_HYDROGENS_ATTRIBUTE == att->name()) {
      if (!_parse_recompute_implicit_hydrogens(att)) {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    } else if (NAME_OF_WEDGE_BOND_ATTRIBUTE == att->name()) {
      if (!_parse_wedge_bonds(att)) {
        cerr << "Cannot parse '" << (*att) << "'\n";
        return 0;
      }
    } else if (NAME_OF_NOOP_REACTION == att->name()) {
      _noop_reaction = 1;
    } else if (NAME_OF_INTER_PARTICLE_BOND_ATTRIBUTE ==
               att->name())  // ignore these, handled by derived object only
      ;
    else if (NAME_OF_INTER_PARTICLE_SUBSTITUTE_ATOM_ATTRIBUTE ==
             att->name())  // ignore these, handled by derived object only
      ;
    else if (NAME_OF_SUBSTITUTE_ATOM_ATTRIBUTE ==
             att->name())  // ignore these, handled by derived object only
      ;
    else if (NAME_OF_SIDECHAIN_SMILES == att->name())
      ;
    else if (NAME_OF_SIDECHAIN_MOLECULE == att->name())
      ;
    else if (NAME_OF_SIDECHAIN_FILE == att->name())
      ;
    else if (NAME_OF_CHANGE_ELEMENT_ATRRIBUTE == att->name()) {
      if (!process_attribute(att, _elements_to_change)) {
        return 0;
      }
    } else if (NAME_OF_FORMAL_CHARGE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _formal_charges_to_assign)) {
        return 0;
      }
    } else if (NAME_OF_CHANGE_FORMAL_CHARGE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _formal_charges_to_change)) {
        return 0;
      }
    } else if (NAME_OF_ISOTOPE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _isotopes_to_assign)) {
        return 0;
      }
    } else if (NAME_OF_INCREMENT_ISOTOPE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _isotopes_to_increment)) {
        return 0;
      }
    } else if (NAME_OF_INVERT_ISOTOPE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _isotopes_to_invert)) {
        return 0;
      }
    } else if (NAME_OF_DIHEDRAL_ANGLE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _reaction_dihedral_angle)) {
        return 0;
      }

      _3d = 1;
    } else if (NAME_OF_BOND_LENGTH_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _reaction_bond_length)) {
        return 0;
      }

      _3d = 1;
    } else if (NAME_OF_BOND_ANGLE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _reaction_bond_angle)) {
        return 0;
      }

      _3d = 1;
    } else if (NAME_OF_3D_REPLACE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _reaction_3d_replace)) {
        return 0;
      }

      _3d = 1;
    } else if (NAME_OF_INVERT_STEREO_CENTRE_ATTRIBUTE == att->name()) {
      int tmp;
      if (!att->value(tmp) || tmp < 0) {
        cerr << "Reaction_Site::construct_from_msi_object: invalid "
             << NAME_OF_INVERT_STEREO_CENTRE_ATTRIBUTE << " value\n";
        return 0;
      }

      _stereo_centres_to_invert.add(tmp);
    } else if (NAME_OF_REACTION_COMMENT == att->name()) {
      att->value(_comment);
    } else if (NAME_OF_REMOVE_CHIRAL_CENTRE_ATTRIBUTE == att->name()) {
      for (int i = 0; i < att->number_int_values(); i++) {
        _chiral_centres_to_remove.add(att->int_multi_value(i));
      }
    } else if (NAME_OF_TOGGLE_KEKULE_FORM_ATTRIBUTE == att->name()) {
      if (!_toggle_kekule_form.add_bond_from_msi_attribute(*att)) {
        cerr << "Reaction_Site::construct_from_msi_object:invalid toggke kekule form "
                "specification '"
             << (*att) << "'\n";
        return 0;
      }
    } else if (NAME_OF_UNLINK_UNMATCHED_ATOMS == att->name()) {
      if (!parse_atom_values(att, _unlink_unmatched_atoms)) {
        cerr << "Reaction_Site::construct_from_msi_object:invalid attribute " << (*att)
             << '\n';
        return 0;
      }
    } else {
      cerr << "Reaction_Site::construct_from_msi_object: unrecognised attribute '"
           << att->name() << "'\n";
      return 0;
    }
  }

  /*
    We don't want things to be dependent on the order in the msi file. Otherwise
    these directives would only work if they preceeded the smarts...
  */

  if (found_query) {
    if (unique_embeddings) {
      set_find_unique_embeddings_only(1);
    }
    if (ignore_symmetry) {
      set_do_not_perceive_symmetry_equivalent_matches(1);
    }
    if (one_embedding) {
      set_find_one_embedding_per_atom(1);
    }
    if (max_matches_allowed) {
      set_max_matches_to_find(max_matches_allowed);
    }
    if (embeddings_do_not_overlap) {
      set_embeddings_do_not_overlap(embeddings_do_not_overlap);
    }

    return 1;
  }

  if (_noop_reaction) {
    return 1;
  }

  if (warn_no_query_in_reaction_object) {
    cerr << "Reaction_Site::construct_from_msi_object: no query found in object\n";
    cerr << "Using default '[*]', which will match the first atom regardless of type\n";
  }

  Substructure_Query::create_from_smarts("[*]");
  Substructure_Query::set_max_matches_to_find(1);

  return 1;
}

int
Scaffold_Reaction_Site::construct_from_msi_object(const msi_object& msi,
                                                  int query_files_in_current_directory,
                                                  const IWString& reaction_directory) {
  if (!Reaction_Site::construct_from_msi_object(msi, query_files_in_current_directory,
                                                reaction_directory)) {
    return 0;
  }

  // for (int i = 0; i < _reaction_stereo_centre.number_elements(); i++)
  //{
  //   _reaction_stereo_centre[i]->all_atoms_in_scaffold();
  // }

  for (int i = 0; i < _reaction_dihedral_angle.number_elements(); i++) {
    _reaction_dihedral_angle[i]->all_atoms_in_scaffold();
  }

  for (int i = 0; i < _reaction_bond_length.number_elements(); i++) {
    _reaction_bond_length[i]->all_atoms_in_scaffold();
  }

  for (int i = 0; i < _reaction_bond_angle.number_elements(); i++) {
    _reaction_bond_angle[i]->all_atoms_in_scaffold();
  }

  // Jul 2025. Need to transfer some of the info to the _match_conditions object.
  // This should not alter anything.
  int nat = msi.attribute_count();
  for (int i = 0; i < nat; i++) {
    const msi_attribute* att = msi.attribute(i);

    if (NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE == att->name()) {
      int s;
      att->value(s);
      _match_conditions.set_embeddings_can_overlap(!s);
    } else if (NAME_OF_REACTION_MAX_MATCHES_ATTRIBUTE == att->name()) {
      int s;
      att->value(s);
      _match_conditions.set_max_matches_to_find(s);
    }
  }

  return 1;
}

/*
  Inter particle bonds will have two or three components.

    (A I join (0 0))          means make a single bond between atoms 0 and 0
    (A I join (0 0 1))        same
    (A I join (1 2 2))        make a double bond between atoms 1 and 2
    (A C join "S1.1 2")       make single bond between atom 1 in sidechain 1 and matched
  atom 2 (A C join "1 S.0 2")
*/

static Inter_Particle_Bond*
parse_inter_particle_bond_attribute(const IWString& s) {
  Inter_Particle_Bond* rc = new Inter_Particle_Bond;
  if (!rc->build(s)) {
    delete rc;
    return nullptr;
  }

  return rc;
}

static Replace_Atom*
parse_replace_atom_attribute(const msi_attribute* att, int default_component) {
  Replace_Atom* r = new Replace_Atom;
  if (!r->construct_from_msi_object(att, default_component)) {
    delete r;
    return nullptr;
  }

  return r;
}

int
Replace_Atom::write_msi(std::ostream& os, const const_IWSubstring& ind,
                        const const_IWSubstring& attribute_name) const {
  os << ind << "  (A C " << attribute_name << " \"" << _a1 << ' ' << _a2 << '"' << ")\n";

  return os.good();
}

/*
  A replace_atom directive must specify two matched atoms
  We get the default component to use from the object that owns us

    (A C replace_atom "S.0 2")      replace scaffold atom 0 with sidechain atom 2
    (A C replace_atom "0 4")        replace sidechain atom 0 with sidechain atom 4
*/

int
Replace_Atom::construct_from_msi_object(const msi_attribute* att, int default_component) {
  const IWString& s = att->stringval();

  if (2 != s.nwords()) {
    cerr << "Replace_Atom::construct_from_msi_object: must have two tokens\n";
    return 0;
  }

  const_IWSubstring s1, s2;
  (void)s.split(s1, ' ', s2);

  if (!_a1.construct(s1, default_component)) {
    cerr << "Replace_Atom::construct_from_msi_object: invalid a1 specification '" << s1
         << "'\n";
    return 0;
  }

  if (!_a2.construct(s2, default_component)) {
    cerr << "Replace_Atom::construct_from_msi_object: invalid a2 specification '" << s2
         << "'\n";
    return 0;
  }

  return 1;
}

int
Sidechain_Reaction_Site::construct_from_msi_object(
    const msi_object& msi, int query_files_in_current_directory,
    const IWString& reaction_directory, const Sidechain_Match_Conditions& smc) {
  _match_conditions = smc;

  if (!Reaction_Site::construct_from_msi_object(msi, query_files_in_current_directory,
                                                reaction_directory)) {
    return 0;
  }

  _copy_match_conditions_to_query();

  const msi_attribute* att;
  int i = 0;
  while (nullptr != (att = msi.attribute(NAME_OF_INTER_PARTICLE_BOND_ATTRIBUTE, i++))) {
    Inter_Particle_Bond* b = parse_inter_particle_bond_attribute(att->stringval());
    if (nullptr == b) {
      return 0;
    }

    b->is_part_of_component(_unique_id);

    _inter_particle_bonds.add(b);
  }

  i = 0;
  while (nullptr != (att = msi.attribute(NAME_OF_SUBSTITUTE_ATOM_ATTRIBUTE, i++))) {
    Replace_Atom* b = parse_replace_atom_attribute(att, _sidechain_number);
    if (nullptr == b) {
      return 0;
    }

    _replace_atom.add(b);
  }

  if (_inter_particle_bonds.empty()) {
    if (!_noop_reaction) {
      cerr << "Sidechain_Reaction_Site::construct_from_msi_object: no joins specified\n";
    }
  }

  const msi_object* nr;
  i = 0;
  while (nullptr != (nr = msi.component(NAME_OF_NO_REACTION_OBJECT, i++))) {
    No_Reaction* tmp = new No_Reaction;
    if (!tmp->construct_from_msi_object(*nr)) {
      delete tmp;
      return 0;
    }

    _no_reaction.add(tmp);
  }

  // Jul 2025. Need to transfer some of the info to the _match_conditions object.
  // This should not alter anything.
  int nat = msi.attribute_count();
  for (int i = 0; i < nat; i++) {
    const msi_attribute* att = msi.attribute(i);

    if (NAME_OF_EMBEDDINGS_DO_NOT_OVERLAP_ATTRIBUTE == att->name()) {
      int s;
      att->value(s);
      _match_conditions.set_embeddings_can_overlap(!s);
    } else if (NAME_OF_REACTION_MAX_MATCHES_ATTRIBUTE == att->name()) {
      int s;
      att->value(s);
      set_max_matches_to_find(s);
    }
  }

  att = msi.attribute(NAME_OF_MAKE_IMPLICIT_HYDROGEN_EXPLICIT_ATTRIBUTE);
  if (nullptr != att) {
    att->value(_make_implicit_hydrogens_explicit);
  }

  // Is the fragment to add present?

  att = msi.attribute(NAME_OF_SIDECHAIN_SMILES);
  if (nullptr != att) {
    return _add_reagents_from_smiles(*att);
  }

  att = msi.attribute(NAME_OF_SIDECHAIN_FILE);
  if (nullptr != att) {
    return _add_reagents_from_file(*att, smc);
  }

  return 1;
}

int
Sidechain_Reaction_Site::_add_reagents_from_smiles(const msi_attribute& att) {
  const_IWSubstring smiles;
  att.value(smiles);

  // ASAN reports this as a memory leak, but I don't see how that could be the case.
  std::unique_ptr<Molecule_and_Embedding> reagent =
      std::make_unique<Molecule_and_Embedding>();

  if (!reagent->build_from_smiles(smiles)) {
    cerr << "Sidechain_Reaction_Site::construct_from_msi_object: cannot parse '" << att
         << '\n';
    return 0;
  }

  // Now see if the query matches

  Substructure_Results sresults;

  const int nhits = _determine_matched_atoms_checking_inactives(*reagent, sresults);

  if (1 == nhits)
    ;
  else if (0 == nhits) {
    if (_noop_reaction)
      ;
    else {
      cerr << "Sidechain_Reaction_Site::construct_from_msi_object: yipes, no match to "
              "query\n";
      cerr << att;
      return 0;
    }
  } else if (nhits > 1) {
    cerr << "Sidechain_Reaction_Site::construct_from_msi_object: yipes, " << nhits
         << " matches to query\n";
    return 0;
  }

  if (!_noop_reaction) {
    reagent->collect_matched_atoms(sresults);

    if (_toggle_kekule_form.active()) {
      if (!reagent->do_toggle_kekule_form(_toggle_kekule_form)) {
        return 0;
      }
    }
  }

  _reagents.add(reagent.release());

  return 1;
}

int
Sidechain_Reaction_Site::_add_reagents_from_file(const msi_attribute& att,
                                                 const Sidechain_Match_Conditions& smc) {
  IWString fname;
  att.value(fname);

  data_source_and_type<Molecule_and_Embedding> input(FILE_TYPE_SMI, fname);

  if (!input.good()) {
    cerr << "Sidechain_Reaction_Site::_add_reagents_from_file:cannot open '" << fname
         << "'\n";
    return 0;
  }

  return add_reagents(fname.null_terminated_chars(), FILE_TYPE_SMI, smc);
}

static int
write_array(std::ostream& os, const const_IWSubstring& ind,
            const resizable_array<int>& things, const const_IWSubstring& attribute_name) {
  int nt = things.number_elements();

  os << ind << "  (A I " << attribute_name << ' ';

  if (nt > 1) {
    os << '(';
  }

  for (int i = 0; i < nt; i++) {
    if (i > 0) {
      os << ' ';
    }
    os << things[i];
  }

  if (nt > 1) {
    os << ')';
  }

  os << ")\n";

  return os.good();
}

template <typename T>
int
write_set(std::ostream& os, const const_IWSubstring& ind,
          const resizable_array_p<T>& things, const const_IWSubstring& attribute_name) {
  for (int i = 0; i < things.number_elements(); i++) {
    (void)things[i]->write_msi(os, ind, attribute_name);
  }

  return os.good();
}

template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Change_Element>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Formal_Charge>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Change_Formal_Charge>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Place_Isotope>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Increment_Isotope>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Invert_Isotope>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Stereo_Centre>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Wedge_Bond>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Dihedral_Angle>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Bond_Length>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_Bond_Angle>&,
                       const const_IWSubstring&);
template int write_set(std::ostream&, const const_IWSubstring&,
                       const resizable_array_p<Reaction_3D_Replace>&,
                       const const_IWSubstring&);

static int
write_bond(std::ostream& os, const IWString& ind, const char* attribute_name,
           const Bond* b) {
  os << ind << "  (A I " << attribute_name << " (" << b->a1() << ' ' << b->a2() << "))\n";

  return os.good();
}

static int
write_atoms(std::ostream& os, const IWString& ind, const char* attribute_name,
            const resizable_array<int>& atoms) {
  os << ind << "  (A I " << attribute_name << ' ';

  int na = atoms.number_elements();
  if (1 == na) {
    os << atoms[0] << ")\n";
    return os.good();
  }

  os << '(';
  for (int i = 0; i < na; i++) {
    if (i > 0) {
      os << ' ';
    }
    os << atoms[i];
  }

  os << "))\n";

  return os.good();
}

int
Inter_Particle_Bond::write_msi(std::ostream& os, const IWString& ind) const {
  os << ind << "  (A C " << NAME_OF_INTER_PARTICLE_BOND_ATTRIBUTE << " \"" << _a1 << ' '
     << _a2 << ' ' << character_representation_of_bond(_bt) << "\")\n";

  return os.good();
}

/*
  writing a reaction_side object is complicated by the fact that the
  underlying object will never be written as a complete object
*/

int
Reaction_Site::_write_msi(std::ostream& os, int& object_id, const IWString& ind) {
  if (_comment.length()) {
    os << ind << "  (A C " << NAME_OF_REACTION_COMMENT << " \"" << _comment << "\")\n";
  }

  if (!Substructure_Query::write_msi(os, object_id, ind.nchars() + 2)) {
    return 0;
  }

  if (_bonds_to_be_made.number_elements()) {
    int nb = _bonds_to_be_made.number_elements();
    for (int i = 0; i < nb; i++) {
      const Bond* b = _bonds_to_be_made[i];
      if (b->is_single_bond()) {
        write_bond(os, ind, NAME_OF_SINGLE_BOND_ATTRIBUTE, b);
      } else if (b->is_double_bond()) {
        write_bond(os, ind, NAME_OF_DOUBLE_BOND_ATTRIBUTE, b);
      } else if (b->is_triple_bond()) {
        write_bond(os, ind, NAME_OF_TRIPLE_BOND_ATTRIBUTE, b);
      } else {
        cerr << "Reaction_Site::write_msi: What kind of bond is this " << (*b) << '\n';
        return 0;
      }
    }
  }

  if (_bonds_to_be_broken.number_elements()) {
    int nb = _bonds_to_be_broken.number_elements();
    for (int i = 0; i < nb; i++) {
      const Bond* b = _bonds_to_be_broken[i];
      write_bond(os, ind, NAME_OF_BOND_TO_BREAK_ATTRIBUTE, b);
    }
  }

  if (_atoms_to_be_removed.number_elements()) {
    write_atoms(os, ind, NAME_OF_REMOVE_ATOM_ATTRIBUTE, _atoms_to_be_removed);
  }

  if (_fragments_to_be_removed.number_elements()) {
    write_atoms(os, ind, NAME_OF_REMOVE_FRAGMENT_ATTRIBUTE, _fragments_to_be_removed);
  }

  if (_fragments_to_be_kept.number_elements()) {
    write_atoms(os, ind, NAME_OF_KEEP_FRAGMENT_ATTRIBUTE, _fragments_to_be_kept);
  }

  if (_elements_to_change.number_elements()) {
    write_set(os, ind, _elements_to_change, NAME_OF_CHANGE_ELEMENT_ATRRIBUTE);
  }

  if (_formal_charges_to_assign.number_elements()) {
    write_set(os, ind, _formal_charges_to_assign, NAME_OF_FORMAL_CHARGE_ATTRIBUTE);
  }

  if (_formal_charges_to_change.number_elements()) {
    write_set(os, ind, _formal_charges_to_change, NAME_OF_CHANGE_FORMAL_CHARGE_ATTRIBUTE);
  }

  if (_isotopes_to_assign.number_elements()) {
    write_set(os, ind, _isotopes_to_assign, NAME_OF_ISOTOPE_ATTRIBUTE);
  }

  if (_isotopes_to_increment.number_elements()) {
    write_set(os, ind, _isotopes_to_increment, NAME_OF_INCREMENT_ISOTOPE_ATTRIBUTE);
  }

  if (_isotopes_to_invert.number_elements()) {
    write_set(os, ind, _isotopes_to_invert, NAME_OF_INVERT_ISOTOPE_ATTRIBUTE);
  }

  if (_unfix_implicit_hydrogens.number_elements()) {
    write_array(os, ind, _unfix_implicit_hydrogens,
                NAME_OF_RECOMPUTE_IMPLICIT_HYDROGENS_ATTRIBUTE);
  }

  if (_wedge_bonds_to_place.number_elements()) {
    write_set(os, ind, _wedge_bonds_to_place, NAME_OF_WEDGE_BOND_ATTRIBUTE);
  }

  if (_reaction_dihedral_angle.number_elements()) {
    write_set(os, ind, _reaction_dihedral_angle, NAME_OF_DIHEDRAL_ANGLE_ATTRIBUTE);
  }

  if (_reaction_bond_length.number_elements()) {
    write_set(os, ind, _reaction_bond_length, NAME_OF_BOND_LENGTH_ATTRIBUTE);
  }

  if (_reaction_bond_angle.number_elements()) {
    write_set(os, ind, _reaction_bond_angle, NAME_OF_BOND_ANGLE_ATTRIBUTE);
  }

  if (_reaction_3d_replace.number_elements()) {
    write_set(os, ind, _reaction_3d_replace, NAME_OF_3D_REPLACE_ATTRIBUTE);
  }

  if (_chiral_centres_to_remove.number_elements()) {
    write_array(os, ind, _chiral_centres_to_remove,
                NAME_OF_REMOVE_CHIRAL_CENTRE_ATTRIBUTE);
  }

  if (_unlink_unmatched_atoms.number_elements()) {
    write_array(os, ind, _unlink_unmatched_atoms, NAME_OF_UNLINK_UNMATCHED_ATOMS);
  }

  for (int i = 0; i < _stereo_centres_to_invert.number_elements(); i++) {
    os << ind << "  (A I " << NAME_OF_INVERT_STEREO_CENTRE_ATTRIBUTE << ' '
       << _stereo_centres_to_invert[i] << ")\n";
  }

  for (int i = 0; i < _replace_atom.number_elements(); i++) {
    _replace_atom[i]->write_msi(os, ind, NAME_OF_SUBSTITUTE_ATOM_ATTRIBUTE);
  }

  return os.good();
}

int
Sidechain_Reaction_Site::write_msi(std::ostream& os, int& object_id, int indentation) {
  assert(ok());

  IWString ind;
  ind.extend(indentation, ' ');

  os << ind << "(" << object_id++ << ' ' << NAME_OF_SIDECHAIN_OBJECT << '\n';

  Reaction_Site::_write_msi(os, object_id, ind);

  if (_make_implicit_hydrogens_explicit) {
    os << ind << "  (" << NAME_OF_MAKE_IMPLICIT_HYDROGEN_EXPLICIT_ATTRIBUTE << ' '
       << _make_implicit_hydrogens_explicit << ")\n";
  }

  for (int i = 0; i < _inter_particle_bonds.number_elements(); i++) {
    _inter_particle_bonds[i]->write_msi(os, ind);
  }

  for (int i = 0; i < _no_reaction.number_elements(); i++) {
    _no_reaction[i]->write_msi(os, object_id, indentation + 2);
  }

  os << ind << ")\n";

  return os.good();
}

int
Scaffold_Reaction_Site::write_msi(std::ostream& os, int& object_id, int indentation) {
  assert(ok());

  IWString ind;
  ind.extend(indentation, ' ');

  os << ind << "(" << object_id++ << ' ' << NAME_OF_SCAFFOLD_OBJECT << '\n';

  Reaction_Site::_write_msi(os, object_id, ind);

  os << ind << ")\n";

  return os.good();
}

int
IWReaction::_construct_from_msi_object(
    const msi_object& msi, const Sidechain_Match_Conditions& sidechain_match_conditions) {
  assert(_sidechains.empty());

  int nmsi = msi.number_elements();
  // if (nmsi < 2)
  //{
  //   cerr << "IWReaction::construct_from_msi_object: a reaction object must have at
  //   least two components\n"; return 0;
  // }

  const msi_attribute* att = msi.attribute(NAME_OF_MATCH_VIA_ATOM_MAP);
  if (nullptr != att) {  // should probably check the value
    _match_via_atom_map = 1;
  }

  extending_resizable_array<int> unique_id_encountered;  // must have unique object ID's

  for (int i = 0; i < nmsi; i++) {
    const msi_object* obj = msi[i];

    int oid = obj->object_id();
    if (0 == unique_id_encountered[oid]) {
      unique_id_encountered[oid] = 1;
    } else {
      cerr << "IWReaction::_construct_from_msi_object:duplicate object id " << oid
           << '\n';
      return 0;
    }

    if (NAME_OF_SCAFFOLD_OBJECT == obj->name()) {
      Reaction_Site::set_match_via_atom_map(_match_via_atom_map);

      if (!Scaffold_Reaction_Site::construct_from_msi_object(
              *obj, _query_files_in_current_directory, _reaction_directory)) {
        return 0;
      }
    } else if (NAME_OF_SIDECHAIN_OBJECT == obj->name() ||
               OLD_NAME_OF_SIDECHAIN_OBJECT == obj->name()) {
      Sidechain_Reaction_Site* s = new Sidechain_Reaction_Site;

      s->set_match_via_atom_map(_match_via_atom_map);

      s->set_sidechain_number(_sidechains.number_elements());

      if (_make_implicit_hydrogens_explicit) {
        s->set_make_implicit_hydrogens_explicit(_make_implicit_hydrogens_explicit);
      }

      if (!s->construct_from_msi_object(*obj, _query_files_in_current_directory,
                                        _reaction_directory,
                                        sidechain_match_conditions)) {
        return 0;
      }

      _sidechains.add(s);
    } else {
      cerr << "IWReaction::_construct_from_msi_object: what kind of object is this?\n";
      cerr << (*obj);
      return 0;
    }
  }

  int nat = msi.attribute_count();
  for (int i = 0; i < nat; i++) {
    const msi_attribute* att = msi.attribute(i);

    if (NAME_OF_REACTION_COMMENT == att->name()) {
      att->value(_comment);
    } else if (NAME_OF_APPEND_TO_NAME == att->name()) {
      att->value(_append_to_name);
    } else if (NAME_OF_FIND_KEKULE_FORM_FOR_INVALID_VALENCES_ATTRIBUTE == att->name()) {
      att->value(_find_kekule_forms_for_bad_valence);
    } else if (NAME_OF_STEREO_CENTRE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _reaction_stereo_centre)) {
        return 0;
      }
    } else if (NAME_OF_OPTIONAL_STEREO_CENTRE_ATTRIBUTE == att->name()) {
      if (!process_attribute(att, _reaction_stereo_centre)) {
        return 0;
      }

      _reaction_stereo_centre.last_item()->set_optional(1);
    } else if (NAME_OF_NOOP_REACTION == att->name()) {
      _noop_reaction = 1;
    } else {
      cerr << "IWReaction::_construct_from_msi_object: unrecognised attribute '" << (*att)
           << "'\n";
      return 0;
    }
  }

  return 1;
}

int
IWReaction::construct_from_msi_object(
    const msi_object& msi, const Sidechain_Match_Conditions& sidechain_match_conditions) {
  int rc = _construct_from_msi_object(msi, sidechain_match_conditions);

  _convert_from_msi_object_numbers_to_internal_indices();

  if (rc) {
    rc = check_internal_consistency();
  }

  if (rc) {
    return rc;
  }

  cerr << msi;

  return 0;
}

/*
  We have read the msi file. Many of the matched atom directives will be in terms of
  msi object numbers. We must convert those to our internal numbering
*/

int
IWReaction::_convert_from_msi_object_numbers_to_internal_indices() {
  extending_resizable_array<int> xref(-9);

  int u = Reaction_Site::unique_id();

  xref[u] = -1;  // special index for scaffold

  for (int i = 0; i < _sidechains.number_elements(); i++) {
    u = _sidechains[i]->unique_id();

    xref[u] = i;
  }

// #define DEBUG_CONVERT_FROM_MSI_OBJECT_NUMBERS_TO_INTERNAL_INDICES
#ifdef DEBUG_CONVERT_FROM_MSI_OBJECT_NUMBERS_TO_INTERNAL_INDICES
  for (int i = 0; i < xref.number_elements(); i++) {
    if (-9 != xref[i]) {
      cerr << "Object ID " << i << " becomes " << xref[i] << '\n';
    }
  }
#endif

  Reaction_Site::adjust_matched_atoms_in_component(xref);

  Reaction_Site::set_unique_id(0);

  for (int i = 0; i < _sidechains.number_elements(); i++) {
    _sidechains[i]->adjust_matched_atoms_in_component(xref);

    _sidechains[i]->set_unique_id(i + 1);
  }

  return 1;
}

int
IWReaction::write_msi(const char* fname) {
  std::ofstream output(fname, std::ios::out);
  if (!output.good()) {
    cerr << "IWReaction::write_msi:cannot open '" << fname << "'\n";
    return 0;
  }

  return write_msi(output);
}

int
IWReaction::write_msi(std::ostream& os) {
  int object_id = 0;

  os << "(" << object_id++ << ' ' << NAME_OF_REACTION_OBJECT << '\n';

  if (_comment.length()) {
    os << "  (A C " << NAME_OF_REACTION_COMMENT << " \"" << _comment << "\")\n";
  }

  if (_append_to_name.length()) {
    os << "  (A C " << NAME_OF_APPEND_TO_NAME << " \"" << _append_to_name << "\")\n";
  }

  if (_find_kekule_forms_for_bad_valence) {
    os << "  (A I " << NAME_OF_FIND_KEKULE_FORM_FOR_INVALID_VALENCES_ATTRIBUTE << ' '
       << _find_kekule_forms_for_bad_valence << ")\n";
  }

  if (_toggle_kekule_form.active()) {
    _toggle_kekule_form.write_msi(os, "  ", NAME_OF_TOGGLE_KEKULE_FORM_ATTRIBUTE);
  }

  // Because the reaction stereo centres may have different attribute names, we do them
  // explicitly Hmmm, they should know their attribute name...

  IWString ind("  ");

  for (int i = 0; i < _reaction_stereo_centre.number_elements(); i++) {
    Reaction_Stereo_Centre* r = _reaction_stereo_centre[i];
    if (r->optional()) {
      r->write_msi(os, ind, NAME_OF_OPTIONAL_STEREO_CENTRE_ATTRIBUTE);
    } else {
      r->write_msi(os, ind, NAME_OF_STEREO_CENTRE_ATTRIBUTE);
    }
  }

  (void)Scaffold_Reaction_Site::write_msi(os, object_id, 2);

  int nsd = _sidechains.number_elements();

  for (int i = 0; i < nsd; i++) {
    _sidechains[i]->write_msi(os, object_id, 2);
  }

  os << ")\n";

  return os.good();
}

/*
  A No_Reaction object looks like

    (0 No_Reaction
      (A C smarts "xxx");
      (A C smarts "yyy")
    )

  where the first query specifies the scaffold and the 2nd specifies the sidechain
*/

int
No_Reaction::construct_from_msi_object(const msi_object& msi) {
  if (0 != msi.number_elements() && 2 != msi.number_elements()) {
    cerr << "No_Reaction::construct_from_msi_object: must specify 0 or 2 query objects\n";
    return 0;
  }

  if (2 == msi.number_elements())  // queries written in expanded form in the object
  {
    if (!_scaffold_no_reaction.construct_from_msi_object(*(msi[0]))) {
      cerr << "No_Reaction::construct_from_msi_object: cannot read first query\n";
      return 0;
    }
    if (!_sidechain_no_reaction.construct_from_msi_object(*(msi[1]))) {
      cerr << "No_Reaction::construct_from_msi_object: cannot read 2nd query\n";
      return 0;
    }
  }

  // Note that we don't do any checking against the case of _scaffold_no_reaction
  // being specified twice, once as an expanded object and once as a smarts

  int na = msi.number_attributes();

  for (int i = 0; i < na; i++) {
    const msi_attribute* att = msi.attribute(i);

    if (NAME_OF_REACTION_COMMENT == att->name()) {
      _comment = att->stringval();
      continue;
    }

    int j = 0;  // processing scaffold or sidechain?
    if (_scaffold_no_reaction.empty()) {
      j = 0;
    } else if (_sidechain_no_reaction.empty()) {
      j = 1;
    }

    int success = 1;

    if (NAME_OF_QUERY_SMARTS_ATTRIBUTE == att->name()) {
      if (0 == j) {
        success = _scaffold_no_reaction.create_from_smarts(att->stringval());
      } else if (1 == j) {
        success = _sidechain_no_reaction.create_from_smarts(att->stringval());
      }
    } else if (NAME_OF_QUERY_FILE_ATTRIBUTE == att->name()) {
      if (0 == j) {
        success = _scaffold_no_reaction.read(att->stringval());
      } else if (1 == j) {
        success = _sidechain_no_reaction.read(att->stringval());
      }
    } else {
      cerr << "No_Reaction::construct_from_msi_object: unrecognised attribute '"
           << att->name() << "'\n";
      return 0;
    }

    if (0 == success) {
      cerr << "No_Reaction::construct_from_msi_object: cannot parse '" << (*att) << "'\n";
      return 0;
    }
  }

  if (_scaffold_no_reaction.empty() || _sidechain_no_reaction.empty()) {
    cerr << "No_Reaction::construct_from_msi_object: one or both queries not specified\n";
    return 0;
  }

  return 1;
}

int
No_Reaction::write_msi(std::ostream& os, int& object_id, int indentation) {
  IWString ind;
  ind.extend(indentation, ' ');

  os << ind << '(' << object_id++ << ' ' << NAME_OF_NO_REACTION_OBJECT << '\n';

  if (_comment.length()) {
    os << ind << "  (A C " << NAME_OF_REACTION_COMMENT << " \"" << _comment << "\")\n";
  }

  _scaffold_no_reaction.write_msi(os, object_id, indentation + 2);
  _sidechain_no_reaction.write_msi(os, object_id, indentation + 2);

  os << ind << ")\n";

  return os.good();
}

int
Reaction_Stereo_Centre::construct_from_msi_attribute(const msi_attribute* msi) {
  const_IWSubstring m;
  msi->value(m);

  int nw = m.nwords();

  if (nw < 4 || nw > 5) {
    cerr << "Reaction_Stereo_Centre::construct_from_msi_attribute: attribute must have 4 "
            "or 5 tokens\n";
    cerr << m << '\n';
    return 0;
  }

  for (int i = 0; i < nw; i++) {
    const_IWSubstring token;
    (void)m.word(i, token);

    if (!_ssc[i].construct(token)) {
      cerr << "Reaction_Stereo_Centre::construct_from_msi_attribute: cannot parse "
              "attribute "
           << i << " '" << token << "'\n";
      cerr << m << '\n';
      return 0;
    }

    //  cerr << "Built stereo component '" << _ssc[i] << "' from '" << token << "'\n";
  }

  if (4 == nw) {
    _ssc[4].set_whatever();
  }

  return 1;
}

int
Reaction_Stereo_Centre::write_msi(std::ostream& os, const const_IWSubstring& ind,
                                  const const_IWSubstring& attribute_name) const {
  os << ind << "  (A C " << attribute_name << " \"";
  for (int i = 0; i < 5; i++) {
    if (i > 0) {
      os << ' ';
    }
    os << _ssc[i];
  }

  os << "\")\n";

  return os.good();
}

/*
  A wedge bond must be of the form

    (A I wedge_bond a1 a2 dir)

  where dir is +1 or -1
*/

int
Reaction_Wedge_Bond::construct_from_msi_attribute(const msi_attribute* att) {
  if (3 != att->number_int_values()) {
    cerr << "Reaction_Wedge_Bond::construct_from_msi_attribute: attribute must have 3 "
            "integer values\n";
    return 0;
  }

  int a = att->int_multi_value(0);
  Pair_of_Atoms::set_a1(a);

  a = att->int_multi_value(1);
  Pair_of_Atoms::set_a2(a);

  if (Pair_of_Atoms::a1() == Pair_of_Atoms::a2()) {
    cerr << "Reaction_Wedge_Bond::construct_from_msi_attribute: atoms the same - "
            "impossible\n";
    return 0;
  }

  int d = att->int_multi_value(2);

  if (1 == d) {
    _direction = 1;
  } else if (-1 == d) {
    _direction = -1;
  } else if (0 == d) {
    _direction = 0;
  } else {
    cerr << "Reaction_Wedge_Bond::construct_from_msi_attribute: invalid direction " << d
         << '\n';
    return 0;
  }

  return 1;
}

int
Reaction_Wedge_Bond::write_msi(std::ostream& os, const IWString& ind,
                               const const_IWSubstring& attribute_name) const {
  os << ind << "  (A I " << attribute_name << " (" << _a1 << ' ' << _a2 << ' '
     << _direction << "))\n";

  return os.good();
}
