#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

namespace {

bool
LooksLikeSdfTagRecord(const std::string& line) {
  const size_t open = line.find('<');
  const size_t close = line.rfind('>');
  return !line.empty() && line[0] == '>' && open != std::string::npos &&
         close != std::string::npos && open < close;
}

std::string
SdfTagName(const std::string& line) {
  const size_t open = line.find('<');
  const size_t close = line.rfind('>');
  if (open == std::string::npos || close == std::string::npos || open >= close) {
    return std::string();
  }

  std::string result = line.substr(open + 1, close - open - 1);
  std::replace(result.begin(), result.end(), ' ', '_');
  return result;
}

nb::dict
SdfTags(const Molecule& mol) {
  nb::dict result;

  std::string current_tag;
  std::string current_value;

  auto flush_current_tag = [&]() {
    if (current_tag.empty()) {
      return;
    }
    result[nb::str(current_tag.c_str())] = nb::str(current_value.c_str());
    current_tag.clear();
    current_value.clear();
  };

  for (int i = 0; i < mol.number_records_text_info(); ++i) {
    const std::string line = mol.text_info(i).AsString();
    if (LooksLikeSdfTagRecord(line)) {
      flush_current_tag();
      current_tag = SdfTagName(line);
      continue;
    }

    if (current_tag.empty()) {
      continue;
    }

    if (line.empty()) {
      flush_current_tag();
      continue;
    }

    if (!current_value.empty()) {
      current_value.append("\n");
    }
    current_value.append(line);
  }

  flush_current_tag();
  return result;
}

}  // namespace

void
BindMolecule(nb::module_& m) {
  nb::class_<Molecule>(m, "Molecule")
      .def(nb::init<>())
      .def("build_from_smiles",
           [](Molecule& mol, const std::string& smiles) {
             return mol.build_from_smiles(smiles);
           },
           nb::arg("smiles"))
      .def("ok", [](const Molecule& mol) { return static_cast<bool>(mol.ok()); })
      .def("natoms", nb::overload_cast<>(&Molecule::natoms, nb::const_),
           "Return the number of atoms")
      .def("natoms",
           [](const Molecule& mol, atomic_number_t atomic_number) {
             return mol.natoms(atomic_number);
           },
           nb::arg("atomic_number"))
      .def("natoms",
           [](const Molecule& mol, const std::string& symbol) {
             const_IWSubstring tmp(symbol);
             const Element* element = get_element_from_symbol_no_case_conversion(tmp);
             if (element == nullptr) {
               throw std::invalid_argument("Unrecognised element");
             }
             return mol.natoms(element);
           },
           nb::arg("symbol"))
      .def("GetNumAtoms", [](const Molecule& mol) { return mol.natoms(); })
      .def("empty", [](const Molecule& mol) { return static_cast<bool>(mol.empty()); })
      .def("resize", &Molecule::resize, nb::arg("natoms"))
      .def("add_atom",
           [](Molecule& mol, int atomic_number) {
             const Element* element = get_element_from_atomic_number(atomic_number);
             if (element == nullptr) {
               throw std::invalid_argument("Invalid atomic number");
             }
             mol.add(element);
             return mol.natoms() - 1;
           },
           nb::arg("atomic_number"))
      .def("atom", &Molecule::atomi, nb::arg("atom"), nb::rv_policy::reference_internal)
      .def("nedges", &Molecule::nedges,
           "Return the number of bonds")
      .def("bond", nb::overload_cast<int>(&Molecule::bondi, nb::const_), nb::arg("bond"),
           nb::rv_policy::reference_internal)
      .def("bonds", &Bonds, nb::rv_policy::reference_internal)
      .def("nrings", [](Molecule& mol) { return mol.nrings(); },
           "Return the number of SSSR rings")
      .def("get_ring_membership", &GetRingMembership,
           "Return SSSR ring membership counts for all atoms")
      .def("ring_membership", [](Molecule& mol) { mol.ring_membership(); },
           "Force ring membership perception")
      .def("is_ring_atom", [](Molecule& mol, atom_number_t atom) {
             return static_cast<bool>(mol.is_ring_atom(atom));
           },
           nb::arg("atom"))
      .def("IsInRing", [](Molecule& mol, atom_number_t atom) {
             return mol.ring_bond_count(atom) > 0;
           },
           nb::arg("atom"))
      .def("in_ring_of_size", &Molecule::in_ring_of_given_size,
           nb::arg("atom"), nb::arg("ring_size"))
      .def("IsAtomInRingOfSize", &Molecule::in_ring_of_given_size,
           nb::arg("atom"), nb::arg("ring_size"))
      .def("NumAtomRings", nb::overload_cast<atom_number_t>(&Molecule::nrings),
           nb::arg("atom"))
      .def("ring_bond_count", nb::overload_cast<atom_number_t>(&Molecule::ring_bond_count),
           nb::arg("atom"))
      .def("ring_bond_counts", &RingBondCount)
      .def("fused_system_identifier", &Molecule::fused_system_identifier,
           nb::arg("atom"))
      .def("fused_system_size", &Molecule::fused_system_size, nb::arg("atom"))
      .def("number_ring_systems", nb::overload_cast<>(&Molecule::number_ring_systems))
      .def("ring", nb::overload_cast<int>(&Molecule::ringi), nb::arg("ring"),
           nb::rv_policy::reference_internal)
      .def("rings", &Rings)
      .def("in_same_ring",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(mol.in_same_ring(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("in_same_ring_system",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) {
             return static_cast<bool>(mol.in_same_ring_system(a1, a2));
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("largest_ring_size", &Molecule::LargestRingSize)
      .def("number_fragments", [](Molecule& mol) { return mol.number_fragments(); },
           "Return the number of fragments")
      .def("atoms_in_largest_fragment", [](Molecule& mol) { return mol.atoms_in_largest_fragment(); },
           "Return the number of atoms in the largest fragment")
      .def("remove_non_periodic_table_elements", &Molecule::remove_all_non_natural_elements,
           "Remove non periodic table elements")
      .def("organic_only", nb::overload_cast<>(&Molecule::organic_only, nb::const_),
           "True if only organic elements are present")
      .def("non_organic_atom_count", &Molecule::non_organic_atom_count,
           "Return the number of non-organic atoms")
      .def("is_organic",
           [](const Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.is_organic(atom)); },
           nb::arg("atom"), "True if atom is organic")
      .def("remove_explicit_hydrogens", nb::overload_cast<>(&Molecule::remove_explicit_hydrogens),
           "Remove explicit hydrogens")
      .def("atomic_number", &Molecule::atomic_number, nb::arg("atom"),
           "Return the atomic number for atom")
      .def("set_atomic_number", &Molecule::set_atomic_number, nb::arg("atom"),
           nb::arg("atomic_number"), "Set the atomic number for atom")
      .def("add_bond",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2, BondType bond_type) {
             return mol.add_bond(a1, a2, ToLillyMolBondType(bond_type));
           },
           nb::arg("a1"), nb::arg("a2"), nb::arg("bond_type"), "Add bond")
      .def("set_bond_type_between_atoms",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2, BondType bond_type) {
             return mol.set_bond_type_between_atoms(a1, a2, ToLillyMolBondType(bond_type));
           },
           nb::arg("a1"), nb::arg("a2"), nb::arg("bond_type"), "Set bond type")
      .def("assign_bond_numbers_to_bonds", &Molecule::assign_bond_numbers_to_bonds)
      .def("bond_between_atoms",
           [](const Molecule& mol, atom_number_t a1, atom_number_t a2) -> const Bond* {
             return mol.bond_between_atoms_if_present(a1, a2);
           },
           nb::arg("a1"), nb::arg("a2"), nb::rv_policy::reference_internal)
      .def("bond_type_between_atoms",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<BondType> {
             mol.compute_aromaticity_if_needed();
             const Bond* bond = mol.bond_between_atoms_if_present(a1, a2);
             if (bond == nullptr) {
               return std::nullopt;
             }
             return ToBondType(*bond);
           },
           nb::arg("a1"), nb::arg("a2"))
      .def("remove_atom", nb::overload_cast<atom_number_t>(&Molecule::remove_atom),
           nb::arg("atom"))
      .def("remove_atoms",
           [](Molecule& mol, const std::vector<int>& to_remove, int flag) {
             return mol.remove_atoms(to_remove.data(), flag);
           },
           nb::arg("to_remove"), nb::arg("flag"))
      .def("remove_atoms",
           [](Molecule& mol, Set_of_Atoms atoms) {
             return mol.remove_atoms(atoms);
           },
           nb::arg("atoms"))
      .def("remove_bond_between_atoms", &Molecule::remove_bond_between_atoms,
           nb::arg("a1"), nb::arg("a2"))
      .def("remove_all_bonds", &Molecule::remove_all_bonds)
      .def("remove_bonds_to_atom",
           [](Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.remove_bonds_to_atom(atom, 0)); },
           nb::arg("atom"), "Remove all bonds to an atom")
      .def("remove_edge", nb::overload_cast<int>(&Molecule::remove_bond),
           nb::arg("bond"), "Remove a bond by bond number")
      .def("chop", &Molecule::chop, nb::arg("natoms"), "Remove the last atoms")
      .def("remove_all", nb::overload_cast<atomic_number_t>(&Molecule::remove_all),
           nb::arg("atomic_number"))
      .def("atomic_numbers", &AtomicNumbers,
           "Return atomic numbers for all atoms")
      .def("ncon", nb::overload_cast<atom_number_t>(&Molecule::ncon, nb::const_),
           nb::arg("atom"), "Return the connection count for atom")
      .def("connections",
           [](const Molecule& mol, atom_number_t atom) {
             Set_of_Atoms connections = mol.connections(atom);
             return std::vector<int>(connections.rawdata(), connections.rawdata() + connections.size());
           },
           nb::arg("atom"))
      .def("other_atom",
           [](const Molecule& mol, atom_number_t atom, int connection) {
             return mol.other(atom, connection);
           },
           nb::arg("atom"), nb::arg("connection"))
      .def("nbonds", nb::overload_cast<atom_number_t>(&Molecule::nbonds, nb::const_),
           nb::arg("atom"), "Return the bond order sum for atom")
      .def("atomic_symbol",
           [](const Molecule& mol, atom_number_t atom) { return mol.atomic_symbol(atom).AsString(); },
           nb::arg("atom"))
      .def("is_aromatic",
           [](Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.IsAromatic(atom)); },
           nb::arg("atom"))
      .def("aromatic_atom_count", nb::overload_cast<>(&Molecule::aromatic_atom_count))
      .def("aromatic_ring_count", nb::overload_cast<>(&Molecule::aromatic_ring_count))
      .def("number_chiral_centres", nb::overload_cast<>(&Molecule::chiral_centres, nb::const_))
      .def("remove_all_chiral_centres", &Molecule::remove_all_chiral_centres)
      .def("chiral_centre_at_atom",
           [](const Molecule& mol, atom_number_t atom) -> std::optional<Chiral_Centre> {
             Chiral_Centre* centre = mol.chiral_centre_at_atom(atom);
             if (centre == nullptr) {
               return std::nullopt;
             }
             return *centre;
           },
           nb::arg("atom"))
      .def("invert_chirality_on_atom", &Molecule::invert_chirality_on_atom,
           nb::arg("atom"))
      .def("remove_chiral_centre_at_atom", &Molecule::remove_chiral_centre_at_atom,
           nb::arg("atom"))
      .def("chiral_centres",
           [](const Molecule& mol) {
             std::vector<Chiral_Centre> result;
             result.reserve(mol.chiral_centres());
             for (const Chiral_Centre* centre : mol.ChiralCentres()) {
               result.push_back(*centre);
             }
             return result;
           })
      .def("isotope", &Molecule::isotope, nb::arg("atom"),
           "Return the isotope value for atom")
      .def("isotopes", &Isotopes,
           "Return isotope values for all atoms")
      .def("set_isotope",
           [](Molecule& mol, atom_number_t atom, isotope_t isotope) {
             return mol.set_isotope(atom, isotope);
           },
           nb::arg("atom"), nb::arg("isotope"))
      .def("remove_isotopes",
           [](Molecule& mol, int unset_implicit_h) {
             return mol.unset_isotopes(unset_implicit_h);
           },
           nb::arg("unset_implicit_h") = 1)
      .def("transform_to_non_isotopic_form",
           [](Molecule& mol, int unset_implicit_h) {
             return mol.transform_to_non_isotopic_form(unset_implicit_h);
           },
           nb::arg("unset_implicit_h") = 1)
      .def("number_isotopic_atoms", nb::overload_cast<>(&Molecule::number_isotopic_atoms, nb::const_))
      .def("first_atom_with_isotope",
           [](const Molecule& mol, isotope_t isotope) { return mol.atom_with_isotope(isotope); },
           nb::arg("isotope"))
      .def("formal_charge", nb::overload_cast<atom_number_t>(&Molecule::formal_charge, nb::const_),
           nb::arg("atom"))
      .def("set_formal_charge", nb::overload_cast<atom_number_t, formal_charge_t>(&Molecule::set_formal_charge),
           nb::arg("atom"), nb::arg("formal_charge"))
      .def("has_formal_charges", nb::overload_cast<>(&Molecule::has_formal_charges, nb::const_))
      .def("number_formal_charges", nb::overload_cast<>(&Molecule::number_formally_charged_atoms, nb::const_))
      .def("net_formal_charge", nb::overload_cast<>(&Molecule::net_formal_charge, nb::const_))
      .def("molecular_formula",
           [](Molecule& mol) { return mol.molecular_formula().AsString(); },
           "Return the molecular formula")
      .def("amw",
           [](const Molecule& mol) { return lillymol::MolecularWeightIsotopesAsLabels(mol); },
           "Return average molecular weight, treating isotope labels as labels")
      .def("amw_ignore_isotopes", &Molecule::molecular_weight_ignore_isotopes,
           "Return average molecular weight with isotope labels ignored strictly")
      .def("exact_mass", [](Molecule& mol) { return mol.exact_mass(); },
           "Return exact mass")
      .def("implicit_hydrogens", nb::overload_cast<atom_number_t>(&Molecule::implicit_hydrogens),
           nb::arg("atom"))
      .def("explicit_hydrogens", nb::overload_cast<atom_number_t>(&Molecule::explicit_hydrogens, nb::const_),
           nb::arg("atom"))
      .def("hcount",
           [](Molecule& mol, atom_number_t atom) { return mol.hcount(atom); },
           nb::arg("atom"))
      .def("implicit_hydrogens_known", &Molecule::implicit_hydrogens_known,
           nb::arg("atom"))
      .def("make_implicit_hydrogens_explicit", nb::overload_cast<>(&Molecule::make_implicit_hydrogens_explicit))
      .def("unset_all_implicit_hydrogen_information", &Molecule::unset_all_implicit_hydrogen_information,
           "Discard all implicit hydrogen known flags")
      .def("AddHs", [](Molecule& mol) { return mol.make_implicit_hydrogens_explicit(); })
      .def("RemoveHs", [](Molecule& mol) { return mol.remove_all(1); })
      .def("valence_ok", nb::overload_cast<>(&Molecule::valence_ok))
      .def("valence_ok", nb::overload_cast<atom_number_t>(&Molecule::valence_ok), nb::arg("atom"))
      .def("lipinski_num_h_donors", &Molecule::LipinskiNumHDonors,
           "Lipinski hydrogen bond donor count")
      .def("lipinski_num_h_acceptors", &Molecule::LipinskiNumHAcceptors,
           "Lipinski hydrogen bond acceptor count")
      .def("rdkit_num_h_donors", &Molecule::RDKitNumHDonors,
           "RDKit compatible hydrogen bond donor count")
      .def("rdkit_num_h_acceptors", &Molecule::RDKitNumHAcceptors,
           "RDKit compatible hydrogen bond acceptor count")
      .def("saturated",
           [](Molecule& mol, atom_number_t atom) { return static_cast<bool>(mol.saturated(atom)); },
           nb::arg("atom"))
      .def("unsaturation", &Molecule::unsaturation, nb::arg("atom"))
      .def("smiles",
           [](Molecule& mol) { return mol.smiles().AsString(); },
           "Return a non-unique SMILES")
      .def("aromatic_smiles",
           [](Molecule& mol) { return mol.aromatic_smiles().AsString(); },
           "Return an aromatic SMILES")
      .def("unique_smiles",
           [](Molecule& mol) { return mol.unique_smiles().AsString(); },
           "Return the unique SMILES")
      .def("random_smiles",
           [](Molecule& mol) { return mol.random_smiles().AsString(); },
           "Return a random SMILES")
      .def("isotopically_labelled_smiles",
           [](Molecule& mol) { return mol.isotopically_labelled_smiles().AsString(); },
           "Return SMILES with isotopes as atom numbers")
      .def("unique_kekule_smiles",
           [](Molecule& mol) { return mol.UniqueKekuleSmiles().AsString(); },
           "Return unique Kekule SMILES")
      .def("smiles_atom_order",
           [](Molecule& mol) {
             std::vector<int> result(mol.natoms());
             mol.smiles_atom_order(result.data());
             return result;
           },
           "Return atom order from the most recent SMILES generation")
      .def("renumber_atoms",
           [](Molecule& mol, const std::vector<int>& new_number) {
             const int matoms = mol.natoms();
             if (static_cast<int>(new_number.size()) != matoms) {
               throw std::invalid_argument("renumber_atoms requires one entry for each atom");
             }
             std::vector<int> seen(matoms, 0);
             for (int i = 0; i < matoms; ++i) {
               const int destination = new_number[i];
               if (destination < 0 || destination >= matoms) {
                 throw std::invalid_argument("renumber_atoms mapping contains an atom number outside [0, natoms)");
               }
               if (seen[destination]) {
                 throw std::invalid_argument("renumber_atoms mapping contains duplicate atom numbers");
               }
               seen[destination] = 1;
             }
             return mol.renumber_atoms(new_number.data());
           },
           nb::arg("new_number"), "Renumber atoms")
      .def("smiles_starting_with_atom",
           [](Molecule& mol, atom_number_t atom) { return mol.smiles_starting_with_atom(atom).AsString(); },
           nb::arg("atom"), "Return SMILES starting at atom")
      .def("name",
           [](const Molecule& mol) { return mol.name().AsString(); },
           "Return the molecule name")
      .def("set_name",
           [](Molecule& mol, const std::string& name) { mol.set_name(name); },
           nb::arg("name"))
      .def("number_records_text_info", &Molecule::number_records_text_info,
           "Return number of retained SDF/text records")
      .def("text_info",
           [](const Molecule& mol) {
             std::vector<std::string> result;
             result.reserve(mol.number_records_text_info());
             for (int i = 0; i < mol.number_records_text_info(); ++i) {
               result.push_back(mol.text_info(i).AsString());
             }
             return result;
           },
           "Return retained SDF/text records")
      .def("sdf_tags", &SdfTags, "Return retained SDF tags as a dict")
      .def("are_bonded", nb::overload_cast<atom_number_t, atom_number_t>(&Molecule::are_bonded, nb::const_),
           nb::arg("a1"), nb::arg("a2"), "True if atoms are bonded")
      .def("bonds_between", nb::overload_cast<atom_number_t, atom_number_t>(&Molecule::bonds_between),
           nb::arg("a1"), nb::arg("a2"), "Return topological distance between atoms")
      .def("longest_path", nb::overload_cast<>(&Molecule::longest_path),
           "Return the longest topological path")
      .def("most_distant_pair",
           [](Molecule& mol) {
             atom_number_t a1 = INVALID_ATOM_NUMBER;
             atom_number_t a2 = INVALID_ATOM_NUMBER;
             int longest_distance = 0;
             const int matoms = mol.natoms();
             for (int i = 0; i < matoms; ++i) {
               for (int j = i + 1; j < matoms; ++j) {
                 if (mol.fragment_membership(i) != mol.fragment_membership(j)) {
                   continue;
                 }
                 const int distance = mol.bonds_between(i, j);
                 if (distance > longest_distance) {
                   a1 = i;
                   a2 = j;
                   longest_distance = distance;
                 }
               }
             }
             return std::make_tuple(a1, a2);
           },
           "Return the most separated atom pair")
      .def("atoms_on_shortest_path",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<std::vector<int>> {
             Set_of_Atoms atoms;
             if (!mol.atoms_between(a1, a2, atoms) || atoms.empty()) {
               return std::nullopt;
             }
             return std::vector<int>(atoms.rawdata(), atoms.rawdata() + atoms.size());
           },
           nb::arg("a1"), nb::arg("a2"), "Return atoms on the shortest path between a1 and a2")
      .def("all_atoms_between",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<std::vector<int>> {
             Set_of_Atoms atoms;
             if (!mol.AllAtomsBetween(a1, a2, atoms) || atoms.empty()) {
               return std::nullopt;
             }
             return std::vector<int>(atoms.rawdata(), atoms.rawdata() + atoms.size());
           },
           nb::arg("a1"), nb::arg("a2"), "Return all atoms on shortest paths between a1 and a2")
      .def("down_the_bond",
           [](Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<std::vector<int>> {
             const int matoms = mol.natoms();
             std::unique_ptr<int[]> dtb = std::make_unique<int[]>(matoms);
             std::optional<int> maybe_n = mol.DownTheBond(a1, a2, dtb.get());
             if (!maybe_n) {
               return std::nullopt;
             }
             std::vector<int> result;
             result.reserve(*maybe_n);
             for (int i = 0; i < matoms; ++i) {
               if (i != a2 && dtb[i]) {
                 result.push_back(i);
               }
             }
             return result;
           },
           nb::arg("a1"), nb::arg("a2"), "Return atoms down the bond from a1 to a2")
      .def("reset_atom_map_numbers", &Molecule::reset_all_atom_map_numbers,
           "Reset atom map numbers")
      .def("set_atom_map_number",
           static_cast<void (Molecule::*)(atom_number_t, int)>(&Molecule::set_atom_map_number),
           nb::arg("atom"), nb::arg("map_number"), "Set atom map number")
      .def("atom_map_number",
           nb::overload_cast<atom_number_t>(&Molecule::atom_map_number, nb::const_),
           nb::arg("atom"), "Return atom map number")
      .def("atom_with_atom_map_number", &Molecule::atom_with_atom_map_number,
           nb::arg("map_number"), "Return atom with atom map number")
      .def("canonical_rank", nb::overload_cast<atom_number_t>(&Molecule::canonical_rank),
           nb::arg("atom"), "Return atom canonical rank")
      .def("canonical_ranks",
           [](Molecule& mol) {
             std::vector<int> result(mol.natoms());
             mol.canonical_ranks(result.data());
             return result;
           },
           "Return canonical ranks for all atoms")
      .def("symmetry_class", nb::overload_cast<atom_number_t>(&Molecule::symmetry_class),
           nb::arg("atom"), "Return atom symmetry class")
      .def("number_symmetry_classes", nb::overload_cast<>(&Molecule::number_symmetry_classes),
           "Return number of symmetry classes")
      .def("symmetry_equivalents",
           [](Molecule& mol, atom_number_t atom) {
             Set_of_Atoms tmp;
             mol.symmetry_equivalents(atom, tmp);
             return std::vector<int>(tmp.rawdata(), tmp.rawdata() + tmp.size());
           },
           nb::arg("atom"), "Return atoms symmetry-equivalent to atom")
      .def("partial_charge_type",
           [](const Molecule& mol) { return mol.partial_charge_type().AsString(); },
           "Return the type of partial charges stored")
      .def("invalidate_partial_charges",
           [](Molecule& mol) { mol.invalidate_charges(); },
           "Discard partial charge information")
      .def("partial_charge",
           [](const Molecule& mol, atom_number_t atom) { return mol.partial_charge(atom); },
           nb::arg("atom"), "Return partial charge on atom")
      .def("compute_Abraham_partial_charges",
           [](Molecule& mol) { return mol.compute_Abraham_partial_charges(); },
           "Compute Abraham partial charges")
      .def("compute_Gasteiger_partial_charges",
           [](Molecule& mol) { return mol.compute_Gasteiger_partial_charges(); },
           "Compute Gasteiger partial charges")
      .def("compute_Huckel_partial_charges",
           [](Molecule& mol) { return mol.compute_Huckel_partial_charges(); },
           "Compute Huckel partial charges")
      .def("compute_Gasteiger_Huckel_partial_charges",
           [](Molecule& mol) { return mol.compute_Gasteiger_Huckel_partial_charges(); },
           "Compute Gasteiger-Huckel partial charges")
      .def("gasteiger_partial_charges",
           [](Molecule& mol) -> std::vector<float> {
             mol.compute_Gasteiger_partial_charges();
             std::vector<float> result;
             result.reserve(mol.natoms());
             for (int i = 0; i < mol.natoms(); ++i) {
               result.push_back(mol.partial_charge(i));
             }
             return result;
           },
           "Return list of Gasteiger partial charges")
      .def("x", nb::overload_cast<atom_number_t>(&Molecule::x, nb::const_),
           nb::arg("atom"), "Return x coordinate")
      .def("y", nb::overload_cast<atom_number_t>(&Molecule::y, nb::const_),
           nb::arg("atom"), "Return y coordinate")
      .def("z", nb::overload_cast<atom_number_t>(&Molecule::z, nb::const_),
           nb::arg("atom"), "Return z coordinate")
      .def("setx", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::setx),
           nb::arg("atom"), nb::arg("x"), "Set x coordinate")
      .def("sety", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::sety),
           nb::arg("atom"), nb::arg("y"), "Set y coordinate")
      .def("setz", static_cast<void (Molecule::*)(atom_number_t, coord_t)>(&Molecule::setz),
           nb::arg("atom"), nb::arg("z"), "Set z coordinate")
      .def("setxyz",
           static_cast<void (Molecule::*)(atom_number_t, coord_t, coord_t, coord_t)>(&Molecule::setxyz),
           nb::arg("atom"), nb::arg("x"), nb::arg("y"), nb::arg("z"),
           "Set atom coordinates")
      .def("get_coordinates",
           [](const Molecule& mol) {
             const int matoms = mol.natoms();
             std::vector<float> result;
             result.reserve(matoms * 3);
             for (int i = 0; i < matoms; ++i) {
               result.push_back(mol.x(i));
               result.push_back(mol.y(i));
               result.push_back(mol.z(i));
             }
             return result;
           },
           "Return coordinates as [x0, y0, z0, x1, y1, z1, ...]")
      .def("set_coordinates",
           [](Molecule& mol, const std::vector<float>& coords) {
             const int expected = 3 * mol.natoms();
             if (static_cast<int>(coords.size()) != expected) {
               throw std::invalid_argument("set_coordinates requires 3 values per atom");
             }
             mol.SetCoordinates(coords.data());
           },
           nb::arg("coords"), "Set coordinates from [x0, y0, z0, x1, y1, z1, ...]")
      .def("bond_length",
           [](const Molecule& mol, atom_number_t a1, atom_number_t a2) -> std::optional<float> {
             if (!mol.are_bonded(a1, a2)) {
               return std::nullopt;
             }
             return mol.bond_length(a1, a2);
           },
           nb::arg("a1"), nb::arg("a2"), "Return distance between bonded atoms")
      .def("bond_angle",
           [](const Molecule& mol, atom_number_t centre, atom_number_t a1, atom_number_t a2) {
             return mol.bond_angle(centre, a1, a2, BondedStatus::kOkNotBonded);
           },
           nb::arg("centre"), nb::arg("a1"), nb::arg("a2"),
           "Return angle defined by three atoms")
      .def("dihedral_angle",
           [](const Molecule& mol, atom_number_t a1, atom_number_t a2,
              atom_number_t a3, atom_number_t a4) {
             return mol.dihedral_angle(a1, a2, a3, a4, BondedStatus::kOkNotBonded);
           },
           nb::arg("a1"), nb::arg("a2"), nb::arg("a3"), nb::arg("a4"),
           "Return dihedral angle defined by four atoms")
      .def("signed_dihedral_angle",
           nb::overload_cast<atom_number_t, atom_number_t, atom_number_t, atom_number_t>(
               &Molecule::signed_dihedral_angle, nb::const_),
           nb::arg("a1"), nb::arg("a2"), nb::arg("a3"), nb::arg("a4"),
           "Return signed dihedral angle defined by four atoms")
      .def("distance_between_atoms",
           nb::overload_cast<atom_number_t, atom_number_t>(
               &Molecule::distance_between_atoms, nb::const_),
           nb::arg("a1"), nb::arg("a2"), "Return spatial distance between atoms")
      .def("longest_intra_molecular_distance", &Molecule::longest_intra_molecular_distance,
           "Return longest spatial distance within the molecule")
      .def("bump_check",
           [](const Molecule& mol, distance_t dist) { return mol.bump_check(dist); },
           nb::arg("dist"), "Return number of non-bonded atom pairs closer than dist")
      .def("highest_coordinate_dimensionality", &Molecule::highest_coordinate_dimensionality,
           "Return highest coordinate dimensionality present")
      .def("discern_chirality_from_3d_structure", &Molecule::discern_chirality_from_3d_structure,
           "Perceive chiral centres from 3D coordinates")
      .def("debug_string", &Molecule::debug_string,
           "Return a dump of internal data structures")
      .def("__repr__",
           [](Molecule& mol) {
             IWString formula;
             mol.isis_like_molecular_formula(formula);
             IWString result;
             result << '<' << mol.name() << " with " << mol.natoms() << " atoms " << formula << '>';
             return result.AsString();
           })
      .def("__str__",
           [](Molecule& mol) {
             IWString result;
             result << mol.smiles() << ' ' << mol.name();
             return result.AsString();
           })
      .def("__iadd__",
           [](Molecule& mol, const Molecule& rhs) -> Molecule& {
             mol += rhs;
             return mol;
           },
           nb::rv_policy::reference_internal)
      .def("__add__",
           [](const Molecule& lhs, const Molecule& rhs) { return lhs + rhs; })
      .def("__len__", [](const Molecule& mol) { return mol.natoms(); })
      .def("__getitem__",
           [](const Molecule& mol, int index) { return mol[index]; },
           nb::rv_policy::reference_internal)
      .def("__iter__",
           [](const Molecule& mol) {
             return nb::make_iterator<nb::rv_policy::reference_internal>(
                 nb::type<Molecule>(), "AtomIterator", mol.begin(), mol.end());
           },
           nb::keep_alive<0, 1>())
      .def("__copy__", [](const Molecule& mol) { return Molecule(mol); })
      .def("__contains__", &MoleculeContainsAtomicNumber)
      .def("__contains__", &MoleculeContainsElementSymbol)
      .def("__contains__", &MoleculeContainsSubstructureQuery);


}

}  // namespace lillymol_nb
