#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

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
      .def("AddHs", [](Molecule& mol) { return mol.make_implicit_hydrogens_explicit(); })
      .def("RemoveHs", [](Molecule& mol) { return mol.remove_all(1); })
      .def("valence_ok", nb::overload_cast<>(&Molecule::valence_ok))
      .def("valence_ok", nb::overload_cast<atom_number_t>(&Molecule::valence_ok), nb::arg("atom"))
      .def("smiles",
           [](Molecule& mol) { return mol.smiles().AsString(); },
           "Return a non-unique SMILES")
      .def("aromatic_smiles",
           [](Molecule& mol) { return mol.aromatic_smiles().AsString(); },
           "Return an aromatic SMILES")
      .def("unique_smiles",
           [](Molecule& mol) { return mol.unique_smiles().AsString(); },
           "Return the unique SMILES")
      .def("name",
           [](const Molecule& mol) { return mol.name().AsString(); },
           "Return the molecule name")
      .def("set_name",
           [](Molecule& mol, const std::string& name) { mol.set_name(name); },
           nb::arg("name"))
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
