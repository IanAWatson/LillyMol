#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

void
BindDescriptors(nb::module_& m) {
  nb::enum_<Hybridization>(m, "Hybridization")
      .value("UNSPECIFIED", Hybridization::kUnspecified)
      .value("S", Hybridization::kS)
      .value("SP", Hybridization::kSp)
      .value("SP2", Hybridization::kSp2)
      .value("SP3", Hybridization::kSp3)
      .value("SP2D", Hybridization::kSp2d)
      .value("SP3D", Hybridization::kSp3d)
      .value("SP3D2", Hybridization::kSp3d2)
      .value("OTHER", Hybridization::kOther);

  nb::enum_<quick_rotbond::QuickRotatableBonds::RotBond>(m, "RotBond")
      .value("UNDEFINED", quick_rotbond::QuickRotatableBonds::RotBond::kUndefined)
      .value("QUICK", quick_rotbond::QuickRotatableBonds::RotBond::kQuick)
      .value("EXPENSIVE", quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);
  m.attr("UNDEFINED") = nb::cast(quick_rotbond::QuickRotatableBonds::RotBond::kUndefined);
  m.attr("QUICK") = nb::cast(quick_rotbond::QuickRotatableBonds::RotBond::kQuick);
  m.attr("EXPENSIVE") = nb::cast(quick_rotbond::QuickRotatableBonds::RotBond::kExpensive);

  nb::class_<alogp::ALogP>(m, "ALogP")
      .def(nb::init<>())
      .def("set_rdkit_phoshoric_acid_hydrogen",
           &alogp::ALogP::set_rdkit_phoshoric_acid_hydrogen,
           "Mimic RDKit handling of hydrogens on phosphoric acids")
      .def("set_use_alcohol_for_acid", &alogp::ALogP::set_use_alcohol_for_acid,
           "Mimic RDKit handling of oxygen atoms in acids")
      .def("logp",
           [](alogp::ALogP& calc, Molecule& mol) { return calc.LogP(mol); },
           nb::arg("mol"), "Compute AlogP");

  nb::class_<quick_rotbond::QuickRotatableBonds>(m, "RotatableBonds")
      .def(nb::init<>())
      .def("rotatable_bonds",
           [](quick_rotbond::QuickRotatableBonds& rotbond, Molecule& mol) {
             return rotbond.Process(mol, nullptr);
           },
           nb::arg("mol"), "Return number of rotatable bonds")
      .def("rotatable_bond_atoms",
           [](quick_rotbond::QuickRotatableBonds& rotbond, Molecule& mol) {
             std::vector<int> bond_rotatable(mol.nedges(), 0);
             rotbond.Process(mol, bond_rotatable.data());
             std::vector<std::tuple<atom_number_t, atom_number_t>> result;
             result.reserve(mol.nedges());
             const int nedges = mol.nedges();
             for (int bond_number = 0; bond_number < nedges; ++bond_number) {
               if (!bond_rotatable[bond_number]) {
                 continue;
               }
               const Bond* bond = mol.bondi(bond_number);
               result.emplace_back(bond->a1(), bond->a2());
             }
             return result;
           },
           nb::arg("mol"), "Return atom pairs defining rotatable bonds")
      .def("set_calculation_type", &quick_rotbond::QuickRotatableBonds::set_calculation_type,
           nb::arg("calculation_type"));

  m.def("MolFromSmiles", &MolFromSmiles, nb::arg("smiles"),
        "Build a Molecule from SMILES, returning None on parse failure");
  m.def("LillyMolFromSmiles", &MolFromSmiles, nb::arg("smiles"),
        "Build a Molecule from SMILES, returning None on parse failure");
  m.def("MolFromSmiles",
        [](const std::vector<std::string>& smiles) {
          std::vector<Molecule> result(smiles.size());
          for (uint32_t i = 0; i < smiles.size(); ++i) {
            if (!result[i].build_from_smiles(smiles[i])) {
              result[i].resize(0);
            }
          }
          return result;
        },
        nb::arg("smiles"),
        "Build molecules from a list of SMILES strings; invalid entries are empty molecules");
  m.def("set_auto_create_new_elements", &set_auto_create_new_elements,
        nb::arg("value"), "Allow arbitrary two-letter elements");
  m.def("set_atomic_symbols_can_have_arbitrary_length",
        &set_atomic_symbols_can_have_arbitrary_length, nb::arg("value"),
        "Allow atomic symbols with arbitrary length");
  m.def("interpret_D_as_deuterium", &element::interpret_d_as_deuterium,
        "Return whether D is interpreted as deuterium");
  m.def("interpret_T_as_deuterium", &element::interpret_t_as_tritium,
        "Return whether T is interpreted as tritium");
  m.def("set_display_strange_chemistry_messages", &set_display_strange_chemistry_messages,
        nb::arg("value"), "Control strange chemistry messages");
  m.def("set_display_smiles_interpretation_error_messages",
        &set_display_smiles_interpretation_error_messages, nb::arg("value"),
        "Control SMILES interpretation error messages");
  m.def("count_atoms_in_smiles",
        [](const std::string& smiles) {
          const const_IWSubstring tmp(smiles);
          return lillymol::count_atoms_in_smiles(tmp);
        },
        nb::arg("smiles"));
  m.def("hybridization_name",
        [](Hybridization hybridization) { return std::string(ToString(hybridization)); },
        nb::arg("hybridization"));
  m.def("hybridization",
        [](Molecule& mol, atom_number_t atom) {
          if (!mol.ok_atom_number(atom)) {
            throw std::invalid_argument("hybridization atom number outside [0, natoms)");
          }
          return HybridizationState(mol, atom);
        },
        nb::arg("mol"), nb::arg("atom"),
        "RDKit-like hybridization of atom, computed on demand");

  m.def("QueryFromSmarts",
        [](const std::string& smarts) -> std::unique_ptr<Substructure_Query> {
          auto result = std::make_unique<Substructure_Query>();
          if (!result->CreateFromSmarts(smarts)) {
            return nullptr;
          }
          return result;
        },
        nb::arg("smarts"),
        "Return a SubstructureQuery built from smarts, or None if parsing fails");
  m.def("HasSubstructMatch",
        [](Molecule& molecule, const std::string& smarts,
           std::optional<int> max_matches_to_find,
           std::optional<bool> unique_embeddings_only,
           std::optional<bool> one_embedding_per_start_atom,
           std::optional<bool> perceive_symmetry_equivalent_matches) {
          std::unique_ptr<Substructure_Query> query = BuildQueryFromSmarts(
              smarts, MakeSmartsSearchOptions(max_matches_to_find, unique_embeddings_only,
                                              one_embedding_per_start_atom,
                                              perceive_symmetry_equivalent_matches));
          nb::gil_scoped_release release;
          return static_cast<bool>(query->substructure_search(&molecule));
        },
        nb::arg("molecule"), nb::arg("smarts"), nb::kw_only(),
        nb::arg("max_matches_to_find") = nb::none(),
        nb::arg("unique_embeddings_only") = nb::none(),
        nb::arg("one_embedding_per_start_atom") = nb::none(),
        nb::arg("perceive_symmetry_equivalent_matches") = nb::none(),
        "Return True if smarts matches molecule");
  m.def("CountSubstructMatches",
        [](Molecule& molecule, const std::string& smarts,
           std::optional<int> max_matches_to_find,
           std::optional<bool> unique_embeddings_only,
           std::optional<bool> one_embedding_per_start_atom,
           std::optional<bool> perceive_symmetry_equivalent_matches) {
          std::unique_ptr<Substructure_Query> query = BuildQueryFromSmarts(
              smarts, MakeSmartsSearchOptions(max_matches_to_find, unique_embeddings_only,
                                              one_embedding_per_start_atom,
                                              perceive_symmetry_equivalent_matches));
          nb::gil_scoped_release release;
          return query->substructure_search(&molecule);
        },
        nb::arg("molecule"), nb::arg("smarts"), nb::kw_only(),
        nb::arg("max_matches_to_find") = nb::none(),
        nb::arg("unique_embeddings_only") = nb::none(),
        nb::arg("one_embedding_per_start_atom") = nb::none(),
        nb::arg("perceive_symmetry_equivalent_matches") = nb::none(),
        "Return number of embeddings for smarts in molecule");
  m.def("GetSubstructMatches",
        [](Molecule& molecule, const std::string& smarts,
           std::optional<int> max_matches_to_find,
           std::optional<bool> unique_embeddings_only,
           std::optional<bool> one_embedding_per_start_atom,
           std::optional<bool> perceive_symmetry_equivalent_matches) {
          std::unique_ptr<Substructure_Query> query = BuildQueryFromSmarts(
              smarts, MakeSmartsSearchOptions(max_matches_to_find, unique_embeddings_only,
                                              one_embedding_per_start_atom,
                                              perceive_symmetry_equivalent_matches));
          return SubstructureSearchMatches(*query, molecule);
        },
        nb::arg("molecule"), nb::arg("smarts"), nb::kw_only(),
        nb::arg("max_matches_to_find") = nb::none(),
        nb::arg("unique_embeddings_only") = nb::none(),
        nb::arg("one_embedding_per_start_atom") = nb::none(),
        nb::arg("perceive_symmetry_equivalent_matches") = nb::none(),
        "Return embeddings as lists of atom numbers, in query atom order");
  m.def("molecular_weight",
        [](const Molecule& mol) { return lillymol::MolecularWeightIsotopesAsLabels(mol); },
        nb::arg("mol"));
  m.def("HbaHbd", &LipinskiHbaHbd, nb::arg("mol"),
        "Return Lipinski HBA/HBD counts as (hba, hbd)");
  m.def("NumHAcceptors", [](const Molecule& mol) { return mol.LipinskiNumHAcceptors(); },
        nb::arg("mol"), "Lipinski hydrogen bond acceptor count");
  m.def("NumHDonors", [](Molecule& mol) { return mol.LipinskiNumHDonors(); },
        nb::arg("mol"), "Lipinski hydrogen bond donor count");
  m.def("RDKitNumHAcceptors", [](Molecule& mol) { return mol.RDKitNumHAcceptors(); },
        nb::arg("mol"), "RDKit-compatible hydrogen bond acceptor count");
  m.def("RDKitNumHDonors", [](Molecule& mol) { return mol.RDKitNumHDonors(); },
        nb::arg("mol"), "RDKit-compatible hydrogen bond donor count");
  m.def("fraction_csp3",
        [](Molecule& mol) {
          int carbon = 0;
          int csp3 = 0;
          const int matoms = mol.natoms();
          for (int i = 0; i < matoms; ++i) {
            if (mol.atomic_number(i) != 6) {
              continue;
            }
            ++carbon;
            if (mol.saturated(i)) {
              ++csp3;
            }
          }
          if (carbon == 0) {
            return 0.0;
          }
          return static_cast<double>(csp3) / static_cast<double>(carbon);
        },
        nb::arg("mol"),
        "Fraction of carbon atoms that are fully saturated");
  m.def("alogp", &ALogP, nb::arg("mol"));
  m.def("xlogp", &XLogPValue, nb::arg("mol"));
  m.def("tpsa",
        [](Molecule& mol) -> std::optional<double> {
          nvrtspsa::NovartisPolarSurfaceArea calc;
          return calc.PolarSurfaceArea(mol);
        },
        nb::arg("mol"));
}

}  // namespace lillymol_nb
