#include "nanobind/lillymol_nb_internal.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Tools/jwcats_lib.h"

namespace lillymol_nb {
namespace {

void
BuildAtomTyping(Atom_Typing_Specification& atom_typing, const std::string& specification) {
  const const_IWSubstring tmp(specification);
  if (!atom_typing.build(tmp)) {
    throw std::runtime_error("AtomTypingSpecification:invalid atom type '" + specification + "'");
  }
}

std::vector<uint32_t>
AssignAtomTypes(Atom_Typing_Specification& atom_typing, Molecule& mol) {
  std::vector<uint32_t> result(mol.natoms());
  if (!atom_typing.assign_atom_types(mol, result.data())) {
    throw std::runtime_error("AtomTypingSpecification:cannot assign atom types");
  }

  return result;
}

std::vector<uint32_t>
AssignAtomTypesFromSpecification(Molecule& mol, const std::string& specification) {
  Atom_Typing_Specification atom_typing;
  BuildAtomTyping(atom_typing, specification);
  return AssignAtomTypes(atom_typing, mol);
}

std::string
AtomTypingStringRepresentation(const Atom_Typing_Specification& atom_typing) {
  IWString result;
  if (!atom_typing.string_representation(result)) {
    throw std::runtime_error("AtomTypingSpecification:cannot form string representation");
  }
  return result.AsString();
}

std::string
AtomTypingTag(const Atom_Typing_Specification& atom_typing, const std::string& stem) {
  IWString result(stem);
  if (!atom_typing.append_to_tag(result)) {
    throw std::runtime_error("AtomTypingSpecification:cannot append atom type tag");
  }
  return result.AsString();
}

void
ThrowForJWCatsStatus(jwcats::ComputeStatus status) {
  switch (status) {
    case jwcats::ComputeStatus::kOk:
      return;
    case jwcats::ComputeStatus::kMissingChargeData:
      throw std::runtime_error("JWCats calculation failed: missing charge data");
    case jwcats::ComputeStatus::kNotInitialised:
      throw std::runtime_error("JWCats calculation failed: object is not initialised");
    case jwcats::ComputeStatus::kError:
      throw std::runtime_error("JWCats calculation failed");
  }

  throw std::runtime_error("JWCats calculation failed: unknown status");
}

bool
BuildJWCatsAssigners(jwcats::JWCats& jwcats, const IWString& charges_dir,
                     const IWString& hbonds_dir) {
  static constexpr int kVerbose = 0;
  return jwcats.charge_assigner().BuildFromDir(charges_dir) &&
         jwcats.donor_acceptor_assigner().BuildFromDir(hbonds_dir, kVerbose);
}

void
BuildJWCatsDefaultAssigners(jwcats::JWCats& jwcats) {
  IWString default_dir("DEF");
  if (!BuildJWCatsAssigners(jwcats, default_dir, default_dir)) {
    throw std::runtime_error(
        "Cannot initialise JWCats assigners from C3TK_DATA_PERSISTENT or LILLYMOL_HOME");
  }
}

void
InitialiseJWCats(jwcats::JWCats& jwcats) {
  if (!jwcats.Initialise()) {
    throw std::runtime_error("Cannot initialise JWCats");
  }
}

std::vector<double>
JWCatsResultToVector(const jwcats::JWCats& jwcats, const jwcats::Result& result) {
  const std::vector<int>& write_array_value = jwcats.write_array_value();
  std::vector<double> values;
  values.reserve(jwcats.FeatureNames().size());

  for (int i = 0; i < static_cast<int>(write_array_value.size()); ++i) {
    if (write_array_value[i]) {
      values.push_back(result.scaled_counts[i]);
    }
  }

  return values;
}

std::vector<double>
ComputeJWCats(jwcats::JWCats& jwcats, Molecule& mol) {
  jwcats::Result result;
  const jwcats::ComputeStatus status = jwcats.Compute(mol, result);
  ThrowForJWCatsStatus(status);
  return JWCatsResultToVector(jwcats, result);
}

}  // namespace

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

  nb::class_<Atom_Typing_Specification>(m, "AtomTypingSpecification")
      .def("__init__",
           [](Atom_Typing_Specification* atom_typing) {
             new (atom_typing) Atom_Typing_Specification();
           })
      .def("__init__",
           [](Atom_Typing_Specification* atom_typing, const std::string& specification) {
             new (atom_typing) Atom_Typing_Specification();
             BuildAtomTyping(*atom_typing, specification);
           },
           nb::arg("specification"))
      .def("build",
           [](Atom_Typing_Specification& atom_typing, const std::string& specification) {
             BuildAtomTyping(atom_typing, specification);
             return true;
           },
           nb::arg("specification"), "Build from an atom typing specification string")
      .def("active",
           [](const Atom_Typing_Specification& atom_typing) {
             return static_cast<bool>(atom_typing.active());
           },
           "True if an atom typing specification has been configured")
      .def("atom_type", &Atom_Typing_Specification::atom_type,
           "Return the configured atom typing mode")
      .def("user_specified_type", &Atom_Typing_Specification::user_specified_type,
           "Return the UST component bit mask")
      .def("string_representation", &AtomTypingStringRepresentation,
           "Return the atom typing string representation")
      .def("append_to_tag", &AtomTypingTag, nb::arg("stem"),
           "Append the atom typing suffix to a tag stem")
      .def("assign_atom_types", &AssignAtomTypes, nb::arg("mol"),
           "Assign atom types and return one integer per atom");

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

  nb::class_<jwcats::JWCats>(m, "JWCats")
      .def("__init__",
           [](jwcats::JWCats* jwcats, bool initialise_default_assigners) {
             new (jwcats) jwcats::JWCats();
             if (initialise_default_assigners) {
               BuildJWCatsDefaultAssigners(*jwcats);
             }
             InitialiseJWCats(*jwcats);
           },
           nb::arg("initialise_default_assigners") = true,
           "Create a JWCats descriptor calculator")
      .def("initialise",
           [](jwcats::JWCats& jwcats) {
             InitialiseJWCats(jwcats);
             return true;
           },
           "Initialise after changing settings")
      .def("build_default_assigners",
           [](jwcats::JWCats& jwcats) {
             BuildJWCatsDefaultAssigners(jwcats);
             return true;
           },
           "Initialise charge and donor/acceptor assigners from default envs")
      .def("build_assigners",
           [](jwcats::JWCats& jwcats, const std::string& charges_dir,
              const std::string& hbonds_dir) {
             if (!BuildJWCatsAssigners(jwcats, IWString(charges_dir), IWString(hbonds_dir))) {
               throw std::runtime_error("Cannot initialise JWCats assigners");
             }
             return true;
           },
           nb::arg("charges_dir"), nb::arg("hbonds_dir"),
           "Initialise assigners from explicit charge and hbonds query directories")
      .def("feature_names", &jwcats::JWCats::FeatureNames,
           "Return descriptor names in the same order as process() values")
      .def("process", &ComputeJWCats, nb::arg("mol"),
           "Compute JWCats descriptors for one molecule as a list of floats")
      .def("set_include_hydrophobic_pairs", &jwcats::JWCats::SetIncludeHydrophobicPairs,
           nb::arg("value"), "Control whether hydrophobe-hydrophobe pair columns are emitted")
      .def("set_min_bond_separation", &jwcats::JWCats::SetMinBondSeparation,
           nb::arg("value"), "Set the minimum bond separation")
      .def("set_max_bond_separation", &jwcats::JWCats::SetMaxBondSeparation,
           nb::arg("value"), "Set the maximum bond separation")
      .def("set_scaling_type", &jwcats::JWCats::SetScalingType, nb::arg("value"),
           "Set scaling type: 0 none, 1 heavy atoms, 2 feature counts, 3 heavy atoms / feature counts")
      .def("set_make_implicit_hydrogens_explicit",
           &jwcats::JWCats::SetMakeImplicitHydrogensExplicit, nb::arg("value"),
           "Make implicit hydrogens explicit during calculation")
      .def("initialised", &jwcats::JWCats::initialised,
           "Return whether the object has been initialised");

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

  m.def("assign_atom_types", &AssignAtomTypesFromSpecification, nb::arg("mol"),
        nb::arg("specification"),
        "Assign atom types from a transient specification and return one integer per atom");

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
