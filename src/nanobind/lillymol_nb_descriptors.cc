#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

void
BindDescriptors(nb::module_& m) {
  m.def("MolFromSmiles", &MolFromSmiles, nb::arg("smiles"),
        "Build a Molecule from SMILES, returning None on parse failure");

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
