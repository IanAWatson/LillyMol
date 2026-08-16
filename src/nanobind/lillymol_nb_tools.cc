// Tool-style bindings for the nanobind LillyMol pilot.

#include "nanobind/lillymol_nb_internal.h"

#include "Molecule_Lib/mol2graph.h"
#include "Molecule_Tools/unique_molecules_api.h"

namespace lillymol_nb {

void
BindTools(nb::module_& m) {
  nb::class_<Mol2Graph>(m, "Mol2Graph")
      .def(nb::init<>())
      .def("set_exclude_triple_bonds_from_graph_reduction",
           &Mol2Graph::set_exclude_triple_bonds_from_graph_reduction,
           nb::arg("value"))
      .def("set_revert_all_directional_bonds_to_non_directional",
           &Mol2Graph::set_revert_all_directional_bonds_to_non_directional,
           nb::arg("value"))
      .def("set_preserve_cc_double_bonds_no_heteroatoms",
           &Mol2Graph::set_preserve_cc_double_bonds_no_heteroatoms,
           nb::arg("value"))
      .def("set_preserve_cc_double_bonds_saturated",
           &Mol2Graph::set_preserve_cc_double_bonds_saturated,
           nb::arg("value"))
      .def("set_append_molecular_formula", &Mol2Graph::set_append_molecular_formula,
           nb::arg("value"))
      .def("set_aromatic_distinguishing_formula",
           &Mol2Graph::set_aromatic_distinguishing_formula, nb::arg("value"))
      .def("set_remove_chiral_centres", &Mol2Graph::set_remove_chiral_centres,
           nb::arg("value"))
      .def("turn_on_most_useful_options", &Mol2Graph::TurnOnMostUsefulOptions)
      .def("set_active", &Mol2Graph::set_active, nb::arg("value"))
      .def("active", &Mol2Graph::active);

  nb::class_<unique_molecules::UniqueMolecules>(m, "UniqueMolecules")
      .def(nb::init<>())
      .def("set_include_chiral_info",
           &unique_molecules::UniqueMolecules::set_include_chiral_info,
           nb::arg("value"), "Control whether chirality is considered")
      .def("set_include_cis_trans_bonding_info",
           &unique_molecules::UniqueMolecules::set_include_cis_trans_bonding_info,
           nb::arg("value"), "Control whether cis/trans bonds are considered")
      .def("set_strip_to_largest_fragment",
           &unique_molecules::UniqueMolecules::set_strip_to_largest_fragment,
           nb::arg("value"), "Strip to largest fragment before comparison")
      .def("set_consider_isotopes",
           &unique_molecules::UniqueMolecules::set_consider_isotopes,
           nb::arg("value"), "Control whether isotopes are considered")
      .def("set_constant_isotope",
           &unique_molecules::UniqueMolecules::set_constant_isotope,
           nb::arg("value"), "Convert all isotope labels to a constant before comparison")
      .def("set_standardize_molecules",
           &unique_molecules::UniqueMolecules::set_standardize_molecules,
           nb::arg("value"), "Control whether molecules are standardised")
      .def("element_transformations",
           [](unique_molecules::UniqueMolecules& unique) -> Element_Transformations& {
             return unique.element_transformations();
           },
           nb::rv_policy::reference_internal, "Element transformations")
      .def("graph_specifications",
           [](unique_molecules::UniqueMolecules& unique) -> Mol2Graph& {
             return unique.graph_specifications();
           },
           nb::rv_policy::reference_internal, "Tautomer matching graph specifications")
      .def("add_to_hash", &unique_molecules::UniqueMolecules::AddToHashes,
           nb::arg("mol"), "Add molecule to internal uniqueness hashes")
      .def("is_unique", &unique_molecules::UniqueMolecules::IsUnique, nb::arg("mol"),
           "Return True if the molecule has not been seen before")
      .def("report", [](const unique_molecules::UniqueMolecules& unique) {
        return unique.Report(std::cerr);
      }, "Write a report to stderr");
}

}  // namespace lillymol_nb
