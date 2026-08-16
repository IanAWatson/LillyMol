// Tool-style bindings for the nanobind LillyMol pilot.

#include "nanobind/lillymol_nb_internal.h"

#include <nanobind/stl/unordered_map.h>

#include "Molecule_Lib/mol2graph.h"
#include "Molecule_Tools/dicer_api.h"
#include "Molecule_Tools/unique_molecules_api.h"
#include "Utilities/GFP_Tools/truncated_distance_matrix.h"

namespace lillymol_nb {
namespace {

void
CheckDistanceMatrixIndex(const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
                         uint32_t ndx) {
  if (ndx >= dm.size()) {
    throw std::out_of_range("TruncatedDistanceMatrix index out of range");
  }
}

std::string
TruncatedDistanceMatrixName(const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
                            uint32_t ndx) {
  CheckDistanceMatrixIndex(dm, ndx);
  return dm.Name(ndx);
}

std::optional<float>
TruncatedDistanceMatrixDistance(const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
                                uint32_t i, uint32_t j) {
  CheckDistanceMatrixIndex(dm, i);
  CheckDistanceMatrixIndex(dm, j);
  return dm.Distance(i, j);
}

float
TruncatedDistanceMatrixDistanceOrDefault(
    const truncated_distance_matrix::TruncatedDistanceMatrix& dm, uint32_t i, uint32_t j) {
  CheckDistanceMatrixIndex(dm, i);
  CheckDistanceMatrixIndex(dm, j);
  return dm.DistanceOrDefault(i, j);
}

std::vector<float>
TruncatedDistanceMatrixDistancesOrDefault(
    const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
    const std::vector<uint32_t>& i, const std::vector<uint32_t>& j) {
  if (i.size() != j.size()) {
    throw std::invalid_argument("index arrays must have the same length");
  }
  for (uint32_t ndx : i) {
    CheckDistanceMatrixIndex(dm, ndx);
  }
  for (uint32_t ndx : j) {
    CheckDistanceMatrixIndex(dm, ndx);
  }
  return dm.DistancesOrDefault(i, j);
}

void
SetTruncatedDistanceMatrixDefaultDistance(
    truncated_distance_matrix::TruncatedDistanceMatrix& dm, float distance) {
  if (!dm.SetDefaultDistance(distance)) {
    throw std::invalid_argument("default distance cannot be below max stored distance");
  }
}

void
SetTruncatedDistanceMatrixDefaultDistanceByte(
    truncated_distance_matrix::TruncatedDistanceMatrix& dm, uint8_t distance) {
  if (!dm.SetDefaultDistanceByte(distance)) {
    throw std::invalid_argument("default distance cannot be below max stored distance");
  }
}

std::unordered_map<std::string, uint32_t>
DiceMolecule(dicer_api::Dicer& dicer, Molecule& mol) {
  std::unordered_map<std::string, uint32_t> result;
  dicer.Dice(mol, result);
  return result;
}

}  // namespace

void
BindTools(nb::module_& m) {
  nb::enum_<truncated_distance_matrix::Storage>(m, "TruncatedDistanceMatrixStorage")
      .value("ROW_SPARSE", truncated_distance_matrix::Storage::kRowSparse)
      .value("ROW_HASH", truncated_distance_matrix::Storage::kRowHash);

  nb::enum_<truncated_distance_matrix::ProtoType>(m, "TruncatedDistanceMatrixProto")
      .value("NEARNEIGHBOURS", truncated_distance_matrix::ProtoType::kNearNeighbours)
      .value("NEARNEIGHBOURS_INDICES",
             truncated_distance_matrix::ProtoType::kNearNeighboursIndices);

  nb::class_<truncated_distance_matrix::TruncatedDistanceMatrix>(
      m, "TruncatedDistanceMatrix")
      .def("__init__",
           [](truncated_distance_matrix::TruncatedDistanceMatrix* dm,
              const std::string& fname, truncated_distance_matrix::Storage storage,
              truncated_distance_matrix::ProtoType proto_type) {
             new (dm) truncated_distance_matrix::TruncatedDistanceMatrix();
             if (!dm->Build(fname, storage, proto_type)) {
               throw std::runtime_error("Cannot build TruncatedDistanceMatrix from " + fname);
             }
           },
           nb::arg("fname"),
           nb::arg("storage") = truncated_distance_matrix::Storage::kRowSparse,
           nb::arg("proto_type") = truncated_distance_matrix::ProtoType::kNearNeighbours)
      .def("size", &truncated_distance_matrix::TruncatedDistanceMatrix::size)
      .def("__len__", &truncated_distance_matrix::TruncatedDistanceMatrix::size)
      .def("number_distances",
           &truncated_distance_matrix::TruncatedDistanceMatrix::number_distances)
      .def("duplicate_distances_differing_by_one",
           &truncated_distance_matrix::TruncatedDistanceMatrix::
               duplicate_distances_differing_by_one)
      .def("max_stored_distance",
           &truncated_distance_matrix::TruncatedDistanceMatrix::MaxStoredDistance)
      .def("default_distance",
           &truncated_distance_matrix::TruncatedDistanceMatrix::DefaultDistance)
      .def("max_stored_distance_byte",
           &truncated_distance_matrix::TruncatedDistanceMatrix::MaxStoredDistanceByte)
      .def("default_distance_byte",
           &truncated_distance_matrix::TruncatedDistanceMatrix::DefaultDistanceByte)
      .def("set_default_distance", &SetTruncatedDistanceMatrixDefaultDistance,
           nb::arg("distance"))
      .def("set_default_distance_byte", &SetTruncatedDistanceMatrixDefaultDistanceByte,
           nb::arg("distance"))
      .def("index",
           [](const truncated_distance_matrix::TruncatedDistanceMatrix& dm,
              const std::string& name) { return dm.Index(name); },
           nb::arg("name"))
      .def("name", &TruncatedDistanceMatrixName, nb::arg("index"))
      .def("distance", &TruncatedDistanceMatrixDistance, nb::arg("i"), nb::arg("j"))
      .def("distance_or_default", &TruncatedDistanceMatrixDistanceOrDefault,
           nb::arg("i"), nb::arg("j"))
      .def("distances_or_default", &TruncatedDistanceMatrixDistancesOrDefault,
           nb::arg("i"), nb::arg("j"));

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

  nb::class_<dicer_api::Dicer>(m, "Dicer")
      .def(nb::init<>())
      .def("set_max_bonds_to_break", &dicer_api::Dicer::set_max_bonds_to_break,
           nb::arg("value"), "Set maximum number of bonds to break")
      .def("set_min_fragment_size", &dicer_api::Dicer::set_min_fragment_size,
           nb::arg("value"), "Set minimum fragment size")
      .def("set_max_fragment_size", &dicer_api::Dicer::set_max_fragment_size,
           nb::arg("value"), "Set maximum fragment size")
      .def("set_break_cc_bonds", &dicer_api::Dicer::set_break_cc_bonds,
           nb::arg("value"), "Control whether C-C bonds are broken")
      .def("set_break_cc_bonds_at_highly_connected",
           &dicer_api::Dicer::set_break_cc_bonds_at_highly_connected,
           nb::arg("value"), "Control whether C-[CD>2] bonds are broken")
      .def("set_break_amide_bonds", &dicer_api::Dicer::set_break_amide_bonds,
           nb::arg("value"), "Control whether amide bonds are broken")
      .def("set_label_join_points", &dicer_api::Dicer::set_label_join_points,
           nb::arg("value"), "Set isotope label for join points")
      .def("set_increment_isotope_for_join_points",
           &dicer_api::Dicer::set_increment_isotope_for_join_points,
           nb::arg("value"), "Increment isotopes at join points")
      .def("set_accumulate_global_fragment_count",
           &dicer_api::Dicer::set_accumulate_global_fragment_count,
           nb::arg("value"), "Accumulate fragments across dice() calls")
      .def("get_global_fragment_count", &dicer_api::Dicer::global_fragment_count,
           "Return accumulated global fragment counts")
      .def("set_perceive_symmetry_equivalent_matches",
           &dicer_api::Dicer::set_perceive_symmetry_equivalent_matches,
           nb::arg("value"), "Control symmetry-equivalent query matches")
      .def("set_determine_fragment_counts",
           &dicer_api::Dicer::set_determine_fragment_counts, nb::arg("value"),
           "Control whether per-molecule fragment counts are determined")
      .def("set_work_like_recap", &dicer_api::Dicer::set_work_like_recap,
           nb::arg("value"), "Work like Recap without recursion")
      .def("set_atom_type", &dicer_api::Dicer::set_atom_type, nb::arg("specification"),
           "Set atom typing used by dicer")
      .def("add_bond_break_smarts", &dicer_api::Dicer::AddBreakBondSmarts,
           nb::arg("smarts"), "Add bond-breaking SMARTS")
      .def("add_bond_break_query", &dicer_api::Dicer::AddBreakBondQuery,
           nb::arg("fname"), "Add bond-breaking query from a file directive")
      .def("add_fragment_requirement_smarts",
           &dicer_api::Dicer::AddFragmentRequirementSmarts, nb::arg("smarts"),
           "Require generated fragments to match SMARTS")
      .def("add_fragment_disqualifier_smarts",
           &dicer_api::Dicer::AddFragmentDisqualifierSmarts, nb::arg("smarts"),
           "Discard generated fragments matching SMARTS")
      .def("dice", &DiceMolecule, nb::arg("mol"),
           "Dice a molecule and return fragment unique smiles mapped to counts");

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
