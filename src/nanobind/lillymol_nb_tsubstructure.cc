#include "nanobind/lillymol_nb_internal.h"

namespace lillymol_nb {

using pybind_substructure::TSubstructure;

void
BindTSubstructure(nb::module_& m) {
  nb::class_<TSubstructure>(m, "TSubstructure")
      .def(nb::init<>())
      .def("read_queries", &TSubstructure::ReadQueries, nb::arg("directive"))
      .def("add_query_from_smarts", &TSubstructure::AddQueryFromSmarts,
           nb::arg("smarts"))
      .def("add_queries_from_smarts", &TSubstructure::AddQueriesFromSmarts,
           nb::arg("smarts"))
      .def("substructure_search",
           [](TSubstructure& ts, const std::string& smiles) {
             nb::gil_scoped_release release;
             return ts.SubstructureSearch(smiles);
           },
           nb::arg("smiles"))
      .def("substructure_search",
           [](TSubstructure& ts, Molecule& mol) {
             nb::gil_scoped_release release;
             return ts.SubstructureSearch(mol);
           },
           nb::arg("mol"))
      .def("substructure_search",
           [](TSubstructure& ts, std::vector<Molecule>& mols) {
             const uint32_t number_molecules = mols.size();
             std::vector<bool> results(number_molecules);
             {
               nb::gil_scoped_release release;
               for (uint32_t i = 0; i < number_molecules; ++i) {
                 results[i] = ts.SubstructureSearch(mols[i]);
               }
             }
             return results;
           },
           nb::arg("mols"))
      .def("substructure_search",
           [](TSubstructure& ts, const std::vector<std::string>& smiles) {
             const uint32_t number_molecules = smiles.size();
             std::vector<bool> results(number_molecules);
             {
               nb::gil_scoped_release release;
               Molecule mol;
               for (uint32_t i = 0; i < number_molecules; ++i) {
                 if (!mol.build_from_smiles(smiles[i])) {
                   std::cerr << "Invalid smiles '" << smiles[i] << "' ignored\n";
                   results[i] = false;
                 } else {
                   results[i] = ts.SubstructureSearch(mol);
                 }
               }
             }
             return results;
           },
           nb::arg("smiles"))
      .def("num_matches",
           [](TSubstructure& ts, const std::string& smiles) {
             nb::gil_scoped_release release;
             return ts.NumofMatches(smiles);
           },
           nb::arg("smiles"))
      .def("num_matches",
           [](TSubstructure& ts, Molecule& mol) {
             nb::gil_scoped_release release;
             return ts.NumofMatches(mol);
           },
           nb::arg("mol"))
      .def("num_matches",
           [](TSubstructure& ts, std::vector<Molecule>& mols) {
             nb::gil_scoped_release release;
             return ts.NumberMatches(mols);
           },
           nb::arg("mols"))
      .def("number_matches",
           [](TSubstructure& ts, const std::vector<std::string>& smiles) {
             nb::gil_scoped_release release;
             return ts.NumberMatches(smiles);
           },
           nb::arg("smiles"))
      .def("label_matched_atoms",
           [](TSubstructure& ts, const std::string& smiles) {
             return ts.LabelMatchedAtoms(smiles);
           },
           nb::arg("smiles"))
      .def("label_matched_atoms",
           [](TSubstructure& ts, Molecule& mol) {
             return ts.LabelMatchedAtoms(mol);
           },
           nb::arg("mol"))
      .def("matched_atoms",
           [](TSubstructure& ts, const std::string& smiles) {
             return ts.MatchedAtoms(smiles);
           },
           nb::arg("smiles"))
      .def("matched_atoms",
           [](TSubstructure& ts, Molecule& mol) {
             return ts.MatchedAtoms(mol);
           },
           nb::arg("mol"))
      .def("all_queries_match",
           [](TSubstructure& ts, Molecule& mol) {
             ts.must_match_all_queries = true;
             nb::gil_scoped_release release;
             return ts.SubstructureSearch(mol);
           },
           nb::arg("mol"))
      .def("number_queries", &TSubstructure::number_queries)
      .def("set_reduce_to_largest_fragment", &TSubstructure::set_reduce_to_largest_fragment,
           nb::arg("value"))
      .def("set_make_implicit_hydrogens_explicit",
           &TSubstructure::set_make_implicit_hydrogens_explicit, nb::arg("value"))
      .def("set_label_by_query_number", &TSubstructure::set_label_by_query_number,
           nb::arg("value"))
      .def("set_unique_embeddings_only", &TSubstructure::set_unique_embeddings_only,
           nb::arg("value"))
      .def("set_find_one_embedding_per_root_atom",
           &TSubstructure::set_find_one_embedding_per_root_atom, nb::arg("value"))
      .def("set_perceive_symmetry_equivalent_matches",
           &TSubstructure::set_perceive_symmetry_equivalent_matches, nb::arg("value"))
      .def("set_max_matches_to_find", &TSubstructure::set_max_matches_to_find,
           nb::arg("value"))
      .def("set_labeled_smiles_are_unique", &TSubstructure::set_labeled_smiles_are_unique,
           nb::arg("value"))
      .def_rw("query", &TSubstructure::query)
      .def_rw("must_match_all_queries", &TSubstructure::must_match_all_queries)
      .def_rw("isotope", &TSubstructure::isotope)
      .def("query_matched", &TSubstructure::QueryMatched)
      .def("query_names", &TSubstructure::query_names);
}

}  // namespace lillymol_nb
