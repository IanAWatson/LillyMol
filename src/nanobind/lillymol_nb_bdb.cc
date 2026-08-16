// BerkeleyDB-backed tool bindings for the nanobind LillyMol pilot.

#include <iostream>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "Molecule_Lib/molecule.h"
#include "Molecule_Tools_Bdb/iwecfp_database_lookup_lib.h"
#include "Molecule_Tools_Bdb/selimsteg.h"

namespace nb = nanobind;

namespace {

std::vector<int>
SyntheticPrecedentPerShellData(iwecfp_database_lookup::SP_Set_of_Databases& databases,
                               Molecule& mol) {
  std::vector<int> result(3);
  if (!databases.PerShellData(mol, result)) {
    std::cerr << "SyntheticPrecedentDatabases::per_shell_data failed\n";
  }
  return result;
}

std::string
SyntheticPrecedentDatabasesRepr(
    const iwecfp_database_lookup::SP_Set_of_Databases&) {
  return "Synthetic Precedent databases";
}

}  // namespace

NB_MODULE(lillymol_nb_bdb, m) {
  nb::class_<iwecfp_database_lookup::SP_Set_of_Databases>(
      m, "SyntheticPrecedentDatabases")
      .def(nb::init<>())
      .def("add_database", &iwecfp_database_lookup::SP_Set_of_Databases::AddDatabase,
           nb::arg("dbname"), "Add an existing synthetic precedent database")
      .def("set_max_radius",
           &iwecfp_database_lookup::SP_Set_of_Databases::set_max_radius,
           nb::arg("max_radius"), "Set maximum shell radius")
      .def("per_shell_data", &SyntheticPrecedentPerShellData, nb::arg("mol"),
           "Report lowest bit prevalence at each radius")
      .def("slurp", &iwecfp_database_lookup::SP_Set_of_Databases::slurp,
           nb::arg("min_examples"),
           "Slurp database entries with at least min_examples examples to memory")
      .def("__repr__", &SyntheticPrecedentDatabasesRepr);

  nb::class_<selimsteg::Selimsteg>(m, "Selimsteg")
      .def(nb::init<>())
      .def("open_database", &selimsteg::Selimsteg::OpenDatabase,
           nb::arg("dbname"), "Open a BerkeleyDB database with id-to-smiles mappings")
      .def("get_smiles", &selimsteg::Selimsteg::Lookup, nb::arg("identifier"),
           "Fetch the smiles for an identifier")
      .def("get_molecule", &selimsteg::Selimsteg::GetMolecule,
           nb::arg("identifier"), "Fetch a Molecule for an identifier")
      .def("get_molecules", &selimsteg::Selimsteg::GetMolecules,
           nb::arg("identifiers"), "Fetch Molecules for a list of identifiers");
}
