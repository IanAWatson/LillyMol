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
#include "Molecule_Tools_Bdb/structure_database.h"

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


std::optional<std::string>
StructureDatabaseLookup(structure_database::StructureDatabase& db, Molecule& mol,
                        uint32_t params) {
  IWString ids_matched;
  if (!db.Lookup(mol, params, ids_matched)) {
    return std::nullopt;
  }
  return ids_matched.AsString();
}

bool
StructureDatabaseOpenRead(structure_database::StructureDatabase& db,
                          const std::string& dbname) {
  IWString tmp(dbname);
  return db.OpenForReading(tmp);
}

}  // namespace

NB_MODULE(lillymol_nb_bdb, m) {
  nb::enum_<structure_database::Lookup>(m, "LookupParams")
      .value("EXACT", structure_database::kExact)
      .value("STRIP", structure_database::kStrip)
      .value("NOCHIRAL", structure_database::kNoChiral)
      .value("GRAPH", structure_database::kGraph)
      .value("NOSTD", structure_database::kNoStandardise);

  m.def("value", [](structure_database::Lookup e) { return static_cast<int>(e); },
        nb::arg("lookup"));

  nb::class_<structure_database::StructureDatabase>(m, "StructureDatabase")
      .def(nb::init<>())
      .def("open_read", &StructureDatabaseOpenRead, nb::arg("dbname"),
           "Open paired structure databases for reading from a stem")
      .def("lookup",
           [](structure_database::StructureDatabase& db, Molecule& mol) {
             return StructureDatabaseLookup(db, mol, structure_database::kExact);
           },
           nb::arg("mol"), "Lookup exact molecule matches")
      .def("lookup", &StructureDatabaseLookup, nb::arg("mol"), nb::arg("params"),
           "Lookup molecule with a bit mask of LookupParams values");

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
           nb::arg("identifiers"), "Fetch Molecules for a list of identifiers")
      .def("close", &selimsteg::Selimsteg::Close, "Close the database. Idempotent")
      // A database holds a file handle, so give it the same shape as Reader - use
      // it in a with block, or call close(). Relying on the object being
      // collected works but is not something a caller can see or control.
      // Returns the python object, not a Selimsteg&. Returning a reference lets
      // the default return value policy copy it, and then the with block operates
      // on a copy while __exit__ closes the original.
      .def("__enter__", [](nb::object self) { return self; })
      // nanobind rejects None for an argument unless .none() says otherwise, and
      // python passes three Nones when the block exits without an exception.
      .def("__exit__",
           [](selimsteg::Selimsteg& s, nb::object, nb::object, nb::object) -> bool {
             s.Close();
             return false;   // do not suppress an exception from the block
           },
           nb::arg("exc_type").none(), nb::arg("exc_value").none(),
           nb::arg("traceback").none());
}
