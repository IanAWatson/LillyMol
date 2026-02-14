// loading a BerkeleyDB of conformers. Use conform_database_lookup for retrievals.

#include "sys/stat.h"
#include "sys/types.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>

#include "db_cxx.h"

#include "Foundational/iwaray/iwaray.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools_Bdb/conformer_database.pb.h"
#else
#include "conformer_database.pb.h"
#endif

namespace conformer_database {

using std::cerr;

void 
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format on
  // clang-format off

  ::exit(rc);
}

class ConformerDB{
  private:
    std::unique_ptr<Db> _db;

    uint64_t _molecules_read;
    uint64_t _duplicates;

    extending_resizable_array<uint32_t> _number_conformers;

    Chemical_Standardisation _chemical_standardisation;

    Report_Progress _report_progress;

    int _verbose;

    // Private functions

  public:
    ConformerDB();
    ~ConformerDB();

    int WriteConformers(Molecule& m, const conformer_database::Conformers& conformers);
    int AddConformer(Molecule& m,
        const std::unique_ptr<int[]>& atom_order_in_smiles,
        conformer_database::Conformers& proto);

    int Initialise(Command_Line& cl);

    int Process(data_source_and_type<Molecule>& input);

    int Report(std::ostream& output) const;
};

ConformerDB::ConformerDB() {
  _molecules_read = 0;
  _duplicates = 0;
  _verbose = 0;
}

ConformerDB::~ConformerDB() {
}

int
ConformerDB::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (! cl.option_present('d')) {
    cerr << "Must specify database via the -d option\n";
    Usage(1);
  }

  if (cl.option_present('d')) {
    const char *dbname = cl.option_value('d');

    int flags;
    DBTYPE dbtype;
    int mode;

    dbtype = DB_BTREE;
    flags = DB_CREATE;
    mode = S_IREAD | S_IWRITE | S_IRGRP | S_IROTH;

    _db.reset(new Db(NULL, DB_CXX_NO_EXCEPTIONS));

    if (int rc = _db->open(NULL, dbname, NULL, dbtype, flags, mode); rc != 0) {
      cerr << "ConformerDB::Initialise:cannot open database '" << dbname << "'\n";
      _db->err(rc, "");
      _db.reset();
      return 0;
    }

    if (_verbose) {
      cerr << "Smiles will be written to database '" << dbname << "'\n";
    }
  }

  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "ConformerDB::Initialise:cannot initialise chemical standardisation (-g)\n";
      return 0;
    }
  }

  if (cl.option_present('r')) {
    if (! _report_progress.initialise(cl, 'r', _verbose)) {
      cerr << "ConformerDB::Initialise:cannot initialise progress reporting (-r)\n";
      return 0;
    }
  }

  return 1;
}

int
ConformerDB::AddConformer(Molecule& m,
        const std::unique_ptr<int[]>& atom_order_in_smiles,
        conformer_database::Conformers& proto) {
  conformer_database::Conformer* c  = proto.mutable_conformers()->Add();
  const int matoms = m.natoms();
  c->mutable_xyz()->Resize(matoms * 3, 0.0f);
  for (int i = 0; i < matoms; ++i) {
    int j = atom_order_in_smiles[i];
    const Atom& a = m[j];
    c->set_xyz(3 * i, a.x());
    c->set_xyz(3 * i + 1, a.y());
    c->set_xyz(3 * i + 2, a.z());
  }

  return 1;
}

int
Reset(Molecule& m, std::unique_ptr<int[]>& atom_order_in_smiles,
      conformer_database::Conformers& proto) {
  proto.Clear();
  const IWString& n = m.name();
  proto.set_id(n.rawdata(), n.length());
  const IWString& s = m.smiles();
  proto.set_smiles(s.rawdata(), s.length());

  const resizable_array<atom_number_t>& aois = m.atom_order_in_smiles();
  const int matoms = m.natoms();
  //for (int i = 0; i < matoms; ++i) {
  //  cerr << i << ' ' << m.smarts_equivalent_for_atom(aois[i]) << " aois " << aois[i] << '\n';
  //}
  atom_order_in_smiles.reset(new int[matoms]);

  std::copy_n(aois.rawdata(), matoms, atom_order_in_smiles.get());

  return 1;
}

int
ConformerDB::WriteConformers(Molecule& m, const conformer_database::Conformers& conformers) {
  // This will happen on the first call.
  if (conformers.conformers_size() == 0) [[ unlikely]] {
    return 1;
  }

  _number_conformers[conformers.conformers_size()]++;

  m.remove_all(1);

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  const IWString& usmi = m.unique_smiles();

  Dbt dkey((void *) usmi.rawdata(), usmi.length());
  Dbt zdata;
  int rc = _db->get(NULL, &dkey, &zdata, 0);
  if (rc == 0) {
    ++_duplicates;
    if (_verbose > 1) {
      cerr << "ConformerDB::WriteConformers:existing data found\n";
      cerr << m.smiles() << ' ' << m.name() << '\n';
    }
    return 0;
  }

  std::string serialized;
  conformers.SerializeToString(&serialized);

  zdata.set_data(serialized.data());
  zdata.set_size(serialized.size());

  rc = _db->put(NULL, &dkey, &zdata, 0);
  if (rc == 0) {
    if (_report_progress()) {
      cerr << "Stored conformers for " << _molecules_read << " molecules\n";
    }

    return 1;
  }

  cerr << "Error storing " << m.name() << '\n';
  _db->err(rc, "");

  return 0;
}

int
ConformerDB::Report(std::ostream& output) const {
  output << "Read " << _molecules_read << " molecules\n";
  output << "Skipped " << _duplicates << " duplicate data\n";
  for (int i = 0; i < _number_conformers; ++i) {
    if (_number_conformers[i]) {
      output << _number_conformers[i] << " molecules had " << i << " conformers\n";
    }
  }

  return 1;
}

int
NextSetOfConformers(data_source_and_type<Molecule>& input,
                    resizable_array_p<Molecule>& molecules) {
  molecules.resize_keep_storage(0);

  static Molecule* next_molecule = nullptr;

  IWString mname;
  // all molecules in a set of conformers must have the same number of atoms.
  int natoms = 0;

  if (next_molecule != nullptr) [[ likely ]] {
    molecules << next_molecule;
    mname = next_molecule->name();
    natoms = next_molecule->natoms();
    next_molecule = nullptr;
  }

  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    if (mname.empty()) [[ unlikely ]] {
      mname = m->name();
      natoms = m->natoms();
      molecules << m;
    } else if (mname == m->name()) {
      if (m->natoms() != natoms) {
        cerr << "NextSetOfConformers:atom count mismatch: expected " << natoms << 
                " got " << m->natoms() << ' ' << m->name() << '\n';
        return 0;
      }

      molecules << m;
    } else {
      next_molecule = m;
      return molecules.number_elements();
    }
  }

  return molecules.number_elements();
}

int
ConformerDB::Process(data_source_and_type<Molecule>& input) {
  resizable_array_p<Molecule> conformers;

  while (NextSetOfConformers(input, conformers)) {
    _molecules_read += conformers.size();

    std::unique_ptr<int[]> atom_order_in_smiles;
    conformer_database::Conformers proto;
    Reset(*conformers[0], atom_order_in_smiles, proto);

    for (Molecule* m : conformers) {
      AddConformer(*m, atom_order_in_smiles, proto);
    }
    WriteConformers(*conformers[0], proto);
  }

  return 1;
}

int
ConformerDatabaseStore(ConformerDB& cdb, const char* fname, FileType input_type) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    assert(0 != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return cdb.Process(input);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:g:d:i:r:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_present('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(5);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(0);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (!cl.option_present('i') && 1 == cl.number_elements() && 0 == strcmp("-", cl[0])) {
    input_type = FILE_TYPE_SMI;
  } else if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  }

  if (0 == input_type && !all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  ConformerDB cdb;
  if (! cdb.Initialise(cl)) {
    cerr << "Cannot initialise database\n";
    return 1;
  }

  for (const char* fname : cl) {
    if (! ConformerDatabaseStore(cdb, fname, input_type)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    cdb.Report(cerr);
  }

  return 0;
}


}  // namespace conformer_database

int
main(int argc, char **argv) {
  int rc = conformer_database::Main(argc, argv);

  return rc;
}
