/*
  Extract the coordinates of atoms matched by a query
*/

#include <stdlib.h>

#include <iostream>
#include <memory>
#include <string>

#include "google/protobuf/text_format.h"

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/charge_assigner.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/qry_wstats.h"
#include "Molecule_Lib/rmele.h"
#include "Molecule_Lib/standardise.h"

#ifdef BUILD_BAZEL
#include "Molecule_Tools/get_coordinates.pb.h"
#else
#include "get_coordinates.pb.h"
#endif

namespace get_coordinates {

using std::cerr;

static void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off

  cerr << R"(Extracts coordinates of the atoms matching query atoms.
get_coordinates -s '[OH]c1ccccc1' file.sdf > file.get_coordinates
The first coordinates written will be those of the Oxygen atom.
Note the symmetry associated with the attached phenyl ring. It is not unique in 2D, but in 3d it is
so both matches may need to be processed.

 -q <query>         specify one or more queries
 -s <smarts>        specify query as smarts
 -z i               ignore queries not hitting
 -X <symbol>        extract/remove all atoms of type <symbol>. Useful for inputs with explicit H atoms.
 -C <symbol>        add a centroid atom to the output.
 -u                 set find unique embeddings only - beware unique embeddings in 2D will not be unique in 3D.
 -o ...             output type. 'def', 'csv', 'textproto'.
 -g ...             chemical standardisation.
 -A ...             standard aromaticity options.
 -v                 verbose output
)";
  // clang-format on

  exit(rc);
}


class Options {
  private:
    int _verbose;

    uint64_t _molecules_read;

    resizable_array_p<Substructure_Hit_Statistics> _queries;

    extending_resizable_array<int> _queries_hit;

    int _ignore_queries_not_matching = 0;

    // If requested with the -C option.
    const Element* _centroid;

    int _csv_output;
    int _write_proto;

    IWString _output_separator;

    Chemical_Standardisation _chemical_standardisation;

    Elements_to_Remove _rmele;

    Charge_Assigner _charge_assigner;

  // private functions.

    int Preprocess(Molecule& m);
    int WriteCsvHeader(IWString_and_File_Descriptor& output) const;
    int WriteProto(const get_coordinates::MoleculeData& proto,
                IWString_and_File_Descriptor& output) const;
    void WriteCoordindatesSymbol(const Atom& a, IWString_and_File_Descriptor& output);
    int WriteCsv(Molecule& m, int query_number, int embedding_number,
                  const Set_of_Atoms& embedding, IWString_and_File_Descriptor& output);
    int WriteCentroid(const Molecule& m, const Set_of_Atoms& embedding,
               IWString_and_File_Descriptor& output) const;
    int Process(Molecule& m, int query_number, int embedding_number,
                 const Set_of_Atoms& embedding, 
                 std::unique_ptr<get_coordinates::MoleculeData>& proto,
                 IWString_and_File_Descriptor& output);
    int Process(Molecule& m, int query_number, const Substructure_Results& sresults,
               std::unique_ptr<get_coordinates::MoleculeData>& maybe_proto,
               IWString_and_File_Descriptor& output);
    int AddToProto(const Molecule& m, uint32_t query_number, int embedding_number,
                const Set_of_Atoms& embedding, get_coordinates::MoleculeData& proto) const;

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _molecules_read = 0;
  _csv_output = 0;
  _write_proto = 0;
  _centroid = nullptr;
  _output_separator = ' ';
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');


  if (cl.option_present('g')) {
    if (!_chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      Usage(6);
    }
  }

  if (!cl.option_present('q') && !cl.option_present('s')) {
    cerr << "Must specify a substructure query via the -q or -s option\n";
    Usage(3);
  }

  if (cl.option_present('q')) {
    if (!process_queries(cl, _queries, _verbose)) {
      cerr << "cannot process queries from -q option(s)\n";
      return 6;
    }
  }

  if (cl.option_present('s')) {
    const_IWSubstring smarts;
    for (int i = 0; cl.value('s', smarts, i); ++i) {
      std::unique_ptr<Substructure_Hit_Statistics> q =
                std::make_unique<Substructure_Hit_Statistics>();
      if (!q->create_from_smarts(smarts)) {
        cerr << "Cannot parse smarts '" << smarts << "'\n";
        return 0;
      }

      _queries << q.release();
    }
  }

  if (cl.option_present('u')) {
    for (Substructure_Hit_Statistics* q : _queries) {
      q->set_find_unique_embeddings_only(1);
    }
  }

  if (cl.option_present('o')) {
    IWString o = cl.string_value('o');
    if (o == "def" || o == "default") {
    } else if (o == "csv") {
      _csv_output = 1;
      _output_separator = ',';
    } else if (o == "textproto") {
      _write_proto = 1;
    } else {
      cerr << "Unrecognised -o qualifier '" << o << "'\n";
      return 0;
    }
  }

  if (cl.option_present('N')) {
    if (!_charge_assigner.construct_from_command_line(cl, _verbose, 'N')) {
      cerr << "Cannot initialise charge assigner (-N option)\n";
      Usage(5);
    }
  }

  if (cl.option_present('z')) {
    _ignore_queries_not_matching = 1;
    if (_verbose) {
      cerr << "Will skip queries not matching\n";
    }
  }

  if (cl.option_present('X')) {
    if (! _rmele.construct_from_command_line(cl, _verbose, 'X')) {
      cerr << "Cannot initialise remove elements (-X)\n";
      return 0;
    }
  }

  if (cl.option_present('C')) {
    IWString s = cl.string_value('C');
    _centroid = get_element_from_symbol_no_case_conversion(s);
    // Should check the return, maybe autocreate...
    if (_verbose) {
      cerr << "Will add a centroid atom of type " << _centroid->symbol() << '\n';
    }
  }

  return 1;
}

int
Options::WriteCsvHeader(IWString_and_File_Descriptor& output) const {
  output << "Smiles";
  output << _output_separator << "Id";
  output << _output_separator << "Qry";
  output << _output_separator << "Match";

  const int n = _queries.front()->max_atoms_in_query();

  for (int i = 0; i < n; ++i) {
    output << _output_separator << "atom" << i;
    output << _output_separator << 'x' << i;
    output << _output_separator << 'y' << i;
    output << _output_separator << 'z' << i;
    output << _output_separator << "Sym" << i;
  }

  if (_centroid != nullptr) {
    output << _output_separator << "centroid";
    output << _output_separator << "centroid_x";
    output << _output_separator << "centroid_y";
    output << _output_separator << "centroid_z";
  }

  output << '\n';

  return 1;
}

#ifdef LEGACY_CODE_NO_LONGER_USED
static int
set_atom_names(const IWString& query_name, const Set_of_Atoms& embedding,
               iwaray<IWString>& atom_label) {
  for (int i = 0; i < embedding.number_elements(); i++) {
    atom_number_t j = embedding[i];

    IWString& aname = atom_label[j];
    if (aname.length()) {
      aname += ':';
    }

    aname += query_name;
  }

  return 1;
}

static int
set_atom_names(const IWString& query_name, const Substructure_Results& sresults,
               iwaray<IWString>& atom_label) {
  int nhits = sresults.hits_found();

  for (int i = 0; i < nhits; i++) {
    const Set_of_Atoms* embedding = sresults.embedding(i);
    (void)set_atom_names(query_name, *embedding, atom_label);
  }

  return 1;
}
#endif

int
Options::WriteCentroid(const Molecule& m, const Set_of_Atoms& embedding,
               IWString_and_File_Descriptor& output) const {
  double xsum = 0;
  double ysum = 0;
  double zsum = 0;
  int n = 0;
  for (atom_number_t atom : embedding) {
    if (atom == kInvalidAtomNumber) {
      continue;
    }
    ++n;
    const Atom& a = m[atom];
    xsum += a.x();
    ysum += a.y();
    zsum += a.z();
  }

  // Make the order the same as WriteCoordindatesSymbol.
  output << _output_separator << static_cast<float>(xsum / n);
  output << _output_separator << static_cast<float>(ysum / n);
  output << _output_separator << static_cast<float>(zsum / n);
  output << _output_separator << _centroid->symbol();

  return 1;
}

void
Options::WriteCoordindatesSymbol(const Atom& a, IWString_and_File_Descriptor& output) {
  output << _output_separator << a.x() <<
            _output_separator << a.y() << _output_separator << a.z() <<
            _output_separator << a.atomic_symbol();
}

int
Options::WriteCsv(Molecule& m, int query_number, int embedding_number,
                  const Set_of_Atoms& embedding, IWString_and_File_Descriptor& output) {

  output << m.smiles();
  output << _output_separator << m.name();
  output << _output_separator << query_number;
  output <<  _output_separator << embedding_number;
  for (atom_number_t atom : embedding) {
    if (atom == kInvalidAtomNumber) {
      continue;
    }
    output << _output_separator << atom;
    WriteCoordindatesSymbol(m[atom], output);
  }

  if (_centroid != nullptr) {
    WriteCentroid(m, embedding, output);
  }

  output << '\n';

  return 1;
}

// First determine if we have a PerQueryResult for this query number.
int
Options::AddToProto(const Molecule& m, uint32_t query_number, int embedding_number,
                const Set_of_Atoms& embedding, get_coordinates::MoleculeData& proto) const {

  get_coordinates::PerQueryResult* result = nullptr;
  for (get_coordinates::PerQueryResult& pqr : *proto.mutable_per_query_result()) {
    if (pqr.query_number() == query_number) {
      result = &pqr;  // this seems a little dangerous...
      break;
    }
  }
  if (result == nullptr) {
    result = proto.add_per_query_result();
    result->set_query_number(query_number);
  }

  get_coordinates::MatchedAtoms* matched_atoms = result->add_atoms();

  for (atom_number_t atom : embedding) {
    if (atom == kInvalidAtomNumber) {
      continue;
    }
    const Atom& a = m[atom];
    get_coordinates::MatchedAtom* extra = matched_atoms->add_atom();
    extra->set_element(a.atomic_symbol().data(), a.atomic_symbol().length());
    extra->set_atom_number(atom);
    extra->set_x(a.x());
    extra->set_y(a.y());
    extra->set_z(a.z());
  }

  return 1;
}

int
Options::Process(Molecule& m, int query_number, int embedding_number,
                 const Set_of_Atoms& embedding, 
                 std::unique_ptr<get_coordinates::MoleculeData>& proto,
                 IWString_and_File_Descriptor& output) {
  if (_csv_output) {
    return WriteCsv(m, query_number, embedding_number, embedding, output);
  } else if (proto) {
    return AddToProto(m, query_number, embedding_number, embedding, *proto);
  }

  output << m.smiles() << _output_separator << m.name() << _output_separator << query_number 
         << _output_separator << embedding_number << '\n';
  const int n = embedding.number_elements();
  for (int i = 0; i < n; ++i) {
    atom_number_t j = embedding[i];
    if (j == kInvalidAtomNumber) {
      continue;
    }
    output << i << _output_separator << j;
    WriteCoordindatesSymbol(m[j], output);
    output << _output_separator << m.smarts_equivalent_for_atom(j);
    output << '\n';
  }
  if (_centroid != nullptr) {
    output << '.';
    output << _output_separator << '.';
    WriteCentroid(m, embedding, output);
    output << _output_separator << '.';
    output << '\n';
  }

  output << "-----------------\n";

  output.write_if_buffer_holds_more_than(4096);

  return output.good();
}

int
Options::Process(Molecule& m, int query_number, const Substructure_Results& sresults,
               std::unique_ptr<get_coordinates::MoleculeData>& maybe_proto,
               IWString_and_File_Descriptor& output) {
  const int n = sresults.number_embeddings();
  for (int i = 0; i < n; ++i) {
    Process(m, query_number, i, *sresults.embedding(i), maybe_proto, output);
  }

  output.write_if_buffer_holds_more_than(4096);

  return output.good();
}

int
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {

  ++_molecules_read;
  if (_molecules_read == 1) [[unlikely]] {
    if (_csv_output) {
      WriteCsvHeader(output);
    }
  }

  Preprocess(m);

  std::unique_ptr<get_coordinates::MoleculeData> maybe_proto;

  if (_write_proto) {
    maybe_proto = std::make_unique<get_coordinates::MoleculeData>();
    maybe_proto->set_name(m.name().data(), m.name().length());
  }

  int nq = _queries.number_elements();

  int rc = 0;  // The number of queries that hit

  for (int i = 0; i < nq; i++) {
    Substructure_Results sresults;

    int tmp = _queries[i]->substructure_search(m, sresults);
    if (_verbose > 1) {
      cerr << ' ' << tmp << " hits to query " << i << " in " << m.name() << '\n';
    }

    if (0 == tmp) {
      if (_ignore_queries_not_matching) {
        continue;
      }

      cerr << "Yipes, zero hits to query " << i << " '" << _queries[i]->comment() << "'\n";
      return 0;
    }

    Process(m, i, sresults, maybe_proto, output);

    rc++;
  }

  _queries_hit[rc]++;

  if (maybe_proto) {
    return WriteProto(*maybe_proto, output);
  }

  return output.good();
}

int
Options::WriteProto(const get_coordinates::MoleculeData& proto,
                IWString_and_File_Descriptor& output) const {
  if (proto.per_query_result_size() == 0) {
    return 1;
  }

  static google::protobuf::TextFormat::Printer printer;
  printer.SetSingleLineMode(true);

  std::string buffer;
  printer.PrintToString(proto, &buffer);
  output << buffer << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (_rmele.active()) {
    _rmele.process(m);
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_charge_assigner.active()) {
    _charge_assigner.process(m);
  }


  return 1;
}

int
Options::Report(std::ostream& output) const {
  cerr << "Read " << _molecules_read << " molecules\n";
  for (Substructure_Hit_Statistics* q : _queries) {
    q->report(output, _verbose);
  }

  for (int i = 0; i < _queries_hit.number_elements(); i++) {
    if (_queries_hit[i]) {
      cerr << _queries_hit[i] << " queries matched " << i << " of the queries\n";
    }
  }

  return 1;
}

int
GetCoordinates(Options& options, Molecule& m, IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
GetCoordinates(Options& options, data_source_and_type<Molecule>& input,
               IWString_and_File_Descriptor& output) {
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (!GetCoordinates(options, *m, output)) {
      return 0;
    }
  }

  return output.good();
}

int
GetCoordinates(Options& options, const char* fname, FileType input_type,
               IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return GetCoordinates(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:q:z:s:N:g:uo:X:C:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(6);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  }

  if (FILE_TYPE_INVALID == input_type && !all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    if (!GetCoordinates(options, fname, input_type, output)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace get_coordinates

int
main(int argc, char** argv) {
  int rc = get_coordinates::Main(argc, argv);

  return rc;
}
