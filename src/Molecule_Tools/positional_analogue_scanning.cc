// Implementation of Positional Analogue Scanning

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "absl/container/flat_hash_set.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/combinations.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/element.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/target.h"

namespace positional_analogue_scanning {

using std::cerr;

class SetOfQueries {
  private:
    resizable_array_p<Substructure_Query> _query;

  public:
    // `fname` is passed to queries_from_file to populate `_query`.
    int Build(IWString& fname);

    int number_queries() const {
      return _query.number_elements();
    }

    int LabelMatchedAtoms(Molecule_to_Match& target, int* matched,
                          int flag, int first_match_only,
                          Set_of_Atoms& destination);
};

int
SetOfQueries::Build(IWString& fname) {
  static constexpr int kInheritDirectoryPath = 1;
  static constexpr int kVerbose = 0;
  if (!  queries_from_file(fname, _query, kInheritDirectoryPath, kVerbose)) {
    cerr << "SetOfQueries::Build:cannot read queries from '" << fname << "'\n";
    return 0;
  }

  return _query.number_elements();
}

// For each query match, fetch the first matched atom, and set the corresponding value
// in `matched` to `flag`.
// Make sure that only one set of queries matches any particular atom.
int
SetOfQueries::LabelMatchedAtoms(Molecule_to_Match& target, int* matched, int flag, int first_match_only,
                Set_of_Atoms& destination) {
  int rc = 0;
  for (Substructure_Query* q : _query) {
    Substructure_Results sresults;
    if (! q->substructure_search(target, sresults)) {
      continue;
    }
    for (const Set_of_Atoms* e : sresults.embeddings()) {
      atom_number_t a = e->first();
      destination << a;

//    if (matched[a] == 0) {
//      matched[a] = flag;
//      destination << a;
//    } 
    }

    if (first_match_only) {
      return 1;
    }
    ++rc;
  }

  return rc;
}

class SetOfMolecules {
  private:
    resizable_array_p<Molecule> _molecule;

  public:
    int Build(IWString& fname);

    int number_molecules() const {
      return _molecule.number_elements();
    }

    Molecule** begin() const {
      return _molecule.begin();
    }

    Molecule** end() const {
      return _molecule.end();
    }
};

int
SetOfMolecules::Build(IWString& fname) {
  data_source_and_type<Molecule> input(FileType::FILE_TYPE_SMI, fname);
  if (! input.good()) {
    cerr << "SetOfMolecules::Build:cannot open '" << fname << "'\n";
    return 0;
  }

  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    cerr << m->smiles() << ' ' << m->name() << " read\n";
    _molecule << m;
  }

  return _molecule.number_elements();
}

void
Usage() {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  cerr << R"(Positional analogue implementation.
 -c <n>         number of combinations to generate for each set of substituents, usually 1, 2 or 3.
 -q <file>      file with queries specifing matche datoms. 
 -v             verbose output.
)";
  // clang-format off

  ::exit(0);
}

class Options {
  private:
    int _verbose;

    int _combo;

    int _nq;
    std::unique_ptr<SetOfQueries[]> _query;

    std::unique_ptr<SetOfMolecules[]> _molecule;

    // For each query, the atoms matched;
    std::unique_ptr<Set_of_Atoms[]> _atoms;

    uint64_t _molecules_processed;
    uint64_t _molecules_written;

    absl::flat_hash_set<IWString> _seen;

    // By default, we do NOT write the starting molecule.
    int _write_starting_molecule;

    // We need to assign a unique name to each product molecule.
    // This is sequential with each starting molecule.
    uint32_t _made_this_molecule;

    extending_resizable_array<uint32_t> _generated;

    int _ignore_no_query_matches;

    uint32_t _max_products_per_starting_molecule;

  // Private functions
    int MaybeOutput(Molecule& m, IWString_and_File_Descriptor& output);
    int Process(Molecule& m,
                 const Set_of_Atoms& matched_atoms,
                 int istart,
                 IWString_and_File_Descriptor& output);
    int Process(Molecule& m,
                 const std::vector<uint32_t>& indices,
                 IWString_and_File_Descriptor& output);

  public:
    Options();

    int Initialise(Command_Line& cl);

    int Process(Molecule& m, IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _combo = 2;
  _nq = 0;
  _molecules_processed = 0;
  _molecules_written = 0;
  _made_this_molecule = 0;
  _ignore_no_query_matches = 0;
  _max_products_per_starting_molecule = std::numeric_limits<uint32_t>::max();
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (cl.option_present('c')) {
    if (! cl.value('c', _combo) || _combo < 1) {
      cerr << "The number of simultaneous attachments (-c) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will generate " << _combo << " way changes\n";
    }
  }

  if (! cl.option_present('q')) {
    cerr << "Options::Initialise:must specify one or more smarts via the -s option\n";
    return 0;
  }

  _nq = cl.option_count('q');
  _query.reset(new SetOfQueries[_nq]);
  _atoms.reset(new Set_of_Atoms[_nq]);

  if (cl.option_present('q')) {
    IWString fname;
    for (int i = 0; cl.value('q', fname, i); ++i) {
      if (! _query[i].Build(fname)) {
        cerr << "Options::Initialise:cannot read queries from '" << fname << "'\n";
        return 0;
      }
    }

    if (_verbose) {
      int totalq = 0;
      for (int i = 0; i < _nq; ++i) {
        totalq += _query[i].number_queries();
      }

      cerr << "Read " << _nq << " query sets containing " << totalq << " substructure queries\n";
    }
  }

  if ((cl.number_elements() - 1) != _nq) {
    cerr << "Options::Initialise:read " << _nq << " sets of queries but " << (cl.size() - 1) << 
                " input files on command line, impossible, must be the same\n";
//  return 0;
  }

  cerr << _nq << " nq\n";
  _molecule.reset(new SetOfMolecules[_nq]);

  for (int i = 1; i < cl.number_elements(); ++i) {
    IWString fname(cl[i]);
    cerr << "i = " << i << " file " << fname << "\n";
    if (! _molecule[i - 1].Build(fname)) {
      cerr << "Cannot read molecules from '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('z')) {
    _ignore_no_query_matches = 1;
    if (_verbose) {
      cerr << "Will ignore molecules not matching any query\n";
    }
  }

  if (cl.option_present('x')) {
    if (! cl.value('x', _max_products_per_starting_molecule) || _max_products_per_starting_molecule < 1) {
      cerr << "The max products per starting molecule (-x) option must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will generate a max of " << _max_products_per_starting_molecule <<
              " products per starting molecule\n";
    }
  }

  if (cl.option_present('p')) {
    _write_starting_molecule = 1;
    if (_verbose) {
      cerr << "Will write the starting molecule\n";
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Options:Read " << _molecules_processed << " wrote " << _molecules_written << '\n';

  for (int i = 0; i < _generated.number_elements() ; ++i) {
    if (_generated[i]) {
      output << _generated[i] << " starting molecules made " << i << " products\n";
    }
  }

  return 1;
}

int
Options::Process(Molecule& m,
                 IWString_and_File_Descriptor& output) {
  ++_molecules_processed;

  _made_this_molecule = 0;

  if (_write_starting_molecule) {
    output << m.smiles() << ' ' << m.name() << '\n';
  }

  const int matoms = m.natoms();
  std::unique_ptr<int[]> label = std::make_unique<int[]>(matoms);
  std::fill_n(label.get(), matoms, 0);

  Molecule_to_Match target(&m);

  static constexpr int kFirstMatchOnly = 1;
  for (int i = 0; i < _nq; ++i) {
    if (! _query[i].LabelMatchedAtoms(target, label.get(), i + 1, kFirstMatchOnly, _atoms[i])) {
      // cerr << "Options::Process:no query matches " << m.name() << " query set " << i << '\n';
      return _ignore_no_query_matches;
    }
    // cerr << "query " << i << " " << _atoms[i] << '\n';
  }

  std::vector<uint32_t> count(_nq);
  for (int i = 0; i < _nq; ++i) {
    count[i] = _atoms[i].size();
  }
  combinations::Combinations<uint32_t> indices(count);

  std::vector<uint32_t> state;
  state.resize(_nq, 0);

  do {
    Process(m, state, output);
    if (_made_this_molecule > _max_products_per_starting_molecule) {
      break;
    }
  } while (indices.Next(state));

  ++_generated[_made_this_molecule];

  return 1;
}
#ifdef N_CHOOSE_K
from https://stackoverflow.com/questions/28711797/generating-n-choose-k-permutations-in-c
void PermGenerator(int n, int k)
{
    std::vector<int> d(n);
    std::iota(d.begin(),d.end(),1);
    cout << "These are the Possible Permutations: " << endl;
    do
    {
        for (int i = 0; i < k; i++)
        {
            cout << d[i] << " ";
        }
        cout << endl;
        std::reverse(d.begin()+k,d.end());
    } while (next_permutation(d.begin(),d.end()));
}

#endif

int
Options::Process(Molecule& m,
                 const std::vector<uint32_t>& indices,
                 IWString_and_File_Descriptor& output) {
  Set_of_Atoms matched_atoms;
  for (int i = 0; i < _nq; ++i) {
    int j = indices[i];
    // cerr << "Query " << i << " atom " << j << '\n';
    if (! matched_atoms.add_if_not_already_present(_atoms[i][j])) {
      return 0;
    }
  }

  return Process(m, matched_atoms, 0, output);
}

int
Options::Process(Molecule& m,
                 const Set_of_Atoms& matched_atoms,
                 int istart,
                 IWString_and_File_Descriptor& output) {
  atom_number_t zatom = matched_atoms[istart];

  if (m.hcount(zatom) == 0) {
    cerr << m.smiles() << " no available H on atom " << zatom <<
            m.smarts_equivalent_for_atom(zatom) << '\n';
    return 0;
  }

  const int initial_natoms = m.natoms();
  // cerr << "Starting moleculehas " << initial_natoms << " atoms\n";

  for (int i = istart; i < _nq; ++i) {
    for (Molecule* frag : _molecule[i]) {
      m.add_molecule(frag);
      m.add_bond(zatom, m.natoms() - 1, SINGLE_BOND);
      MaybeOutput(m, output);
      if (istart < (_nq - 1)) {
        Process(m, matched_atoms, istart + 1, output);
      }
      m.resize(initial_natoms);
    }
    if (_made_this_molecule > _max_products_per_starting_molecule) {
      return 1;
    }
  }
  // cerr << "Finished iteration\n";

  return 1;
}

int
Options::MaybeOutput(Molecule& m, IWString_and_File_Descriptor& output) {
  if (const auto iter = _seen.find(m.unique_smiles()); iter != _seen.end()) {
    return 0;
  }

  _seen.insert(m.unique_smiles());

  ++_made_this_molecule;

  output << m.smiles() << ' ' << m.name() << '.' << _made_this_molecule << '\n';

  ++_molecules_written;

  return 1;
}

int
PositionalAnalogueScanning(Molecule& m,
                Options& options,
                IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int PositionalAnalogueScanning(data_source_and_type<Molecule>& input,
                Options& options,
                IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    if (! PositionalAnalogueScanning(*m, options, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

int
PositionalAnalogueScanning(const char* fname, Options& options,
                        IWString_and_File_Descriptor& output) {
  data_source_and_type<Molecule> input(fname);
  if (! input.good()) {
    cerr << "PositionalAnalogueScanning:cannot open '" << fname << "'\n";
    return 0;
  }

  return PositionalAnalogueScanning(input, options, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:Eq:c:z:x:p");

  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage();
  }

  const int verbose = cl.option_present('v');

  Options options;

  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage();
  }

  if (cl.empty()) {
    cerr << "Insufficent arguments\n";
    Usage();
  }

  IWString_and_File_Descriptor output(1);

  if (! PositionalAnalogueScanning(cl[0], options, output)) {
    return 1;
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace pos

int 
main(int argc, char ** argv) {
  int rc = positional_analogue_scanning::Main(argc, argv);

  return rc;
}
