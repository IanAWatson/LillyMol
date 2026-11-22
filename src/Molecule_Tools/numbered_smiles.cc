/*
  Frequently I need to know the atom numbers of the atoms in a molecule.
  Produce isotopic smiles
*/

#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <numeric>
#include <vector>
#include <iostream>
#include <memory>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/output.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/smiles.h"

using std::cerr;

static int verbose = 0;
const char* prog_name = nullptr;

static IWString output_file_name;

/*
  The -n switch governs the numbering assigned.
*/

enum class NumberingType {
  kAtomNumber,
  kSymmetryClass,
  kCanonicalRank,
  kRingSystem,
  kAtomType,
  kAtomicNumber
};

// #define NUMBER_BY_ATOM_NUMBER 1
// #define NUMBER_BY_SYMMETRY_CLASS 2
// #define NUMBER_BY_CANONICAL_RANK 3

static NumberingType numbering_type = NumberingType::kAtomNumber;

static int atom_number_offset = 0;

static int number_step = 1;

static int overwrite_existing_isotopic_information = 0;

static int sort_molecule_by_numbering = 0;

static int reduce_to_largest_fragment = 0;

static int make_implicit_hydrogens_explicit = 0;

static int compute_bond_symmetry_small = 0;
static int compute_bond_symmetry_large = 0;

static int ring_systems_span_spiro_fusions = 0;

static int apply_isotopic_labels = 1;

static Atom_Typing_Specification atom_typing;

static int remove_chiral_centres = 0;

static void
usage(int rc = 0)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "  -n atom        number by atom number (the default)\n";
  cerr << "  -n symm        number by symmetry class\n";
  cerr << "  -n canon       number by canonical order\n";
  cerr << "  -n ring        number by ring system\n";
  cerr << "  -n z           number by atomic number\n";
  cerr << "  -n <number>    number atoms starting at a non-negative <number>\n";
  cerr << "  -e <step>      atom-number increment; multiplier for other numbering modes (default 1)\n";
  cerr << "  -I             overwrite any existing isotopic information\n";
  cerr << "  -P <atype>     number by atom typing (no -n option needed)\n";
  cerr << "  -r             stably sort atoms by their assigned labels\n";
  cerr << "  -H             make all implicit Hydrogens explicit\n";
  cerr << "  -b ...         write bond symmetry information to stderr\n";
  cerr << "  -m             overwrite atom maps with numbering; preserve isotopes\n";
  cerr << "  -c             remove all chiral centres\n";
  cerr << "  -l             reduce to the largest fragment before processing\n";
  cerr << "  -Y ...         more options, enter '-Y help' for info\n";
  cerr << "  -S <stem>      specify output file name stem\n";
  cerr << "  -K ...         standard smiles control options\n";
  cerr << "  -i <type>      specify input type\n";
  display_standard_aromaticity_options(cerr);
  cerr << "  -v             verbose output\n";
  // clang-format on

  exit(rc);
}

static int
write_bond_symmetries(Molecule& m, const int* bs, std::ostream& output)
{
  output << m.name() << '\n';

  int nb = m.nedges();

  for (int i = 0; i < nb; i++) {
    const Bond* b = m.bondi(i);

    output << b->a1() << " (" << m.symmetry_class(b->a1()) << ") " << b->a2() << " ("
           << m.symmetry_class(b->a2()) << ") " << bs[i] << '\n';
  }

  return output.good();
}

static int
do_compute_bond_symmetry_small(Molecule& m, std::ostream& output)
{
  int* bs = new int[m.nedges()];
  std::unique_ptr<int[]> free_bs(bs);

  m.bond_symmetry_class_small_memory(bs);

  return write_bond_symmetries(m, bs, output);
}

static int
do_compute_bond_symmetry_large(Molecule& m, std::ostream& output)
{
  int* bs = new int[m.nedges()];
  std::unique_ptr<int[]> free_bs(bs);

  m.bond_symmetry_class_large_memory(bs);

  return write_bond_symmetries(m, bs, output);
}

// Build an old-to-new permutation from labels, without changing the labels.
// Stable ordering retains the input order of atoms with equal labels.
static int
SortByNumbering(Molecule& m) {
  const int matoms = m.natoms();
  std::vector<atom_number_t> order(matoms);
  std::iota(order.begin(), order.end(), 0);
  if (apply_isotopic_labels) {
    std::stable_sort(order.begin(), order.end(), [&m](atom_number_t a, atom_number_t b) {
      return m.isotope(a) < m.isotope(b);
    });
  } else {
    std::stable_sort(order.begin(), order.end(), [&m](atom_number_t a, atom_number_t b) {
      return m.atom_map_number(a) < m.atom_map_number(b);
    });
  }
  auto destination = std::make_unique<int[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    destination[order[i]] = i;
  }
  return m.renumber_atoms(destination.get());
}

static int
ApplyNumbering(Molecule& m, const isotope_t* labels) {
  // Atom maps are independent of isotopes and always replace existing maps.
  if (!apply_isotopic_labels) {
    for (int i = 0; i < m.natoms(); ++i) {
      m.set_atom_map_number(i, static_cast<int>(labels[i]));
    }
    return 1;
  }
  if (overwrite_existing_isotopic_information) {
    return m.set_isotopes(labels);
  }
  for (int i = 0; i < m.natoms(); ++i) {
    if (m.isotope(i) == 0) {
      m.set_isotope(i, labels[i]);
    }
  }
  return 1;
}

// Compute in a wider type so large starts, steps, and atom types cannot wrap.
static int
ScaleNumbering(const Molecule& m, isotope_t* labels) {
  const uint64_t max_label = apply_isotopic_labels
      ? std::numeric_limits<isotope_t>::max() : std::numeric_limits<int>::max();
  for (int i = 0; i < m.natoms(); ++i) {
    const uint64_t value = numbering_type == NumberingType::kAtomNumber
        ? static_cast<uint64_t>(atom_number_offset) + static_cast<uint64_t>(i) * number_step
        : static_cast<uint64_t>(labels[i]) * number_step;
    if (value > max_label) {
      cerr << "Numbering exceeds the supported label range for '" << m.name() << "'\n";
      return 0;
    }
    labels[i] = static_cast<isotope_t>(value);
  }
  return 1;
}

static int
fetch_symmetry_class(Molecule& m, isotope_t* isotope)
{
  std::copy_n(m.symmetry_classes(), m.natoms(), isotope);

  return 1;
}

// The canonical rank is an integer, so a temporary is needed.
static int
fetch_canonical_rank(Molecule& m, isotope_t* isotope)
{
  const int matoms = m.natoms();
  std::unique_ptr<int[]> as_int = std::make_unique<int[]>(matoms);

  m.canonical_ranks(as_int.get());

  std::copy_n(as_int.get(), matoms, isotope);

  return 1;
}

static int
AddDoublyBondedExocyclic(Molecule& m,
                         isotope_t* isotope) {

  // Add doubly bonded exocyclic bonds
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    if (isotope[i] == 0) {
      continue;
    }
    const Atom* a = m.atomi(i);
    if (a->ncon() == 2 || a->nbonds() == a->ncon()) {
      continue;
    }
    for (const Bond* b : *a) {
      if (!b->is_double_bond()) {
        continue;
      }
      const atom_number_t j = b->other(i);
      if (isotope[j] == isotope[i] || m.atomic_number(j) == 6) {
        continue;
      }
      // We have an exocyclic double bond.
      // If we want to restrict to just =O or =S, uncomment this.
      // if (m.ncon(j) != 1) {
      //   continue;
      // }
      isotope[j] = isotope[i];
    }
  }

  return 1;
}

static int
fetch_ring_system_number(Molecule& m, isotope_t* isotope)
{
  const int nrings = m.nrings();
  if (nrings == 0) {
    return 1;
  }

  const int matoms = m.natoms();
  std::unique_ptr<int[]> label(new_int(matoms));

  if (ring_systems_span_spiro_fusions) {
    m.label_atoms_by_ring_system_including_spiro_fused(label.get());
  } else {
    m.label_atoms_by_ring_system(label.get());
  }
  std::copy_n(label.get(), matoms, isotope);

  return AddDoublyBondedExocyclic(m, isotope);
}

int
fetch_atomic_number(const Molecule& m, isotope_t* isotope)
{
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    isotope[i] = m.atomic_number(i);
  }

  return matoms;
}

// Needs a temporary array to hold the right int type for
// atom typing.
int
fetch_atom_types(Molecule& m, isotope_t* isotope) {
  const int matoms = m.natoms();
  std::unique_ptr<uint32_t[]> atype = std::make_unique<uint32_t[]>(matoms);

  atom_typing.assign_atom_types(m, atype.get());

  std::copy_n(atype.get(), matoms, isotope);

  return 1;
}

static int
numbered_smiles(Molecule& m,
                isotope_t* isotope,
                Molecule_Output_Object& output)
{
  switch (numbering_type) {
    case NumberingType::kAtomNumber:
      // ScaleNumbering constructs the arithmetic progression.
      break;
    case NumberingType::kSymmetryClass:
      fetch_symmetry_class(m, isotope);
      break;
    case NumberingType::kCanonicalRank:
      fetch_canonical_rank(m, isotope);
      break;
    case NumberingType::kRingSystem:
      fetch_ring_system_number(m, isotope);
      break;
    case NumberingType::kAtomType:
      fetch_atom_types(m, isotope);
      break;
    case NumberingType::kAtomicNumber:
      fetch_atomic_number(m, isotope);
      break;
  }

  if (!ScaleNumbering(m, isotope) || !ApplyNumbering(m, isotope)) {
    return 0;
  }
  if (sort_molecule_by_numbering && !SortByNumbering(m)) {
    return 0;
  }

  if (compute_bond_symmetry_small) {
    do_compute_bond_symmetry_small(m, cerr);
  }
  if (compute_bond_symmetry_large) {
    do_compute_bond_symmetry_large(m, cerr);
  }

  return output.write(m);
}

static int
numbered_smiles(Molecule& m, Molecule_Output_Object& output)
{
  assert(m.ok());

  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (remove_chiral_centres) {
    m.remove_all_chiral_centres();
  }

  int matoms = m.natoms();

  if (0 == matoms) {
    cerr << "numbered_smiles: nostructure encountered '" << m.name() << "'\n";
    return 1;
  }

  if (make_implicit_hydrogens_explicit) {
    m.make_implicit_hydrogens_explicit();
    matoms = m.natoms();
  }

  std::unique_ptr<isotope_t[]> isotope = std::make_unique<isotope_t[]>(matoms);

  return numbered_smiles(m, isotope.get(), output);
}

static int
numbered_smiles(data_source_and_type<Molecule>& input, Molecule_Output_Object& output)
{
  Molecule* m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    if (verbose) {
      cerr << "Processing '" << input.molecules_read() << " " << m->name() << "'\n";
    }

    if (!numbered_smiles(*m, output)) {
      return 0;
    }
  }

  if (verbose) {
    cerr << "Processed " << input.molecules_read() << " molecules\n";
  }

  return !input.stopped_because_of_error() && output.good();
}

// Check buffered writes while failures can still affect the exit status.
static int
FinishOutputs(Molecule_Output_Object& output) {
  int ok = output.empty() || output.number_elements() == output.number_output_types();
  for (auto* stream : output) {
    stream->flush();
    if (!stream->good()) {
      ok = 0;
    }
  }
  std::cout.flush();
  if (!std::cout.good()) {
    ok = 0;
  }
  if (!ok) {
    cerr << "numbered_smiles: error writing output\n";
  }
  return ok;
}

// Remove only the filename's final extension, retaining directory and other dots.
static IWString
DefaultOutputStem(const char* fname) {
  std::filesystem::path path(fname);
  path.replace_extension();
  path += "_num";
  return IWString(path.string());
}


// Check all ordinary output names before any output stream is opened.
static int
CheckOutputNames(const Command_Line& cl, const Molecule_Output_Object& output) {
  if (output_file_name == "-") {
    return 1;
  }
  std::vector<IWString> destinations;
  for (const char* input : cl) {
    const IWString stem = output_file_name.empty() ? DefaultOutputStem(input) : output_file_name;
    const int keep_suffix = output_file_name.empty();
    for (FileType type : output.output_types()) {
      IWString fname;
      if (!create_file_with_appropriate_name(stem, fname, type, keep_suffix)) {
        return 0;
      }
      for (const char* source : cl) {
        if (std::string_view(source) != "-") {
          const auto same = iwmisc::SameFile(fname, IWString(source));
          if (!same) {
            cerr << "Cannot compare file paths\n";
            return 0;
          }
          if (*same) {
            cerr << "Output '" << fname << "' would overwrite input '" << source << "'\n";
            return 0;
          }
        }
      }
      for (const IWString& previous : destinations) {
        const auto same = iwmisc::SameFile(fname, previous);
        if (!same) {
          cerr << "Cannot compare file paths\n";
          return 0;
        }
        if (*same) {
          cerr << "Output filenames '" << fname << "' and '" << previous << "' coincide\n";
          return 0;
        }
      }
      destinations.push_back(fname);
    }
    if (!output_file_name.empty()) {
      break;  // All inputs share the explicitly specified output.
    }
  }
  return 1;
}

static int
numbered_smiles(const char* fname, FileType input_type, Molecule_Output_Object& output)
{
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
    if (input_type == FILE_TYPE_INVALID) {
      cerr << "Cannot determine input type for '" << fname << "'\n";
      return 0;
    }
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  static int first_call = 1;

  if (output_file_name.empty()) {
    // A new stem closes the preceding file; check its writes first.
    if (!FinishOutputs(output)) {
      return 0;
    }
    const IWString ofname = DefaultOutputStem(fname);
    // The stem may contain periods which are part of the original filename.
    if (!output.new_stem(ofname, 1)) {
      cerr << "Cannot set new stem '" << ofname << "'\n";
      return 0;
    }
  } else if (first_call) {
    if (!output.new_stem(output_file_name)) {
      cerr << "Cannot set stem for '" << output_file_name << "'\n";
      return 0;
    }
    first_call = 0;
  }

  return numbered_smiles(input, output);
}

static void
DisplayDashYOptions(int rc) {
  std::cerr << " -Y spiro          ring systems span spiro fusions\n";

  std::exit(rc);
}

static int
numbered_smiles(int argc, char** argv)
{
  Command_Line cl(argc, argv, "IvS:i:o:n:E:A:rHb:e:mK:P:Y:cl");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(6);
  }

  DisplayDashHelpIfRequested(cl, 'Y', DisplayDashYOptions);

  verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;
    if (verbose) {
      cerr << "Will reduce molecules to their largest fragment\n";
    }
  }

  if (!process_elements(cl, verbose)) {
    usage(1);
  }

  if (!process_standard_smiles_options(cl, verbose, 'K')) {
    cerr << "Cannot process smiles options\n";
    return 1;
  }

  if (!process_standard_aromaticity_options(cl, verbose)) {
    cerr << "Cannot process standard aromaticity options\n";
    usage(3);
  }

  if (cl.option_present('I')) {
    overwrite_existing_isotopic_information = 1;
    if (verbose) {
      cerr << "Existing isotopic information will be overwritten\n";
    }
  }

  if (cl.option_present('P')) {
    const IWString p = cl.string_value('P');
    if (! atom_typing.build(p)) {
      cerr << "Cannot build atom typing '" << p << "'\n";
      return 1;
    }
    numbering_type = NumberingType::kAtomType;
  }

  if (cl.option_present('c')) {
    remove_chiral_centres = 1;
    if (verbose) {
      cerr << "Will remove chiral centres from molecules\n";
    }
  }

  if (cl.option_present('Y')) {
    const_IWSubstring y;
    for (int i = 0; cl.value('Y', y, i); ++i) {
      if (y == "spiro") {
        ring_systems_span_spiro_fusions = 1;
      } else if (y == "help") {
        DisplayDashYOptions(0);
      } else {
        cerr << "Unrecognised -Y qualifier '" << y << "'\n";
        DisplayDashYOptions(1);
      }
    }
  }

  if (cl.option_present('n')) {
    IWString tmp;
    int i = 0;
    while (cl.value('n', tmp, i++)) {
      if (tmp.starts_with("atom")) {
        numbering_type = NumberingType::kAtomNumber;
        if (verbose) {
          cerr << "Atoms will be numbered according to their atom number\n";
        }
      } else if (tmp.starts_with("symm")) {
        numbering_type = NumberingType::kSymmetryClass;
        if (verbose) {
          cerr << "Atoms will be numbered by symmetry class\n";
        }
      } else if (tmp.starts_with("canon")) {
        numbering_type = NumberingType::kCanonicalRank;
        if (verbose) {
          cerr << "Atoms will be numbered by canonical number\n";
        }
      } else if (tmp.starts_with("ring")) {
        numbering_type = NumberingType::kRingSystem;
        if (verbose) {
          cerr << "Atoms numbered by ring system\n";
        }
      } else if (tmp.starts_with("z")) {
        numbering_type = NumberingType::kAtomicNumber;
        if (verbose) {
          cerr << "Atoms numbered by atomic number\n";
        }
      } else if (tmp.is_int(atom_number_offset)) {
        if (atom_number_offset < 0) {
          cerr << "The starting atom number (-n) must be non-negative\n";
          return 1;
        }
        numbering_type = NumberingType::kAtomNumber;
        if (verbose) {
          cerr << "Atoms numbered by atom number, offset = " << atom_number_offset
               << '\n';
        }
      } else {
        cerr << "Unrecognised numbering type '" << tmp << "'\n";
        usage(9);
      }
    }
  }

  if (cl.option_present('e')) {
    if (!cl.value('e', number_step) || number_step < 1) {
      cerr << "The number step value (-e) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Gap between numbers set to " << number_step << '\n';
    }
  }

  if (cl.option_present('b')) {
    const_IWSubstring b = cl.string_value('b');

    if ("small" == b) {
      compute_bond_symmetry_small = 1;
      if (verbose) {
        cerr << "Will compute bond symmetry with the small memory model\n";
      }
    } else if ("large" == b) {
      compute_bond_symmetry_large = 1;
      if (verbose) {
        cerr << "Will compute bond symmetry with the large memory model\n";
      }
    } else {
      cerr << "Unrecognised -b qualifier '" << b << "'\n";
      usage(5);
    }
  }

  if (cl.option_present('m')) {
    apply_isotopic_labels = 0;
    if (verbose) {
      cerr << "Will apply atom map numbers\n";
    }
  }

  // By default we echo any chiral info present in the input

  set_include_chiral_info_in_smiles(1);

  if (cl.option_present('H')) {
    make_implicit_hydrogens_explicit = 1;

    if (verbose) {
      cerr << "Implicit hydrogens will be made explicit\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot determine input types\n";
    return 32;
  } else  // good
  {
  }

  Molecule_Output_Object output;
  if (cl.option_present('o')) {
    if (!output.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output types (-o)\n";
      return 77;
    }
  } else {
    output.add_output_type(FILE_TYPE_SMI);  // our default
  }

  if (cl.option_present('S')) {
    cl.value('S', output_file_name);
    if (0 == output_file_name.length()) {
      cerr << "The -S option must be followed by a stem\n";
      usage(10);
    }

    if (verbose) {
      cerr << "Stem for output is '" << output_file_name << "'\n";
    }
  }

  if (cl.option_present('r')) {
    sort_molecule_by_numbering = 1;

    if (verbose) {
      cerr << "Molecule will be sorted by numbering scheme\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(9);
  }

  if (!CheckOutputNames(cl, output)) {
    return 1;
  }

  int rc = 0;
  for (const char* fname : cl) {
    if (!numbered_smiles(fname, input_type, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      rc = 1;
      break;
    }
  }

  if (!FinishOutputs(output)) {
    rc = 1;
  }
  return rc;
}

int
main(int argc, char** argv)
{
  prog_name = argv[0];

  int rc = numbered_smiles(argc, argv);

  return rc;
}
