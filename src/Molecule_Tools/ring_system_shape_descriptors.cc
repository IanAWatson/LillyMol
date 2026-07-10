#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>

#include "Foundational/cmdline/cmdline.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/molecule_preprocessing.h"

#include "Molecule_Tools/ring_system_shape.h"

namespace ring_system_shape_descriptors {

using std::cerr;

using molecule_processing::MoleculePreprocessing;
using ring_system_shape::AnalyseRingSystemShape;
using ring_system_shape::RingSystemShape;
using ring_system_shape::RingSystemShapeClass;

void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on

  // clang-format off
  cerr << R"(Computes simple ring-system rod-like shape descriptors.

For each molecule, ring systems are labelled and each ring system is classified
by the number and placement of meaningful exit points.

Output columns are:
  Name nringsys nringsys_terminal nringsys_applicable nringsys_rodlike nringsys_not_rodlike nringsys_multisub nringsys_invalid rodlike_molecule rod_deficit_min rod_deficit_max rod_deficit_ave

Options:
 -x           ignore terminal single-atom ring substituents when larger substituents exist
 -o <sep>     output separator, recognised names include tab and space
 -i <type>    input type
 -g ...       chemical standardisation
 -l           reduce to largest fragment
 -c           remove chirality
 -A ...       aromaticity options
 -E ...       element options
 -v           verbose output
)";
  // clang-format on

  ::exit(rc);
}

struct RingSystemShapeDescriptors {
  int number_ring_systems = 0;
  int terminal = 0;
  int applicable = 0;
  int rod_like = 0;
  int not_rod_like = 0;
  int multi_substituted = 0;
  int invalid = 0;

  int rod_deficit_count = 0;
  int rod_deficit_min = 0;
  int rod_deficit_max = 0;
  int rod_deficit_sum = 0;

  void
  AddRodDeficit(int rod_deficit) {
    if (rod_deficit < 0) {
      return;
    }

    if (rod_deficit_count == 0) {
      rod_deficit_min = rod_deficit;
      rod_deficit_max = rod_deficit;
    } else {
      if (rod_deficit < rod_deficit_min) {
        rod_deficit_min = rod_deficit;
      }
      if (rod_deficit > rod_deficit_max) {
        rod_deficit_max = rod_deficit;
      }
    }

    ++rod_deficit_count;
    rod_deficit_sum += rod_deficit;
  }

  float
  AverageRodDeficit() const {
    if (rod_deficit_count == 0) {
      return 0.0f;
    }

    return static_cast<float>(rod_deficit_sum) / static_cast<float>(rod_deficit_count);
  }

  int
  RodLikeMolecule() const {
    if (number_ring_systems == 0) {
      return 0;
    }

    return not_rod_like == 0 && multi_substituted == 0 && invalid == 0;
  }
};

class Options {
 private:
  int _verbose = 0;

  MoleculePreprocessing _preprocessing;

  char _output_separator = ' ';

  int _remove_terminal_single_atom_substituents = 0;

  uint64_t _molecules_read = 0;
  uint64_t _atoms_removed = 0;

 public:
  Options();

  int Initialise(Command_Line& cl);

  int
  verbose() const {
    return _verbose;
  }

  int Preprocess(Molecule& m);

  int RemoveTerminalSingleAtomSubstituents(Molecule& m);

  int WriteHeader(IWString_and_File_Descriptor& output) const;

  int Process(Molecule& m, IWString_and_File_Descriptor& output);

  int Report(std::ostream& output) const;
};

Options::Options() {
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (!_preprocessing.Initialise(cl)) {
    cerr << "Options::Initialise:cannot initialise preprocessing\n";
    return 0;
  }

  if (cl.option_present('o')) {
    IWString sep = cl.string_value('o');
    if (!char_name_to_char(sep)) {
      cerr << "Invalid output separator '" << sep << "'\n";
      return 0;
    }
    _output_separator = sep[0];
  }

  if (cl.option_present('x')) {
    _remove_terminal_single_atom_substituents = 1;
  }

  return 1;
}

int
Options::RemoveTerminalSingleAtomSubstituents(Molecule& m) {
  const int matoms = m.natoms();
  if (matoms == 0) {
    return 1;
  }

  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(matoms);
  const int number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());
  if (number_ring_systems == 0) {
    return 1;
  }

  std::vector<Set_of_Atoms> single_atom_substituents(number_ring_systems + 1);
  std::vector<int> has_larger_attachment(number_ring_systems + 1, 0);

  for (int i = 0; i < matoms; ++i) {
    const int ring_system_id = ring_system_membership[i];
    if (ring_system_id <= 0) {
      continue;
    }

    const Atom& atom = m[i];
    for (const Bond* bond : atom) {
      const atom_number_t other = bond->other(i);
      if (ring_system_membership[other] == ring_system_id) {
        continue;
      }

      if (bond->is_single_bond() && m.ncon(other) == 1) {
        single_atom_substituents[ring_system_id].add_if_not_already_present(other);
      } else if (bond->is_single_bond() || m.ncon(other) > 1) {
        has_larger_attachment[ring_system_id] = 1;
      }
    }
  }

  Set_of_Atoms atoms_to_remove;
  for (int ring_system_id = 1; ring_system_id <= number_ring_systems; ++ring_system_id) {
    if (!has_larger_attachment[ring_system_id]) {
      continue;
    }

    for (atom_number_t atom : single_atom_substituents[ring_system_id]) {
      atoms_to_remove.add_if_not_already_present(atom);
    }
  }

  if (atoms_to_remove.empty()) {
    return 1;
  }

  const int atoms_removed = m.remove_atoms(atoms_to_remove);
  _atoms_removed += atoms_removed;

  return 1;
}

int
Options::Preprocess(Molecule& m) {
  if (m.empty()) {
    return 0;
  }

  if (_preprocessing.active()) {
    _preprocessing.Process(m);
  }

  if (_remove_terminal_single_atom_substituents) {
    RemoveTerminalSingleAtomSubstituents(m);
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Processed " << _molecules_read << " molecules\n";
  if (_remove_terminal_single_atom_substituents) {
    output << "Removed " << _atoms_removed << " terminal single atom substituents\n";
  }
  return 1;
}

RingSystemShapeDescriptors
ComputeRingSystemShapeDescriptors(Molecule& m) {
  RingSystemShapeDescriptors result;

  const int matoms = m.natoms();
  if (matoms == 0) {
    return result;
  }

  std::unique_ptr<int[]> ring_system_membership = std::make_unique<int[]>(matoms);
  result.number_ring_systems =
      m.label_atoms_by_ring_system_including_spiro_fused(ring_system_membership.get());

  if (result.number_ring_systems == 0) {
    return result;
  }

  for (int ring_system_id = 1; ring_system_id <= result.number_ring_systems;
       ++ring_system_id) {
    RingSystemShape ring_system_shape;
    if (!AnalyseRingSystemShape(m, ring_system_membership.get(), ring_system_id,
                                ring_system_shape)) {
      ++result.invalid;
      continue;
    }

    switch (ring_system_shape.shape_class) {
      case RingSystemShapeClass::kNotApplicable:
        ++result.terminal;
        break;
      case RingSystemShapeClass::kRodLike:
        ++result.applicable;
        ++result.rod_like;
        result.AddRodDeficit(ring_system_shape.rod_deficit);
        break;
      case RingSystemShapeClass::kNotRodLike:
        ++result.applicable;
        ++result.not_rod_like;
        result.AddRodDeficit(ring_system_shape.rod_deficit);
        break;
      case RingSystemShapeClass::kMultiSubstituted:
        ++result.multi_substituted;
        break;
      case RingSystemShapeClass::kInvalid:
        ++result.invalid;
        break;
    }
  }

  return result;
}

int
Options::WriteHeader(IWString_and_File_Descriptor& output) const {
  output << "Name" << _output_separator << "nringsys" << _output_separator
         << "nringsys_terminal" << _output_separator << "nringsys_applicable"
         << _output_separator << "nringsys_rodlike" << _output_separator
         << "nringsys_not_rodlike" << _output_separator << "nringsys_multisub"
         << _output_separator << "nringsys_invalid" << _output_separator
         << "rodlike_molecule" << _output_separator << "rod_deficit_min"
         << _output_separator << "rod_deficit_max" << _output_separator
         << "rod_deficit_ave\n";
  return 1;
}

int
Options::Process(Molecule& m, IWString_and_File_Descriptor& output) {
  ++_molecules_read;

  const RingSystemShapeDescriptors descriptors = ComputeRingSystemShapeDescriptors(m);

  output << m.name() << _output_separator << descriptors.number_ring_systems
         << _output_separator << descriptors.terminal << _output_separator
         << descriptors.applicable << _output_separator << descriptors.rod_like
         << _output_separator << descriptors.not_rod_like << _output_separator
         << descriptors.multi_substituted << _output_separator << descriptors.invalid
         << _output_separator << descriptors.RodLikeMolecule() << _output_separator
         << descriptors.rod_deficit_min << _output_separator
         << descriptors.rod_deficit_max << _output_separator
         << descriptors.AverageRodDeficit() << '\n';

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

int
RingSystemShapeDescriptors(Options& options, Molecule& m,
                           IWString_and_File_Descriptor& output) {
  return options.Process(m, output);
}

int
RingSystemShapeDescriptors(Options& options, data_source_and_type<Molecule>& input,
                           IWString_and_File_Descriptor& output) {
  Molecule* m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (!options.Preprocess(*m)) {
      continue;
    }

    if (!RingSystemShapeDescriptors(options, *m, output)) {
      return 0;
    }
  }

  return 1;
}

int
RingSystemShapeDescriptors(Options& options, const char* fname, FileType input_type,
                           IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "RingSystemShapeDescriptors:cannot open '" << fname << "'\n";
    return 0;
  }

  if (options.verbose() > 1) {
    input.set_verbose(1);
  }

  return RingSystemShapeDescriptors(options, input, output);
}

int
RingSystemShapeDescriptors(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:lcg:i:o:x");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose)) {
    Usage(1);
  }
  if (!process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process elements\n";
    Usage(1);
  }

  Options options;
  if (!options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    return 1;
  }

  FileType input_type = FILE_TYPE_INVALID;

  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (cl.number_elements() == 1 && strcmp(cl[0], "-") == 0) {
    input_type = FILE_TYPE_SMI;
  } else if (!all_files_recognised_by_suffix(cl)) {
    return 1;
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  options.WriteHeader(output);

  for (const char* fname : cl) {
    if (!RingSystemShapeDescriptors(options, fname, input_type, output)) {
      cerr << "RingSystemShapeDescriptors::fatal error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace ring_system_shape_descriptors

int
main(int argc, char** argv) {
  return ring_system_shape_descriptors::RingSystemShapeDescriptors(argc, argv);
}
