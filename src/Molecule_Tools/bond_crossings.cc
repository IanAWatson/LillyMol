// Score crossed bonds in existing 2D depictions. Coordinate generation is not needed.
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include "Depict/bond_crossings.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"

namespace {
using std::cerr;

void Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "Score crossed bonds in existing 2D molecular depictions\n"
       << "bond_crossings [options] molecules.sdf\n"
       << " -l          retain the largest fragment, preserving its coordinates\n"
       << " -z i        skip molecules whose highest coordinate dimensionality is not 2\n"
       << " -d          include degenerate/overlapping bonds, near misses, atoms on ring\n"
       << "             bonds, coincident atoms, and closest nonbonded distance\n"
       << " -t <f>      near-miss distance as a fraction of median bond length (default 0.1; -d)\n"
       << " -i <type>   input type; normally inferred from suffix\n"
       << " -A ...      standard aromaticity options\n"
       << " -v          verbose statistics\n";
  // clang-format on
  ::exit(rc);
}

class Options {
 private:
  bool _largest_fragment = false;
  bool _skip_not_2d = false;
  int _verbose = 0;
  int _read = 0;
  int _written = 0;
  int _skipped = 0;
  bond_crossings::Options _analysis;
  Fraction_as_String _numbers;

  void Number(double value, IWString_and_File_Descriptor& output) {
    output << ' ';
    // Match other descriptor tools: integers stay integers, other values use
    // the stable fixed-decimal formatter (no large lookup table).
    if (std::isfinite(value) && value == std::round(value) &&
        value >= 0 && value <= std::numeric_limits<int>::max()) {
      output << static_cast<int>(value);
    } else {
      _numbers.append_number(output, static_cast<float>(value));
    }
  }

 public:
  bool Initialise(Command_Line& cl) {
    _verbose = cl.option_count('v');
    _largest_fragment = cl.option_present('l');
    _analysis.set_diagnostics(cl.option_present('d'));
    if (cl.option_present('z')) {
      const_IWSubstring z;
      for (int i = 0; cl.value('z', z, i); ++i) {
        if (z != "i") {
          cerr << "Unrecognised -z qualifier '" << z << "'\n";
          return false;
        }
      }
      _skip_not_2d = true;
    }
    if (cl.option_present('t')) {
      double fraction;
      if (!cl.option_present('d') || !cl.value('t', fraction) ||
          !std::isfinite(fraction) || fraction < 0) {
        cerr << "-t requires -d and a finite non-negative fraction\n";
        return false;
      }
      _analysis.set_near_miss_fraction(fraction);
    }
    return _numbers.initialise_without_lookup(4);
  }

  void Header(IWString_and_File_Descriptor& output) const {
    output << "Name bc_natoms bc_nrs bc_nratoms bc_crossings bc_ring_ring bc_ring_chain"
           << " bc_inter_system bc_crossed_bonds bc_nrs_crossed bc_max_crossings bc_score"
           << " bc_score_per_atom bc_score_per_rs bc_score_per_ring_atom bc_worst_density"
           << " bc_min_2d bc_min_2d_per_atom bc_max_2d bc_max_2d_per_atom"
           << " bc_topological_diameter bc_max_2d_per_topology";
    if (_analysis.diagnostics()) {
      output << " bc_degenerate_bonds bc_overlaps bc_near_misses bc_atoms_on_ring_bonds"
             << " bc_coincident_atoms bc_min_nonbonded_2d";
    }
    output << '\n';
  }

  bool Process(Molecule& mol, IWString_and_File_Descriptor& output) {
    ++_read;
    if (_largest_fragment && mol.number_fragments() > 1) {
      mol.reduce_to_largest_fragment_carefully();
    }
    bond_crossings::Result result;
    const auto status = bond_crossings::Analyse(mol, _analysis, result);
    if (status == bond_crossings::Status::kNot2D && _skip_not_2d) {
      ++_skipped;
      return true;
    }
    if (status != bond_crossings::Status::kOk) {
      cerr << "bond_crossings: '" << mol.name() << "': "
           << bond_crossings::StatusMessage(status) << '\n';
      return false;
    }
    IWString id;
    if (!mol.name().word(0, id)) {
      id << "molecule_" << _read;
    }
    output << id;
    Number(result.atoms(), output);
    Number(result.systems().size(), output);
    Number(result.ring_atoms(), output);
    Number(result.crossings(), output);
    Number(result.ring_ring(), output);
    Number(result.ring_chain(), output);
    Number(result.inter_system(), output);
    Number(result.crossed_bonds(), output);
    Number(result.systems_crossed(), output);
    Number(result.max_crossings(), output);
    Number(result.score(), output);
    Number(result.score_per_atom(), output);
    Number(result.score_per_system(), output);
    Number(result.score_per_ring_atom(), output);
    Number(result.worst_density(), output);
    Number(result.min_distance(), output);
    Number(result.min_distance_per_atom(), output);
    Number(result.max_distance(), output);
    Number(result.max_distance_per_atom(), output);
    Number(result.topological_diameter(), output);
    Number(result.extent_per_topology(), output);
    if (_analysis.diagnostics()) {
      Number(result.degenerate_bonds(), output);
      Number(result.overlaps(), output);
      Number(result.near_misses(), output);
      Number(result.atoms_on_ring_bonds(), output);
      Number(result.coincident_atoms(), output);
      Number(result.min_nonbonded(), output);
    }
    output << '\n';
    ++_written;
    return output.write_if_buffer_holds_more_than(32768);
  }

  void Report() const {
    if (_verbose || _skipped) {
      cerr << "Read " << _read << " molecules, wrote " << _written
           << ", skipped " << _skipped << " with coordinate dimensionality other than 2\n";
    }
  }
};

int Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vli:A:z:dt:");
  if (cl.unrecognised_options_encountered() || cl.empty()) {
    Usage(1);
  }
  if (cl.option_present('A')) {
    if (!process_standard_aromaticity_options(cl, cl.option_count('v'), 'A')) {
      return 1;
    }
  } else {
    set_global_aromaticity_type(Daylight);
  }
  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      return 1;
    }
  } else if (cl.size() == 1 && cl[0] == IWString("-")) {
    input_type = FILE_TYPE_SDF;
  } else if (!all_files_recognised_by_suffix(cl)) {
    cerr << "Cannot infer input type; use -i\n";
    return 1;
  }
  Options options;
  if (!options.Initialise(cl)) {
    return 1;
  }
  IWString_and_File_Descriptor output(1);
  options.Header(output);
  for (const char* fname : cl) {
    const FileType type = input_type == FILE_TYPE_INVALID ? discern_file_type_from_name(fname) : input_type;
    data_source_and_type<Molecule> input(type, fname);
    if (!input.good()) {
      cerr << "Cannot open '" << fname << "'\n";
      return 1;
    }
    while (Molecule* molecule = input.next_molecule()) {
      std::unique_ptr<Molecule> mol(molecule);
      if (!options.Process(*mol, output)) {
        return 1;
      }
    }
    if (input.stopped_because_of_error()) {
      return 1;
    }
  }
  options.Report();
  return output.flush() ? 0 : 1;
}
}  // namespace

int main(int argc, char** argv) { return Main(argc, argv); }
