/*
  Computes ring substitution fingerprint
*/

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <memory>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/ring_substitution.h"

namespace ring_substitution_main {

using std::cerr;

using ring_substitution::RingSubstitutionGenerator;

const char* prog_name = nullptr;

static int molecules_read = 0;

static int molecules_with_no_rings = 0;

static IWString fingerprint_tag("NCRS<");

static Chemical_Standardisation chemical_standardisation;

static int reduce_to_largest_fragment = 0;

static int verbose = 0;

static int flatten_fingerprints_to_01 = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");

static int work_as_filter = 0;

static Accumulator_Int<int> acc_nset;

static void
usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << "  -J <tag>       tag for fingerprints\n";
  cerr << "  -P <bonds>    maximum path length around edge of ring system\n";
  cerr << "  -M ...        options for specifying just what kind of information to generate\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -g ...        chemical standardisation options\n";
  cerr << "  -i <type>     input specification\n";
  cerr << "  -A ...        aromaticity specifications, enter '-A help' for info\n";
  cerr << "  -E ...        element specifications, enter '-E help' for info\n";
  cerr << "  -v            verbose output\n";
  // clang-format on

  exit(rc);
}

static void
preprocess(Molecule& m) {
  if (reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment();
  }

  if (chemical_standardisation.active()) {
    chemical_standardisation.process(m);
  }

  return;
}

static int
RingSubstitution(RingSubstitutionGenerator& gen,
                  Molecule& m, IWString_and_File_Descriptor& output) {
  preprocess(m);

#ifdef DEBUG_RING_SUBSTITUTION
  for (int i = 0; i < m.natoms(); i++) {
    cerr << "Atom " << i << " assigned type " << atype[i] << '\n';
  }
#endif

  Sparse_Fingerprint_Creator sfpc;

  if (m.nrings() > 0) {
    gen.Generate(m, sfpc);
  }

  if (!work_as_filter) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  if (flatten_fingerprints_to_01) {
    sfpc.flatten_to_01();
  }

  IWString tmp;
  sfpc.daylight_ascii_form_with_counts_encoded(fingerprint_tag, tmp);

  output << tmp << '\n';

  if (verbose > 2) {
    sfpc.debug_print(cerr);
  }

  if (!work_as_filter) {
    output << "|\n";
  }

  if (verbose) {
    acc_nset.extra(sfpc.nset());
  }

  if (verbose > 2) {
    cerr << molecules_read << " '" << m.name() << "' hits " << sfpc.nbits() << " bits\n";
  }

  return 1;
}

static int
RingSubstitution(RingSubstitutionGenerator& gen, data_source_and_type<Molecule>& input,
                  IWString_and_File_Descriptor& output) {
  Molecule* m;

  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);

    molecules_read++;

    if (!RingSubstitution(gen, *m, output)) {
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static int
RingSubstitution(RingSubstitutionGenerator& gen,
                  const char* fname, FileType input_type,
                  IWString_and_File_Descriptor& output) {
  if (FILE_TYPE_INVALID == input_type) {
    input_type = discern_file_type_from_name(fname);
    assert(FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose > 2) {
    input.set_verbose(1);
  }

  return RingSubstitution(gen, input, output);
}

static int
RingSubstitutionFilterRecord(RingSubstitutionGenerator& gen,
                                const const_IWSubstring& buffer,
                                IWString_and_File_Descriptor& output) {
  const_IWSubstring tmp(buffer);
  tmp.remove_leading_chars(smiles_tag.length());
  tmp.chop(1);

  Molecule m;
  if (!m.build_from_smiles(tmp)) {
    cerr << "Invalid smiles '" << tmp << "'\n";
    return 0;
  }

  return RingSubstitution(gen, m, output);
}

static int
RingSubstitutionFilter(RingSubstitutionGenerator& gen,
                         iwstring_data_source& input,
                         IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    output << buffer << "\n";

    if (!buffer.starts_with(smiles_tag)) {
      continue;
    }

    if (!RingSubstitutionFilterRecord(gen, buffer, output)) {
      cerr << "Fatal error, line " << input.lines_read() << '\n';
      return 0;
    }

    output.write_if_buffer_holds_more_than(4096);
  }

  return 1;
}

static void
display_dash_M_options(std::ostream& os) {
  os << "  -M posi    set bits to indicate positions on rings only\n";
  os << "  -M simp    set bits to indicate positions and simple atom types\n";
  os << "  -M full    set bits to indicate positions and complete atom types\n";
  os << "  -M sfb     set bits to indicate presence or absence of single features\n";
  os << "  -M 01      remove count information from fingerprints, just presence of "
        "absence\n";

  exit(2);
}

static int
RingSubstitutionFilter(RingSubstitutionGenerator& gen,
                        const char* fname, IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open filter input '" << fname << "'\n";
    return 0;
  }

  return RingSubstitutionFilter(gen, input, output);
}

static int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:E:i:J:lg:M:fP:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!process_standard_aromaticity_options(cl, verbose > 1)) {
    cerr << "Cannot process -A option\n";
    usage(11);
  }

  if (!process_elements(cl, verbose > 1, 'E')) {
    cerr << "Cannot initialise elements\n";
    usage(8);
  }

  if (cl.option_present('g')) {
    if (!chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g')) {
      cerr << "Cannot process chemical standardisation options (-g)\n";
      usage(32);
    }
  }

  if (cl.option_present('l')) {
    reduce_to_largest_fragment = 1;

    if (verbose) {
      cerr << "Will reduce to largest fragment\n";
    }
  }

  RingSubstitutionGenerator gen;

  if (cl.option_present('b')) {
    //  break_at_subsitutents = 1;
    if (verbose) {
      cerr << "Will break after encountering closest substituent\n";
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (!process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      usage(6);
    }
  } else if (cl.option_present('f'))
    ;
  else if (!all_files_recognised_by_suffix(cl)) {
    return 4;
  }

  if (cl.option_present('P')) {
    int max_path_length;
    if (!cl.value('P', max_path_length) || max_path_length < 2) {
      cerr << "The maximum path length (-m) must be a whole +ve number\n";
      usage(7);
    }

    if (verbose) {
      cerr << "Max path length " << max_path_length << '\n';
    }
    gen.set_max_path_length(max_path_length);
  }

  if (cl.option_present('M')) {
    int nset = 0;

    int i = 0;
    const_IWSubstring m;
    while (cl.value('M', m, i++)) {
      if (m.starts_with("posi")) {
        gen.set_positional_information_only(1);
        nset++;
        if (verbose) {
          cerr << "Only positional information will be recorded\n";
        }
      } else if (m.starts_with("simp")) {
        gen.set_simple_atom_types(1);
        nset++;
        if (verbose) {
          cerr << "Simple atom types for substitutent atoms\n";
        }
//    } else if ("full" == m) {
//      gen.set_full_atom_types(1);
//      nset++;
//      if (verbose) {
//        cerr << "Full atom typing\n";
//      }
      } else if ("sfb" == m) {
        gen.set_place_single_feature_bits(1);
        if (verbose) {
          cerr << "Will set bits according to presence or absence of features\n";
        }
      } else if ("01" == m) {
        flatten_fingerprints_to_01 = 1;
        if (verbose) {
          cerr << "Fingerprints will be flattened to just 01 forms\n";
        }
      } else {
        cerr << "Unrecognised -M qualifier '" << m << "'\n";
        display_dash_M_options(cerr);
        usage(4);
      }
    }

    if (1 != nset) {
      cerr << "Sorry, must specify one, and exactly one of 'posi', 'simple' or 'full' "
              "for substituent atom types\n";
      usage(6);
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  IWString_and_File_Descriptor output(1);

  int rc = 0;

  if (cl.option_present('f')) {
    work_as_filter = 1;
    RingSubstitutionFilter(gen, cl[0], output);
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      if (!RingSubstitution(gen, cl[i], input_type, output)) {
        rc = i + 1;
        break;
      }
    }
  }

  output.flush();

  if (verbose) {
    cerr << "Processed " << molecules_read << " molecules, " << molecules_with_no_rings
         << " had no rings\n";
    if (acc_nset.n() > 1) {
      cerr << "Fingerprints have between " << acc_nset.minval() << " and "
           << acc_nset.maxval() << " bits set, ave "
           << static_cast<float>(acc_nset.average()) << '\n';
    }
  }

  return rc;
}

}  // namespace ring_substitution_main

int
main(int argc, char** argv) {
  ring_substitution_main::prog_name = argv[0];

  int rc = ring_substitution_main::Main(argc, argv);

  return rc;
}
