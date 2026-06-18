
#include <iostream>

#include "Foundational/accumulator/accumulator.h"

#include "Molecule_Tools/iwdescr_internal.h"

using std::cerr;
void
DescriptorsToCompute::SetAll(int s) {
  adjacent_ring_fusion_descriptors = s;
  bonds_between_rings = s;
  charge_descriptors = s;
  complexity_descriptors = s;
  crowding_descriptors = s;
  distance_matrix_descriptors = s;
  donor_acceptor = s;
  hbond_descriptors = s;
  ncon_descriptors = s;
  perform_expensive_chirality_perception = s;
  partial_symmetry_descriptors = s;
  polar_bond_descriptors = s;
  ring_chain_descriptors = s;
  ring_substitution_descriptors = s;
  compute_alogp = s;
  compute_xlogp = s;
  ring_substitution_ratio_descriptors = s;
  simple_hbond_descriptors = s;
  specific_groups = s;
  spinach_descriptors = s;
  symmetry_descriptors = s;
  psa = s;
  long_carbon_chains = s;
  saturated_chains = s;
  planarity = s;
  ring_fusion_descriptors = s;
  ramey_descriptors = s;
}

int
DescriptorsToCompute::DisplayUsage(std::ostream& output) const {
  // clang-format off
  output << R"(Specify descriptor(s) to be disabled or enabled.
If the feature name starts with '-', the feature is disabled.
 -O none          turn off all optional descriptors
 -O all           turn on  all optional descriptors
 -O adjring       enable adjacent to ring fusion descriptors
 -O bbr           enable bonds between rings descriptors
 -O charge        enable all formal charge descriptors
 -O chiral        enable all expensive chirality perception descriptors
 -O complex       enable planar fused ring descriptors
 -O crowd         enable atomic crowding descriptors
 -O lcc           enable atomic crowding descriptors
 -O dm            enable descriptors based on the distance matrix
 -O donacc        enable donor acceptor derived descriptors
 -O hbond         enable donor/acceptor derived descriptors
 -O shbond        enable simplistic donor/acceptor derived descriptors
 -O ncon          enable ncon and fncon descriptors
 -O pbond         enable polar bond derived descriptors
 -O psa           enable Novartis polar surface area descriptors
 -O planarity     enable convex planarity descriptor
 -O psymm         enable partial symmetry derived descriptors
 -O ramey         enable Ramey (element count) descriptors
 -O rcj           enable ring chain join descriptors
 -O rfuse         enable ring fusion descriptors
 -O rss           enable ring substitution descriptors
 -O rssr          enable ring substitution ratio descriptors
 -O spch          enable spinach related descriptors
 -O spcgrp        enable specific group descriptors
 -O satchain      enable computation of saturated chain descriptors
 -O symm          enable symmetry related descriptors
 -O xlogp         enable xlogp
)";
  // clang-format on

  return 0;
}

int
DescriptorsToCompute::Initialise(Command_Line& cl) {
  static constexpr char kFlag = 'O';
  const int verbose = cl.option_present('v');
  const_IWSubstring o;
  for (int i = 0; cl.value(kFlag, o, i); ++i) {
    int s = 1;
    if (o[0] == '-') {
      s = 0;
      o++;
    }
    if (o == "none") {
      SetAll(0);
      if (verbose) {
        cerr << "All optional descriptors turned off\n";
      }
    } else if (o == "all") {
      SetAll(1);
      if (verbose) {
        cerr << "All optional descriptors enabled\n";
      }
    } else if (o == "adjring") {
      adjacent_ring_fusion_descriptors = s;
      if (verbose) {
        cerr << "bonding adjacent to ring fusion descriptors enabled\n";
      }
    } else if (o == "bbr") {
      bonds_between_rings = s;
      if (verbose) {
        cerr << "bonds between rings descriptors enabled\n";
      }
    } else if (o == "charge") {
      charge_descriptors = s;
      if (verbose) {
        cerr << "charge derived descriptors enabled\n";
      }
    } else if (o == "complex") {
      complexity_descriptors = s;
      if (verbose) {
        cerr << "complexity derived descriptors enabled\n";
      }
    } else if (o == "crowd") {
      crowding_descriptors = s;
      if (verbose) {
        cerr << "atomic crowding descriptors enabled\n";
      }
    } else if (o == "lcc") {
      long_carbon_chains = s;
      if (verbose) {
        cerr << "long carbon chain enabled\n";
      }
    } else if (o == "dm") {
      distance_matrix_descriptors = s;
      if (verbose) {
        cerr << "distance matrix derived descriptors enabled\n";
      }
    } else if (o == "donacc") {
      donor_acceptor = s;
      if (verbose) {
        cerr << "donor acceptor derived descriptors enabled\n";
      }
    } else if (o == "hbond") {
      hbond_descriptors = s;
      if (verbose) {
        cerr << "hydrogen bonding derived descriptors enabled\n";
      }
    } else if (o == "shbond") {
      simple_hbond_descriptors = s;
      if (verbose) {
        cerr << "simplistic hydrogen bonding derived descriptors enabled\n";
      }
    } else if (o == "chiral" || o == "expchiral") {
      perform_expensive_chirality_perception = s;
      if (verbose) {
        cerr << "Will check all unlabelled atoms for possible chirality\n";
      }
    } else if (o == "ncon") {
      ncon_descriptors = s;
      if (verbose) {
        cerr << "ncon and fncon descriptors enabled\n";
      }
    } else if (o == "pbond") {
      polar_bond_descriptors = s;
      if (verbose) {
        cerr << "polar bond derived descriptors enabled\n";
      }
    } else if (o == "psa") {
      psa = s;
      if (verbose) {
        cerr << "Novartis polar surface area descriptor enabled\n";
      }
    } else if (o == "psymm") {
      partial_symmetry_descriptors = s;
      if (verbose) {
        cerr << "partial symmetry derived descriptors enabled\n";
      }
    } else if (o == "symm") {
      symmetry_descriptors = s;
      if (verbose) {
        cerr << "symmetry derived descriptors enabled\n";
      }
    } else if (o == "ramey") {
      ramey_descriptors = s;
      if (verbose) {
        cerr << "element count descriptors enabled\n";
      }
    } else if (o == "rcj") {
      ring_chain_descriptors = s;
      if (verbose) {
        cerr << "ring chain join descriptors enabled\n";
      }
    } else if (o == "rfuse") {
      ring_fusion_descriptors = s;
      if (verbose) {
        cerr << "ring fusion descriptors enabled\n";
      }
    } else if (o == "rss") {
      ring_substitution_descriptors = s;
      if (verbose) {
        cerr << "ring substitution descriptors enabled\n";
      }
    } else if (o == "rssr") {
      ring_substitution_ratio_descriptors = s;
      if (verbose) {
        cerr << "ring substitution ratio descriptors enabled\n";
      }
    } else if (o == "spch") {
      spinach_descriptors = s;
      if (verbose) {
        cerr << "spinach derived descriptors enabled\n";
      }
    } else if (o == "spcgrp") {
      specific_groups = s;
      if (verbose) {
        cerr << "specific substructure descriptors enabled\n";
      }
    } else if (o == "alogp") {
      compute_alogp = s;
      if (verbose) {
        cerr << "alogp will be computed\n";
      }
    } else if (o == "xlogp") {
      compute_xlogp = s;
      if (verbose) {
        cerr << "xlogp will be computed\n";
      }
    } else if (o == "satchain") {
      saturated_chains = s;
      if (verbose) {
        cerr << "Saturated chain descriptors computed\n";
      }
    } else if (o == "planarity") {
      planarity = s;
      if (verbose) {
        cerr << "Will compute planarity descriptors\n";
      }
    } else if (o.starts_with("F:")) {
      o.remove_leading_chars(2);
      if (!ReadDescriptorsToCompute(o)) {
        cerr << "DescriptorsToCompute::Initialise:cannot process '" << o << "'\n";
        return 0;
      }
    } else if (o == "help") {
      return DisplayUsage(cerr);
    } else {
      cerr << "DescriptorsToCompute::Initialise:unrecognised '" << o << "'\n";
      return DisplayUsage(cerr);
    }
  }

  return 1;
}

int
DescriptorsToCompute::ReadDescriptorsToCompute(const const_IWSubstring& fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "DescriptorsToCompute::ReadDescriptorsToCompute:cannot open '" << fname
         << "'\n";
    return 0;
  }

  return DescriptorsToCompute::ReadDescriptorsToCompute(input);
}

int
DescriptorsToCompute::ReadDescriptorsToCompute(iwstring_data_source& input) {
  input.set_strip_trailing_blanks(1);

  // Turn off everything to start
  SetAll(0);

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (buffer.starts_with("#")) {
      continue;
    }

    // This code is copied from above, so keep the same variable name.
    // Should set up an Activate() function...
    const const_IWSubstring o(buffer);

    if (o == "none") {
      SetAll(0);
    } else if (o == "all") {
      SetAll(1);
    } else if (o == "adjring") {
      adjacent_ring_fusion_descriptors = 1;
    } else if (o == "bbr") {
      bonds_between_rings = 1;
    } else if (o == "charge") {
      charge_descriptors = 1;
    } else if (o == "complex") {
      complexity_descriptors = 1;
    } else if (o == "crowd") {
      crowding_descriptors = 1;
    } else if (o == "dm") {
      distance_matrix_descriptors = 1;
    } else if (o == "donacc") {
      donor_acceptor = 1;
    } else if (o == "hbond") {
      hbond_descriptors = 1;
    } else if (o == "shbond") {
      simple_hbond_descriptors = 1;
    } else if (o == "chiral" || o == "expchiral") {
      perform_expensive_chirality_perception = 1;
    } else if (o == "ncon") {
      ncon_descriptors = 1;
    } else if (o == "pbond") {
      polar_bond_descriptors = 1;
    } else if (o == "psa") {
      psa = 1;
    } else if (o == "psymm") {
      partial_symmetry_descriptors = 1;
    } else if (o == "symm") {
      symmetry_descriptors = 1;
    } else if (o == "ramey") {
      ramey_descriptors = 1;
    } else if (o == "rcj") {
      ring_chain_descriptors = 1;
    } else if (o == "rfuse") {
      ring_fusion_descriptors = 1;
    } else if (o == "rss") {
      ring_substitution_descriptors = 1;
    } else if (o == "rssr") {
      ring_substitution_ratio_descriptors = 1;
    } else if (o == "spch") {
      spinach_descriptors = 1;
    } else if (o == "spcgrp") {
      specific_groups = 1;
    } else if (o == "alogp") {
      compute_alogp = 1;
    } else if (o == "xlogp") {
      compute_xlogp = 1;
    } else {
      cerr << "DescriptorsToCompute::ReadDescriptorsToCompute:unrecognised '" << o
           << "'\n";
      return 0;
    }
  }

  return 1;
}
