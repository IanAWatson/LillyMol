// Convert molecules to SAFE form.
// Gotta be SAFE: A New Framework for Molecular Design
// Emanuel Noutahi, Christian Gabellini, Michael Craig
// Jonathan Lim, Prudencio Tosou. Valence Labs.
// https://arxiv.org/pdf/2310.10773.pdf

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"

#include "Molecule_Lib/molecule.h"

#include "highest_ring_number.h"

namespace mol2safe {

using std::cerr;

void
Usage(int rc) {
  ::exit(rc);
}

class Options {
  private:
    int _verbose = 0;
    int _molecules_read = 0;
    int _bad_smiles = 0;
    int _too_many_ring_numbers = 0;

  public:
};

int
Mol2SAFE(Molecule& m,
         const int hring,
         IWString_and_File_Descriptor& output) {
  return 1;
}

int
Mol2SAFEInner(const IWString& buffer,
              IWString_and_File_Descriptor& output) {
  Molecule m;
  if (! m.build_from_smiles(buffer)) {
    cerr << "Mol2SAFEInner:invalid smiles\n";
    return 0;
  }

  std::optional<int> hring = lillymol::HighestRingNumber(buffer);
  if (! hring) {
    return 1;
  }

  return Mol2SAFE(m, *hring, output);
}

int
Mol2SAFE(iwstring_data_source& input,
         IWString_and_File_Descriptor& output) {
  IWString buffer;
  while (input.next_record(buffer)) {
    if (! Mol2SAFEInner(buffer, output)) {
      cerr << "Mol2SAFEInner:error processing '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Mol2SAFE(const char* fname,
         IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Mol2SAFE:cannot open '" << fname << "'\n";
    return 0;
  }

  return Mol2SAFE(input, output);
}

int
Mol2SAFE(int argc, char** argv) {
  Command_Line cl(argc, argv, "v");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }
  
  IWString_and_File_Descriptor output(1);
  for (const char* fname : cl) {
    if (! Mol2SAFE(fname, output)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  return 0;
}


}  // namespace mol2safe


int
main(int argc, char ** argv) {

  int rc = mol2safe::Mol2SAFE(argc, argv);

  return rc;
}
