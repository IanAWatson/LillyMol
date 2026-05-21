// Tester for writing zstd compressed files.

#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iwstring_and_file_descriptor.h"

namespace test_zstd {

using std::cerr;

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "v");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    return 1;
  }

  if (cl.empty()) {
    cerr << "Must specify input file as command line argument\n";
    return 1;
  }

  IWString ifname(cl[0]);
  IWString ofname(cl[0]);
  ofname << ".zst";

  iwstring_data_source input(ifname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Cannot open '" << ifname << "'\n";
    return 1;
  }

  iwstring::IWString_and_File_Descriptor output;
  if (! output.open(ofname)) {
    cerr << "Cannot open compressed output file '" << ofname << "'\n";
    return 1;
  }

  IWString file_contents;

  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(1024);
    file_contents << buffer << '\n';
  }

  output << file_contents;

//cerr << "Calling close\n";
//output.close();

  return 0;
}

}  // namespace test_zstd
int
main(int argc, char** argv) {
  return test_zstd::Main(argc, argv);
}
