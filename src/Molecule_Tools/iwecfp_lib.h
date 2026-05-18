#ifndef MOLECULE_TOOLS_IWECFP_LIB_H_
#define MOLECULE_TOOLS_IWECFP_LIB_H_

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <unordered_map>

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/substructure.h"

namespace iwecfp {

using atype_t = unsigned int;

class Iwecfp {
 private:
  int _verbose = 0;

  int _add_tails = 0;
  int _reduce_to_largest_fragment = 0;

  IWString _smiles_tag = "$SMI<";
  IWString _identifier_tag = "PCN<";
  IWString _tag;

  int _min_shell_radius = 0;
  int _max_shell_radius = 0;
  int _additive = 1;
  int _each_shell_gets_own_fingerprint = 0;
  int _all_bonds_same_type = 0;

  std::unordered_map<unsigned int, unsigned int> _bits_to_investigate;
  int _looking_for_bit_meanings = 0;
  int _write_smiles_with_bit_meanings = 0;
  int _bits_found = 0;

  IWString_and_File_Descriptor _stream_for_bit_meanings;
  IWString_and_File_Descriptor _stream_for_all_bits;

  int _equalise_atom_coverage = 0;
  int _label_by_visited = 0;

  uint32_t _fixed_width_fingerprint = 0;

  resizable_array_p<Substructure_Query> _start_atom_query;

  int _central_atom_possible_chiral = 0;
  int _bit_replicates = 0;
  int _bit_replicate_offset = 741;

  IWDigits _iwdigits_center;
  IWDigits _iwdigits;

  Molecule* _current_molecule = nullptr;
  atom_number_t _centre_of_shell = INVALID_ATOM_NUMBER;
  IWString _smarts_for_centre_of_shell;
  isotope_t _centre_atom_isotope = 0;

 public:
  Iwecfp();
  ~Iwecfp();

  int Initialise(Command_Line& cl);

  int Fingerprint(Molecule& m, const atype_t* atom_constant,
                  Sparse_Fingerprint_Creator* sfc);

  int Finalise();
  int Report(std::ostream& output) const;

  const IWString& tag() const { return _tag; }
  int verbose() const { return _verbose; }
  int max_radius() const { return _max_shell_radius;}
  int each_shell_gets_own_fingerprint() const { return _each_shell_gets_own_fingerprint;}

 private:
  int ReadBitsToInvestigate(iwstring_data_source& input);
  int ReadBitsToInvestigate(const char* fname);

  int CheckAgainstList(Molecule& m, const IWString& smarts,
                       atom_number_t centre_of_shell, unsigned int sum_so_far,
                       int radius);
  void WriteLabelledSmiles(const Molecule& m, int centre_of_shell, int radius,
                           IWString_and_File_Descriptor& output);
  void WriteBit(int centre_of_shell,
                const IWString& smarts_for_centre_of_shell,
                int radius, unsigned int bit, Molecule& m,
                IWString_and_File_Descriptor& output);

  int BondConstant(const Bond* bond) const;
  void Increment(unsigned int& sum_so_far, int bc, atype_t atom_constant) const;

  int GenerateShells(int matoms, int radius, int max_radius,
                     const Atom* const* atoms, const atype_t* atom_constant,
                     int* processing_status, unsigned int sum_so_far,
                     Molecule& m, Sparse_Fingerprint_Creator* sfc);

  int AddFingerprint(const Sparse_Fingerprint_Creator& from,
                     Sparse_Fingerprint_Creator& to);
  void FormBit(Molecule& m, const atype_t* atom_constant,
               const Atom* const* atoms, atom_number_t zatom, int max_r,
               int* processing_status, Sparse_Fingerprint_Creator* sfc);

  void IdentifyAtomsWithinRange(Molecule& m, Set_of_Atoms* atoms_within_range);
  int DoEqualiseAtomCoverage(Molecule& m, const Atom* const* atoms,
                             const atype_t* atom_constant,
                             int* processing_status,
                             Sparse_Fingerprint_Creator* sfc);

  int IdentifyStartAtoms(Molecule& m, int* processing_status, int matched_flag);
  int HandleNoStartAtoms(Molecule& m, IWString_and_File_Descriptor& output);
};

void Usage(int rc);
void DisplayDashGOptions(std::ostream& output);
void DisplayDashYOptions(std::ostream& output);

}  // namespace iwecfp

#endif  // MOLECULE_TOOLS_IWECFP_LIB_H_
