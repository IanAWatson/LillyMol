#include "absl/container/flat_hash_map.h"

#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwaray/iwaray.h"
#include "Foundational/iwstring/absl_hash.h"

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/gfp_standard.h"

#include "Molecule_Lib/iwreaction.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#ifdef BUILD_BAZEL
#include "Retrosynthesis/retrosynthesis.pb.h"
#else
#include "retrosynthesis.pb.h"
#endif

namespace retrosynthesis {

// The overall class which parses command lines may have some structure
// altering directives. These need to be passed down to all the objects
// that get created by it.
class Preprocessor {
  private:
    int _remove_chirality;

    Chemical_Standardisation _chemical_standardisation;
  public:
    Preprocessor();

    // Only the -c and -g options are recognised.
    int Initialise(Command_Line& cl);

    int Process(Molecule& m);
};

// Molecules need to be decomposed.
class MoleculeAndFragmentation : public Molecule {
  private:
    MoleculeAndFragmentation* _parent;

    // THe _found molecules will have the ID of the reagent.
    resizable_array_p<MoleculeAndFragmentation> _found;

    resizable_array_p<MoleculeAndFragmentation> _notfound;

  public:
    MoleculeAndFragmentation();

    void set_parent_molecule(const Molecule& m) {
      Molecule::operator=(m);
    }

    int DebugPrint(std::ostream& output) const;
    int DebugPrint(const IWString& indentation, std::ostream& output);

    resizable_array_p<MoleculeAndFragmentation>& found() {
      return _found;
    }
    resizable_array_p<MoleculeAndFragmentation>& notfound() {
      return _notfound;
    }

    uint32_t number_found() const {
      return _found.size();
    }
    uint32_t number_not_found() const {
      return _notfound.size();
    }

    // Assume ownership of `m` in the _found array.
    void BeenFound(MoleculeAndFragmentation* m) {
      _found << m;
    }
    void NotFound(MoleculeAndFragmentation* m) {
      _notfound << m;
    }
};

std::ostream&
operator<<(std::ostream& output, MoleculeAndFragmentation& mf);

class SetOfFingerprints {
  private:
    IW_General_Fingerprint* _fp;
    int _number_fingerprints;
    isotope_t _isotope;

    // a mapping from the unique smiles of a building block to its name.
    absl::flat_hash_map<IWString, IWString> _usmi_to_id;

  public:
    SetOfFingerprints();
    ~SetOfFingerprints();

    int Build(IWString& fname, Preprocessor& preprocess);
    int Build(iwstring_data_source& input, Preprocessor& preprocess);

    isotope_t isotope() const {
      return _isotope;
    }

    // If the m.unique_smiles() is in _usmi_to_id, set m.name();
    int InDatabase(Molecule& m);
};

// A reversed reaction generates 
class Reaction {
  private:
    // Perhaps it will be necessary to have multiple reactions to
    // handle a specific set of building blocks - boronic...
    IWReaction* _reaction;
    int _number_reactions;

    // The name of the first _reaction. It is an error to NOT have a name.
    IWString _name;

    // All molecules in the set will have an atom with this isotope.
    isotope_t _isotope;

    // A mapping from isotope number to a set of fingerprints.
    absl::flat_hash_map<int, SetOfFingerprints*> _iso_to_fp;

  // private functions
    int ReadFingerprints(const IWString& dirname, const std::string& fname, Preprocessor& preprocess);
    int ReadFingerprints(iwstring_data_source & input, Preprocessor& preprocess);

    int ProcessMatch(Molecule& m, MoleculeAndFragmentation& results);
    int InDatabase(Molecule& m);

  public:
    Reaction();
    ~Reaction();

    int Build(const ReactionData& proto, Preprocessor& preprocess, const IWString& dirname);

    const IWString& name() const {
      return _name;
    }

    int Process(Molecule& m, MoleculeAndFragmentation& result);
};

class Retrosynthesis {
  private:
    Reaction* _reaction;
    int _number_reactions;

    int _verbose;

  // Private functions
    int Build(const RetrosynthesisData& proto, Preprocessor& preprocess, const IWString& dirname);

    int ReadFingerprints(IWString& fname);
    int ReadFingerprints(iwstring_data_source& input);

    int Process2(Molecule& m, MoleculeAndFragmentation& results);

  public:
    Retrosynthesis();
    ~Retrosynthesis();

    int Initialise(Command_Line& cl, Preprocessor& preprocess);

    int Process(Molecule& m, MoleculeAndFragmentation& results);
};  // class Retrosynthesis

}  // namespace retrosynthesis
