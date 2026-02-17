#include <iostream>
#include <filesystem>
#include <memory>
#include <optional>

namespace fs = std::filesystem;

#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "Utilities/GFP_Tools/gfp.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/create_components.h"

#include "retrosynthesis.h"

namespace retrosynthesis {

using std::cerr;

IWString smiles_tag("$SMI<");
IWString identifier_tag("PCN<");

Preprocessor::Preprocessor() {
  _remove_chirality = 0;
}

int
Preprocessor::Initialise(Command_Line& cl) {
  const int verbose = cl.option_present('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (verbose) {
      cerr << "Chirality removed from all molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, 0, 'g')) {
      cerr << "Preprocessor::Initialise:cannot initialise standardisation\n";
      return 0;
    }
  }

  return 1;
}

int
Preprocessor::Process(Molecule& m) {
  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  _chemical_standardisation.process(m);

  return 1;
}

MoleculeAndFragmentation::MoleculeAndFragmentation() {
  _parent = nullptr;
}

std::ostream&
operator<<(std::ostream& output, MoleculeAndFragmentation& mf) {
  output << "MoleculeAndFragmentation for " << mf.smiles() << " found " << mf.number_found();

  return output;
}

int
MoleculeAndFragmentation::DebugPrint(std::ostream& output) const {
  output << "MoleculeAndFragmentation " << name() << " found " << _found.size() << " and " << _notfound.size() << " not found\n";
  return 1;
}

int
MoleculeAndFragmentation::DebugPrint(const IWString& indentation, std::ostream& output) {
  output << indentation << "MoleculeAndFragmentation " << aromatic_smiles() << ' ' << name();

  if (_found.empty() && _notfound.empty()) {
    output << '\n';
    return 1;
  }

  output << " found " << _found.size() << " and " << _notfound.size() << " not found\n";
  IWString next_indentation(indentation);
  next_indentation << "  ";

  if (! _found.empty()) {
    output << indentation << "Found\n";
    for (MoleculeAndFragmentation* f : _found) {
      f->DebugPrint(next_indentation, output);
    }
  }

  if (! _notfound.empty()) {
    output << indentation << "NotFound\n";
    for (MoleculeAndFragmentation* f : _notfound) {
      f->DebugPrint(next_indentation, output);
    }
  }

  return 1;
}


Reaction::Reaction() {
  _reaction = nullptr;
  _number_reactions = 0;

  _isotope = 0;
}

Reaction::~Reaction() {
  if (_reaction != nullptr) {
    delete [] _reaction;
  }
}

// if `fname` is a valid file, return that.
// If we get a valid file by prepending `dirname` return that.
// Otherwise return `fname` and the read will fail.
IWString
MaybeWithDirectoryPrepended(const IWString& dirname, const std::string& fname) {
  cerr << "What about '" << fname << "'\n";
  if (fs::exists(fname)) {
    return IWString(fname);
  }

  // If a full path name, cannot do anything.
  if (fname[0] == '/') {
    return IWString(fname);
  }

  const std::string sdirname(dirname.rawdata(), dirname.length());
  cerr << "Did not exist, try " << sdirname << '\n';
  fs::path fullpath = sdirname;
  fullpath /= fname;

  cerr << "fullpath " <<fullpath << '\n';

  if (fs::exists(fullpath)) {
    return IWString(fullpath);
  }

  return IWString(fname);
}

int
ReadReaction(IWString& fname, Sidechain_Match_Conditions& smc,
             IWReaction& rxn) {

  std::optional<ReactionProto::Reaction> maybe_proto =
      iwmisc::ReadTextProtoCommentsOK<ReactionProto::Reaction>(fname);
  if (!maybe_proto) {
    cerr << "ReadReaction:cannot read reaction proto from '" << fname << "'\n";
    return 0;
  }

  if (!rxn.ConstructFromProto(*maybe_proto, fname, smc)) {
    cerr << "ReadReaction:cannot parse reaction proto\n";
    cerr << maybe_proto->ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

int
Reaction::Build(const ReactionData& proto, Preprocessor& preprocess,
                const IWString& dirname) {
  if (proto.reaction_file_size() == 0) {
    cerr << "Reaction::Build:no reaction file\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  _reaction = new IWReaction[proto.reaction_file_size()];

  Sidechain_Match_Conditions smc;

  _number_reactions = proto.reaction_file_size();
  for (int i = 0; i < _number_reactions; ++i) {
    IWString fname = MaybeWithDirectoryPrepended(dirname, proto.reaction_file(i));
    cerr << "Reaction being read from '" << fname << "'\n";
    if (! ReadReaction(fname, smc, _reaction[i])) {
      cerr << "Retrosynthesis::Build:cannot read reaction '" << fname << "'\n";
      return 0;
    }
  }

  if (proto.fingerprint_file_size() == 0) {
    cerr << "Retrosynthesis::Build:no fingerprints\n";
    return 0;
  }

  for (const std::string& fname : proto.fingerprint_file()) {
    if (! ReadFingerprints(dirname, fname, preprocess)) {
      cerr << "Retrosynthesis::Build:cannot read fingerprints '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

Retrosynthesis::Retrosynthesis() {
  _reaction = nullptr;
  _number_reactions = 0;

  _verbose = 0;

  set_copy_user_specified_atom_void_ptrs_during_create_subset(1);
}

Retrosynthesis::~Retrosynthesis() {
  if (_reaction != nullptr) {
    delete [] _reaction;
  }
}

int
Retrosynthesis::Initialise(Command_Line& cl, Preprocessor& preprocess) {
  _verbose = cl.option_present('v');

  if (! cl.option_present('C')) {
    cerr << "Retrosynthesis::Initialise:must specify data proto via the -C option\n";
    return 0;
  }

  IWString fname = cl.string_value('C');
  std::optional<RetrosynthesisData> maybe_proto = iwmisc::ReadTextProto<RetrosynthesisData>(fname);
  if (! maybe_proto) {
    cerr << "RetrosynthesisData::Initialise:invalid proto '" << fname << "'\n";
    return 0;
  }

  IWString dirname;
  cerr << "Config file is '" << fname << "'\n";
  fs::path pname(fname.null_terminated_chars());
  dirname = pname.parent_path();
  cerr << "From '" << fname << " get dirname " << dirname << '\n';

  if (! Build(*maybe_proto, preprocess, dirname)) {
    cerr << "RetrosynthesisData::Initialise:invalid proto\n";
    cerr << maybe_proto->ShortDebugString() << '\n';
    return 0;
  }

  return 1;
}

int
Retrosynthesis::Build(const RetrosynthesisData& proto, Preprocessor& preprocess,
                const IWString& dirname) {
  _number_reactions = proto.reaction_size();
  if (_number_reactions == 0) {
    cerr << "Retrosynthesis::Build:no building blocks in proto\n";
    cerr << proto.ShortDebugString() << '\n';
    return 0;
  }

  _reaction = new Reaction[_number_reactions];

  for (int i = 0; i < _number_reactions; ++i) {
    if (! _reaction[i].Build(proto.reaction(i), preprocess, dirname)) {
      cerr << "Retrosynthesis::Build:cannot build reaction " << i << '\n';
      cerr << proto.reaction(i).ShortDebugString() << '\n';
      return 0;
    }
  }

  return _number_reactions;
}

int
Reaction::ReadFingerprints(const IWString& dirname, const std::string& fname, Preprocessor& preprocess) {
  IWString fullpath = MaybeWithDirectoryPrepended(dirname, fname);
  iwstring_data_source input;
  if (! input.open(fullpath)) {
    cerr << "Retrosynthesis::ReadFingerprints:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadFingerprints(input, preprocess);
}

int
Reaction::ReadFingerprints(iwstring_data_source& input,
        Preprocessor& preprocess) {
  std::unique_ptr<SetOfFingerprints> fp = std::make_unique<SetOfFingerprints>();

  if(! fp->Build(input, preprocess)) {
    cerr << "Retrosynthesis::ReadFingerprints:error reading fingerprints\n";
    return 0;
  }

  const isotope_t iso = fp->isotope();
  if (iso == 0) {
    cerr << "return::ReadFingerprints:no isotope or multiple isotopes\n";
    return 0;
  }

  auto iter = _iso_to_fp.find(iso);
  if (iter != _iso_to_fp.end()) {
    cerr << "Retrosynthesis::ReadFingerprints:duplicate isotope " << iso << '\n';
    return 0;
  }

  _iso_to_fp[iso] = fp.release();

  return 1;
}

// Hopefully we are dealing with fragments that have a single isotopic atom.
// If we find `m` has a single isotopic atom, return that atom number
atom_number_t
SingleIsotope(const Molecule& m) {
  isotope_t single_isotope = 0;
  atom_number_t rc = kInvalidAtomNumber;

  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    const isotope_t iso = m.isotope(i);

    if (iso == 0) [[likely]] {
      continue;
    }

    // We must have just one instance of that isotope.
    if (single_isotope > 0) [[unlikely]] {
      return kInvalidAtomNumber;
    }

    // There can be just one atom with this isotope.
    if (rc != kInvalidAtomNumber) [[unlikely]] {
      return kInvalidAtomNumber;
    }

    rc = i;
    single_isotope = iso;
  }

  return rc;
}

SetOfFingerprints::SetOfFingerprints() {
  _fp = nullptr;
  _number_fingerprints = 0;
  _isotope = 0;
}

SetOfFingerprints::~SetOfFingerprints() {
  if (_fp != nullptr) {
    delete [] _fp;
  }
}

int
SetOfFingerprints::Build(iwstring_data_source& input,
        Preprocessor& preprocess) {
  _number_fingerprints = input.count_records_starting_with(identifier_tag);
  if (_number_fingerprints == 0) {
    cerr << "Retrosynthesis::ReadFingerprints:no fingerprints in file\n";
    return 0;
  }

  _fp = new IW_General_Fingerprint[_number_fingerprints];

  IW_TDT tdt;
  for (int ndx = 0; tdt.next(input); ++ndx) {
    int fatal = 0;
    if (! _fp[ndx].construct_from_tdt(tdt, fatal)) {
      cerr << "SetOfFingerprints::Build:bad data\n";
      return 0;
    }

    const IWString& id = _fp[ndx].id();
    if (id.empty()) {
      cerr << "Retrosynthesis::ReadFingerprints:empty id\n";
      return 0;
    }

    IWString smiles;
    if (! tdt.dataitem_value(smiles_tag, smiles)) {
      cerr << "Retrosynthesis::ReadFingerprints:no smiles in tdt\n";
      return 0;
    }

    Molecule m;
    if (! m.build_from_smiles(smiles)) {
      cerr << "Retrosynthesis::ReadFingerprints:invalid smiles\n";
      cerr << smiles << '\n';
      return 0;
    }

    preprocess.Process(m);

    if (atom_number_t a = SingleIsotope(m); a == kInvalidAtomNumber) {
      cerr << smiles << " no isotope\n";
      return 0;
    } else if (_isotope == 0) {
      _isotope = m.isotope(a);
    } else if (m.isotope(a) == _isotope) {
      // Same as what was seen before, good
    } else {
      cerr << "SetOfFingerprints::ReadFingerprints:more than one isotopic value\n";
      cerr << smiles << '\n';
      cerr << "Previously " << _isotope << " now read " << m.isotope(a) << '\n';
      return 0;
    }

    // Should check for duplicates.
    _usmi_to_id[m.unique_smiles()] = id;
  }

  return _number_fingerprints;
}

int
SetOfFingerprints::InDatabase(Molecule& m) {
  const IWString& usmi = m.unique_smiles();
  cerr << _isotope << " looking for " << m.smiles() << '\n';

  if (auto iter = _usmi_to_id.find(usmi); iter == _usmi_to_id.end()) {
    return 0;
  } else {
    m.set_name(iter->second);
    return 1;
  }
}

int
Retrosynthesis::Process(Molecule& m, MoleculeAndFragmentation& results) {
  const int matoms = m.natoms();
  std::unique_ptr<int[]> atom_number = std::make_unique<int[]>(matoms);
  for (int i = 0; i < matoms; ++i) {
    atom_number[i] = i;
    m.set_user_specified_atom_void_ptr(i, atom_number.get() + i);
  }

  cerr << "Retrosynthesis has " << _number_reactions << " reactions\n";

  return Process2(m, results);
}

int
Retrosynthesis::Process2(Molecule& m, MoleculeAndFragmentation& results) {
  int rc = 0;

  for (int i = 0; i < _number_reactions; ++i) {
    rc += _reaction[i].Process(m, results);
  }

  if (rc == 0) {
    return 0;
  }

  // We found matches, and the number not found has not increased, we are done.
  if (results.number_not_found() == 0) {
    return rc;
  }

  for (MoleculeAndFragmentation* subset : results.notfound()) {
    rc += Process2(*subset, *subset);
  }

  return rc;
}

int
Reaction::Process(Molecule& m, MoleculeAndFragmentation& results) {
  int rc = 0;

  Substructure_Results sresults;
  for (int i = 0; i < _number_reactions; ++i) {
    if (_reaction[i].substructure_search(m, sresults) == 0) {
      continue;
    }
    cerr << "match to qeruy\n";

    for (const Set_of_Atoms* s : sresults.embeddings()) {
      Molecule tmp(m);
      if (_reaction[i].in_place_transformation(tmp, s)) [[unlikely]] {
        cerr << " After reaction " << tmp.smiles() << '\n';
        rc += ProcessMatch(tmp, results);
      }
    }
  }

  return rc;
}

// `m` has been processed by one of our reactions and should be multi-fragment
int
Reaction::ProcessMatch(Molecule& m, MoleculeAndFragmentation& results) {
  if (m.number_fragments() < 2) {
    return 0;
  }

  resizable_array_p<MoleculeAndFragmentation> frags;
  m.create_components(frags);

  cerr << "Created " << frags.size() << " components\n";
  for (auto* q : frags) {
    cerr << q->smiles() << '\n';
  }

  for (MoleculeAndFragmentation* f : frags) {
    atom_number_t single_isotope = SingleIsotope(*f);
    if (single_isotope == kInvalidAtomNumber) [[unlikely]] {
      continue;
    }

    // f->set_isotope(single_isotope, 1);
    cerr << "Looking for " << f->smiles() << '\n';
    if (InDatabase(*f)) {
      cerr << "Found\n";
      results.BeenFound(f);
    } else {
      cerr << " Not found\n";
      results.NotFound(f);
    }
  }

  frags.resize_no_delete (0);

  return 1;
}

// `m` is a newly created fragment with isotope 1 at the attachment point.
// If this molecule is in the smiles hash of any reagent, set the name and return 1.
int
Reaction::InDatabase(Molecule& m) {
  for (const auto[k, v] : _iso_to_fp) {
    if (v->InDatabase(m)) {
      return 1;
    }
  }

  return 0;
}


}  // namespace retrosynthesis
