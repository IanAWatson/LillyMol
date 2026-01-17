// Spread variant where molecules must be selected in groups

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"


#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#define RESIZABLE_ARRAY_IMPLEMENTATION 1
#include "Foundational/iwaray/iwaray.h"
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION 1
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iwstring.h"
#include "Foundational/iw_tdt/iw_tdt.h"

#include "gfp_standard.h"

namespace gfp_group_spread {

using std::cerr;

// clang-format off
static void
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(Spread variant that selects molecules in groups.
)";

  ::exit(rc);
}

class Fingerprint : public GFP_Standard {
  private:
    IWString _name;

    IWString _group;
    float _distance;
    float _desirability;

  public:
    Fingerprint();

    void set_group(const IWString& s) {
      _group = s;
    }
    const IWString& group() const {
      return _group;
    }

    void set_name(const IWString& s) {
      _name = s;
    }
    const IWString& name() const {
      return _name;
    }

    void set_distance(float s) {
      _distance = s;
    }
    float distance() const {
      return _distance;
    }

    void set_desirability(float s) {
      _desirability = s;
    }
    float desirability() const {
      return _desirability;
    }
};

Fingerprint::Fingerprint() {
  _distance = 0.0f;
}

// A group is desirable because it contains molecules that are
// distinct from what has already been selected.
// A group is desirable because it contains desirable molecules.
// This defines how those completely separate objectives are combined.
struct Objective {
  public:
    double _distance_weight;
    double _desirability_weight;
  public:
    Objective();
};

Objective::Objective() {
  _distance_weight = 0.5;
  _desirability_weight = 0.5;
}

class Group {
  private:
    // Will be set once this group has been selected
    // It will be the order in which our group was found.
    uint32_t _selected;

    IWString _group_name;

    resizable_array<Fingerprint*> _fp;

    // Each fingerprint can have an associated desirability based on the
    // desirability of each member of the group.
    Accumulator<double> _desirability;

    // The distances to previously selected for all our fingerprints.
    Accumulator<double> _distance;

    // Maybe it would be interesting to keep track of all distances.
    // Convert to a number in the range [0,100)
    // We allocate 101 values just in case a distance of 1.0 ever happens
    // which is hard to imagine.
    uint32_t* _distance_histogram;

    // We don't want to give too much weight to long distances
    // that are likely to be non meaningful.
    static inline float _truncate_long_distances = std::numeric_limits<float>::max();
    // If _truncate_long_distances is set, we can do a one-time computation of the
    // index into _distance_histogram.
    uint32_t _truncate_index;

    static inline char _output_separator = ' ';

    // We may have a penalty on exact matches.
    uint32_t _zero_distances;

    // private functions
    void AnotherDistance(float d);

  public:
    Group(const IWString& id);
    ~Group();

    const IWString& name() const {
      return _group_name;
    }

    void set_selected(uint32_t s) {
      _selected = s;
    }
    uint32_t selected() const {
      return _selected;
    }

    void SetTruncateLongDistances(float s);

    int EstablishInitialDistances(Group& rhs);

    int AnotherMember(Fingerprint* fp);

    double Score(const Objective& objective) const;

    void GroupHasBeenSelected(const Group* sel);

    double AverageDistance() const {
      return _distance.average();
    }
    double AverageDesirability() const {
      return _desirability.average();
    }

    int Write(IWString_and_File_Descriptor& output) const;
};

Group::Group(const IWString& id) : _group_name(id) {
  _selected = 0;
  _zero_distances = 0;
  _distance_histogram = new uint32_t[101];
  std::fill_n(_distance_histogram, 101, 0);
  _truncate_index = 0;
}

Group::~Group() {
  if (_distance_histogram != nullptr) {
    delete [] _distance_histogram;
  }
}

void
Group::SetTruncateLongDistances(float s) {
  _truncate_long_distances = s;
  _truncate_index = static_cast<uint32_t>(s * 100.0f);
}

int
Group::AnotherMember(Fingerprint* fp) {
  _fp << fp;
  _desirability.extra(fp->desirability());
  AnotherDistance(fp->distance());


  return 1;
}

int
Group::EstablishInitialDistances(Group& rhs) {
  const uint32_t n1 = _fp.size();
  const uint32_t n2 = rhs._fp.size();
  for (uint32_t i = 0; i < n1; ++i) {
    const Fingerprint* fpi = _fp[i];
    for (uint32_t j = 0; j < n2; ++j) {
      float d = fpi->tanimoto_distance(*rhs._fp[j]);
      AnotherDistance(d);
      rhs.AnotherDistance(d);
    }
  }

  return 1;
}


void
Group::AnotherDistance(float d) {
  if (d > _truncate_long_distances) {
    _distance.extra(_truncate_long_distances);
    if (_distance_histogram != nullptr) {
      ++_distance_histogram[_truncate_index];
    }
  } else {
    _distance.extra(d);
    if (d == 0.0f) {
      ++_zero_distances;
    }

    if (_distance_histogram != nullptr) {
      uint32_t ndx = static_cast<uint32_t>(d * 100.0f);
      cerr << "Updating histogram, ndx " << ndx << "\n'";
      _distance_histogram[ndx] += 1;
    }
  }
}

void
Group::GroupHasBeenSelected(const Group* sel) {
  for (const Fingerprint* rhs : sel->_fp) {
    for (Fingerprint* lhs : _fp) {
      float d = rhs->tanimoto_distance(*lhs);
      AnotherDistance(d);
    }
  }
}

int
Group::Write(IWString_and_File_Descriptor& output) const {
  output << _group_name;
  for (const Fingerprint* f : _fp) {
    output << _output_separator << f->name();
  }
  output << '\n';

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

class Options {
  private:
    // The number of groups to select - the -n option.
    uint32_t _nsel;

    uint32_t _number_fingerprints;
    Fingerprint* _fp;

    // A mapping from names to index in the _fp array.
    absl::flat_hash_map<IWString, uint32_t> _id_to_ndx;

    IWString _identifier_tag;

    // Each unique group name will populate one of these.
    resizable_array_p<Group> _group;

    // The overall desirability of a group of molecules is a composite
    // of how desirable the molecules are and the distance to previously
    // selected items.
    // Note that the desirability should be in the range [0,1] in order
    // to ensure comparability.
    double _distance_weight;
    double _desirability_weight;

    int _verbose;

    // private functions
    int ReadDesirabilty(IWString& fname);
    int ReadDesirabilty(iwstring_data_source& input);
    int ReadDesirabiltyRecord(const const_IWSubstring& buffer);
    int ReadFingerprints(iwstring_data_source& input);
    int ReadGroupMembership(iwstring_data_source& input);
    int ReadPreviousDistance(iwstring_data_source& input);
    int ReadGroupMembershipRecord(const const_IWSubstring& buffer);
    int ReadPreviousDistanceRecord(const const_IWSubstring& buffer);

    int DetermineInitialDistances();

    Group* GetNextGroup() const;
    double Score(double ave_dist, double ave_desirability) const;

  public:
    Options();
    ~Options();

    int Initialise(const Command_Line& cl);

    int ReadFingerprints(const char* fname);
    int ReadGroupMembership(IWString& fname);
    int ReadPreviousDistance(IWString& fname);

    int Spread(IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _number_fingerprints = 0;
  _fp = nullptr;
  _verbose = 0;
  _nsel = 0;

  _distance_weight = 0.5;
  _desirability_weight = 0.5;

  _identifier_tag = "PCN<";
}

Options::~Options() {
  if (_fp != nullptr) {
    delete [] _fp;
  }
}

void
DisplayDashwOptions(int rc) {
  cerr << R"(The -w option controls the relative weight given to average desirability and distance to previously selected.
These two properties need to be combined in order to come up with a unified score for each group of molecules.
The -w option specifies the relative weights assigned to each measure of fitness.
All weights are relative weights and should sum to 1.0, although this is not enforced.

 -w dist=<w>            the relative weight assigned to distance from previously selected.
 -w desr=<w>            the relative weight assigned to the score of the molecules in the group.
)";

  ::exit(rc);
}

int
Options::Initialise(const Command_Line& cl) {
  _verbose = cl.option_present('v');

  if (! ReadFingerprints(cl[0])) {
    cerr << "Options::Initialise:cannot read fingerprints '" << cl[0] << "'\n";
    return 0;
  }

  if (_verbose) {
    cerr << "Read " << _number_fingerprints << " fingerprints\n";
  }

  if (cl.option_present('T')) {
    float tld;
    if (! cl.value('T', tld) || tld < 0.0f || tld > 1.0f) {
      cerr << "Options::Initialise:invalid distance cutoff (-T)\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Will truncate distances at " << tld << '\n';
    }

    Group notused("notused");
    notused.SetTruncateLongDistances(tld);
  }

  DetermineInitialDistances();

  if (cl.option_present('A')) {
    IWString fname = cl.string_value('A');
    if (! ReadPreviousDistance(fname)) {
      cerr << "Options::Initialise:cannot read previously computed distances file '" << fname << "'\n";
      return 0;
    }
  }

  if (! cl.option_present('G')) {
    cerr << "Options::Initialise:Must specify Group memberships with the -G option\n";
    Usage(1);
  }

  if (cl.option_present('G')) {
    IWString fname = cl.string_value('G');
    if (! ReadGroupMembership(fname)) {
      cerr << "Options::Initialise:cannot read previously computed distances '" << fname << "'\n";
      return 1;
    }
  }

  // THis will be an index into the _group array.
  absl::flat_hash_map<IWString, uint32_t> group_to_ndx;
  for (uint32_t i = 0; i < _number_fingerprints; ++i) {
    const IWString& g = _fp[i].group();
    if (g.empty()) {
      cerr << "Options::Initialise:no group " << _fp[i].name() << '\n';
      return 0;
    }

    auto iter = group_to_ndx.find(g);
    if (iter == group_to_ndx.end()) {
      uint32_t n = group_to_ndx.size();
      group_to_ndx[g] = n;
    }
    _group << new Group(g);
  }

  const uint32_t ngroups = group_to_ndx.size();

  if (ngroups < 2) {
    cerr << "Options::Initialise:only " << group_to_ndx.size() << " groups\n";
    return 0;
  }

  _group.resize(ngroups);

  for (uint32_t i = 0; i < _number_fingerprints; ++i) {
    const IWString& g = _fp[i].group();
    auto iter = group_to_ndx.find(g);
    assert(iter != group_to_ndx.end());
    cerr << "group '" << g << "' found at " << iter->second << " index\n";
    _group[iter->second]->AnotherMember(_fp + i);
  }

  if (cl.option_present('D')) {
    IWString fname = cl.string_value('D');
    if (! ReadDesirabilty(fname)) {
      cerr << "Options::Initialise:cannot read desirability data '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('n')) {
    if (! cl.value('n', _nsel) || _nsel > _group.size()) {
      cerr << "Invalid number to select (-n) must be less than ngroups " <<
              _group.size() << '\n';
      return 0;
    }
  } else {
    _nsel = _group.size();
  }

  if (cl.option_present('w')) {
    IWString w;
    for (int i = 0; cl.value('w', w, i); ++i) {
      if (w.starts_with("dist=")) {
        w.remove_leading_chars(5);
        if (! w.numeric_value(_distance_weight) || _distance_weight < 0.0) {
          cerr << "Options::Initialise:invalid distance weight '" << w << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Distance weight " << _distance_weight << '\n';
        }
      } else if (w.starts_with("desr")) {
        w.remove_leading_chars(5);
        if (! w.numeric_value(_desirability_weight) || _desirability_weight < 0.0) {
          cerr << "Options::Initialise:invalid desirability weight '" << w << "'\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Desirability weight " << _desirability_weight << '\n';
        }
      } else if (w == "help") {
        DisplayDashwOptions(0);
      } else {
        cerr << "Options::Initialise:unrecognised -w qualifier '" << w << "'\n";
        DisplayDashwOptions(1);
      }
    }
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  return 1;
}

int
Options::DetermineInitialDistances() {
  const uint32_t n = _group.size();
  for (uint32_t i = 0; i < n; ++i) {
    Group* gi = _group[i];
    for (uint32_t j = i + 1; j < n; ++j) {
      _group[j]->EstablishInitialDistances(*gi);
    }
  }

  return 1;
}

int
Options::ReadFingerprints(const char* fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadFingerprints:cannot open '" << fname << "'\n";
    return 0;
  }

  _number_fingerprints = input.count_records_starting_with(_identifier_tag);

  if (0 == _number_fingerprints) {
    cerr << "Options::ReadFingerprints:no occurrences of " << _identifier_tag << "' in input\n";
    return 0;
  }

  _fp = new Fingerprint[_number_fingerprints];

  return ReadFingerprints(input);
}

int
Options::ReadFingerprints(iwstring_data_source& input) {
  IW_TDT tdt;

  uint32_t ndx = 0;

  for (; tdt.next(input) && ndx < _number_fingerprints; ndx++) {
//  _fp[ndx].set_initial_ndx(ndx);

//  tdt.dataitem_value(smiles_tag, smiles[ndx]);

//  tdt.dataitem_value(identifier_tag, pcn[ndx]);

    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read fingerprint\n";
      return 0;
    }

    if (0 == ndx) {
      if (!standard_fingerprints_present()) {
        return 0;
      }
    }

    _fp[ndx].build_molecular_properties(gfp.molecular_properties_integer());
    _fp[ndx].build_iw(gfp[0]);
    _fp[ndx].build_mk(gfp[1]);
    _fp[ndx].build_mk2(gfp[2]);
    _fp[ndx].set_name(gfp.id());

    _id_to_ndx[gfp.id()] = ndx;
  }

  return _number_fingerprints;
}

int
Options::ReadGroupMembership(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadGroupMembership:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadGroupMembership(input);
}

int
Options::ReadGroupMembership(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadGroupMembership:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadGroupMembershipRecord(buffer)) {
      cerr << "Options::ReadGroupMembership:error " << buffer << '\n';
      return 0;
    }
  }

  return 1;
}

int
Options::ReadGroupMembershipRecord(const const_IWSubstring& buffer) {
  int i = 0;
  IWString id, group;
  if (! buffer.nextword(id, i) || ! buffer.nextword(group, i) ||
      id.empty() || group.empty()) {
    return 0;
  }

  const auto iter = _id_to_ndx.find(id);
  if (iter == _id_to_ndx.end()) {
    cerr << "Options::ReadGroupMembershipRecord:no fingerprint for '" << id << "'\n";
    return 0;
  }

  _fp[iter->second].set_group(group);

  return 1;
}

int
Options::ReadPreviousDistance(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadPreviousDistance:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadPreviousDistance(input);
}

int
Options::ReadPreviousDistance(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadPreviousDistance:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadPreviousDistanceRecord(buffer)) {
      cerr << "Options::ReadPreviousDistance:error reading '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadPreviousDistanceRecord(const const_IWSubstring& buffer) {
  int i = 0;
  IWString id, dist;
  if (! buffer.nextword(id, i) || ! buffer.nextword(dist, i) ||
      id.empty() || dist.empty()) {
    return 0;
  }

  auto iter = _id_to_ndx.find(id);
  if (iter == _id_to_ndx.end()) {
    cerr << "Options::ReadPreviousDistanceRecord:no fingerprint for '" << id << "'\n";
    return 0;
  }

  float d;
  if (! dist.numeric_value(d) || d < 0.0f) {
    cerr << "Options::ReadPreviousDistanceRecord:invalid numeric '" << dist << "'\n";
    return 0;
  }

  _fp[iter->second].set_distance(d);

  return 1;
}

int
Options::ReadDesirabilty(IWString& fname) {
  iwstring_data_source input(fname.null_terminated_chars());
  if (! input.good()) {
    cerr << "Options::ReadDesirabilty:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadDesirabilty(input);
}

int
Options::ReadDesirabilty(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadDesirabilty:cannot read header\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! ReadDesirabiltyRecord(buffer)) {
      cerr << "Options::ReadDesirabilty:error reading '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::ReadDesirabiltyRecord(const const_IWSubstring& buffer) {
  int i = 0;
  IWString id, desirability;
  if (! buffer.nextword(id, i) || ! buffer.nextword(desirability, i) ||
      id.empty() || desirability.empty()) {
    return 0;
  }

  auto iter = _id_to_ndx.find(id);
  if (iter == _id_to_ndx.end()) {
    cerr << "Options::ReadDesirabiltyRecord:no fingerprint for '" << id << "'\n";
    return 0;
  }

  float d;
  if (! desirability.numeric_value(d) || d < 0.0f) {
    cerr << "Options::ReadDesirabiltyRecord:invalid desirability '" << desirability << "\n";
    return 0;
  }
  
  _fp[iter->second].set_desirability(d);

  return 1;
}

int
Options::Spread(IWString_and_File_Descriptor & output) {
  uint32_t nsel = 0;
  cerr << "Will select " << _nsel << " groups\n";
  for ( ; nsel < _nsel; ++nsel) {
    Group* sel = GetNextGroup();
    cerr << sel << " selected\n";
    if (sel == nullptr) {
      return 1;
    }

    for (Group* g : _group) {
      if (g == sel) {
        continue;
      }
      g->GroupHasBeenSelected(sel);
    }

    sel->set_selected(nsel);
    sel->Write(output);
  }

  // This is not necessary, groups were written as they formed.
  _group.iwqsort_lambda([](const Group* g1, const Group* g2) {
      if (g1->selected() < g2->selected()) {
        return -1;
      }
      if (g1->selected() > g2->selected()) {
        return 1;
      }
      return 0;
    });

  return 1;
}

double
Options::Score(double ave_dist, double ave_desirability) const {
  double rc = _distance_weight * ave_dist +
              _desirability_weight * ave_desirability;
  return rc;
}

Group*
Options::GetNextGroup() const {
  double highest_score = -1.0;
  Group* highest_score_group = nullptr;

  for (Group* g : _group) {
    cerr << " group? " << g->name() << " sel " << g->selected() << '\n';
    if (g->selected()) {
      continue;
    }

    double score = Score(g->AverageDistance(), g->AverageDesirability());
    if (score > highest_score) {
      highest_score = score;
      highest_score_group = g;
    }
  }

  if (highest_score_group == nullptr) {
    return nullptr;
  }

  return highest_score_group;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:G:w:n:D:T:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, cannot handle multiple input files, concatenate them and try again\n";
    Usage(2);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);
  if (! options.Spread(output)) {
    cerr << "gfp_group_spread failure\n";
    return 1;
  }

  output.flush();

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_group_spread

int
main(int argc, char** argv) {
  return gfp_group_spread::Main(argc, argv);
}
