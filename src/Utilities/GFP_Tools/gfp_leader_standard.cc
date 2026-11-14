/*
  Leader implementation with fixed fingerprints
*/

#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

using std::cout;

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/report_progress.h"

#include "omp.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#else
#include "nearneighbours.pb.h"
#endif

#include "gfp_standard.h"
#include "nndata.h"
#include "smiles_id_dist.h"
#include "sparse_collection.h"

using std::cerr;
using std::endl;

const char* prog_name = nullptr;

static int verbose = 0;

static int reassign_max_size = 0;
static bool reassign_to_any = false;

static int clusters_to_form = std::numeric_limits<int>::max();

static Report_Progress report_progress;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");
static IWString mk_tag("FPMK<");
static IWString mk2_tag("FPMK2<");
static IWString iw_tag("FPIW<");

static float threshold = 0.0;

class Leader_Item : public Smiles_ID_Dist {
 private:
  const GFP_Standard* _gfp;

 public:
  Leader_Item();

  void
  set_gfp(const GFP_Standard* s) {
    _gfp = s;
  }
};

Leader_Item::Leader_Item() {
  _gfp = nullptr;

  return;
}

static int pool_size = 0;

static Leader_Item* leader_item = nullptr;
static GFP_Standard* fingerprints = nullptr;
static int* selected = nullptr;
static float* distances = nullptr;

static int write_neighbours_as_protos = 0;

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
  cerr << "Leader implementation requiring MPR IW MK MK2\n";
  cerr << " -t <dist>      distance threshold\n";
  cerr << " -R <n>         reassign clusters of at most <n> members to nearest surviving leader\n";
  cerr << " -R <n>:any     search all members of larger clusters; ignores threshold and size cap\n";
  cerr << " -n <number>    number of clusters to form\n";
  cerr << " -A <fname>     file of previously selected fingerprints\n";
  cerr << " -s <size>      size of fingerprint file\n";
  cerr << " -p <n>         1 (default) or 2 leaders per scan; same clustering results\n";
  cerr << " -h <n>         maximum number of OMP threads to use\n";
  cerr << " -r <n>         report progress every <n> clusters formed\n";
  cerr << " -k .           write neighbours as textproto form\n";
  cerr << " -v             verbose output\n";
  // clang-format on

  exit(rc);
}

// Fingerprints stay in their original storage for output and reassignment.
// Only this input-ordered list of candidate indices is squeezed. Selected slots
// are marked -1 during parallel scans, then removed when they occupy a quarter
// of the remaining range. Compaction preserves priority and bounds wasted work.
class ActiveCandidates {
 private:
  std::vector<int> _indices;
  int _first = 0;
  int _holes = 0;

 public:
  ActiveCandidates() {
    _indices.reserve(pool_size);
    for (int i = 0; i < pool_size; ++i) {
      if (!selected[i]) {  // Exclude -A matches before the first scan.
        _indices.push_back(i);
      }
    }
  }

  int first() const { return _first; }
  int size() const { return _indices.size(); }
  int* data() { return _indices.data(); }

  int NextLeaders(int* leaders, int requested) const {
    int found = 0;
    for (int i = _first; i < size() && found < requested; ++i) {
      if (_indices[i] >= 0) {
        leaders[found++] = _indices[i];
      }
    }
    return found;
  }

  void RemoveLeaders(int count) {
    while (count > 0) {
      if (_indices[_first++] < 0) {
        --_holes;
      } else {
        --count;
      }
    }
  }

  // Called after the parallel scan has completed; no worker retains a pointer
  // when storage is resized or indices are moved.
  void RemoveSelected(int count) {
    _holes += count;
    if (_holes == 0 || _holes < (size() - _first + 3) / 4) {
      return;
    }
    int destination = 0;
    for (int i = _first; i < size(); ++i) {
      if (_indices[i] >= 0) {
        _indices[destination++] = _indices[i];
      }
    }
    _indices.resize(destination);
    _first = 0;
    _holes = 0;
  }
};

static void
AssignToLeader(ActiveCandidates& candidates, int centre, int cluster) {
  int* indices = candidates.data();
  const int first = candidates.first();
  const int stop = candidates.size();
  int removed = 0;
#pragma omp parallel for schedule(static) reduction(+:removed) if (stop - first > 500)
  for (int pos = first; pos < stop; ++pos) {
    const int item = indices[pos];
    if (item < 0) {
      continue;
    }
    const float d = fingerprints[centre].tanimoto_distance(fingerprints[item]);
    if (d <= threshold) {
      selected[item] = cluster;
      distances[item] = d;
      indices[pos] = -1;
      ++removed;
    }
  }
  candidates.RemoveSelected(removed);
}

static void
AssignToTwoLeaders(ActiveCandidates& candidates, const int* leaders, int first_cluster) {
  int* indices = candidates.data();
  const int first = candidates.first();
  const int stop = candidates.size();
  int removed = 0;
#pragma omp parallel for schedule(static) reduction(+:removed) if (stop - first > 500)
  for (int pos = first; pos < stop; ++pos) {
    const int item = indices[pos];
    if (item < 0) {
      continue;
    }
    float d[2];
    fingerprints[item].tanimoto_distance_2(fingerprints[leaders[0]],
                                         fingerprints[leaders[1]], d);
    if (d[0] <= threshold) {
      selected[item] = first_cluster;
      distances[item] = d[0];
    } else if (d[1] <= threshold) {
      selected[item] = first_cluster + 1;
      distances[item] = d[1];
    } else {
      continue;
    }
    indices[pos] = -1;
    ++removed;
  }
  candidates.RemoveSelected(removed);
}

static int
leader(IWString_and_File_Descriptor& output) {
  ActiveCandidates candidates;
  int cluster_id = 0;
  int centre;
  while (cluster_id < clusters_to_form && candidates.NextLeaders(&centre, 1)) {
    ++cluster_id;
    if (report_progress()) {
      cerr << "Formed " << cluster_id << " clusters, "
           << count_non_zero_occurrences_in_array(selected, pool_size)
           << " items in clusters\n";
    }
    candidates.RemoveLeaders(1);
    selected[centre] = cluster_id;
    distances[centre] = 0.0f;
    AssignToLeader(candidates, centre, cluster_id);
  }
  return cluster_id;
}

static int
leader_two_at_once() {
  ActiveCandidates candidates;
  int cluster_id = 0;
  int leaders[2];
  int found;
  while (cluster_id < clusters_to_form &&
         (found = candidates.NextLeaders(leaders, std::min(2, clusters_to_form - cluster_id)))) {
    if (found == 2 && fingerprints[leaders[0]].tanimoto_distance(
                         fingerprints[leaders[1]]) <= threshold) {
      found = 1;  // Leave the second candidate available for this leader's scan.
    }
    candidates.RemoveLeaders(found);
    for (int i = 0; i < found; ++i) {
      selected[leaders[i]] = cluster_id + i + 1;
      distances[leaders[i]] = 0.0f;
    }
    if (found == 1) {
      AssignToLeader(candidates, leaders[0], cluster_id + 1);
    } else {
      AssignToTwoLeaders(candidates, leaders, cluster_id + 1);
    }
    cluster_id += found;
  }
  return cluster_id;
}

static int
do_previously_selected(const IW_General_Fingerprint& gfp, float threshold) {
  GFP_Standard sgfp;

  sgfp.build_molecular_properties(gfp.molecular_properties_integer());
  sgfp.build_mk(gfp[kMkFp]);
  sgfp.build_mk2(gfp[kMk2Fp]);
  sgfp.build_iw(gfp[kIwFp]);

  for (int i = 0; i < pool_size; i++) {
    if (selected[i]) {
      continue;
    }

    float d = 1.0f - sgfp.tanimoto(fingerprints[i]);

    if (d > threshold) {
      continue;
    }

    selected[i] = -1;
  }

  return 1;
}

static int
do_previously_selected(iwstring_data_source& input, float threshold) {
  IW_TDT tdt;

  while (tdt.next(input)) {
    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read previously selected fingerprint\n";
      return 0;
    }

    do_previously_selected(gfp, threshold);
  }

  return 1;
}

static int
do_previously_selected(const char* fname, float threshold) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected(input, threshold);
}

static int
read_pool(iwstring_data_source& input) {
  IW_TDT tdt;

  IWString tmp;

  int ndx = 0;

  for (; tdt.next(input) && ndx < pool_size; ndx++) {
    leader_item[ndx].set_gfp(fingerprints + ndx);

    tdt.dataitem_value(smiles_tag, tmp);
    leader_item[ndx].set_smiles(tmp);

    tdt.dataitem_value(identifier_tag, tmp);
    leader_item[ndx].set_id(tmp);

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

    fingerprints[ndx].build_molecular_properties(gfp.molecular_properties_integer());
    fingerprints[ndx].build_mk(gfp[1]);
    fingerprints[ndx].build_mk2(gfp[2]);
    fingerprints[ndx].build_iw(gfp[0]);
  }

  if (ndx < 2) {
    cerr << "Yipes, did not read enough fingerprints\n";
    return 0;
  }

  pool_size = ndx;

  return ndx;
}

static int
WriteNeighboursAsProtos(int clusters_formed, IWString_and_File_Descriptor& output) {
  int start = 0;
  for (int c = 1; c <= clusters_formed; c++) {
    int items_in_cluster = 0;

    for (int j = start; j < pool_size; j++) {
      if (selected[j] != c) {
        continue;
      }

      if (0 == items_in_cluster) {
        start = j;
      }

      items_in_cluster++;
    }

    nnbr::NearNeighbours cluster;
    const IWString& smi = leader_item[start].smiles();
    cluster.set_smiles(smi.data(), smi.length());
    const IWString& id = leader_item[start].id();
    cluster.set_name(id.data(), id.length());
    cluster.set_cluster(c - 1);

    start++;

    for (int j = start; j < pool_size; j++) {
      if (c != selected[j]) {
        continue;
      }

      nnbr::Nbr* nbr = cluster.add_nbr();
      const IWString& smi = leader_item[j].smiles();
      nbr->set_smi(smi.data(), smi.length());
      const IWString& id = leader_item[j].id();
      nbr->set_id(id.data(), id.length());
      nbr->set_dist(distances[j]);
    }

    if (!gfp::WriteNNData(cluster, output)) {
      cerr << "Cannot write cluster " << (c - 1) << '\n';
      return 0;
    }
  }

  return 1;
}

// Membership and leaders describe the initial clustering. Only initial members
// of larger clusters (or just their leaders) are targets; assignments are independent.
static int
ReassignSmallClusters(const std::vector<int>& sizes, const std::vector<int>& leaders) {
  std::vector<int> candidates;
  std::vector<int> targets;
  for (int i = 0; i < pool_size; ++i) {
    if (selected[i] <= 0) {  // Excluded by -A or not clustered due to -n.
      continue;
    }
    if (sizes[selected[i]] <= reassign_max_size) {
      candidates.push_back(i);
    } else if (reassign_to_any || leaders[selected[i]] == i) {
      targets.push_back(i);
    }
  }
  if (candidates.empty()) {
    return 0;
  }
  if (targets.empty()) {
    cerr << "Warning: no clusters larger than " << reassign_max_size
         << "; retaining original clustering\n";
    return 0;
  }
  std::vector<int> assignments(candidates.size());
#pragma omp parallel for schedule(static) if (candidates.size() > 32)
  for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
    const int item = candidates[i];
    float closest = std::numeric_limits<float>::max();
    int destination = -1;
    // Targets are in input order; strict comparison retains the first equal hit.
    for (int j : targets) {
      const float d = fingerprints[item].tanimoto_distance(fingerprints[j]);
      if (d < closest) {
        closest = d;
        destination = selected[j];
      }
    }
    assignments[i] = destination;
    // DIST always means distance to the leader, not to the closest member.
    distances[item] = reassign_to_any ? fingerprints[item].tanimoto_distance(
        fingerprints[leaders[destination]]) : closest;
  }
  for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
    selected[candidates[i]] = assignments[i];
  }
  if (verbose) {
    cerr << "Reassigned " << candidates.size() << " members of clusters of size <= "
         << reassign_max_size << '\n';
  }
  return candidates.size();
}

// After reassignment input order no longer identifies the leader. Build linked
// member lists in one pass, then write each item once, keeping original leaders.
static int
WriteReassignedClusters(int clusters_formed, const std::vector<int>& leaders,
                        IWString_and_File_Descriptor& output) {
  std::vector<int> head(clusters_formed + 1, -1);
  std::vector<int> next(pool_size, -1);
  std::vector<int> sizes(clusters_formed + 1, 0);
  for (int i = pool_size - 1; i >= 0; --i) {
    if (selected[i] <= 0) {
      continue;
    }
    const int c = selected[i];
    next[i] = head[c];
    head[c] = i;
    ++sizes[c];
  }
  int number = 0;
  for (int c = 1; c <= clusters_formed; ++c) {
    if (head[c] < 0) {
      continue;
    }
    const int centre = leaders[c];
    nnbr::NearNeighbours proto;
    if (write_neighbours_as_protos) {
      const IWString& smi = leader_item[centre].smiles();
      const IWString& id = leader_item[centre].id();
      proto.set_smiles(smi.data(), smi.length());
      proto.set_name(id.data(), id.length());
      proto.set_cluster(number++);
    } else {
      output << smiles_tag << leader_item[centre].smiles() << ">\n";
      output << identifier_tag << leader_item[centre].id() << ">\n";
      output << "CLUSTER<" << number++ << ">\nCSIZE<" << sizes[c] << ">\n";
    }
    for (int j = head[c]; j >= 0; j = next[j]) {
      if (j == centre) {
        continue;
      }
      if (write_neighbours_as_protos) {
        nnbr::Nbr* nbr = proto.add_nbr();
        const IWString& smi = leader_item[j].smiles();
        const IWString& id = leader_item[j].id();
        nbr->set_smi(smi.data(), smi.length());
        nbr->set_id(id.data(), id.length());
        nbr->set_dist(distances[j]);
      } else {
        output << smiles_tag << leader_item[j].smiles() << ">\n";
        output << identifier_tag << leader_item[j].id() << ">\n";
        output << distance_tag << distances[j] << ">\n";
        output.write_if_buffer_holds_more_than(8192);
      }
    }
    if (write_neighbours_as_protos) {
      if (!gfp::WriteNNData(proto, output)) {
        return 0;
      }
    } else {
      output << "|\n";
    }
  }
  return 1;
}

static int
allocate_pool(int s) {
  leader_item = new Leader_Item[s];
  selected = new_int(s);
  fingerprints = new GFP_Standard[s];
  distances = new float[s];

  if (nullptr == leader_item || nullptr == selected || nullptr == fingerprints ||
      nullptr == distances) {
    cerr << "Yipes, could not allocate " << pool_size << " fingerprints\n";
    return 0;
  }

  pool_size = s;

  return 1;
}

static int
gfp_leader_standard(int argc, char** argv) {
  Command_Line cl(argc, argv, "vt:A:s:p:r:h:n:k:R:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  set_report_fingerprint_status(0);

  if (!cl.option_present('t')) {
    cerr << "Must specify cluster threshold via the -t option\n";
    usage(2);
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', threshold) || threshold <= 0.0 || threshold >= 1.0) {
      cerr << "The threshold value must be a valid distance\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Threshold set to " << threshold << endl;
    }
  }

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 2) {
      cerr << "The pool size specification (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Problem sized for " << pool_size << " fingerprints\n";
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', clusters_to_form) || clusters_to_form < 1) {
      cerr << "The number of clusters to form (-n) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will stop processing after forming " << clusters_to_form << " clusters\n";
    }
  }

  if (cl.option_present('R')) {
    const_IWSubstring specification = cl.option_value('R');
    reassign_to_any = specification.ends_with(":any");
    if (reassign_to_any) {
      specification.chop(4);
    }
    if (!specification.numeric_value(reassign_max_size) || reassign_max_size < 1) {
      cerr << "The -R option requires a positive size, optionally followed by ':any'\n";
      usage(2);
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise report progress option (-r)\n";
      usage(2);
    }
  }

  if (cl.option_present('k')) {
    write_neighbours_as_protos = 1;
    if (verbose) {
      cerr << "Will generate output in proto form\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(2);
  }

  if (cl.number_elements() > 1) {
    cerr << "Sorry, cannot handle multiple input files, concatenate them and try again\n";
    usage(2);
  }

  const char* fname = cl[0];

  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open fingerprint file '" << fname << "'\n";
    return 2;
  }

  if (0 == pool_size) {
    pool_size = input.count_records_starting_with(identifier_tag);

    if (0 == pool_size) {
      cerr << "No occurrences of " << identifier_tag << "' in input\n";
      return 0;
    }

    if (!allocate_pool(pool_size)) {
      return 3;
    }

    cerr << "Job automatically sized to " << pool_size << " fingerprints\n";
  } else {
    if (!allocate_pool(pool_size)) {
      return 3;
    }
  }

  if (!read_pool(input)) {
    cerr << "Cannot read fingerprints\n";
    return 0;
  }

  if (cl.option_present('h')) {
    int h;
    if (!cl.value('h', h) || h < 0) {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve "
              "number\n";
      usage(2);
    }

    omp_set_num_threads(h);
    if (verbose) {
      cerr << "Running with " << h << " threads\n";
    }
  }

  if (cl.option_present('A')) {
    const char* fname = cl.option_value('A');

    if (!do_previously_selected(fname, threshold)) {
      cerr << "Cannot process previously selected file '" << fname << "'\n";
      return 0;
    }
  }

  int clusters_formed;

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('p')) {
    int p;
    if (!cl.value('p', p) || p <= 0 || p > 2) {
      cerr << "The only valid values for the -p option are 1 and 2\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Will process " << p << " leaders at each step\n";
    }

    clusters_formed = p == 1 ? leader(output) : leader_two_at_once();
  } else {
    clusters_formed = leader(output);
  }

  if (verbose) {
    cerr << "Clustered " << pool_size << " fingerprints into " << clusters_formed
         << " clusters\n";
  }

  if (reassign_max_size > 0) {
    std::vector<int> sizes(clusters_formed + 1, 0);
    std::vector<int> leaders(clusters_formed + 1, -1);
    for (int i = 0; i < pool_size; ++i) {
      if (selected[i] <= 0) {
        continue;
      }
      const int c = selected[i];
      if (sizes[c]++ == 0) {
        leaders[c] = i;
      }
    }
    if (ReassignSmallClusters(sizes, leaders) > 0) {
      return WriteReassignedClusters(clusters_formed, leaders, output) ? 0 : 1;
    }
  }

  if (write_neighbours_as_protos) {
    WriteNeighboursAsProtos(clusters_formed, output);
  }

  int start = 0;
  for (int c = 1; c <= clusters_formed; c++) {
    int items_in_cluster = 0;

    for (int j = start; j < pool_size; j++) {
      if (selected[j] != c) {
        continue;
      }

      if (0 == items_in_cluster) {
        start = j;
      }

      items_in_cluster++;
    }

    output << smiles_tag << leader_item[start].smiles() << ">\n";
    output << identifier_tag << leader_item[start].id() << ">\n";
    output << "CLUSTER<" << (c - 1) << ">\n";
    output << "CSIZE<" << items_in_cluster << ">\n";

    start++;

    if (start == pool_size) {
      output << "|\n";
      break;
    }

    for (int j = start; j < pool_size; j++) {
      if (c != selected[j]) {
        continue;
      }

      output << smiles_tag << leader_item[j].smiles() << ">\n";
      output << identifier_tag << leader_item[j].id() << ">\n";
      output << distance_tag << distances[j] << ">\n";
      output.write_if_buffer_holds_more_than(8192);
    }

    output << "|\n";
  }

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = gfp_leader_standard(argc, argv);

  return rc;
}
