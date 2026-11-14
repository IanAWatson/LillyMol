/*
  Spread implementation with fixed fingerprints
*/

#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <optional>
#include <random>

#define USE_OMP
#ifdef USE_OMP
#include "omp.h"
#endif

#define REPORT_PROGRESS_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/numeric_data_from_file.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "gfp_standard.h"

using std::cerr;

const char* prog_name = nullptr;

static int verbose = 0;

static int choose_first_point_randomly = 0;

static int nsel = 0;

static float stop_once_distance_drops_below = 0.0f;

static int previously_selected_file_specified = 0;

static int previous_computed_distance_column = -1;

static IWString previously_computed_nn_distance_tag;

static int brief_output = 0;

static Numeric_Data_From_File<float> previously_computed_distances;

static Report_Progress report_progress;

static int squeeze = 0;
static int next_squeeze = 0;

static IWString smiles_tag("$SMI<");
static IWString identifier_tag("PCN<");
static IWString distance_tag("DIST<");
static IWString scale_tag("SCALE<");

static Fraction_as_String fraction_as_string;

// Squeezing moves active fingerprints, while SMILES and identifiers retain
// input order. Each item therefore keeps its original metadata index.

class SSpread_Item : public GFP_Standard {
 private:
  int _selected;

  float _distance;

  int _nearest_previously_selected;

  int _initial_ndx;

 public:
  SSpread_Item();

  int selected() const {
    return _selected;
  }

  void set_selected(int s) {
    _selected = s;
  }

  float distance() const {
    return _distance;
  }
  // Initial -A/-N distances need not have a known in-pool neighbour identity.
  void set_distance(float d) {
    _distance = d;
  }

  int nearest_previously_selected() const {
    return _nearest_previously_selected;
  }

  SSpread_Item& operator=(const SSpread_Item& rhs);
  void set_initial_ndx(int s) {
    _initial_ndx = s;
  }

  int initial_ndx() const {
    return _initial_ndx;
  }

  int SelectedOrZeroDistance() const;

  int NotifySelectedUseThreshold(const SSpread_Item& rhs) {
    std::optional<float> d = GFP_Standard::tanimoto_distance_if_less(rhs, _distance);
    if (!d) {
      return 0;
    }
    _distance = *d;
    _nearest_previously_selected = rhs.initial_ndx();
    return 1;
  }

  const IWString& smiles(const IWString* smiles) const {
    return smiles[_initial_ndx];
  }

  const IWString& pcn(const IWString* pcn) const {
    return pcn[_initial_ndx];
  }

};

SSpread_Item::SSpread_Item() {
  _selected = 0;
  _distance = 1.0;
  _initial_ndx = -1;
  _nearest_previously_selected = -1;
}

SSpread_Item&
SSpread_Item::operator=(const SSpread_Item& rhs) {
  _selected = rhs._selected;
  _distance = rhs._distance;
  _nearest_previously_selected = rhs._nearest_previously_selected;
  _initial_ndx = rhs._initial_ndx;
  GFP_Standard* me = this;
  const GFP_Standard* other = &rhs;
  *me = *other;
  return *this;
}

int
SSpread_Item::SelectedOrZeroDistance() const {
  if (_selected || 0.0f == _distance) {
    return 1;
  }

  return 0;
}

// A struct to describe the results of the computation.
struct Selected_Item {
  int _sel;
  int _nsn;
  float _dist;
  float _scale;
};

static int pool_size = 0;

// Each of these will be sized for the job.
static SSpread_Item* fingerprints = nullptr;
static IWString* smiles = nullptr;
static IWString* pcn = nullptr;
static int* selected = nullptr;
static float* distances = nullptr;
static int* nearest_previously_selected = nullptr;
static float* weight = nullptr;
static Selected_Item* selected_item = nullptr;

// clang-format off
static void
usage(int rc)
{
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
// clang-format on
// clang-format off
  cerr << R"(Spread implementation requiring MPR IW MK MK2 fingerprints, the gfp_make default.
 -n <nsel>       number of items to select
 -A <fname>      file of previously selected fingerprints
 -p FILE=<fname> read weight values from <fname>
 -p COL=<col>    weight values are column <col> in the name (1-based, >=2)
 -p FCOL=<col>   weight column in FILE= input (1-based including id; default 2)
 -N <tag>        gfp_nearneighbours has been run and initial distances are in <tag> "
                  (also FILE=, COL=)
 -S rand         choose first item randomly
 -s <size>       size of fingerprint file
 -t <dist>       stop selection once distance drops below <dist>
 -h <n>          maximum number of OMP threads to use
 -r <n>          report progress every <n> items selected
 -q <n>          squeeze out already selected items every <n> selections (def 1000)
 -b              brief output (smiles id dist)
 -v              verbose output
)";
// clang-format on

  exit(rc);
}

// clang-format on

static int
get_previously_computed_nearest_neighbour(
    const IWString& id, const IW_STL_Hash_Map_float& previously_computed_distances,
    float& d) {
  IWString key(id);
  key.truncate_at_first(' ');
  IW_STL_Hash_Map_float::const_iterator f = previously_computed_distances.find(key);

  if (f ==
      previously_computed_distances.end()) {  // OK if no previously computed distance
    return 1;
  }

  d = f->second;
  return 1;
}

static int
get_previously_computed_nearest_neighbour(
    const IW_TDT& tdt, const IWString& previously_computed_nn_distance_tag, float& d) {
  if (!tdt.dataitem_value(previously_computed_nn_distance_tag, d)) {
    return 0;
  }

  return 1;
}

static int
get_previously_computed_from_name_column(const IWString& id,
                                         int previous_computed_distance_column,
                                         float& d) {
  const_IWSubstring c;
  if (!id.word(previous_computed_distance_column, c)) {
    cerr << "get_previously_computed_from_name_column:cannot extract column "
         << (previous_computed_distance_column + 1) << " from '" << id << "'\n";
    return 0;
  }

  if (!c.numeric_value(d) || d < 0.0f) {
    cerr << "get_previously_computed_from_name_column:invalid numeric '" << id << "'\n";
    return 0;
  }

  return 1;
}

static void
set_selected_item(Selected_Item& s, const int isel, const int nsn, const float d) {
  s._sel = isel;
  s._nsn = nsn;
  s._dist = d;

  return;
}

int
RandomPoolMember(int pool_size) {
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution<int> u(0, pool_size - 1);
  return u(generator);
}

// Neighbour identities are original indices, so they remain valid when active
// fingerprints and their per-item state move. Selected output also stores its
// scale explicitly rather than looking it up in the squeezed weight array.
template <typename F>
static int
SqueezeWeighted(F* pool, int& size, int next, float* weights) {
  if (size < 1000) {
    return next;
  }
  int destination = 0;
  int result = -1;
  for (int i = 0; i < size; ++i) {
    if (selected[i]) {
      continue;
    }
    if (i == next) {
      result = destination;
    }
    if (destination != i) {
      pool[destination] = pool[i];
      selected[destination] = 0;
      distances[destination] = distances[i];
      nearest_previously_selected[destination] = nearest_previously_selected[i];
      weights[destination] = weights[i];
    }
    ++destination;
  }
  size = destination;
  return result;
}

template <typename F>
int
spread_with_weights(F* fingerprints, int& pool_size, float* weight,
                    Selected_Item* selected_item) {
  int s;

  // Everywhere else in this function distances[] holds *weighted* values - see
  // the update in the selection loop below. The -A seed distances arrive
  // unweighted, so scale them once here and the seed is then on the same footing
  // as every subsequent update. Done before the branch so that -S rand -A gets
  // weighted seeds too.
  if (previously_selected_file_specified) {
    for (int i = 0; i < pool_size; ++i) {
      distances[i] *= weight[i];
    }
  }

  if (choose_first_point_randomly) {
    s = RandomPoolMember(pool_size);
  } else if (previously_selected_file_specified) {
    s = 0;
    float dmax = distances[0];
    for (int i = 1; i < pool_size; i++) {
      if (distances[i] > dmax) {
        dmax = distances[i];
        s = i;
      }
    }
  } else {
    s = 0;
  }

  if (!previously_selected_file_specified) {
    for (int i = 0; i < pool_size; i++) {
      distances[i] = weight[i] * fingerprints[s].tanimoto_distance(fingerprints[i]);
      assert(std::isfinite(distances[i]) && distances[i] >= 0.0f);
      nearest_previously_selected[i] = fingerprints[s].initial_ndx();
    }

    distances[s] = 1.0f;
    nearest_previously_selected[s] = -1;

  }

  int items_written = 0;

  for (int items_selected = 0; items_selected < nsel; items_selected++) {

    // Neighbour ids are original indices, including after squeezing. With -A
    // or -N, the first selection can have no known in-pool neighbour (-1).
    const int nps = nearest_previously_selected[s];

    set_selected_item(selected_item[items_selected], fingerprints[s].initial_ndx(),
                      nps, distances[s]);

    selected_item[items_selected]._scale = weight[s];
    items_written = items_selected + 1;

    selected[s] = 1;
    if (items_written == nsel) {
      break;
    }

    int nexts_shared = -1;
    float dmax_shared = -1.0f;

#pragma omp parallel shared(nexts_shared, dmax_shared)
    {
      float mydmax = -1.0f;
      int mynexts = -1;

#pragma omp for schedule(dynamic, 256) nowait
      for (int i = 0; i < pool_size; i++) {
        if (selected[i] || 0.0f == distances[i]) {
          continue;
        }

        float d = weight[i] * fingerprints[s].tanimoto_distance(fingerprints[i]);
        if (d < distances[i]) {
          distances[i] = d;
          nearest_previously_selected[i] = fingerprints[s].initial_ndx();
        }

        if (distances[i] > mydmax ||
            (distances[i] == mydmax && (mynexts < 0 || i < mynexts))) {
          mydmax = distances[i];
          mynexts = i;
        }
      }

#pragma omp critical
      {
        if (mydmax > dmax_shared ||
            (mydmax == dmax_shared && mynexts >= 0 &&
             (nexts_shared < 0 || mynexts < nexts_shared))) {
          dmax_shared = mydmax;
          nexts_shared = mynexts;
        }
      }
    }

    if (nexts_shared < 0) {
      // Zero-distance (including zero-weight) items still fill the requested
      // count when no positive stopping threshold has been specified.
      if (stop_once_distance_drops_below == 0.0f) {
        for (int i = 0; i < pool_size && items_written < nsel; ++i) {
          if (selected[i]) {
            continue;
          }
          set_selected_item(selected_item[items_written], fingerprints[i].initial_ndx(),
                            nearest_previously_selected[i], distances[i]);
          selected_item[items_written]._scale = weight[i];
          selected[i] = 1;
          ++items_written;
        }
      }
      break;
    }
    if (dmax_shared < stop_once_distance_drops_below) {
      break;
    }

    s = nexts_shared;
    if (next_squeeze == items_written) {
      s = SqueezeWeighted(fingerprints, pool_size, s, weight);
      next_squeeze += squeeze;
    }

    if (report_progress()) {
      cerr << "Selected " << items_selected << " fingerprints\n";
    }
  }

  // The caller writes selected_item[0 .. return value); the loop above can stop
  // early, and the untouched tail of the array is uninitialised.
  return items_written;
}

template <typename T>
int
do_squeeze(T* pool, int& pool_size, int next_selected) {
  if (pool_size < 1000) {
    return next_selected;
  }

  int ndx = 0;
  int rc = -1;
  for (int i = 0; i < pool_size; ++i) {
    if (pool[i].selected()) {
      continue;
    }

    if (i == next_selected) {
      rc = ndx;
    }

    if (i == ndx) {
      ++ndx;
      continue;
    }

    pool[ndx] = pool[i];
    ++ndx;
  }

  pool_size = ndx;

  return rc;
}

template <typename T>
int
WriteSelectedBrief(const T* pool, int selected, IWString_and_File_Descriptor& output) {
  const T& sel = pool[selected];
  output << sel.smiles(smiles) << ' ' << sel.pcn(pcn) << ' ';
  fraction_as_string.append_number(output, sel.distance());
  output << '\n';
  output.write_if_buffer_holds_more_than(4096);
  return 1;
}

template <typename T>
int
WriteSelected(T* pool, int selected, IWString_and_File_Descriptor& output) {
  if (brief_output) {
    return WriteSelectedBrief(pool, selected, output);
  }

  const T& sel = pool[selected];

  output << smiles_tag << sel.smiles(smiles) << ">\n";
  output << identifier_tag << sel.pcn(pcn) << ">\n";

  if (sel.nearest_previously_selected() >= 0) {
    int ndx = sel.nearest_previously_selected();
    output << smiles_tag << smiles[ndx] << ">\n";
    output << identifier_tag << pcn[ndx] << ">\n";
    output << distance_tag;
    fraction_as_string.append_number(output, sel.distance());
    output << ">\n";
  } else {
    output << smiles_tag << "*>\n";
    output << identifier_tag << "*>\n";
    output << distance_tag;
    fraction_as_string.append_number(output, sel.distance());
    output << ">\n";
  }

  output << "|\n";

  output.write_if_buffer_holds_more_than(4096);

  return 1;
}

// Return the index of the member of `pool` with the longest distance.
int
MaxInitialDistance(const SSpread_Item* pool, int pool_size) {
  float max_dist = -1.0f;
  int rc = -1;

  for (int i = 0; i < pool_size; ++i) {
    if (pool[i].distance() > max_dist) {
      max_dist = pool[i].distance();
      rc = i;
    }
  }

  return rc;
}

template <typename F>
int
spread_unweighted(F* pool, int& pool_size, IWString_and_File_Descriptor& output) {
  int s;

  if (choose_first_point_randomly) {
    s = RandomPoolMember(pool_size);
  } else if (previously_selected_file_specified) {
    s = MaxInitialDistance(pool, pool_size);
  } else {
    s = 0;
  }

  for (int items_selected = 0; items_selected < nsel; items_selected++) {
    WriteSelected(pool, s, output);
    pool[s].set_selected(1);
    if (items_selected + 1 == nsel) {
      return nsel;
    }

    // Variables shared between threads.
    int nexts_shared = -1;
    float dmax_shared = -1.0f;

#pragma omp parallel shared(nexts_shared, dmax_shared)
    {
      float mydmax = -1.0f;
      int mynexts = -1;

#pragma omp for schedule(static) nowait
      for (int i = 0; i < pool_size; i++) {
        if (pool[i].SelectedOrZeroDistance()) {
          continue;
        }

        pool[i].NotifySelectedUseThreshold(pool[s]);
        if (pool[i].distance() > mydmax ||
            (pool[i].distance() == mydmax && (mynexts < 0 || i < mynexts))) {
          mydmax = pool[i].distance();
          mynexts = i;
        }
      }

#pragma omp critical
      {
        if (mydmax > dmax_shared ||
            (mydmax == dmax_shared && mynexts >= 0 &&
             (nexts_shared < 0 || mynexts < nexts_shared))) {
          dmax_shared = mydmax;
          nexts_shared = mynexts;
        }
      }
    }

    if (nexts_shared < 0) {  // we are ignoring nsel, unlikely to be a problem...
      int written = items_selected + 1;
      if (stop_once_distance_drops_below == 0.0f) {
        for (int i = 0; i < pool_size && written < nsel; ++i) {
          if (!pool[i].selected()) {
            WriteSelected(pool, i, output);
            pool[i].set_selected(1);
            ++written;
          }
        }
      }
      return written;
    }

    s = nexts_shared;

    if (dmax_shared < stop_once_distance_drops_below) {
      return items_selected + 1;
    }

    if (report_progress()) {
      cerr << "Selected " << items_selected << " fingerprints, dist " << dmax_shared
           << '\n';
    }

    if (next_squeeze == items_selected + 1) {
      s = do_squeeze(pool, pool_size, s);
      next_squeeze += squeeze;
    }
  }

  return nsel;
}

static int
assign_scaling_factors_from_token(const IWString* pcn, const int pool_size,
                                  int scaling_factor_column, float* weight) {
  for (int i = 0; i < pool_size; i++) {
    const_IWSubstring s;

    if (!pcn[i].word(scaling_factor_column, s)) {
      cerr << "Cannot extract column " << (scaling_factor_column + 1) << " from '"
           << pcn[i] << "'\n";
      return 0;
    }

    float v;

    if (!s.numeric_value(v) || !std::isfinite(v) || v < 0.0f) {
      cerr << "Invalid scaling factor '" << pcn[i] << "'\n";
      return 0;
    }

    weight[i] = v;
  }

  return 1;
}

static int
assign_scaling_factors_from_file(iwstring_data_source& input, int col,
                                 const IWString* pcn, int verbose) {
  input.set_dos(1);
  input.set_translate_tabs(1);
  IW_STL_Hash_Map_int id_to_ndx;
  for (int i = 0; i < pool_size; ++i) {
    IWString id = pcn[i];
    id.truncate_at_first(' ');
    id_to_ndx[id] = i;
  }
  std::vector<unsigned char> seen(pool_size, 0);
  const_IWSubstring buffer;
  int count = 0;
  while (input.next_record(buffer)) {
    IWString id;
    const_IWSubstring token;
    if (!buffer.word(0, id) || !buffer.word(col, token)) {
      cerr << "Cannot extract weight column " << col + 1 << " from '" << buffer << "'\n";
      return 0;
    }
    const auto found = id_to_ndx.find(id);
    float value;
    if (!token.numeric_value(value)) {
      if (input.lines_read() == 1 && found == id_to_ndx.end()) {
        continue;  // Optional header, not a malformed weight for a known id.
      }
      cerr << "Invalid fingerprint weight '" << buffer << "'\n";
      return 0;
    }
    if (!std::isfinite(value) || value < 0.0f) {
      cerr << "Invalid fingerprint weight '" << buffer << "'\n";
      return 0;
    }
    if (found == id_to_ndx.end()) {
      if (verbose > 2) {
        cerr << "Weight specified for '" << id << "' but no fingerprint\n";
      }
      continue;
    }
    const int ndx = found->second;
    if (seen[ndx]) {
      cerr << "Duplicate fingerprint weight for '" << id << "'\n";
      return 0;
    }
    seen[ndx] = 1;
    weight[ndx] = value;
    ++count;
  }
  if (count != pool_size) {
    cerr << "Only got " << count << " weight values for " << pool_size << " fingerprints\n";
    return 0;
  }
  return 1;
}

static int
assign_scaling_factors_from_file(const char* fname, int col, const IWString* pcn,
                                 int verbose) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open scaling data file '" << fname << "'\n";
    return 0;
  }

  return assign_scaling_factors_from_file(input, col, pcn, verbose);
}

static void
do_previously_selected(const IW_General_Fingerprint& gfp) {
  GFP_Standard sgfp;

  sgfp.build_molecular_properties(gfp.molecular_properties_integer());
  sgfp.build_iw(gfp[0]);
  sgfp.build_mk(gfp[1]);
  sgfp.build_mk2(gfp[2]);

#pragma omp parallel for schedule(static)
  for (int i = 0; i < pool_size; i++) {
    float d = sgfp.tanimoto_distance(fingerprints[i]);

    assert(d >= 0.0f);

    if (d < distances[i]) {
      distances[i] = d;
    }
  }

  return;
}

static int
do_previously_selected(iwstring_data_source& input) {
  IW_TDT tdt;

  Report_Progress my_report_progress;
  if (report_progress.active()) {
    my_report_progress.set_report_every(report_progress.report_every());
  }

  int fingerprints_read = 0;

  while (tdt.next(input)) {
    IW_General_Fingerprint gfp;

    int fatal;
    if (!gfp.construct_from_tdt(tdt, fatal)) {
      cerr << "Cannot read previously selected fingerprint\n";
      return 0;
    }

    do_previously_selected(gfp);

    ++fingerprints_read;

    if (my_report_progress()) {
      cerr << "Processed " << my_report_progress.times_called() << " previously selected fingerprints\n";
    }
  }

  // A file that yields no fingerprints - empty, or not TDT at all - leaves every
  // initial distance at 1.0, which is indistinguishable from -A never having been
  // specified. Fail rather than silently ignore the option.
  if (fingerprints_read == 0) {
    cerr << "do_previously_selected:no fingerprints read\n";
    return 0;
  }

  if (verbose) {
    cerr << "Read " << fingerprints_read << " previously selected fingerprints\n";
  }

  return 1;
}

static int
do_previously_selected(const char* fname) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open previously selected file '" << fname << "'\n";
    return 0;
  }

  return do_previously_selected(input);
}

template <typename F>
int
read_pool(iwstring_data_source& input, F* pool, int& pool_size) {
  IW_TDT tdt;

  int ndx = 0;

  for (; tdt.next(input) && ndx < pool_size; ndx++) {
    pool[ndx].set_initial_ndx(ndx);

    tdt.dataitem_value(smiles_tag, smiles[ndx]);

    tdt.dataitem_value(identifier_tag, pcn[ndx]);

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
    fingerprints[ndx].build_iw(gfp[0]);
    fingerprints[ndx].build_mk(gfp[1]);
    fingerprints[ndx].build_mk2(gfp[2]);

    if (previously_computed_nn_distance_tag.length()) {
      if (!get_previously_computed_nearest_neighbour(
              tdt, previously_computed_nn_distance_tag, distances[ndx])) {
        return 0;
      }
    } else if (previously_computed_distances.size()) {
      if (!get_previously_computed_nearest_neighbour(
              pcn[ndx], previously_computed_distances, distances[ndx])) {
        return 0;
      }
    } else if (previous_computed_distance_column >= 0) {
      if (!get_previously_computed_from_name_column(
              pcn[ndx], previous_computed_distance_column, distances[ndx])) {
        return 0;
      }
    }
  }

  for (int i = 0; i < ndx; ++i) {
    if (pool[i].initial_ndx() != i) {
      cerr << "Initial index wrong\n";
    }
  }

  if (ndx < 2) {
    cerr << "Yipes, did not read enough fingerprints\n";
    return 0;
  }

  pool_size = ndx;

  return pool_size;
}

static int
allocate_pool(int s, int scaling_factor_specified) {
  selected = new_int(s);
  fingerprints = new SSpread_Item[s];
  distances = new float[s];
  std::fill_n(distances, s, 1.0f);
  nearest_previously_selected = new_int(s, -1);
  smiles = new IWString[s];
  pcn = new IWString[s];
  selected_item = new Selected_Item[s];

  if (nullptr == selected || nullptr == fingerprints || nullptr == distances ||
      nullptr == nearest_previously_selected || nullptr == smiles) {
    cerr << "Yipes, could not allocate " << pool_size << " fingerprints\n";
    return 0;
  }

  if (scaling_factor_specified) {
    weight = new_float(s, 1.0f);
    if (nullptr == weight) {
      cerr << "Cannot allocate weight array\n";
      return 0;
    }
  }

  pool_size = s;

  return 1;
}

static int
gfp_spread_standard(int argc, char** argv) {
  Command_Line cl(argc, argv, "vA:s:r:n:bh:p:S:N:t:q:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  set_report_fingerprint_status(0);

  set_default_iwstring_float_concatenation_precision(4);

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 2) {
      cerr << "The pool size specification (-s) must be a whole +ve number\n";
      usage(2);
    }

    if (verbose) {
      cerr << "Problem sized for " << pool_size << " fingerprints\n";
    }
  }

  if (cl.option_present('r')) {
    if (!report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialise report progress option (-r)\n";
      usage(2);
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

  int scaling_factor_specified = 0;
  IWString scaling_factor_file_name;
  int scaling_factor_column_in_file = 1;
  int scaling_factor_column = -1;

  if (cl.option_present('p')) {
    IWString col;

    int i = 0;
    const_IWSubstring p;

    while (cl.value('p', p, i++)) {
      if (p.starts_with("FILE=") && 0 == scaling_factor_file_name.length()) {
        p.remove_leading_chars(5);
        scaling_factor_file_name = p;
      } else if (p.starts_with("FCOL=")) {
        p.remove_leading_chars(5);
        if (!p.numeric_value(scaling_factor_column_in_file) ||
            scaling_factor_column_in_file < 2) {
          cerr << "Invalid scaling factor column number 'FCOL=" << p << "'\n";
          return 2;
        }
        scaling_factor_column_in_file--;
      } else if (p.starts_with("COL=") && 0 == col.length()) {
        p.remove_leading_chars(4);
        col = p;
      } else {
        cerr << "Unrecognised -p qualifier '" << p << "'\n";
        usage(4);
      }
    }

    if (scaling_factor_file_name.length() && col.length()) {
      cerr << "Cannot specify both scaling factor file name (FILE=) and scaling factor "
              "column "
              "(COL=)\n";
      usage(1);
    }

    if (col.length()) {
      if (!col.numeric_value(scaling_factor_column) || scaling_factor_column < 2) {
        cerr << "Invalid scaling factor column '" << col << "'\n";
        return 4;
      }

      scaling_factor_column--;
    }

    if (!scaling_factor_file_name.length() && !col.length()) {
      cerr << "The -p option requires FILE= or COL=\n";
      return 2;
    }
    scaling_factor_specified = 1;
  }

  if (cl.option_present('N')) {
    int nset = 0;

    IWString fname;

    int i = 0;
    const_IWSubstring n;

    while (cl.value('N', n, i++)) {
      if (n.starts_with("FILE=")) {
        n.remove_leading_chars(5);
        fname = n;
        nset++;
      } else if (n.starts_with("COL=")) {
        n.remove_leading_chars(4);
        if (!n.numeric_value(previous_computed_distance_column) ||
            previous_computed_distance_column < 1) {
          cerr << "The previously computed distance column must be a whole +ve number\n";
          return 3;
        }
        if (verbose) {
          cerr << "Previously computed near neighbour distances in column "
               << previous_computed_distance_column << '\n';
        }

        previous_computed_distance_column--;
        nset++;
      } else if (0 == previously_computed_nn_distance_tag.length()) {
        previously_computed_nn_distance_tag = n;
        if (verbose) {
          cerr << "Previously computed near neighbour distances in the '"
               << previously_computed_nn_distance_tag << "' tag\n";
        }

        if (!previously_computed_nn_distance_tag.ends_with('<')) {
          previously_computed_nn_distance_tag += '<';
        }

        nset++;
      } else {
        cerr << "Unrecognised -N qualifier '" << n << "'\n";
        usage(3);
      }
    }

    if (nset > 1) {
      cerr << "Can specify just one of FILE=, COL= or tag for previously computed "
              "distances\n";
      usage(3);
    }

    if (fname.length()) {
      if (!previously_computed_distances.read_data(fname)) {
        cerr << "Cannot read previously computed distances from '" << fname << "'\n";
        return 4;
      }

      if (verbose) {
        cerr << "Read " << previously_computed_distances.size()
             << " previously computed nn distances from '" << fname << "'\n";
      }
    }

    previously_selected_file_specified = 1;
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
      return 3;
    }

    if (!allocate_pool(pool_size, scaling_factor_specified)) {
      return 3;
    }

    if (verbose) {
      cerr << "Job automatically sized to " << pool_size << " fingerprints\n";
    }
  } else {
    if (!allocate_pool(pool_size, scaling_factor_specified)) {
      return 3;
    }
  }

  if (!read_pool(input, fingerprints, pool_size)) {
    cerr << "Cannot read fingerprints\n";
    return 3;
  }

  if (scaling_factor_column >= 1) {
    if (!assign_scaling_factors_from_token(pcn, pool_size, scaling_factor_column,
                                           weight)) {
      cerr << "Cannot extract scaling factors from column " << scaling_factor_column
           << '\n';
      return 3;
    }
  } else if (scaling_factor_file_name) {
    if (!assign_scaling_factors_from_file(scaling_factor_file_name,
                                          scaling_factor_column_in_file, pcn, verbose)) {
      cerr << "Cannot assign scaling factors from '" << scaling_factor_file_name << "'\n";
      return 3;
    }
  }

#ifdef USE_OMP
  if (cl.option_present('h')) {
    int h;
    if (!cl.value('h', h) || h < 1) {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve "
              "number\n";
      usage(2);
    }

    omp_set_num_threads(h);
  }
#endif

  if (cl.option_present('N')) {
    for (int i = 0; i < pool_size; ++i) {
      if (!std::isfinite(distances[i]) || distances[i] < 0.0f || distances[i] > 1.0f) {
        cerr << "Invalid initial distance for '" << pcn[i] << "'\n";
        return 3;
      }
    }
  }

  if (cl.option_present('A')) {
    const char* fname = cl.option_value('A');

    if (!do_previously_selected(fname)) {
      cerr << "Cannot process previously selected file '" << fname << "'\n";
      return 3;
    }

    if (verbose) {
      cerr << "Finished processing already selected file\n";
    }

    previously_selected_file_specified = 1;
  }

  if (previously_selected_file_specified) {
    for (int i = 0; i < pool_size; ++i) {
      if (!std::isfinite(distances[i]) || distances[i] < 0.0f || distances[i] > 1.0f) {
        cerr << "Invalid initial distance for '" << pcn[i] << "': " << distances[i] << '\n';
        return 3;
      }
      fingerprints[i].set_distance(distances[i]);
    }
  }

  if (cl.option_present('S')) {
    const_IWSubstring s = cl.string_value('S');

    if ("rand" == s) {
      choose_first_point_randomly = 1;

      if (verbose) {
        cerr << "Will choose first item randomly\n";
      }
    } else {
      cerr << "Unrecognised -S qualifier '" << s << "'\n";
      return 2;
    }
  }

  if (cl.option_present('n')) {
    if (!cl.value('n', nsel) || nsel < 1) {
      cerr << "The number if fingerprints to select must be a whole +ve number\n";
      usage(2);
    } else if (nsel > pool_size) {
      cerr << "ONly " << pool_size << " fingerprints read, nsel " << nsel
           << " impossible\n";
      return 3;
    }
  } else {
    nsel = pool_size;
  }

  if (cl.option_present('q')) {
    if (!cl.value('q', squeeze) || squeeze < 0) {
      cerr << "The squeeze every option (-q) must be a whole non negative number\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will squeeze selected items every " << squeeze << " items selected\n";
    }

    if (squeeze == 0) {
      next_squeeze = -1;
    } else {
      next_squeeze = squeeze;
    }
  } else {
    squeeze = 1000;
    next_squeeze = 1000;
  }

  if (cl.option_present('t')) {
    if (!cl.value('t', stop_once_distance_drops_below) ||
        stop_once_distance_drops_below < 0.0f || stop_once_distance_drops_below >= 1.0f) {
      cerr << "The stop selection distance option (-t) must be a valid distance\n";
      usage(1);
    }

    if (verbose) {
      cerr << "Will stop selection once distance drops below "
           << stop_once_distance_drops_below << '\n';
    }
  }

  // Share raw numeric strings across DIST, SCALE, and brief output. Values
  // outside [0, 1] use Fraction_as_String's stable formatting fallback.
  if (!fraction_as_string.initialise(0.0f, 1.0f, 4)) {
    cerr << "Cannot initialise distance and scale formatting\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  // Set before selection: spread_unweighted writes via WriteSelected, which
  // consults brief_output.
  if (cl.option_present('b')) {
    brief_output = 1;

    if (verbose) {
      cerr << "Brief output\n";
    }
  }

  cerr.flush();

  // Ordinary spread writes as it selects. Weighted spread fills selected_item[]
  // and leaves writing to the block below.
  int items_selected = 0;
  if (scaling_factor_specified) {
    items_selected =
        spread_with_weights(fingerprints, pool_size, weight, selected_item);
  } else {
    spread_unweighted(fingerprints, pool_size, output);
  }

  if (verbose) {
    cerr << "Writing results\n";
  }

  // Only the weighted path defers output. Use its actual selection count:
  // an early stop leaves the rest of selected_item[] uninitialised.
  if (!scaling_factor_specified) {
    // spread_unweighted has already written everything.
  } else if (brief_output) {
    for (int i = 0; i < items_selected; ++i) {
      const Selected_Item& s = selected_item[i];

      output << smiles[s._sel] << ' ' << pcn[s._sel] << ' ';
      fraction_as_string.append_number(output, s._dist);
      output << '\n';
      output.write_if_buffer_holds_more_than(4096);
    }
  } else {
    for (int i = 0; i < items_selected; ++i) {
      const Selected_Item& s = selected_item[i];

      const int id = s._sel;
      const int nsn = s._nsn;

      output << smiles_tag << smiles[id] << ">\n";
      output << identifier_tag << pcn[id] << ">\n";
      if (nsn >= 0) {
        output << smiles_tag << smiles[nsn] << ">\n";
        output << identifier_tag << pcn[nsn] << ">\n";
        output << scale_tag;
        fraction_as_string.append_number(output, s._scale);
        output << ">\n";
        output << distance_tag;
        fraction_as_string.append_number(output, s._dist);
        output << ">\n";
      } else {
        output << smiles_tag << "*>\n";
        output << identifier_tag << "*>\n";
        output << scale_tag << "*>\n";
        output << distance_tag;
        fraction_as_string.append_number(output, s._dist);
        output << ">\n";
      }
      output << "|\n";

      output.write_if_buffer_holds_more_than(4096);
    }
  }

  output.flush();

  delete[] fingerprints;
  delete[] smiles;
  delete[] pcn;
  delete[] selected;
  delete[] distances;
  delete[] nearest_previously_selected;
  delete[] weight;
  delete[] selected_item;
  delete_gfp_file_scope_static_objects();

  return 0;
}

int
main(int argc, char** argv) {
  prog_name = argv[0];

  int rc = gfp_spread_standard(argc, argv);

  return rc;
}
