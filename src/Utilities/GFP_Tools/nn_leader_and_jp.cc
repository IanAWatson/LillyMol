/*
  We have run gfp_nearneighbours and can quickly do several clusterings.
*/

#include <stdlib.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/set_or_unset.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#ifdef BUILD_BAZEL
#include "Utilities/GFP_Tools/nearneighbours.pb.h"
#else
#include "nearneighbours.pb.h"
#endif

using std::cerr;

typedef float similarity_type_t;

static int verbose = 0;

static Fraction_as_String fraction_as_string;

static float maximum_distance = 1.0f;

static int write_tdt_output = 1;

/*
  The variables which control the clustering
*/

static int max_clusters_to_find = std::numeric_limits<int>::max();
static int max_cluster_size = 0;

/*
  We can process any number of thresholds
*/

static int number_thresholds = 0;
static similarity_type_t* threshold = nullptr;

/*
  The identifier tag used in each TDT
*/

static IWString smiles_tag("$SMI<");

static IWString identifier_tag("PCN<");
static IWString neighbour_tag("NBR<");

static IWString distance_tag("DIST<");

static int max_cluster_size_column = -1;

static IWString cluster_number_tag("CLUSTER<");
static IWString cluster_size_tag("CSIZE<");

static int taylor_butina = 0;

static int suppress_self_neighbours = 0;

class JarvisPatrickParameters {
 private:
  int _nbr_list_size;
  int _nbrs_in_common;

 public:
  JarvisPatrickParameters();

  int
  active() const {
    return _nbrs_in_common > 0;
  }

  int
  nbr_list_size() const {
    return _nbr_list_size;
  }

  int
  nbrs_in_common() const {
    return _nbrs_in_common;
  }

  int build(const const_IWSubstring&);
};

JarvisPatrickParameters::JarvisPatrickParameters() {
  _nbr_list_size = 0;
  _nbrs_in_common = 0;

  return;
}

int
JarvisPatrickParameters::build(const const_IWSubstring& buffer) {
  const_IWSubstring snl, snb;
  if (!buffer.split(snl, ',', snb) || 0 == snl.length() || 0 == snb.length()) {
    cerr << "JarvisPatrickParameters::build:JP parameters must be of the form 'J,K' "
            "where J is the neighbour list size\n";
    cerr << "  and K is the number of neighbours in common\n";
    return 0;
  }

  if ('.' == snl)  // leave as default
    ;
  else if (!snl.numeric_value(_nbr_list_size) || _nbr_list_size < 1) {
    cerr << "JarvisPatrickParameters::build:invalid neighbour list size '" << buffer
         << "'\n";
    return 0;
  }

  if (!snb.numeric_value(_nbrs_in_common) || _nbrs_in_common < 1) {
    cerr << "JarvisPatrickParameters::build:invalid neighbours in common size '" << buffer
         << "'\n";
    return 0;
  }

  if (0 == _nbr_list_size)
    ;
  else if (_nbr_list_size < _nbrs_in_common) {
    cerr << "JarvisPatrickParameters::build:invalid combination '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

class Leader_Item {
 private:
  int _nbrs;

  int* _nbr;
  float* _dist;

  //  Each item can limit the number of items in its cluster

  Set_or_Unset<int> _max_cluster_size;

  int _ndx;

  //  private functions

  void RemoveExcludedNeighbours();

  int Build(const nnbr::NearNeighbours& proto,
                   const IW_STL_Hash_Map_int& id_to_ndx,
                   int n, bool ignore_missing_neighbours);
 public:
  Leader_Item();
  ~Leader_Item();

  int
  neighbour_count() const {
    return _nbrs;
  }

  void
  set_ndx(int s) {
    _ndx = s;
  }

  void update_max_distance(float& d) const;

  int
  max_cluster_size(int& m) const {
    return _max_cluster_size.value(m);
  }

  int build(iwstring_data_source& input, const IW_STL_Hash_Map_int& id_to_ndx,
            int n, bool ignore_missing_neighbours);
  int Build(iw_tf_data_record::TFDataReader& input,
                   const IW_STL_Hash_Map_int& id_to_ndx,
                   int n, bool ignore_missing_neighbours);

  int form_cluster_threshold(int* selected, const float threshold,
                             const int cluster_number,
                             const resizable_array_p<IWString>& smiles,
                             const resizable_array_p<IWString>& id,
                             IWString_and_File_Descriptor& output);
  int form_cluster_max_cluster_size(int* selected, const float threshold,
                                    const int max_cluster_size, const int cluster_number,
                                    const resizable_array_p<IWString>& smiles,
                                    const resizable_array_p<IWString>& id,
                                    IWString_and_File_Descriptor& output);

  int mark_neighbours(int* n, const int maxn, const int x) const;
  bool neighbours_in_common(const int* n, const JarvisPatrickParameters&,
                            const int x) const;

  int unselected_neighbours(const int* selected) const;
  int unselected_neighbours(const int* selected, const similarity_type_t within) const;
};

Leader_Item::Leader_Item() {
  _nbrs = 0;
  _nbr = nullptr;
  _dist = nullptr;
  _ndx = -1;

  return;
}

Leader_Item::~Leader_Item() {
  if (nullptr != _nbr) {
    delete[] _nbr;
  }

  if (nullptr != _dist) {
    delete[] _dist;
  }

  return;
}

// Normalise a loaded list once, keeping indices and distances aligned. The
// clustering kernels need no checks for omitted neighbours or suppressed self IDs.
void
Leader_Item::RemoveExcludedNeighbours() {
  int retained = 0;
  for (int i = 0; i < _nbrs; ++i) {
    if (_nbr[i] < 0 || (suppress_self_neighbours && _nbr[i] == _ndx)) {
      continue;
    }
    _nbr[retained] = _nbr[i];
    _dist[retained] = _dist[i];
    ++retained;
  }
  _nbrs = retained;
}

void
Leader_Item::update_max_distance(float& d) const {
  if (0 == _nbrs) {
    return;
  }

  if (_dist[_nbrs - 1] > d) {
    d = _dist[_nbrs - 1];
  }

  return;
}

// Per-record readers must distinguish an error from clean EOF. The first
// pass must never turn a damaged suffix into a successfully loaded prefix.
enum class RecordReadResult { kRecord, kEof, kError };

static int
data_from_tdt(const const_IWSubstring buffer, resizable_array_p<IWString>& s) {
  const int openangle = buffer.index('<');
  if (openangle < 0 || !buffer.ends_with('>')) {
    return 0;
  }

  IWString* tmp = new IWString;
  // An empty SMILES is allowed, but from_to requires a nonempty range.
  if (buffer.length() > openangle + 2) {
    buffer.from_to(openangle + 1, buffer.length() - 2, *tmp);
  }
  s.add(tmp);
  return 1;
}

static RecordReadResult
ReadSmilesAndIdsAndNbrCount(const nnbr::NearNeighbours& proto,
                           resizable_array_p<IWString>& smiles,
                           resizable_array_p<IWString>& id,
                           resizable_array<int>& ncount) {
  if (proto.name().empty() ||
      std::isspace(static_cast<unsigned char>(proto.name().front()))) {
    cerr << "ReadSmilesAndIdsAndNbrCount:missing or empty identifier\n";
    return RecordReadResult::kError;
  }
  id << new IWString(proto.name());
  smiles << new IWString(proto.smiles());
  ncount << proto.nbr_size();
  return RecordReadResult::kRecord;
}

static RecordReadResult
ReadSmilesAndIdsAndNbrCount(iw_tf_data_record::TFDataReader& input,
                           resizable_array_p<IWString>& smiles,
                           resizable_array_p<IWString>& id,
                           resizable_array<int>& ncount) {
  // Read the frame separately so a protobuf parse failure cannot be confused
  // with EOF, even if the reader has already buffered the final frame.
  const auto data = input.Next();
  if (!data) {
    if (input.good() && input.eof()) {
      return RecordReadResult::kEof;
    }
    cerr << "ReadSmilesAndIdsAndNbrCount:record read error after "
         << input.items_read() << " records\n";
    return RecordReadResult::kError;
  }
  nnbr::NearNeighbours proto;
  if (!proto.ParseFromString(std::string_view(data->data(), data->length()))) {
    cerr << "ReadSmilesAndIdsAndNbrCount:cannot parse protobuf record "
         << input.items_read() << '\n';
    return RecordReadResult::kError;
  }
  return ReadSmilesAndIdsAndNbrCount(proto, smiles, id, ncount);
}

static RecordReadResult
read_smiles_and_ids_and_nbr_count(iwstring_data_source& input,
                                resizable_array_p<IWString>& smiles,
                                resizable_array_p<IWString>& id,
                                resizable_array<int>& ncount) {
  const_IWSubstring buffer;
  bool got_smiles = false;
  bool got_id = false;
  bool got_record_data = false;
  bool terminated = false;
  int nbrs = 0;

  while (input.next_record(buffer)) {
    if (buffer.empty()) {
      continue;
    }
    got_record_data = true;
    if (buffer == '|') {
      terminated = true;
      break;
    }
    if (buffer.starts_with(smiles_tag)) {
      if (!got_smiles) {
        if (!data_from_tdt(buffer, smiles)) {
          cerr << "Invalid SMILES dataitem at line " << input.lines_read() << '\n';
          return RecordReadResult::kError;
        }
        got_smiles = true;
      }
    } else if (buffer.starts_with(identifier_tag)) {
      // Validate parent and neighbour IDs before substring extraction. In
      // particular, an empty parent must not be mistaken for a missing one.
      if (buffer.length() <= identifier_tag.length() + 1 || !buffer.ends_with('>') ||
          std::isspace(static_cast<unsigned char>(buffer[identifier_tag.length()]))) {
        cerr << "Invalid or empty identifier at line " << input.lines_read() << '\n';
        return RecordReadResult::kError;
      }
      if (!got_id) {
        if (!data_from_tdt(buffer, id)) {
          return RecordReadResult::kError;
        }
        got_id = true;
      } else {
        ++nbrs;
      }
    }
  }

  if (!input.good()) {
    cerr << "TDT read error near line " << input.lines_read() << '\n';
    return RecordReadResult::kError;
  }
  if (!terminated) {
    if (!got_record_data && input.eof()) {
      return RecordReadResult::kEof;
    }
    cerr << "Unterminated neighbour record near line " << input.lines_read() << '\n';
    return RecordReadResult::kError;
  }
  if (!got_smiles || !got_id) {
    cerr << "Neighbour record missing SMILES or identifier near line "
         << input.lines_read() << '\n';
    return RecordReadResult::kError;
  }
  ncount.add(nbrs);
  return RecordReadResult::kRecord;
}

static int
read_smiles_and_ids_and_nbr_count(iwstring_data_source& input, const int max_pool_size,
                                resizable_array_p<IWString>& smiles,
                                resizable_array_p<IWString>& id,
                                resizable_array<int>& ncount) {
  while (smiles.number_elements() < max_pool_size) {
    const auto status = read_smiles_and_ids_and_nbr_count(input, smiles, id, ncount);
    if (status == RecordReadResult::kError) {
      return 0;
    }
    if (status == RecordReadResult::kEof) {
      break;
    }
  }
  if (smiles.empty()) {
    cerr << "read_smiles_and_ids_and_nbr_count:no data\n";
    return 0;
  }
  return smiles.number_elements();
}

static int
ReadSerializedProtos(iw_tf_data_record::TFDataReader& input, const int max_pool_size,
                     resizable_array_p<IWString>& smiles,
                     resizable_array_p<IWString>& id,
                     resizable_array<int>& ncount) {
  while (smiles.number_elements() < max_pool_size) {
    const auto status = ReadSmilesAndIdsAndNbrCount(input, smiles, id, ncount);
    if (status == RecordReadResult::kError) {
      return 0;
    }
    if (status == RecordReadResult::kEof) {
      break;
    }
  }
  if (smiles.empty()) {
    cerr << "ReadSerializedProtos:no data\n";
    return 0;
  }
  return smiles.number_elements();
}

static int
read_smiles_and_ids_and_nbr_count(const char* fname, const int max_pool_size,
                                  resizable_array_p<IWString>& smiles,
                                  resizable_array_p<IWString>& id,
                                  resizable_array<int>& ncount) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "read_smiles_and_ids_and_nbr_count:cannot open '" << fname << "'\n";
    return 0;
  }

  return read_smiles_and_ids_and_nbr_count(input, max_pool_size, smiles, id, ncount);
}

static int
ReadSerializedProtos(const char* fname, const int max_pool_size,
                                  resizable_array_p<IWString>& smiles,
                                  resizable_array_p<IWString>& id,
                                  resizable_array<int>& ncount) {
  iw_tf_data_record::TFDataReader input(fname);
  if (! input.good()) {
    cerr << "ReadSerializedProtos:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSerializedProtos(input, max_pool_size, smiles, id, ncount);
}

static int
get_nbr_id(const const_IWSubstring& buffer, const IW_STL_Hash_Map_int& id_to_ndx,
           int& nbr, bool ignore_missing_neighbours) {
  int openangle = buffer.index('<');
  if (openangle < 0) {
    return 0;
  }

  IWString id;
  buffer.from_to(openangle + 1, buffer.length() - 2, id);

  id.truncate_at_first(' ');

  const auto f = id_to_ndx.find(id);

  if (f == id_to_ndx.end()) {
    if (ignore_missing_neighbours) {
      nbr = -1;
      return 1;
    }
    cerr << "get_nbr_id:no number for '" << id << "'\n";
    return 0;
  }

  nbr = f->second;

  // cerr << "From '" << buffer << "' id '" << id << "' nbr " << nbr << '\n';

  return 1;
}

static int
GetNbrId(const std::string& id, const IW_STL_Hash_Map_int& id_to_ndx,
         int& nbr, bool ignore_missing_neighbours) {
  IWString s(id);

  s.truncate_at_first(' ');

  const auto f = id_to_ndx.find(s);

  if (f == id_to_ndx.end()) {
    if (ignore_missing_neighbours) {
      nbr = -1;
      return 1;
    }
    cerr << "GetNbrId:no number for '" << id << "'\n";
    return 0;
  }

  nbr = f->second;

  return 1;
}

template <typename T>
int
convert_to_float(const const_IWSubstring& buffer, T& d, const T minv, const T maxv) {
  const int openangle = buffer.index('<');
  if (openangle < 0 || buffer.length() <= openangle + 2 || !buffer.ends_with('>')) {
    cerr << "convert_to_float:invalid numeric dataitem '" << buffer << "'\n";
    return 0;
  }

  const_IWSubstring s;
  buffer.from_to(openangle + 1, buffer.length() - 2, s);

  if (!s.numeric_value(d) || !std::isfinite(d) || d < minv || d > maxv) {
    cerr << "convert_to_float:invalid distance '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
Leader_Item::build(iwstring_data_source& input, const IW_STL_Hash_Map_int& id_to_ndx,
                   int n, bool ignore_missing_neighbours) {
  _nbrs = n;
  _nbr = new int[_nbrs];
  _dist = new float[_nbrs];

  IWString id;

  int ndx_nbrs = 0;  // index into _nbrs
  int ndx_dist = 0;  // index into _dist array

  bool terminated = false;
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    //  cerr << "Leader_Item::build:reading '" << buffer << "'\n";

    if ('|' == buffer) {
      terminated = true;
      break;
    }

    if (buffer.starts_with(identifier_tag)) {
      if (0 == id.length()) {
        // identifier_tag already includes '<'; the value starts immediately after it.
        buffer.from_to(identifier_tag.length(), buffer.length() - 2, id);
      } else if (ndx_nbrs >= _nbrs) {
        cerr << "Leader_Item::build:too many neighbour identifiers for '" << id
             << "', expected " << _nbrs << ", line " << input.lines_read() << '\n';
        return 0;
      } else if (!get_nbr_id(buffer, id_to_ndx, _nbr[ndx_nbrs], ignore_missing_neighbours)) {
        cerr << "Leader_Item::build:invalid neighbour specification '" << buffer << "'\n";
        return 0;
      } else {
        ndx_nbrs++;
      }
    } else if (buffer.starts_with(distance_tag)) {
      if (ndx_dist >= _nbrs) {
        cerr << "Leader_Item::build:too many distances for '" << id
             << "', expected " << _nbrs << ", line " << input.lines_read() << '\n';
        return 0;
      }
      if (!convert_to_float(buffer, _dist[ndx_dist], 0.0f, maximum_distance)) {
        cerr << "Leader_Item::build:invalid distance '" << buffer << "'\n";
        return 0;
      }

      ndx_dist++;
    }
  }

  if (!terminated || id.empty() || !input.good()) {
    cerr << "Leader_Item::build:incomplete or unreadable neighbour record at line "
         << input.lines_read() << '\n';
    return 0;
  }

  if (ndx_nbrs != _nbrs || ndx_dist != _nbrs) {
    cerr << "Leader_Item::build:expected " << _nbrs << " neighbours for '" << id
         << "', read " << ndx_nbrs << " identifiers and " << ndx_dist << " distances\n";
    return 0;
  }

  if (ignore_missing_neighbours || suppress_self_neighbours) {
    RemoveExcludedNeighbours();
  }

#ifdef DEBUG_BUILD
  cerr << "Item " << _ndx << " has " << _nbrs << " neighbours\n";
  for (int i = 0; i < _nbrs; ++i) {
    cerr << ' ' << _nbr[i] << '\n';
  }
#endif

  if (max_cluster_size_column >= 0) {
    IWString s;
    if (!id.word(max_cluster_size_column, s)) {
      cerr << "Leader_Item::build:cannot extract column " << (max_cluster_size_column + 1)
           << " from '" << id << "'\n";
      return 0;
    }

    int m;
    if (!s.numeric_value(m) || m < 1) {
      cerr << "Leader_Item::build:invalid max cluster size " << id << " column "
           << (max_cluster_size_column + 1) << '\n';
      return 0;
    }

    _max_cluster_size.set(m);
  }

  return 1;
}


int
Leader_Item::Build(iw_tf_data_record::TFDataReader& input,
                   const IW_STL_Hash_Map_int& id_to_ndx,
                   int n, bool ignore_missing_neighbours) {
  std::optional<nnbr::NearNeighbours> maybe_proto = 
        input.ReadProto<nnbr::NearNeighbours>();
  if (! maybe_proto) {
    return 0;
  }

  return Build(*maybe_proto, id_to_ndx, n, ignore_missing_neighbours);
}

int
Leader_Item::Build(const nnbr::NearNeighbours& proto,
                   const IW_STL_Hash_Map_int& id_to_ndx,
                   int n, bool ignore_missing_neighbours) {

  if (proto.nbr_size() != n) {
    cerr << "Leader_Item::Build:size mismatch, expected " << n << " got " <<
            proto.nbr_size() << '\n';
    return 0;
  }

  _nbrs = n;
  _nbr = new int[_nbrs];
  _dist = new float[_nbrs];

  int ndx = 0;  // Index into the arrays.

  for (const nnbr::Nbr& nbr : proto.nbr()) {
    if (!GetNbrId(nbr.id(), id_to_ndx, _nbr[ndx], ignore_missing_neighbours)) {
      cerr << "Leader_Item::build:invalid neighbour specification '" << nbr.id() << "'\n";
      return 0;
    }

    float d = nbr.dist();
    if (!std::isfinite(d) || d < 0.0f || d > maximum_distance) {
      cerr << "Leader_Item::Build:invalid distance " << proto.ShortDebugString() << '\n';
      return 0;
    }
    _dist[ndx] = d;

    ++ndx;
  }

  if (ignore_missing_neighbours || suppress_self_neighbours) {
    RemoveExcludedNeighbours();
  }

#ifdef DEBUG_BUILD
  cerr << "Item " << _ndx << " has " << _nbrs << " neighbours\n";
  for (int i = 0; i < _nbrs; ++i) {
    cerr << ' ' << _nbr[i] << '\n';
  }
#endif

  if (max_cluster_size_column >= 0) {
    const IWString id = proto.name();
    IWString s;
    if (!id.word(max_cluster_size_column, s)) {
      cerr << "Leader_Item::build:cannot extract column " << (max_cluster_size_column + 1)
           << " from '" << id << "'\n";
      return 0;
    }

    int m;
    if (!s.numeric_value(m) || m < 1) {
      cerr << "Leader_Item::build:invalid max cluster size " << id << " column "
           << (max_cluster_size_column + 1) << '\n';
      return 0;
    }

    _max_cluster_size.set(m);
  }

  return 1;
}


static int
write_cluster_header(const int ndx, const resizable_array_p<IWString>& smiles,
                     const resizable_array_p<IWString>& id, const int cluster_number,
                     const int csize, IWString_and_File_Descriptor& output) {
  if (write_tdt_output) {
    output << smiles_tag << *(smiles[ndx]) << ">\n";
    output << identifier_tag << *(id[ndx]) << ">\n";
    output << cluster_number_tag << cluster_number << ">\n";
    output << cluster_size_tag << csize << ">\n";
  } else {
    output << *(smiles[ndx]) << ' ' << *(id[ndx]) << " CLUSTER " << cluster_number << " ("
           << csize << " members)\n";
  }

  return 1;
}

static void
write_smiles_id_dist(const IWString* smiles, const IWString* id, const float d,
                     IWString_and_File_Descriptor& output) {
  if (write_tdt_output) {
    output << smiles_tag << *(smiles) << ">\n";
    output << identifier_tag << *(id) << ">\n";
    if (d >= 0.0f) {
      output << distance_tag;
      fraction_as_string.append_number(output, d);
      output << ">\n";
    }
  } else {
    output << *smiles << ' ' << *id;
    if (d >= 0.0f) {
      output << ' ';
      fraction_as_string.append_number(output, d);
    }
    output << '\n';
  }

  output.write_if_buffer_holds_more_than(32768);

  return;
}

static int
finish_cluster(const int csize, IWString_and_File_Descriptor& output) {
  if (write_tdt_output) {
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(32768);

  return csize;
}

int
Leader_Item::form_cluster_threshold(int* selected, const float threshold,
                                    const int cluster_number,
                                    const resizable_array_p<IWString>& smiles,
                                    const resizable_array_p<IWString>& id,
                                    IWString_and_File_Descriptor& output) {
  // First need to find the cluster size so we can write the CSIZE tag.
  int csize = 1;

  // This will be the range of the second scan through the neighbour list.
  int istop = _nbrs;
  for (int i = 0; i < _nbrs; ++i) {
    if (_dist[i] > threshold) {
      istop = i;
      break;
    }

    const int j = _nbr[i];

    if (selected[j]) {
      continue;
    }

    csize++;
  }

  write_cluster_header(_ndx, smiles, id, cluster_number, csize, output);

  if (1 == csize) {
    selected[_ndx] = 1;
    return finish_cluster(1, output);
  }

  for (int i = 0; i < istop; ++i) {
    const int j = _nbr[i];

    if (selected[j]) {
      continue;
    }

    write_smiles_id_dist(smiles[j], id[j], _dist[i], output);

    selected[j] = 1;
  }

  selected[_ndx] = 1;

  return finish_cluster(csize, output);
}

int
Leader_Item::form_cluster_max_cluster_size(int* selected, const float threshold,
                                           const int max_cluster_size,
                                           const int cluster_number,
                                           const resizable_array_p<IWString>& smiles,
                                           const resizable_array_p<IWString>& id,
                                           IWString_and_File_Descriptor& output) {
  // The limit includes the leader. Any self neighbours suppressed by -h
  // have already been removed during ingestion.
  int csize = 1;
  int istop = 0;
  for (; istop < _nbrs && csize < max_cluster_size; ++istop) {
    if (_dist[istop] > threshold) {
      break;
    }
    if (!selected[_nbr[istop]]) {
      ++csize;
    }
  }

  write_cluster_header(_ndx, smiles, id, cluster_number, csize, output);

  // Reuse the same neighbour prefix so counting and writing apply the same
  // threshold, limit and selection rules. A limit of one leaves it empty.
  for (int i = 0; i < istop; ++i) {
    const int j = _nbr[i];
    if (selected[j]) {
      continue;
    }

    write_smiles_id_dist(smiles[j], id[j], _dist[i], output);
    selected[j] = 1;
  }

  selected[_ndx] = 1;
  return finish_cluster(csize, output);
}

int
Leader_Item::mark_neighbours(int* n, const int maxn, const int x) const {
  int istop = _nbrs;
  if (0 == maxn)  // use all neighbours
    ;
  else if (istop > maxn) {
    istop = maxn;
  }

  for (int i = 0; i < istop; ++i) {
    n[_nbr[i]] = x;
  }

  return _nbrs;
}

bool
Leader_Item::neighbours_in_common(const int* n, const JarvisPatrickParameters& jp,
                                  const int x) const {
  int istop = _nbrs;
  if (0 == jp.nbr_list_size())
    ;
  else if (istop > jp.nbr_list_size()) {
    istop = jp.nbr_list_size();
  }

  int nic = 0;

  for (int i = 0; i < istop; ++i) {
    if (x != n[_nbr[i]]) {
      continue;
    }

    nic++;

    if (nic >= jp.nbrs_in_common()) {
      return true;
    }
  }

  return false;
}

int
Leader_Item::unselected_neighbours(const int* sel) const {
  int rc = 0;

  for (int i = 0; i < _nbrs; ++i) {
    if (!sel[_nbr[i]]) {
      rc++;
    }
  }

  return rc;
}

int
Leader_Item::unselected_neighbours(const int* sel, const similarity_type_t within) const {
  if (0 == _nbrs) {
    return 0;
  }

  if (_dist[_nbrs - 1] <=
      within) {  // all nbrs within distance, no need to check individually
    return unselected_neighbours(sel);
  }

  int rc = 0;

  for (int i = 0; i < _nbrs; ++i) {
    if (sel[_nbr[i]]) {
      continue;
    }

    if (_dist[i] > within) {
      break;
    }

    rc++;
  }

  return rc;
}

/*
  Items c1 and c2 belong in the same cluster. But they had already been assigned cluster
  id's. Merge those two clusters
*/

static void
merge_clusters(int* cluster, const int n, const int c1, const int c2) {
  int cfrom = cluster[c2];
  int cto = cluster[c1];

  for (int i = 0; i < n; ++i) {
    if (cfrom == cluster[i]) {
      cluster[i] = cto;
    }
  }

  return;
}

static int
Jarvis_Patrick(Leader_Item* pool, const int pool_size, const JarvisPatrickParameters& jp,
               int* cluster) {
  if (verbose) {
    cerr << "Doing JP clustering, nbrs ";
    if (0 == jp.nbr_list_size()) {
      cerr << '.';
    } else {
      cerr << jp.nbr_list_size();
    }
    cerr << " must have " << jp.nbrs_in_common() << " in common to merge\n";
  }

  int* tmp = new_int(pool_size, -1);
  std::unique_ptr<int[]> free_tmp(tmp);

  for (int i = 0; i < pool_size; ++i) {
    pool[i].mark_neighbours(tmp, jp.nbr_list_size(), i);

    for (int j = i + 1; j < pool_size; ++j) {
      if (pool[j].neighbours_in_common(tmp, jp, i)) {
        //      cerr << " i = " << i << " j = " << j << " common " <<
        //      pool[j].neighbours_in_common(tmp, jp, i) << " merge \n";
        merge_clusters(cluster, pool_size, i, j);
      }
    }
  }

  return 1;
}

static int
Jarvis_Patrick(Leader_Item* pool, const int pool_size, const JarvisPatrickParameters& jp,
               const resizable_array_p<IWString>& smiles,
               const resizable_array_p<IWString>& id,
               IWString_and_File_Descriptor& output) {
  int* cluster = new int[pool_size];
  std::unique_ptr<int[]> free_cluster(cluster);

  for (int i = 0; i < pool_size; ++i) {
    cluster[i] = i + 1;
  }

  if (!Jarvis_Patrick(pool, pool_size, jp, cluster)) {
    return 0;
  }

  extending_resizable_array<int> cluster_size;

  int cluster_number = 0;

  for (int i = 0; i < pool_size; ++i) {
    if (0 == cluster[i]) {
      continue;
    }

    int csize = std::count(cluster + i, cluster + pool_size, cluster[i]);

    write_cluster_header(i, smiles, id, cluster_number, csize, output);

    cluster_number++;

    for (int j = i + 1; j < pool_size; ++j) {
      if (cluster[i] != cluster[j]) {
        continue;
      }

      cluster[j] = 0;

      write_smiles_id_dist(smiles[j], id[j], -1.0, output);
    }

    if (verbose) {
      cluster_size[csize]++;
    }

    finish_cluster(csize, output);
  }

  if (verbose) {
    cerr << "Clustered " << smiles.size() << " items into " << cluster_number
         << " clusters\n";
    int in_clusters = 0;
    for (int i = 0; i < cluster_size.number_elements(); ++i) {
      if (cluster_size[i] > 0) {
        cerr << cluster_size[i] << " clusters of size " << i << '\n';
        in_clusters += (cluster_size[i] * i);
      }
    }
    cerr << in_clusters << " items in clusters\n";
  }

  return output.flush();
}

static int
report_clustering(const int* selected, const resizable_array<int>& cluster_size,
                  int pool_size, int clusters_found) {
  cerr << "Clustered " << pool_size << " fingerprints into " << clusters_found
       << " clusters\n";

  int isum = 0;
  for (int i = 0; i < cluster_size.number_elements(); i++) {
    int j = cluster_size[i];
    if (0 == j) {
      continue;
    }

    cerr << j << " clusters of size " << i << " members\n";

    isum += j * i;
  }

  if (isum != pool_size) {
    cerr << "In clusters " << isum << ", started with " << pool_size << " fingerprints\n";
  }

  return cerr.good();
}

static int
identify_next_tb_leader(const Leader_Item* pool, const int pool_size, const int* selected,
                        const similarity_type_t global_threshold, int& leader) {
  leader = -1;
  int nmax = -1;

  for (int i = 0; i < pool_size; ++i) {
    if (selected[i]) {
      continue;
    }

    const int n = pool[i].unselected_neighbours(selected, global_threshold);

    // Scan in input order and replace only for a strictly larger count.
    // Earlier molecules have higher priority when neighbour counts tie.
    if (n > nmax) {
      nmax = n;
      leader = i;
    }
  }

  return leader >= 0;
}

static int
choose_next_centre(const int* selected, const Leader_Item* pool, const int pool_size,
                   const similarity_type_t global_threshold,
                   int& istart, int& icentre) {
  if (taylor_butina) {
    return identify_next_tb_leader(pool, pool_size, selected, global_threshold, icentre);
  }

  icentre = -1;

  // Input order defines priority. Resume after the previous leader since
  // earlier items have already been assigned to clusters.
  for (int i = istart; i < pool_size; ++i) {
    if (!selected[i]) {
      icentre = i;
      istart = i + 1;
      return 1;
    }
  }

  return 0;
}

int
leader(Leader_Item* pool, const int pool_size, int* selected,
       const resizable_array_p<IWString>& smiles, const resizable_array_p<IWString>& id,
       similarity_type_t global_threshold, IWString_and_File_Descriptor& output) {
  if (verbose > 1) {
    cerr << "Beginning clustering at " << global_threshold << '\n';
  }

  int clusters_found = 0;
  std::fill_n(selected, pool_size, 0);

  extending_resizable_array<int> cluster_size;

  // Each threshold starts a fresh scan, while selections within a run advance it.
  int istart = 0;
  int icentre;
  if (!choose_next_centre(selected, pool, pool_size, global_threshold,
                          istart, icentre)) {
    cerr << "Yipes, cannot find initial leader\n";
    return 0;
  }

  resizable_array<int> cluster;

  while (1) {
    Leader_Item& ldr = pool[icentre];

    int max_cluster_size_this_molecule;

    if (ldr.max_cluster_size(max_cluster_size_this_molecule)) {
      ;
    } else {
      max_cluster_size_this_molecule = max_cluster_size;
    }

    if (verbose > 1) {
      cerr << "Start cluster " << clusters_found << ". Threshold = " << global_threshold
           << ", max size " << max_cluster_size_this_molecule << '\n';
    }

    int csize;

    if (max_cluster_size_this_molecule) {
      csize = pool[icentre].form_cluster_max_cluster_size(
          selected, global_threshold, max_cluster_size_this_molecule, clusters_found,
          smiles, id, output);
    } else {
      csize = pool[icentre].form_cluster_threshold(selected, global_threshold, clusters_found,
                                                   smiles, id, output);
    }

    clusters_found++;

    if (verbose) {
      cluster_size[csize]++;
    }

    if (clusters_found >= max_clusters_to_find) {
      break;
    }

    if (!choose_next_centre(selected, pool, pool_size, global_threshold,
                            istart, icentre)) {
      break;
    }
  }

  if (verbose) {
    report_clustering(selected, cluster_size, pool_size, clusters_found);
  }

  return output.flush();
}

static int
leader(Leader_Item* pool, const int pool_size, int* selected,
       const resizable_array_p<IWString>& smiles, const resizable_array_p<IWString>& id,
       similarity_type_t global_threshold, const const_IWSubstring& stem, int ndx) {
  IWString fname(stem);
  fname << ndx;
  fname << ".ldr";

  IWString_and_File_Descriptor output;
  if (!output.open(fname)) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  if (verbose) {
    cerr << "Begin clustering at " << global_threshold << '\n';
  }

  return leader(pool, pool_size, selected, smiles, id, global_threshold, output);
}

static int
build_pool(iwstring_data_source& input, Leader_Item* pool, const int istop, int& ndx,
           const IW_STL_Hash_Map_int& id_to_ndx, resizable_array<int>& ncount,
           bool ignore_missing_neighbours) {
  for (; ndx < istop; ++ndx) {
    if (!pool[ndx].build(input, id_to_ndx, ncount[ndx], ignore_missing_neighbours)) {
      cerr << "Cannot read pool item " << ndx << ", near line " << input.lines_read()
           << '\n';
      return 0;
    }
  }

  return ndx;
}

static int
BuildPool(iw_tf_data_record::TFDataReader& input, Leader_Item* pool,
           const int istop, int& ndx,
           const IW_STL_Hash_Map_int& id_to_ndx, resizable_array<int>& ncount,
           bool ignore_missing_neighbours) {
  for (; ndx < istop; ++ndx) {
    if (!pool[ndx].Build(input, id_to_ndx, ncount[ndx], ignore_missing_neighbours)) {
      cerr << "Cannot read pool item " << ndx << ", near line " << input.items_read()
           << '\n';
      return 0;
    }
  }

  return ndx;
}

static int
build_pool(const char* fname, Leader_Item* pool, const int istop, int& ndx,
           const IW_STL_Hash_Map_int& id_to_ndx, resizable_array<int>& ncount,
           bool ignore_missing_neighbours) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "build_pool:cannot open '" << fname << "'\n";
    return 0;
  }

  return build_pool(input, pool, istop, ndx, id_to_ndx, ncount, ignore_missing_neighbours);
}

static int
BuildPool(const char* fname, Leader_Item* pool, const int istop, int& ndx,
           const IW_STL_Hash_Map_int& id_to_ndx, resizable_array<int>& ncount,
           bool ignore_missing_neighbours) {
  iw_tf_data_record::TFDataReader input(fname);

  if (!input.good()) {
    cerr << "build_pool:cannot open '" << fname << "'\n";
    return 0;
  }

  return BuildPool(input, pool, istop, ndx, id_to_ndx, ncount, ignore_missing_neighbours);
}

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
  cerr << R"(Performs leader clustering, OR Jarvis Patrick clustering on a near neighbour file,
from gfp_nearneighbours_single_file for example.

Options applicable to Leader have an L prepended, J for Jarvis Patrick functionality,
and T for Taylor Butina.
Input order defines molecule priority (highest priority first). Leader clustering
uses this order. Taylor Butina selects the molecule with
the most unselected neighbours, breaking ties by earliest input position.

Usage:
L -t <dis>         distance threshold(s) in [0, 1], inclusive; CSV '-t 0.15,0.20,0.25'
L -C <number>      maximum number of clusters to find
L -m <number>      maximum cluster size
L -M col=<col>     max cluster size for each molecule in column <col> of name

T -T                 do Taylor Butina clustering

J -J <n,c>         Jarvis Patrick: size of nbr list (. for no limit), neighbours in common

 -r               input file is TFDataRecord serialized protos (nn2proto -S ...)
 -B <stem>        file name stem for multiple outputs
 -s <number>      max items across all files; required for multiple input files
 -a               write smiles as output rather than tdt form
 -h               suppress self neighbours during Taylor Butina and Leader
 -v               verbose output
)";
  // clang-format on

  exit(rc);
}

static int
do_jarvis_patrick(const JarvisPatrickParameters* jp, const int njp, Leader_Item* pool,
                  const int pool_size, const resizable_array_p<IWString>& smiles,
                  const resizable_array_p<IWString>& id, const IWString& stem) {
  if (1 == njp) {
    if (0 == stem.length()) {
      IWString_and_File_Descriptor output(1);
      return Jarvis_Patrick(pool, pool_size, jp[0], smiles, id, output);
    }

    IWString fname;
    fname << stem << ".jpc";

    IWString_and_File_Descriptor output;
    if (!output.open(fname)) {
      cerr << "Cannot open output file '" << fname << "'\n";
      return 0;
    }

    return Jarvis_Patrick(pool, pool_size, jp[0], smiles, id, output);
  }

  for (int i = 0; i < njp; ++i) {
    IWString fname;

    fname << stem << i << ".jpc";

    IWString_and_File_Descriptor output;
    if (!output.open(fname)) {
      cerr << "Cannot open output file '" << fname << "'\n";
      return 0;
    }

    if (!Jarvis_Patrick(pool, pool_size, jp[i], smiles, id, output)) {
      return 0;
    }
  }

  return 1;
}

static int
do_leader_clustering(Leader_Item* pool, const int pool_size,
                     const resizable_array_p<IWString>& smiles,
                     const resizable_array_p<IWString>& id, const IWString& stem,
                     const float* threshold, const int number_thresholds) {
  int* selected = new int[pool_size];
  std::unique_ptr<int[]> free_selected(selected);

  IWString_and_File_Descriptor output(1);

  if (number_thresholds == 1) {
    return leader(pool, pool_size, selected, smiles, id, threshold[0], output);
  }

  for (int i = 0; i < number_thresholds; ++i) {
    if (!leader(pool, pool_size, selected, smiles, id, threshold[i], stem, i)) {
      cerr << "Clustering at " << threshold[i] << " failed\n";
      return 0;
    }
  }

  return 1;
}

static int
ReportClusterSizes(const Leader_Item* pool, int pool_size, std::ostream& output) {
  extending_resizable_array<uint32_t> nbrs;
  for (int i = 0; i < pool_size; ++i) {
    ++nbrs[pool[i].neighbour_count()];
  }

  for (int i = 0; i < nbrs.number_elements(); ++i) {
    if (nbrs[i] > 0) {
      output << nbrs[i] << " items had " << i << " neighbours\n";
    }
  }

  return 1;
}


static int
leader(int argc, char** argv) {
  Command_Line cl(argc, argv, "vs:I:t:C:m:M:B:J:aThr");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (!fraction_as_string.initialise(0.0f, 1.0f, 4)) {
    cerr << "Cannot initialise distance formatting\n";
    return 1;
  }

  if (cl.option_present('a')) {
    write_tdt_output = 0;

    if (verbose) {
      cerr << "Smiles output\n";
    }
  }

  int pool_size = 0;
  Leader_Item* pool = nullptr;

  if (cl.option_present('s')) {
    if (!cl.value('s', pool_size) || pool_size < 1) {
      cerr << "The -s option must be followed by a whole positive number\n";
      usage(3);
    }

    if (verbose) {
      cerr << "system sized to " << pool_size << '\n';
    }
  }

  if (cl.option_present('C')) {
    if (!cl.value('C', max_clusters_to_find) || max_clusters_to_find < 1) {
      cerr << "The -C option (max clusters to find) must be followed by a positive "
              "integer\n";
      usage(41);
    }

    if (verbose) {
      cerr << "Will find a max of " << max_clusters_to_find << " clusters\n";
    }
  }

  if (cl.option_present('I')) {
    (void)cl.value('I', identifier_tag);

    if (verbose) {
      cerr << "Identifiers in dataitem '" << identifier_tag << "'\n";
    }
  }

  if (cl.option_present('J'))
    ;
  else if (!cl.option_present('t')) {
    cerr << "Threshold distance must be specified via -t\n";
    usage(28);
  }

  number_thresholds = cl.option_count('t');

  if (1 == number_thresholds)  // maybe it is of the form 0.1,0.2,
  {
    const_IWSubstring t = cl.string_value('t');

    number_thresholds = t.nwords(',');
    if (number_thresholds == 0) {
      cerr << "The -t option requires a distance in [0, 1]\n";
      return 1;
    }
    threshold = new similarity_type_t[number_thresholds];

    int i = 0;
    const_IWSubstring token;
    for (int ndx = 0; t.nextword(token, i, ','); ndx++) {
      if (!token.numeric_value(threshold[ndx]) || threshold[ndx] < 0.0f ||
          threshold[ndx] > 1.0f) {
        cerr << "Invalid threshold '" << token << "' in '" << t
             << "': expected a distance in [0, 1]\n";
        return 1;
      }

      if (verbose) {
        cerr << "Threshold " << ndx << " set to " << threshold[ndx] << '\n';
      }
    }
  } else if (number_thresholds > 1) {
    threshold = new similarity_type_t[number_thresholds];
    for (int i = 0; i < number_thresholds; i++) {
      if (!cl.value('t', threshold[i], i) || threshold[i] < 0.0 || threshold[i] > 1.0) {
        cerr << "Invalid threshold: expected a distance in [0, 1]\n";
        return 6;
      }

      if (verbose) {
        cerr << "Threshold " << i << " set to " << threshold[i] << '\n';
      }
    }
  }

  if (cl.option_present('m')) {
    if (!cl.value('m', max_cluster_size) || max_cluster_size < 2) {
      cerr << "The -m (max cluster size) option must be followed by a whole number > 1\n";
      usage(43);
    }

    if (verbose) {
      cerr << "Max cluster size " << max_cluster_size << '\n';
    }
  }

  if (cl.option_present('M')) {
    const_IWSubstring m = cl.string_value('M');
    if (!m.starts_with("col=")) {
      cerr << "The -M option only supports col=<col> (a column in the molecule name)\n";
      return 1;
    }
    m.remove_leading_chars(4);
    if (!m.numeric_value(max_cluster_size_column) || max_cluster_size_column < 1) {
      cerr << "The column for the per molecule maximum cluster size must be a whole "
              "positive number\n";
      usage(11);
    }

    if (verbose) {
      cerr << "The maximum cluster size per molecule will be in column "
           << max_cluster_size_column << '\n';
    }
    --max_cluster_size_column;
  }

  JarvisPatrickParameters* jp = nullptr;
  int njp = 0;

  if (cl.option_present('J')) {
    njp = cl.option_count('J');
    jp = new JarvisPatrickParameters[njp];

    const_IWSubstring j;
    for (int i = 0; cl.value('J', j, i); ++i) {
      if (!jp[i].build(j)) {
        cerr << "Cannot initialise Jarvis Patrick parameters '" << j << "'\n";
        return 2;
      }
    }

    if (verbose) {
      cerr << "Recognised " << njp << " Jarvis Patrick clustering specifications\n";
    }
  }

  if (cl.option_present('T')) {
    if (nullptr != jp) {
      cerr << "Cannot do Jarvis Patrick and Taylor Butina\n";
      return 1;
    }

    taylor_butina = 1;

    if (verbose) {
      cerr << "Will do Taylor Butina leader clustering\n";
    }
  }

  if (cl.option_present('h') && njp == 0) {
    suppress_self_neighbours = 1;

    if (verbose) {
      cerr << "Will suppress self neighbours\n";
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    usage(1);
  }

  if (cl.number_elements() > 1 && !cl.option_present('s')) {
    cerr << "Multiple input files require an explicit pool size via -s\n";
    return 1;
  }

  int input_is_serialized_protos = 0;
  if (cl.option_present('r')) {
    input_is_serialized_protos = 1;
  }

  int max_pool_size;
  if (pool_size > 0) {
    max_pool_size = pool_size;
  } else {
    max_pool_size = std::numeric_limits<int>::max();
  }

  resizable_array_p<IWString> smiles, id;
  resizable_array<int> ncount;
  resizable_array<int> istop;

  for (int i = 0; i < cl.number_elements() && smiles.number_elements() < max_pool_size; ++i) {
    int rc;
    if (input_is_serialized_protos) {
      rc = ReadSerializedProtos(cl[i], max_pool_size, smiles, id, ncount);
    } else {
      rc = read_smiles_and_ids_and_nbr_count(cl[i], max_pool_size, smiles, id, ncount);
    }

    if (! rc) {
      cerr << "Cannot extract identifiers from " << cl[i] << "'\n";
      return i + 1;
    }

    if (verbose > 1) {
      cerr << "After reading " << cl[i] << " have " << smiles.size() << " molecules\n";
    }

    istop.add(smiles.number_elements());
  }

  pool_size = smiles.number_elements();
  // Once the explicit limit is reached, unseen IDs can belong to the unread
  // portion of the input. Cluster the retained subset and omit those neighbours.
  const bool ignore_missing_neighbours =
      cl.option_present('s') && pool_size == max_pool_size;

  IW_STL_Hash_Map_int id_to_ndx;

  for (int i = 0; i < pool_size; ++i) {
    IWString s = *(id[i]);

    s.truncate_at_first(' ');

    if (id_to_ndx.contains(s)) {
      cerr << "Duplicate identifier '" << s << "' index " << i << ", cannot process\n";
      return 1;
    }

    id_to_ndx[s] = i;
  }

  pool = new Leader_Item[pool_size];
  std::unique_ptr<Leader_Item[]> free_pool(pool);
  for (int i = 0; i < pool_size; ++i) {
    pool[i].set_ndx(i);
  }

  int ndx = 0;

  for (int i = 0; i < istop.number_elements(); ++i) {
    int rc;
    if (input_is_serialized_protos) {
      rc = BuildPool(cl[i], pool, istop[i], ndx, id_to_ndx, ncount, ignore_missing_neighbours);
    } else {
      rc = build_pool(cl[i], pool, istop[i], ndx, id_to_ndx, ncount, ignore_missing_neighbours);
    }
    if (! rc) {
      cerr << "Cannot read fingerprints from '" << cl[i] << "'\n";
      return i + 1;
    }
  }

  if (verbose) {
    ReportClusterSizes(pool, pool_size, cerr);
  }

  float max_distance = 0.0f;
  Accumulator_Int<int> nbr_count;

  for (int i = 0; i < pool_size; ++i) {
    if (verbose) {
      pool[i].update_max_distance(max_distance);
      nbr_count.extra(pool[i].neighbour_count());
    }
  }

  if (verbose) {
    cerr << "Read " << pool_size << " NN data sets. Max distance " << max_distance
         << ". ";
    if (nbr_count.minval() == nbr_count.maxval()) {
      cerr << nbr_count.minval() << " neighbours";
    } else {
      cerr << "Neighbours btw " << nbr_count.minval() << " and " << nbr_count.maxval()
           << " ave " << static_cast<float>(nbr_count.average());
    }
    cerr << '\n';
  }

  if (number_thresholds > 1 && !cl.option_present('B')) {
    cerr << "If doing multiple thresholds, must specify file name stem via the -B "
            "option\n";
    usage(12);
  }

  if (njp > 1 && !cl.option_present('B')) {
    cerr << "If doing multiple JP clusters, must specify file name stem via the -B "
            "option\n";
    usage(12);
  }

  IWString stem;
  if (cl.option_present('B')) {
    cl.value('B', stem);

    if (verbose) {
      cerr << "Output file(s) created with stem '" << stem << "'\n";
    }
  }

  int rc = 0;

  if (njp) {
    if (!do_jarvis_patrick(jp, njp, pool, pool_size, smiles, id, stem)) {
      rc = 1;
    }

    delete[] jp;
  } else {
    if (!do_leader_clustering(pool, pool_size, smiles, id, stem, threshold,
                              number_thresholds)) {
      rc = 1;
    }

    if (nullptr != threshold) {
      delete[] threshold;
    }
  }

  return rc;
}

int
main(int argc, char** argv) {
  int rc = leader(argc, argv);

  return rc;
}
