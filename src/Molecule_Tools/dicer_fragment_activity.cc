// Associate dicer fragments with activity, without interpreting molecular structures.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "google/protobuf/struct.pb.h"
#include "google/protobuf/text_format.h"
#include "google/protobuf/util/json_util.h"

#include "Foundational/accumulator/accumulator.h"
#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/data_source/tfdatarecord.h"

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/strings/string_view.h"
#ifdef BUILD_BAZEL
#include "Molecule_Tools/dicer_fragments.pb.h"
#else
#include "dicer_fragments.pb.h"
#endif

namespace dicer_fragment_activity {
using google::protobuf::Struct;
using google::protobuf::Value;
using std::cerr;

void
Usage(int rc) {
  cerr << R"(Associate per-molecule dicer output with an activity file.

  dicer -k 3 -I 1 -B serialized_proto -S diced.dat molecules.smi
  dicer_fragment_activity -A activity.txt -tfdata -support 10 diced.dat
  dicer_fragment_activity -A activity.txt -tfdata -support 10 -json report.json diced.dat

 -A <fname>              id/activity file with a header
 -tfdata                 read TFDataRecord serialized DicedMolecule protos
                         default: one DicedMolecule textproto per line (-B proto)
 -support <n>            minimum number of molecules containing a fragment
 -C <natoms>             maximum fragment atom count
 -json <fname>           write a structured activity report instead of the table
                         use '-' for JSON on stdout
 -unmatched_support <n>  minimum absent-group size for JSON (default 1)
 -bins <n>               common histogram bins for JSON (default 20)
 -examples <n>           parent examples retained in JSON (default 3; 0 allowed)
 -v                      report population and fragment counts

The population is the unique molecules actually read, including those with no
fragments. Repeated molecule ids and duplicate activity ids are errors. Each
fragment counts once per molecule; fragment occurrence counts are not weights.

The default table retains columns Frag Parent Natoms N Fraction Minval Ave
Maxval Diff. Diff is matched mean minus the processed population mean. JSON
also reports matched-minus-unmatched mean and median differences.
)";
  ::exit(rc);
}

struct Activity {
  double _value;
  bool _seen = false;
};

struct Parent {
  std::string _id;
  std::string _smiles;
  double _activity;
};

struct Fragment {
  Accumulator<double> _activity;
  std::string _parent;
  uint32_t _natoms = 0;
  // Populated only for JSON. Initially parent indices, later activity ranks.
  std::vector<uint32_t> _members;
};

void
Number(Struct& object, const std::string& name, double value) {
  Value& field = (*object.mutable_fields())[name];
  if (std::isfinite(value)) {
    field.set_number_value(value);
  } else {
    field.set_null_value(google::protobuf::NULL_VALUE);
  }
}

void
String(Struct& object, const std::string& name, const std::string& value) {
  (*object.mutable_fields())[name].set_string_value(value);
}

Struct&
Object(Struct& object, const std::string& name) {
  return *(*object.mutable_fields())[name].mutable_struct_value();
}

template <typename T>
void
Array(Struct& object, const std::string& name, const std::vector<T>& values) {
  auto* list = (*object.mutable_fields())[name].mutable_list_value();
  for (const auto value : values) {
    list->add_values()->set_number_value(value);
  }
}

// Quantiles use linear interpolation at (n - 1) * percentile, matching the
// dictionary-feature analysis. value_at accesses either a group or its implicit
// complement; no vector containing all absent molecules is needed.
template <typename Getter>
Struct
Summary(uint32_t n, long double sum, Getter value_at) {
  Struct result;
  Number(result, "n", n);
  if (n == 0) {
    for (const char* key : {"min", "q1", "median", "mean", "q3", "max"}) {
      (*result.mutable_fields())[key].set_null_value(google::protobuf::NULL_VALUE);
    }
    return result;
  }
  auto quantile = [&](double fraction) {
    const double pos = (n - 1) * fraction;
    const uint32_t lower = std::floor(pos);
    const uint32_t upper = std::ceil(pos);
    const long double weight = pos - lower;
    return static_cast<double>((1 - weight) * value_at(lower) + weight * value_at(upper));
  };
  Number(result, "min", value_at(0));
  Number(result, "q1", quantile(0.25));
  Number(result, "median", quantile(0.5));
  Number(result, "mean", sum / n);
  Number(result, "q3", quantile(0.75));
  Number(result, "max", value_at(n - 1));
  return result;
}

bool
WriteObject(std::ostream& output, const Struct& object) {
  std::string json;
  const auto status = google::protobuf::util::MessageToJsonString(object, &json);
  if (!status.ok()) {
    cerr << "Cannot encode activity report: " << status.ToString() << '\n';
    return false;
  }
  output << json;
  return output.good();
}

class Analysis {
 private:
  bool _verbose = false;
  bool _tfdata = false;
  uint32_t _support = 0;
  int _max_atoms = 0;
  uint32_t _unmatched_support = 1;
  int _bins = 20;
  int _examples = 3;
  std::string _activity_fname;
  std::string _json_fname;
  absl::flat_hash_map<std::string, Activity> _activities;
  absl::flat_hash_map<std::string, Fragment> _fragments;
  Accumulator<double> _population;
  std::vector<Parent> _parents;
  std::vector<uint32_t> _order;
  std::vector<uint32_t> _ranks;
  std::vector<double> _values;
  std::vector<double> _edges;
  std::vector<uint32_t> _histogram;
  long double _sum = 0;

  using Entry = absl::flat_hash_map<std::string, Fragment>::value_type;
  bool ReadActivity();
  bool AddMolecule(const dicer_data::DicedMolecule& proto);
  bool ReadText(const char* fname);
  bool ReadBinary(const char* fname);
  bool OutputOverwritesInput(const Command_Line_v2& cl) const;
  std::vector<Entry*> OutputFragments();
  bool Supported(const Fragment& fragment) const;
  void PrepareReport();
  uint32_t HistogramBin(double value) const;
  double AbsentValue(const std::vector<uint32_t>& members, uint32_t index) const;
  Struct FragmentReport(const std::string& smiles, Fragment& fragment, uint32_t index);
  bool WriteJson();
  bool WriteTable();

 public:
  bool Initialise(Command_Line_v2& cl);
  bool Read(const char* fname);
  bool Write();
};

bool
Analysis::OutputOverwritesInput(const Command_Line_v2& cl) const {
  if (_json_fname.empty() || _json_fname == "-") {
    return false;
  }
  std::error_code error;
  const auto destination = std::filesystem::weakly_canonical(_json_fname, error);
  if (error) {
    cerr << "Cannot resolve JSON output path: " << error.message() << '\n';
    return true;
  }
  auto same_file = [&](const std::string& name) {
    std::error_code ec;
    if (std::filesystem::equivalent(name, _json_fname, ec) && !ec) {
      return true;
    }
    ec.clear();
    return std::filesystem::weakly_canonical(name, ec) == destination && !ec;
  };
  if (same_file(_activity_fname)) {
    return true;
  }
  for (const char* fname : cl) {
    if (same_file(fname)) {
      return true;
    }
  }
  return false;
}

bool
Analysis::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_present('v');
  _tfdata = cl.option_present("tfdata");
  _activity_fname = cl.string_value('A').AsString();
  if (cl.option_present("support") && (!cl.value("support", _support) || _support < 1)) {
    return false;
  }
  if (cl.option_present('C') && (!cl.value('C', _max_atoms) || _max_atoms < 1)) {
    return false;
  }
  if (cl.option_present("unmatched_support") &&
      (!cl.value("unmatched_support", _unmatched_support) || _unmatched_support < 1)) {
    return false;
  }
  if (cl.option_present("bins") && (!cl.value("bins", _bins) || _bins < 1)) {
    return false;
  }
  if (cl.option_present("examples") &&
      (!cl.value("examples", _examples) || _examples < 0)) {
    cerr << "Examples must be a non-negative integer\n";
    return false;
  }
  if (cl.option_present("json")) {
    _json_fname = cl.string_value("json").AsString();
    if (OutputOverwritesInput(cl)) {
      cerr << "JSON output would overwrite an input file\n";
      return false;
    }
  }
  return ReadActivity();
}

bool
Analysis::ReadActivity() {
  iwstring_data_source input(_activity_fname.c_str());
  if (!input.good()) {
    cerr << "Cannot read activity file '" << _activity_fname << "'\n";
    return false;
  }
  const_IWSubstring line;
  if (!input.next_record(line) || line.empty()) {
    cerr << "Activity file needs a header\n";
    return false;
  }
  while (input.next_record(line)) {
    line.strip_leading_blanks();
    if (line.empty() || line.starts_with('#')) {
      continue;
    }
    int i = 0;
    const_IWSubstring id, token;
    double value;
    if (!line.nextword(id, i) || !line.nextword(token, i) ||
        !token.numeric_value(value) || !std::isfinite(value)) {
      cerr << "Invalid id/activity record: " << line << '\n';
      return false;
    }
    if (!_activities.emplace(id.AsString(), Activity{value}).second) {
      cerr << "Duplicate activity identifier '" << id << "'\n";
      return false;
    }
  }
  if (!input.good() || !input.eof()) {
    cerr << "Error reading activity file '" << _activity_fname << "'\n";
    return false;
  }
  if (_activities.empty()) {
    cerr << "No activity values read\n";
    return false;
  }
  return true;
}

bool
Analysis::AddMolecule(const dicer_data::DicedMolecule& proto) {
  auto activity = _activities.find(proto.name());
  if (activity == _activities.end()) {
    const size_t end = proto.name().find_first_of(" \t\r\n");
    activity = _activities.find(proto.name().substr(0, end));
  }
  if (activity == _activities.end()) {
    cerr << "No activity for '" << proto.name() << "'\n";
    return false;
  }
  if (activity->second._seen) {
    cerr << "Duplicate molecule identifier '" << activity->first << "'\n";
    return false;
  }
  activity->second._seen = true;
  const double value = activity->second._value;
  const uint32_t parent = _population.n();
  _population.extra(value);
  if (!_json_fname.empty()) {
    _parents.push_back({activity->first, proto.smiles(), value});
  }
  absl::flat_hash_set<absl::string_view> seen;
  for (const auto& fragment : proto.fragment()) {
    const std::string& key = fragment.usmi().empty() ? fragment.smi() : fragment.usmi();
    if (key.empty() || fragment.nat() == 0) {
      cerr << "Fragment needs SMILES and a positive atom count: "
           << fragment.ShortDebugString() << '\n';
      return false;
    }
    auto [iter, inserted] = _fragments.try_emplace(key);
    Fragment& result = iter->second;
    if (inserted) {
      result._parent = activity->first;
      result._natoms = fragment.nat();
    } else if (result._natoms != fragment.nat()) {
      cerr << "Inconsistent atom counts for fragment '" << key << "'\n";
      return false;
    }
    if (!seen.insert(key).second) {
      continue;
    }
    result._activity.extra(value);
    if (!_json_fname.empty()) {
      result._members.push_back(parent);
    }
  }
  return true;
}

bool
Analysis::ReadText(const char* fname) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open dicer input '" << fname << "'\n";
    return false;
  }
  const_IWSubstring line;
  while (input.next_record(line)) {
    line.strip_leading_blanks();
    if (line.empty() || line.starts_with('#')) {
      continue;
    }
    dicer_data::DicedMolecule proto;
    if (!google::protobuf::TextFormat::ParseFromString(line.AsString(), &proto) ||
        !AddMolecule(proto)) {
      cerr << "Invalid dicer textproto in '" << fname << "'\n";
      return false;
    }
  }
  if (!input.good() || !input.eof()) {
    cerr << "Error reading dicer input '" << fname << "'\n";
    return false;
  }
  return true;
}

bool
Analysis::ReadBinary(const char* fname) {
  iw_tf_data_record::TFDataReader input(fname);
  if (!input.good()) {
    cerr << "Cannot open dicer input '" << fname << "'\n";
    return false;
  }
  while (true) {
    const auto record = input.Next();
    if (!record) {
      // A failed protobuf parse is not EOF. Check framing separately, too.
      if (input.good() && input.eof()) {
        return true;
      }
      cerr << "Invalid TFDataRecord in '" << fname << "'\n";
      return false;
    }
    dicer_data::DicedMolecule proto;
    if (!proto.ParseFromArray(record->data(), record->length()) || !AddMolecule(proto)) {
      cerr << "Invalid serialized DicedMolecule in '" << fname << "'\n";
      return false;
    }
  }
}

bool
Analysis::Read(const char* fname) {
  return _tfdata ? ReadBinary(fname) : ReadText(fname);
}

std::vector<Analysis::Entry*>
Analysis::OutputFragments() {
  std::vector<Entry*> result;
  // Apply support/size limits before sorting: most generated fragments may
  // have low support, so sorting every discovered SMILES would waste time and
  // a large pointer array on entries that will never be written.
  for (auto& entry : _fragments) {
    if (!Supported(entry.second)) {
      continue;
    }
    if (!_json_fname.empty() &&
        _population.n() - entry.second._activity.n() < _unmatched_support) {
      continue;
    }
    result.push_back(&entry);
  }
  std::sort(result.begin(), result.end(),
            [](const Entry* lhs, const Entry* rhs) { return lhs->first < rhs->first; });
  return result;
}

bool
Analysis::Supported(const Fragment& fragment) const {
  return fragment._activity.n() >= _support &&
         (_max_atoms == 0 || fragment._natoms <= static_cast<uint32_t>(_max_atoms));
}

void
Analysis::PrepareReport() {
  _order.resize(_parents.size());
  std::iota(_order.begin(), _order.end(), uint32_t{0});
  std::sort(_order.begin(), _order.end(), [&](uint32_t lhs, uint32_t rhs) {
    if (_parents[lhs]._activity != _parents[rhs]._activity) {
      return _parents[lhs]._activity < _parents[rhs]._activity;
    }
    return _parents[lhs]._id < _parents[rhs]._id;
  });
  _ranks.resize(_parents.size());
  for (uint32_t rank = 0; rank < _order.size(); ++rank) {
    _ranks[_order[rank]] = rank;
    _values.push_back(_parents[_order[rank]]._activity);
    _sum += _values.back();
  }
  if (_values.front() == _values.back()) {
    _edges = {_values.front(), _values.back()};
    _bins = 1;
  } else {
    const long double low = _values.front();
    const long double range = static_cast<long double>(_values.back()) - low;
    for (int i = 0; i < _bins; ++i) {
      _edges.push_back(low + range * i / _bins);
    }
    _edges.push_back(_values.back());
  }
  _histogram.resize(_bins, 0);
  for (double value : _values) {
    ++_histogram[HistogramBin(value)];
  }
}

uint32_t
Analysis::HistogramBin(double value) const {
  if (_values.front() == _values.back()) {
    return 0;
  }
  const long double fraction =
      (static_cast<long double>(value) - _values.front()) /
      (static_cast<long double>(_values.back()) - _values.front());
  return std::min<uint32_t>(_bins - 1, fraction * _bins);
}

// Select the kth value absent from this fragment. At rank r, the number of
// absent values is (r + 1) minus the number of member ranks <= r. Binary search
// this monotone count. This costs O(log N log support), not a scan over all N
// molecules for every fragment, and handles tied activities by their ranks.
double
Analysis::AbsentValue(const std::vector<uint32_t>& members, uint32_t index) const {
  uint32_t low = 0;
  uint32_t high = _values.size() - 1;
  while (low < high) {
    const uint32_t mid = low + (high - low) / 2;
    const auto hits =
        std::upper_bound(members.begin(), members.end(), mid) - members.begin();
    if (mid + 1 - hits > index) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }
  return _values[low];
}

Struct
Analysis::FragmentReport(const std::string& smiles, Fragment& fragment, uint32_t index) {
  Struct result;
  Number(result, "index", index);
  String(result, "name", "Fragment " + std::to_string(index));
  String(result, "source", smiles);
  String(result, "source_type", "dicer");
  Number(result, "natoms", fragment._natoms);
  auto* examples = (*result.mutable_fields())["examples"].mutable_list_value();
  const uint32_t nexamples = std::min<size_t>(_examples, fragment._members.size());
  for (uint32_t i = 0; i < nexamples; ++i) {
    const Parent& parent = _parents[fragment._members[i]];
    Struct& example = *examples->add_values()->mutable_struct_value();
    String(example, "identifier", parent._id);
    String(example, "smiles", parent._smiles);
    Number(example, "activity", parent._activity);
  }
  for (uint32_t& member : fragment._members) {
    member = _ranks[member];
  }
  std::sort(fragment._members.begin(), fragment._members.end());
  std::vector<uint32_t> matched_histogram(_bins, 0);
  long double matched_sum = 0;
  for (uint32_t rank : fragment._members) {
    matched_sum += _values[rank];
    ++matched_histogram[HistogramBin(_values[rank])];
  }
  Struct& matched = Object(result, "matched");
  Object(matched, "summary") =
      Summary(fragment._members.size(), matched_sum,
              [&](uint32_t i) { return _values[fragment._members[i]]; });
  Array(matched, "histogram", matched_histogram);
  std::vector<uint32_t> unmatched_histogram(_histogram);
  for (int i = 0; i < _bins; ++i) {
    unmatched_histogram[i] -= matched_histogram[i];
  }
  Struct& unmatched = Object(result, "unmatched");
  Object(unmatched, "summary") =
      Summary(_values.size() - fragment._members.size(), _sum - matched_sum,
              [&](uint32_t i) { return AbsentValue(fragment._members, i); });
  Array(unmatched, "histogram", unmatched_histogram);
  for (const char* statistic : {"mean", "median"}) {
    const auto& a = matched.fields().at("summary").struct_value().fields().at(statistic);
    const auto& b =
        unmatched.fields().at("summary").struct_value().fields().at(statistic);
    if (a.kind_case() == Value::kNumberValue && b.kind_case() == Value::kNumberValue) {
      Number(result, std::string("delta_") + statistic,
             a.number_value() - b.number_value());
    } else {
      (*result.mutable_fields())[std::string("delta_") + statistic].set_null_value(
          google::protobuf::NULL_VALUE);
    }
  }
  return result;
}

bool
Analysis::WriteJson() {
  PrepareReport();
  std::ofstream file;
  std::ostream* output = &std::cout;
  if (_json_fname != "-") {
    file.open(_json_fname);
    if (!file.good()) {
      cerr << "Cannot open JSON output '" << _json_fname << "'\n";
      return false;
    }
    output = &file;
  }
  Struct metadata;
  String(metadata, "activity", _activity_fname);
  String(metadata, "feature_kind", "fragment");
  String(metadata, "multiplicity_mode", "presence");
  Number(metadata, "molecules_read_with_activity", _population.n());
  Number(metadata, "activities_unused", _activities.size() - _population.n());
  Number(metadata, "queries", _fragments.size());
  Number(metadata, "min_matches", _support);
  Number(metadata, "min_unmatched", _unmatched_support);
  Number(metadata, "histogram_bins", _bins);
  const auto ordered = OutputFragments();
  Number(metadata, "queries_reported", ordered.size());
  Struct global;
  Object(global, "summary") =
      Summary(_values.size(), _sum, [&](uint32_t i) { return _values[i]; });
  Array(global, "histogram", _histogram);
  Struct& bins = Object(global, "bins");
  Number(bins, "min", _values.front());
  Number(bins, "max", _values.back());
  Array(bins, "edges", _edges);
  *output << "{\"metadata\":";
  if (!WriteObject(*output, metadata)) {
    return false;
  }
  *output << ",\"global\":";
  if (!WriteObject(*output, global)) {
    return false;
  }
  *output << ",\"queries\":[";
  uint32_t index = 0;
  for (Entry* entry : ordered) {
    if (index) {
      *output << ',';
    }
    if (!WriteObject(*output, FragmentReport(entry->first, entry->second, ++index))) {
      return false;
    }
  }
  *output << "]}\n";
  output->flush();
  return output->good();
}

bool
Analysis::WriteTable() {
  std::cout << "Frag Parent Natoms N Fraction Minval Ave Maxval Diff\n"
            << std::setprecision(7);
  for (const Entry* entry : OutputFragments()) {
    const Fragment& fragment = entry->second;
    if (!Supported(fragment)) {
      continue;
    }
    std::cout << entry->first << ' ' << fragment._parent << ' ' << fragment._natoms << ' '
              << fragment._activity.n() << ' '
              << static_cast<double>(fragment._activity.n()) / _population.n() << ' '
              << fragment._activity.minval() << ' ' << fragment._activity.average() << ' '
              << fragment._activity.maxval() << ' '
              << fragment._activity.average() - _population.average() << '\n';
  }
  std::cout.flush();
  return std::cout.good();
}

bool
Analysis::Write() {
  if (_population.n() == 0) {
    cerr << "No diced molecules read\n";
    return false;
  }
  if (_verbose) {
    cerr << "Read " << _population.n() << " molecules and " << _fragments.size()
         << " fragments; " << _activities.size() - _population.n()
         << " activity ids unused\n";
  }
  return _json_fname.empty() ? WriteTable() : WriteJson();
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv,
                     "-v-A=sfile-tfdata-support=ipos-C=ipos-json=s-bins=ipos-examples=i-"
                     "unmatched_support=ipos");
  if (cl.unrecognised_options_encountered() || !cl.option_present('A') || cl.empty()) {
    Usage(1);
  }
  Analysis analysis;
  if (!analysis.Initialise(cl)) {
    return 1;
  }
  for (const char* fname : cl) {
    if (!analysis.Read(fname)) {
      return 1;
    }
  }
  return analysis.Write() ? 0 : 1;
}
}  // namespace dicer_fragment_activity

int
main(int argc, char** argv) {
  return dicer_fragment_activity::Main(argc, argv);
}
