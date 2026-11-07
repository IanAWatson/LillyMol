// Compare two collections via fingerprint distances.
// Stop when the distribution of distances converges.

#include <stdlib.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>

#define USE_OMP
#ifdef USE_OMP
#include "omp.h"
#endif

#define REPORT_PROGRESS_IMPLEMENTATION
#include "re2/re2.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iw_tdt/iw_tdt.h"
#include "Utilities/GFP_Tools/gfp_standard.h"

namespace gfp_compare_collections {

using std::cerr;

void 
Usage(int rc) {
// clang-format off
#if defined(GIT_HASH) && defined(TODAY)
  cerr << __FILE__ << " compiled " << TODAY << " git hash " << GIT_HASH << '\n';
#else
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
#endif
  // clang-format on
  // clang-format off
  cerr << R"(Compares a set of fingerprints against another, sampling the profile of distances until converged.
Requires standard MPR/FPIW/FPMK/FPMK2 fingerprints.
Sampling counts are exact, regardless of batch size.
A final distribution is written at EOF even without convergence.
Empty inputs or no distances retained by -T are errors.

 -p <fname>     fingerprint file of the comparison set - likely Chembl or some other reference set.
 -s <n>         use at most the first <n> reference fingerprints; use all if omitted.
                Randomise the reference order or use a deliberately diverse prefix.
                Larger values than the file size are allowed but reserve extra memory.
 -T <dist>      upper distance threshold - ignore distances larger than <dist>.
 -n <number>    check convergence and write an output file every <number> fingerprints (default 1000).
 -f <number>    first check after this many queries (minimum 10, default 1000).
                Convergence is assessed either by absolute (-a) or relative (-r) tolerance.
                A comparison is made between the previous distribution and the current values.
 -r <tol>       Relative tolerance. Converged if every bucketised part of the distribution is within <tol>.
 -a <tol>       Absolute tolerance. Converged if every fractional difference is less than <tol>.
 -S <stem>      output files are written as <stem>0.csv, <stem>1.csv...
 -b <size>      batch size. Query fingerprints are processed in groups - default 10.
 -h <threads>   number of OMP threads to use.
 -v             verbose output
)";

  exit(rc);
}

constexpr int kNBins = 1001;

class Options {
  private:
    int _verbose;
    uint64_t _fingerprints_read;

    float _max_distance;

    uint64_t _next_sample;
    uint64_t _sample_interval;

    double _absolute_tolerance;
    double _relative_tolerance;

    GFP_Standard* _fp;
    uint32_t _number_fingerprints;
    // Established once from the first reference fingerprint, then reused.
    int _iw_index = -1;
    int _mk_index = -1;
    int _mk2_index = -1;

    // From ChatGPT we find that efficiency improves if we can amortize the
    // overhead of thread creation across multiple fingerprints.
    // The initial version processed fingerprints one at a time and parallel
    // processing happened with each fingerprint.
    uint32_t _batch_size = 10;

    IWString _stem;
    IWString _reference_fname;
    IWString _query_fname;
    int _next_fname_index;
    uint64_t _last_written_fingerprints = 0;

    uint64_t _acc[kNBins];
    double _distribution[kNBins];
    double _previous_distribution[kNBins];

  // Private functions
    int InitialiseFingerprintLayout(const IW_General_Fingerprint& gfp);
    int ReadFingerprint(IW_TDT& tdt, GFP_Standard& destination,
                        bool initialise_layout = false);
    int ReadComparisonSet(IWString& fname);
    int ReadComparisonSet(iwstring_data_source& input, bool limited);
    int Process(iwstring_data_source& input, int& converged);
    int WriteDistribution();

    int UpdateDistribution();
    int Sample(int& converged);
    int AtolConverged() const;
    int RtolConverged() const;
    int ProcessBatch(const GFP_Standard* queries, uint32_t n, int& converged);

  public:
    Options();
    ~Options();

    int Initialise(Command_Line& cl);

    int Process(const char* fname, int& converged);

    int Report(std::ostream& output) const;
};

Options::Options() {
  _verbose = 0;
  _fingerprints_read = 0;
  _next_sample = 1000;
  _sample_interval = 1000;
  _max_distance = 1.0f;
  _absolute_tolerance = 0.0;
  _relative_tolerance = 0.0;
  _fp = nullptr;
  _number_fingerprints = 0;
  _next_fname_index = 0;

  std::fill_n(_acc, kNBins, 0);
  std::fill_n(_distribution, kNBins, 0.0);
  std::fill_n(_previous_distribution, kNBins, 0.0);
}

Options::~Options() {
  if (_fp != nullptr) {
    delete [] _fp;
  }
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('s')) {
    if (! cl.value('s', _number_fingerprints) || _number_fingerprints == 0) {
      cerr << "The number of fingerprints in the -p file (-s option) must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Fingerprint array sized to " << _number_fingerprints << " fingerprints\n";
    }
  }

  if (cl.option_present('p')) {
    IWString fname = cl.string_value('p');
    if (! ReadComparisonSet(fname)) {
      cerr << "Cannot read comparison set '" << fname << "'\n";
      return 0;
    }
  }

  if (cl.option_present('S')) {
    cl.value('S', _stem);
    if (_verbose) {
      cerr << "output files written with stem '" << _stem << "'\n";
    }
  }

  if (cl.option_present('f')) {
    if (! cl.value('f', _next_sample) || _next_sample < 10) {
      cerr << "The first check option (-f) must be a reasonably large +ve integer\n";
      return 0;
    }

    if (_verbose) {
      cerr << "First check will be after " << _next_sample << " fingerprints\n";
    }
  }

  if (cl.option_present('n'))  {
    if (! cl.value('n', _sample_interval) || _sample_interval == 0) {
      cerr << "The sample interval (-n) must be a whole +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will check for convergence every " << _sample_interval << " fingerprints\n";
    }
  }

  if (cl.option_present('a') && cl.option_present('r')) {
    cerr << "Only one of absolute (-a) or relative (-r) tolerance can be specified\n";
    return 0;
  }

  if (cl.option_present('a')) {
    if (! cl.value('a', _absolute_tolerance) || _absolute_tolerance <= 0.0) {
      cerr << "The convergence tolerance (-a) option must be a +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will be converged if distribution remains within " << _absolute_tolerance << '\n';
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _relative_tolerance) || _relative_tolerance <= 0.0) {
      cerr << "The relative tolerance convergence criterion (-r) must be a positive value\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Relative tolerance " << _relative_tolerance << '\n';
    }
  }

  if (_absolute_tolerance == 0.0 && _relative_tolerance == 0.0) {
    _relative_tolerance = 0.01;
  }

  if (cl.option_present('T')) {
    if (! cl.value('T', _max_distance) || _max_distance <= 0.0f || _max_distance > 1.0f) {
      cerr << "The maximum distance value (-T) must be a valid distance\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will only record distances <= " << _max_distance << '\n';
    }
  }

  if (cl.option_present('b')) {
    if (! cl.value('b', _batch_size) || _batch_size < 1) {
      cerr << "The batch size parameter must be a whole +ve number\n";
      return 0;
    }
    if (_verbose) {
      cerr << "Query fingerprints processed in batches of size " << _batch_size << '\n';
    }
  }

  return 1;
}

int
Options::ReadComparisonSet(IWString& fname) {
  _reference_fname = fname;
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadComparisonSet:cannot open '" << fname << "'\n";
    return 0;
  }

  const bool limited = _number_fingerprints != 0;
  if (!limited) {
    RE2 vbar("^\\|$");
    const uint64_t count = input.grep(vbar);
    if (count > std::numeric_limits<uint32_t>::max()) {
      cerr << "Options::ReadComparisonSet:too many fingerprints\n";
      return 0;
    }
    _number_fingerprints = count;
    if (_number_fingerprints == 0) {
      cerr << "Options::ReadComparisonSet:no fingerprints\n";
      return 0;
    }
  }

  _fp = new GFP_Standard[_number_fingerprints];

  return ReadComparisonSet(input, limited);
}

// Validate the complete layout once. As in other standard GFP tools, all
// subsequent reference and query fingerprints are assumed to use this layout.
int
Options::InitialiseFingerprintLayout(const IW_General_Fingerprint& gfp) {
  if (!standard_fingerprints_present()) {
    return 0;
  }
  for (int i = 0; i < number_fingerprints(); ++i) {
    const IWString& tag = fixed_fingerprint_tag(i);
    if (tag.starts_with("FPIW")) {
      _iw_index = i;
    } else if (tag.starts_with("FPMK2")) {
      _mk2_index = i;
    } else if (tag.starts_with("FPMK")) {
      _mk_index = i;
    }
  }
  if (gfp.molecular_properties_integer().nproperties() != 8 ||
      gfp[_iw_index].nbits() != 2048 || gfp[_mk_index].nbits() == 0 ||
      gfp[_mk_index].nbits() > 256 || gfp[_mk_index].nbits() % 32 != 0 ||
      gfp[_mk2_index].nbits() != gfp[_mk_index].nbits()) {
    cerr << "Invalid standard fingerprint layout: expected 8 MPR properties, "
            "2048 FPIW bits and equal FPMK/FPMK2 sizes (multiples of 32, at most 256)\n";
    return 0;
  }
  return 1;
}

int
Options::ReadFingerprint(IW_TDT& tdt, GFP_Standard& destination, bool initialise_layout) {
  IW_General_Fingerprint gfp;
  int fatal = 0;
  if (!gfp.construct_from_tdt(tdt, fatal)) {
    cerr << "Options::ReadFingerprint:cannot read fingerprint\n";
    return 0;
  }
  if (initialise_layout && !InitialiseFingerprintLayout(gfp)) {
    return 0;
  }
  destination.build_molecular_properties(gfp.molecular_properties_integer());
  destination.build_iw(gfp[_iw_index]);
  destination.build_mk(gfp[_mk_index]);
  destination.build_mk2(gfp[_mk2_index]);
  return 1;
}

int
Options::ReadComparisonSet(iwstring_data_source& input, bool limited) {
  IW_TDT tdt;
  uint32_t ndx = 0;
  // An explicit -s limits reading, not just storage. Leave the rest unread.
  while (!limited || ndx < _number_fingerprints) {
    if (!tdt.next(input)) {
      if (!tdt.rawdata().empty()) {
        cerr << "Malformed reference TDT near line " << input.lines_read() << '\n';
        return 0;
      }
      break;
    }
    if (ndx == _number_fingerprints) {
      cerr << "Reference fingerprint count exceeds allocated capacity "
           << _number_fingerprints << '\n';
      return 0;
    }
    if (!ReadFingerprint(tdt, _fp[ndx], ndx == 0)) {
      cerr << "Options::ReadComparisonSet:error reading fingerprint, line "
           << input.lines_read() << '\n';
      return 0;
    }
    ++ndx;
  }

  if (ndx == 0) {
    cerr << "No reference fingerprints read\n";
    return 0;
  }
  _number_fingerprints = ndx;

  if (_verbose) {
    cerr << "Read " << _number_fingerprints << " fingerprints\n";
  }

  return 1;
}

int
Options::Process(const char* fname, int& converged) {
  _query_fname = fname;
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, converged);
}

int
Options::Process(iwstring_data_source& input, int& converged) {
  // Retain the requested batch size after a small first sampling interval,
  // but do not allocate more than any sampling interval can use.
  const uint32_t capacity = std::min<uint64_t>(
      _batch_size, std::max(_next_sample, _sample_interval));
  auto batch = std::make_unique<GFP_Standard[]>(capacity);

  IW_TDT tdt;
  uint32_t ndx = 0;
  while (tdt.next(input)) {
    if (!ReadFingerprint(tdt, batch[ndx])) {
      cerr << "Options::Process:error reading fingerprint, line "
           << input.lines_read() << '\n';
      return 0;
    }
    ++ndx;

    // Split batches at sampling points so -f and -n are independent of -b.
    if (ndx == capacity || _fingerprints_read + ndx == _next_sample) {
      if (!ProcessBatch(batch.get(), ndx, converged)) {
        return 0;
      }
      ndx = 0;
      if (converged) {
        break;
      }
    }
  }

  if (!converged && !tdt.rawdata().empty()) {
    cerr << "Malformed query TDT near line " << input.lines_read() << '\n';
    return 0;
  }
  if (ndx > 0 && !ProcessBatch(batch.get(), ndx, converged)) {
    return 0;
  }

  if (_fingerprints_read == 0) {
    cerr << "No query fingerprints read\n";
    return 0;
  }
  if (!UpdateDistribution()) {
    cerr << "No distances retained; check the reference set and -T threshold\n";
    return 0;
  }
  // A scheduled sample may already have written this exact set of observations.
  if (_last_written_fingerprints != _fingerprints_read && !WriteDistribution()) {
    return 0;
  }
  if (_verbose) {
    cerr << (converged ? "Converged after " : "End of input without convergence after ")
         << _fingerprints_read << " query fingerprints\n";
  }
  return 1;
}

int
Options::ProcessBatch(const GFP_Standard* queries, uint32_t n,
                      int& converged) {
  if (_verbose > 1) {
    cerr << "ProcessBatch, processed " << _fingerprints_read << " fingerprints\n";
  }
  
  // Cache member variables - ChatGPT suggestion.
  const float max_distance = _max_distance;
  const uint32_t nfp = _number_fingerprints;
  GFP_Standard* fp = _fp;

#pragma omp parallel
  {
    uint64_t local_acc[kNBins] = {0};

#pragma omp for schedule(static)
    for (uint32_t q = 0; q < n; ++q) {
      const GFP_Standard& query = queries[q];

      for (uint32_t i = 0; i < nfp; ++i) {
        const float d = query.tanimoto_distance(fp[i]);
        if (d > max_distance) {
          continue;
        }

        const int ndx = static_cast<int>(d * 1000.0f);
        ++local_acc[ndx];
      }
    }

#pragma omp critical
    {
      for (int i = 0; i < kNBins; ++i) {
        _acc[i] += local_acc[i];
      }
    }
  }

  _fingerprints_read += n;

  if (_fingerprints_read < _next_sample) {
    return 1;
  }

  if (!Sample(converged)) {
    return 0;
  }
  if (_sample_interval > std::numeric_limits<uint64_t>::max() - _fingerprints_read) {
    _next_sample = std::numeric_limits<uint64_t>::max();
  } else {
    _next_sample = _fingerprints_read + _sample_interval;
  }
  return 1;
}

// Normalise accumulated counts. No retained distances is not a distribution.
int
Options::UpdateDistribution() {
  uint64_t tot = 0;
  for (int i = 0; i < kNBins; ++i) {
    tot += _acc[i];
  }
  if (tot == 0) {
    return 0;
  }
  const double multiplier = 1.0 / static_cast<double>(tot);
  for (int i = 0; i < kNBins; ++i) {
    _distribution[i] = static_cast<double>(_acc[i]) * multiplier;
  }
  return 1;
}

// Keep convergence status separate from errors opening or writing a sample.
int
Options::Sample(int& converged) {
  if (_next_fname_index > 0) {
    std::copy_n(_distribution, kNBins, _previous_distribution);
  }
  if (!UpdateDistribution()) {
    // A later batch can still contain distances within -T.
    if (_verbose) {
      cerr << "No distances retained after " << _fingerprints_read << " queries\n";
    }
    return 1;
  }
  if (_next_fname_index > 0) {
    converged = _absolute_tolerance > 0.0 ? AtolConverged() : RtolConverged();
  }
  if (_verbose) {
    cerr << "Sample after " << _fingerprints_read << " query fingerprints\n";
  }
  return WriteDistribution();
}

// Return true if the distribution is converted via absolute diff.
// This is slightly inefficient, since we could return false once we
// have found a difference greater than _absolute_tolerance.
int
Options::AtolConverged() const {
  double max_diff = 0.0;

  for (int i = 0; i < kNBins; ++i) {
    const double d = std::abs(_distribution[i] - _previous_distribution[i]);
    if (d > max_diff) {
      max_diff = d;
    }
  }

  if (max_diff <= _absolute_tolerance) {
    return 1;
  }

  if (_verbose) {
    cerr << _fingerprints_read << " fingerprints read, atol max_diff " << max_diff << '\n';
  }

  return 0;
}

// Returns true of the distribution is converged via the relative
// tolerance criteria.
// Slightly inefficient because once a failure in rtol is detected
// we should return 0 immediately, but we keep going to if _verbose
// is set, we can return the max_rtol.
int
Options::RtolConverged() const {
  int failed = 0;
  double max_rtol = 0.0;

  for (int i = 0; i < kNBins; ++i) {
    if (_previous_distribution[i] == 0.0 && _distribution[i] > 0.0) {
      failed = 1;
      continue;
    }

    const double d = std::abs(_distribution[i] - _previous_distribution[i]);
    if (d == 0.0) {
      continue;
    }

    const double ave = (_distribution[i] + _previous_distribution[i]) * 0.5;
    const double rtol = d / ave;

    if (rtol > max_rtol) {
      max_rtol = rtol;
    }

    if (rtol > _relative_tolerance) {
      failed = 1;
    }
  }

  if (! failed) {
    return 1;
  }

  if (_verbose) {
    cerr << _fingerprints_read
         << " fingerprints read, max relative tolerance " << max_rtol << '\n';
  }

  return 0;
}

int
Options::WriteDistribution() {
  IWString fname;
  fname << _stem << _next_fname_index << ".csv";
  for (IWString* input : {&_reference_fname, &_query_fname}) {
    std::error_code ec;
    if (std::filesystem::equivalent(fname.null_terminated_chars(),
                                    input->null_terminated_chars(), ec)) {
      cerr << "Output file '" << fname << "' would overwrite input '" << *input << "'\n";
      return 0;
    }
  }
  if (_verbose) {
    cerr << "WriteDistribution '" << fname << "'\n";
  }
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Options::WriteDistribution:cannot open '" << fname << "'\n";
    return 0;
  }

  static constexpr char kSep = ',';

  double mean_distance = 0.0;

  output << "Dist" << kSep << "Fraction\n";
  for (int i = 0; i < kNBins; ++i) {
    const float dist = static_cast<float>(i / 1000.0f);
    if (dist > _max_distance) {
      break;
    }
    output << dist << kSep << static_cast<float>(_distribution[i]) << '\n';

    // Work with fractions rather than narrowing potentially huge bin counts
    // to the int count accepted by Accumulator::extra.
    mean_distance += static_cast<double>(dist) * _distribution[i];
  }

  if (_verbose) {
    cerr << "Average " << static_cast<float>(mean_distance) << '\n';
  }

  if (!output.close()) {
    cerr << "Options::WriteDistribution:error writing '" << fname << "'\n";
    return 0;
  }
  ++_next_fname_index;
  _last_written_fingerprints = _fingerprints_read;
  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _fingerprints_read << " fingerprints\n";
  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vs:S:p:a:r:n:f:h:T:b:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  int verbose = cl.option_count('v');

  if (! cl.option_present('S')) {
    cerr << "Must specify output file name stem via the -S option\n";
    Usage(1);
  }

  if (! cl.option_present('p')) {
    cerr << "Must specify comparison set via the -p option\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Must specify input file as command line argument\n";
    Usage(1);
  }

  if (cl.number_elements() > 1) {
    cerr << "Only one input query fingerprint file can be specified\n";
    Usage(1);
  }

#ifdef USE_OMP
  if (cl.option_present('h')) {
    int h;
    if (!cl.value('h', h) || h <= 0) {
      cerr << "The maximum number of threads to use (-h) must be a valid whole +ve "
              "number\n";
      Usage(2);
    }

    omp_set_num_threads(h);
  }
#endif

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  int converged = 0;
  if (! options.Process(cl[0], converged)) {
    cerr << "Error processing '" << cl[0] << "'\n";
    return 1;
  }

  if (verbose) {
    options.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_compare_collections

int
main(int argc, char **argv) {
  int rc = gfp_compare_collections::Main(argc, argv);

  return rc;
}
