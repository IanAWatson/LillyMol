// Compare two collections.

#include <stdlib.h>

#include <algorithm>
#include <iostream>

#define USE_OMP
#ifdef USE_OMP
#include "omp.h"
#endif

#define REPORT_PROGRESS_IMPLEMENTATION
#include "re2/re2.h"

#include "Foundational/accumulator/accumulator.h"
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
 -p <fname>     fingerprint file of the comparison set - likely Chembl or some other reference set.
 -s <n>         number of items in the -p file. This tool will count if not specified.
 -n <number>    check convergence and write an output file every <number> fingerprints (default 1000).
 -f <number>    first time to check for convergence. If not specified it will be the -n value.
                Convgence as assessed either by absoluate (-a) or relative (-r) tolerange.
                A comparison is made between the previous distribution and the current values.
 -r <tol>       Relative tolerance. Converged if every bucketised part of the distribution is within <tol>.
 -a <tol>       Absolute tolerance. Converged if every difference is less than <tol>.
 -S <stem>      output files are written as <stem>_n.csv.
 -h <threads>   number of OMP threads to use.
 -v             verbose output
)";

  exit(rc);
}

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

    IWString _stem;
    int _next_fname_index;

    uint64_t _acc[1001];
    double _distribution[1001];
    double _previous_distribution[1001];

  // Private functions
    int ReadComparisonSet(IWString& fname);
    int ReadComparisonSet(iwstring_data_source& input);
    void AddToBinnedDistances(float d);
    int Process(iwstring_data_source& input, int& converged);
    int Process(GFP_Standard& gfp_standard, int& converged);
    int WriteDistribution(uint64_t tot);

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

  std::fill_n(_acc, 1001, 0);
  std::fill_n(_distribution, 1001, 0.0);
  std::fill_n(_previous_distribution, 1001, 0.0);
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
      cerr << "The convergence tolerance (-t) option must be a +ve number\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Will be converged if distribution remains within " << _absolute_tolerance << '\n';
    }
  }

  if (cl.option_present('r')) {
    if (! cl.value('r', _relative_tolerance) || _relative_tolerance <= 0.0) {
      cerr << "The absoluate tolerance convergence criterion (-r) must be a non negative value\n";
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

  return 1;
}

int
Options::ReadComparisonSet(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadComparisonSet:cannot open '" << fname << "'\n";
    return 0;
  }

  if (_number_fingerprints == 0) {
    RE2 vbar("^|");
    _number_fingerprints = input.grep(vbar);
    if (_number_fingerprints == 0) {
      cerr << "Options::ReadComparisonSet:no fingerprints\n";
      return 0;
    }
  }

  _fp = new GFP_Standard[_number_fingerprints];

  return ReadComparisonSet(input);
}

int
FromGfp(const IW_General_Fingerprint& gfp,
        GFP_Standard& destination) {

  destination.build_molecular_properties(gfp.molecular_properties_integer());
  destination.build_iw(gfp[0]);
  destination.build_mk(gfp[1]);
  destination.build_mk2(gfp[2]);

  return 1;
}

int
FromGfp(IW_TDT& tdt,
        GFP_Standard& destination,
        int& fatal) {
  IW_General_Fingerprint gfp;
  fatal = 0;
  if (!gfp.construct_from_tdt(tdt, fatal)) {
    cerr << "Options::ReadComparisonSet:cannot read fingerprint\n";
    return 0;
  }

  return FromGfp(gfp, destination);
}

int
Options::ReadComparisonSet(iwstring_data_source& input) {
  IW_TDT tdt;
  int fatal;
  uint32_t ndx = 0;
  for (; ndx < _number_fingerprints && tdt.next(input); ++ndx) {
    if (! FromGfp(tdt, _fp[ndx], fatal)) {
      if (fatal) {
        cerr << "Options::ReadComparisonSet:error reading fingerprint, line " << input.lines_read() << '\n';
        return 0;
      }
      break;
    }
  }

  _number_fingerprints = ndx;

  if (_verbose) {
    cerr << "Read " << _number_fingerprints << " fingerprints\n";
  }

  return 1;
}

int
Options::Process(const char* fname, int& converged) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, converged);
}

int
Options::Process(iwstring_data_source& input, int& converged) {
  IW_TDT tdt;
  int fatal = 0;
  while (tdt.next(input)) {
    GFP_Standard gfp_standard;
    if (FromGfp(tdt, gfp_standard, fatal)) {
      //  great
    } else if (fatal) {
      cerr << "Options::Process:error reading fingerprint, line " << input.lines_read() << '\n';
      return 0;
    } else {
      break;
    }

    if (Process(gfp_standard, converged)) {
      // great
    } else if (converged) {
      return 0;
    }
  }

  return 1;
}

int
Options::Process(GFP_Standard& gfp_standard, int& converged) {
#pragma omp parallel shared(_acc)
  {
#pragma omp for schedule(dynamic, 256)
    for (uint32_t i = 0; i < _number_fingerprints; ++i) {
      float d = gfp_standard.tanimoto_distance(_fp[i]);
      if (d > _max_distance) {
        continue;
      }
      AddToBinnedDistances(d);
    }
  }

  ++_fingerprints_read;
  if (_fingerprints_read < _next_sample) {
    return 1;
  }

  if (_next_fname_index > 0) {
    std::copy_n(_distribution, 1001, _previous_distribution);
  }

  uint64_t tot = 0;
  for (int i = 0; i < 1001; ++i) {
    tot += _acc[i];
  }

  double multiplier = 1.0 / static_cast<float>(tot);
  for (int i = 0; i < 1001; ++i) {
    _distribution[i] = static_cast<double>(_acc[i]) * multiplier;
  }

  converged = 0;

  if (_next_fname_index == 0) {
  } else if (_absolute_tolerance > 0.0) {
    double max_diff = 0.0;
    for (int i = 0; i < 1001; ++i) {
      double d = std::abs(_distribution[i] - _previous_distribution[i]);
      if (d > max_diff) {
        max_diff = d;
      }
    }

    if (max_diff <= _absolute_tolerance) {
      converged = 1;
    }

    if (_verbose) {
      cerr << _fingerprints_read << " fingerprints read, max_diff " << max_diff;
      if (converged) {
        cerr << " converged\n";
      } else {
        cerr << '\n';
      }
    }
  } else if (_relative_tolerance > 0.0) {
    int failed = 0;
    double max_rtol = 0.0;
    for (int i = 0; i < 1001; ++i) {
      if (_previous_distribution[i] == 0 && _distribution[i] > 0) {
        failed = 1;
        continue;
      }

      double d = std::abs(_distribution[i] - _previous_distribution[i]);
      if (d == 0.0) {
        continue;
      }
      double ave = (_distribution[i] + _previous_distribution[i]) * 0.5;
      double rtol = (d / ave);
      if (rtol > 1.0) {
        cerr << "rtol " << rtol << ' ' <<  _distribution[i] << ' ' << _previous_distribution[i] << " diff " << d << " ave " << ave << '\n';
      }
      if (rtol > max_rtol) {
        max_rtol = rtol;
      }
      if (rtol > _relative_tolerance) {
        failed = 1;
        break;
      }
    }

    if (! failed) {
      converged = 1;
    }

    if (_verbose) {
      cerr << _fingerprints_read << " fingerprints read, max relative tolerance " << max_rtol;
      if (failed) {
        cerr << " failed\n";
      } else {
        cerr << " converged\n";
      }
    }
  }

  WriteDistribution(tot);

  _next_sample += _sample_interval;

  if (converged) {
    return 0;
  }

  return 1;
}

int
Options::WriteDistribution(uint64_t tot) {
  IWString fname;
  fname << _stem << _next_fname_index << ".csv";
  ++_next_fname_index;
  IWString_and_File_Descriptor output;
  if (! output.open(fname)) {
    cerr << "Options::WriteDistribution:cannot open '" << fname << "'\n";
    return 0;
  }

  static constexpr char kSep = ',';

  Accumulator<double> stats;

  output << "Dist" << kSep << "Fraction\n";
  float next_distance = 0.05;
  uint64_t sum = 0;
  for (int i = 0; i < 1001; ++i) {
    float dist = static_cast<float>(i / 1000.0f);
    output << dist << kSep << static_cast<float>(_distribution[i]) << '\n';

    sum += _acc[i];

    if (_acc[i] > 0) {
      stats.extra(dist, _acc[i]);
    }

    if (dist > _max_distance) {
      break;
    }

    if (dist < next_distance) {
      continue;
    }
    if (dist > 0.50) {  // uninteresting.
      continue;
    }
    if (dist >= next_distance) {
      cerr << dist << ' ' << sum << ' ' << static_cast<float>(sum) / static_cast<float>(tot) << '\n';
      next_distance += 0.05;
    }
  }

  if (_verbose) {
    cerr << "Average " << static_cast<float>(stats.average()) << " std " << std::sqrt(stats.variance()) << '\n';
    cerr << stats.minval() << " to " << stats.maxval() << ' ' << stats.sum() << " sum\n";
  }

  return 1;
}

void
Options::AddToBinnedDistances(float d) {
  int ndx = static_cast<int>(d * 1000.0f + 0.49999f);
  ++_acc[ndx];
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << _fingerprints_read << " fingerprints\n";
  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vs:S:p:a:r:n:f:h:T:");

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

#ifdef USE_OMP
  if (cl.option_present('h')) {
    int h;
    if (!cl.value('h', h) || h < 0) {
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

  if (cl.empty()) {
    cerr << "Must specify input file as command line argument\n";
    Usage(1);
  }

  for (const char* fname : cl) {
    int converged = 0;
    if (options.Process(fname, converged)) {
      continue;
    }
    if (converged) {
      break;
    }
    cerr << "Error processing '" << fname << "'\n";
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
