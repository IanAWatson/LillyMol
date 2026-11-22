#include "tsubstructure_fp.h"

#include <algorithm>
#include <iostream>
#include <limits>

#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"

using std::cerr;

static IWString smiles_tag("$SMI<");
static IWString pcn_tag("PCN<");

TSubstructure_FP::TSubstructure_FP() {
  _work_as_filter = 0;
  _bit_replicates = 1;
  _default_fingerprint_nbits = 1;
  _extra_bit_total_hits = 0;
  _each_query_writes_bit = 1;
}

template <typename OUTPUT>
int
TSubstructure_FP::do_fingerprint_output(Molecule& m, int nq, const int* hits,
                                        OUTPUT& output) const {
  if (_tag.empty()) {
    cerr << "TSubstructure_FP::do_fingerprint_output:no tag\n";
    return 0;
  }

  const uint64_t feature_count =
      static_cast<uint64_t>(std::max(nq, _default_fingerprint_nbits)) * _bit_replicates +
      _extra_bit_total_hits;
  if (feature_count > (std::numeric_limits<int>::max() - 63) / 8) {
    cerr << "TSubstructure_FP: too many fingerprint features\n";
    return 0;
  }
  if (!_work_as_filter) {
    write_smiles_and_pcn(m, output);
  }

  if (_tag.starts_with("FP")) {
    return _do_fingerprint_output(nq, hits, output);
  }
  if (_tag.starts_with("NC")) {
    return _do_sparse_fingerprint_output(nq, hits, output);
  }
  if (_tag.starts_with("FC")) {
    return _do_fixed_counted_fingerprint_output(nq, hits, output);
  }

  cerr << "TSubstructure_FP::do_fingerprint_output:unknown tag type " << _tag << '\n';
  return 0;
}

template int TSubstructure_FP::do_fingerprint_output(Molecule& m, int nq, const int* hits,
                                                     std::ostream& output) const;

template <typename OUTPUT>
int
TSubstructure_FP::_do_sparse_fingerprint_output(const int nq, const int* hits,
                                                OUTPUT& output) const {
  Sparse_Fingerprint_Creator sfc;

  uint64_t total_hits = 0;
  uint32_t queries_matching = 0;
  for (int i = 0; i < nq; ++i) {
    if (0 == hits[i]) {
      continue;
    }

    total_hits += hits[i];
    ++queries_matching;

    if (!_each_query_writes_bit) {
      continue;
    }

    sfc.hit_bit(i * _bit_replicates, hits[i]);  // always.
    for (int j = 1; j < _bit_replicates; ++j) {
      sfc.hit_bit(i * _bit_replicates + j, hits[i]);
    }
  }

  if (!_each_query_writes_bit && queries_matching > 0) {
    const int count = std::min<uint32_t>(queries_matching, 255);
    for (int i = 0; i < _bit_replicates; ++i) {
      sfc.hit_bit(i, count);
    }
  } else if (_extra_bit_total_hits && total_hits > 0) {
    const int count = std::min<uint64_t>(total_hits, 255);
    for (int i = 0; i < _extra_bit_total_hits; ++i) {
      sfc.hit_bit(_bit_replicates * nq + i, count);
    }
  }

  sfc.write_fingerprint(_tag, output);
  if (!_work_as_filter) {
    output << "|\n";
  }

  return output.good();
}

template <typename OUTPUT>
int
TSubstructure_FP::_do_fingerprint_output(const int nq, const int* hits,
                                         OUTPUT& output) const {
  int bits_needed = 0;
  if (!_each_query_writes_bit) {
    bits_needed = _bit_replicates;
  } else if (_default_fingerprint_nbits) {
    bits_needed = _default_fingerprint_nbits * _bit_replicates;
  } else if (_bit_replicates > 0) {
    bits_needed = nq * _bit_replicates;
  } else {
    bits_needed = nq;
  }

  const int extra_bit_bstart = bits_needed;

  // In a binary fingerprint this just indicates whether any of the queries matched.
  if (!_each_query_writes_bit) {
  } else if (_extra_bit_total_hits) {
    bits_needed += _extra_bit_total_hits;
  }

  if (0 != bits_needed % 8) {
    bits_needed = (bits_needed / 8 + 1) * 8;
  }

  IW_Bits_Base fp;
  fp.allocate_space_for_bits(bits_needed);

  // If we are writing an extra bit for the total number of hits.
  int total_hits = 0;
  for (int i = 0; i < nq; ++i) {
    if (0 == hits[i]) {
      continue;
    }

    ++total_hits;
    if (!_each_query_writes_bit) {
      continue;
    }

    fp.set(i * _bit_replicates);
    for (int j = 1; j < _bit_replicates; ++j) {
      fp.set(i * _bit_replicates + j);
    }
  }

  if (total_hits > 0) {
    if (!_each_query_writes_bit) {
      for (int i = 0; i < _bit_replicates; ++i) {
        fp.set(i);
      }
    } else {
      for (int i = 0; i < _extra_bit_total_hits; ++i) {
        fp.set(extra_bit_bstart + i);
      }
    }
  }

  IWString tmp;
  fp.daylight_ascii_representation_including_nset_info(tmp);
  output << _tag << tmp << ">\n";
  if (!_work_as_filter) {
    output << "|\n";
  }

  return output.good();
}

static inline uint8_t
ToByte(int v) {
  static constexpr int kMaxuint8 = std::numeric_limits<uint8_t>::max();

  if (v > kMaxuint8) {
    return std::numeric_limits<uint8_t>::max();
  } else {
    return static_cast<uint8_t>(v);
  }
}

template <typename OUTPUT>
int
TSubstructure_FP::_do_fixed_counted_fingerprint_output(int nq, const int* hits,
                                                       OUTPUT& output) const {
  // One byte per feature; pad storage to whole 64-bit words.
  const int query_features = _each_query_writes_bit ? nq * _bit_replicates : 0;
  const int extra_features =
      _each_query_writes_bit ? _extra_bit_total_hits : _bit_replicates;
  const int nfeatures = query_features + extra_features;
  const int nbits = ((nfeatures + 7) / 8) * 64;
  fixed_bit_vector::FixedBitVector fp(nbits);
  uint8_t* counts = reinterpret_cast<uint8_t*>(const_cast<uint64_t*>(fp.bits()));

  int queries_matching = 0;
  uint64_t total_hits = 0;
  int nonzero_features = 0;
  for (int i = 0; i < nq; ++i) {
    if (hits[i] == 0) {
      continue;
    }
    ++queries_matching;
    total_hits += hits[i];
    if (_each_query_writes_bit) {
      for (int j = 0; j < _bit_replicates; ++j) {
        counts[i * _bit_replicates + j] = ToByte(hits[i]);
        ++nonzero_features;
      }
    }
  }
  const uint8_t extra_count =
      _each_query_writes_bit ? static_cast<uint8_t>(std::min<uint64_t>(total_hits, 255))
                             : ToByte(queries_matching);
  if (extra_count > 0) {
    for (int i = 0; i < extra_features; ++i) {
      counts[query_features + i] = extra_count;
      ++nonzero_features;
    }
  }

  IWString as_daylight = fp.DaylightAsciiRepresentation();
  as_daylight << ';' << nfeatures << ';' << nonzero_features;
  as_daylight << ';' << nfeatures << ';' << nonzero_features;

  output << _tag << as_daylight << ">\n";

  if (!_work_as_filter) {
    output << "|\n";
  }

  return output.good();
}

template <typename OUTPUT>
int
write_smiles_and_pcn(Molecule& m, OUTPUT& output) {
  output << smiles_tag << m.smiles() << ">\n";
  output << pcn_tag << m.name() << ">\n";

  return output.good();
}

template int write_smiles_and_pcn(Molecule&, std::ostream&);
