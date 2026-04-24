#include <cstdint>
#include <cstddef>
#include <immintrin.h>

#if __has_include(<bit>)
  #include <bit>
#endif

namespace iwbits {

// ------------------------------------------------------------
// Thanks ChatGPT
// Unfortunately it looks like this does not make any difference between
// the 64 bit version.

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__)

uint64_t
intersection_popcount_avx512(const uint64_t* a,
                             const uint64_t* b,
                             std::size_t n_words) {
    std::size_t i = 0;
    uint64_t total = 0;

    alignas(64) uint64_t lanes[8];

    for (; i + 8 <= n_words; i += 8) {
        __m512i va = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(a + i));
        __m512i vb = _mm512_loadu_si512(reinterpret_cast<const __m512i*>(b + i));
        __m512i vand = _mm512_and_si512(va, vb);

        // Popcount each 64-bit lane
        __m512i vcnt = _mm512_popcnt_epi64(vand);

        // Store and horizontally sum 8 lanes
        _mm512_store_si512(reinterpret_cast<__m512i*>(lanes), vcnt);
        total += lanes[0] + lanes[1] + lanes[2] + lanes[3]
              +  lanes[4] + lanes[5] + lanes[6] + lanes[7];
    }

    // Tail
    for (; i < n_words; ++i) {
        total += _mm_popcnt_u64(a[i] & b[i]);
    }

    return total;
}

uint64_t
IntersectionPopcount(const uint64_t* a, const uint64_t* b, std::size_t nwords) {
  uint64_t rc = 0;

  for (uint64_t i = 0; i < nwords; ++i) {
    rc += _mm_popcnt_u64(a[i] & b[i]);
  }

  return rc;
}

uint64_t
BitsInCommonAvx(const uint64_t* a, const uint64_t* b, std::size_t nwords) {
  uint64_t bic;
  if (nwords < 8) {
    bic = IntersectionPopcount(a, b, nwords);
  } else {
     bic = intersection_popcount_avx512(a, b, nwords);
  }

  return bic;
}

#else
#endif

}  // namespace iwbits
