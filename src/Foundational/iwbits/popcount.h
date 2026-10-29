#ifndef POPCOUNT_H
#define POPCOUNT_H

// The compiler builtin also works when the target does not enable POPCNT.
// The _mm_popcnt_* intrinsics require -mpopcnt and break portable builds.
#if defined(__x86_64__)
#define POPCOUNT __builtin_popcountll
#elif defined(__i386__)
#define POPCOUNT __builtin_popcount
#elif __APPLE__
#define POPCOUNT __builtin_popcountll
#else
#include <bit>
#define POPCOUNT std::popcount
#endif

#endif /* POPCOUNT_H */
