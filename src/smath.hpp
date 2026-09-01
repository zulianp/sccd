#ifndef SCCD_MATH_HPP
#define SCCD_MATH_HPP

#include <stdint.h>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace sccd {

/** \brief Number of set bits in a 32-bit lane mask. */
static inline int popcount32(const uint32_t x) {
#ifdef _MSC_VER
  return (int)__popcnt(x);
#else
  return __builtin_popcount(x);
#endif
}

/** \brief Index of the lowest set bit. Undefined for x == 0. */
static inline int ctz32(const uint32_t x) {
#ifdef _MSC_VER
  unsigned long i;
  _BitScanForward(&i, x);
  return (int)i;
#else
  return __builtin_ctz(x);
#endif
}

template <typename T> static inline T max(const T a, const T b) {
  return (a > b) ? a : b;
}

template <typename T> static inline T min(const T a, const T b) {
  return (a < b) ? a : b;
}

template <typename T> static inline T abs(const T x) { return x < 0 ? -x : x; }

template <typename T> static inline T pow2(const T x) { return x*x; }
template <typename T> static inline T pow3(const T x) { return x*x*x; }

template <typename T> static inline T array_min(const int n, const T* arr) {
  T min = arr[0];
  for (int i = 1; i < n; i++) {
    if (arr[i] < min) {
      min = arr[i];
    }
  }
  return min;
}
template <typename T> static inline T array_max(const int n, const T* arr) {
  T max = arr[0];
  for (int i = 1; i < n; i++) {
    if (arr[i] > max) {
      max = arr[i];
    }
  }
  return max;
}
} // namespace sccd

#endif // SCCD_MATH_HPP
