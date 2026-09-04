#ifndef SCCD_MATH_HPP
#define SCCD_MATH_HPP

#include <stdint.h>
#include <string.h>

#include <type_traits>

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
/**
 * \brief NaN test that survives -ffast-math.
 *
 * std::isnan compiles to a constant false under -ffinite-math-only, because the
 * flag promises the value cannot be NaN. Anything that *validates* a result must
 * not rest on that promise: the benchmark datasets encode "this query has no
 * collision" as a NaN root, so a folded isnan silently turns every correct
 * no-collision answer into a reported missed collision, and inverts the meaning
 * of the acceptance gate. Reading the bits cannot be folded away.
 */
template <typename T> static inline bool is_nan_bits(const T v) {
  static_assert(std::is_floating_point<T>::value, "is_nan_bits expects a floating-point type");
  typedef typename std::conditional<sizeof(T) == 4, uint32_t, uint64_t>::type U;
  static_assert(sizeof(U) == sizeof(T), "unsupported floating-point width");
  U bits;
  memcpy(&bits, &v, sizeof(T));
  const U exponent = (sizeof(T) == 4) ? (U)0x7F800000u : (U)0x7FF0000000000000ull;
  const U mantissa = (sizeof(T) == 4) ? (U)0x007FFFFFu : (U)0x000FFFFFFFFFFFFFull;
  return (bits & exponent) == exponent && (bits & mantissa) != 0;
}

/** \brief std::isfinite without the -ffast-math folding. See is_nan_bits. */
template <typename T> static inline bool is_finite_bits(const T v) {
  static_assert(std::is_floating_point<T>::value, "is_finite_bits expects a floating-point type");
  typedef typename std::conditional<sizeof(T) == 4, uint32_t, uint64_t>::type U;
  U bits;
  memcpy(&bits, &v, sizeof(T));
  const U exponent = (sizeof(T) == 4) ? (U)0x7F800000u : (U)0x7FF0000000000000ull;
  return (bits & exponent) != exponent;
}

} // namespace sccd

#endif // SCCD_MATH_HPP
