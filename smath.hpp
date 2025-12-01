#ifndef SCCD_MATH_HPP
#define SCCD_MATH_HPP

namespace sccd {
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
