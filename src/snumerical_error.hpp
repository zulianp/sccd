#ifndef SCCD_SNUMERICAL_ERROR_HPP
#define SCCD_SNUMERICAL_ERROR_HPP

#include "smath.hpp"

#include <cmath>
#include <limits>

namespace sccd {

    template <bool is_vertex_face, typename T>
    static inline T numerical_error_bound_component(const T s0,
                                                     const T s1,
                                                     const T s2,
                                                     const T s3,
                                                     const T e0,
                                                     const T e1,
                                                     const T e2,
                                                     const T e3) {
        const T max_abs = sccd::max<T>(
            sccd::max<T>(sccd::max<T>(sccd::abs<T>(s0), sccd::abs<T>(s1)),
                          sccd::max<T>(sccd::abs<T>(s2), sccd::abs<T>(s3))),
            sccd::max<T>(sccd::max<T>(sccd::abs<T>(e0), sccd::abs<T>(e1)),
                          sccd::max<T>(sccd::abs<T>(e2), sccd::abs<T>(e3))));
        const T delta = sccd::min<T>(max_abs, T(1));
        const T filter = T(is_vertex_face ? 30 : 28) * std::numeric_limits<T>::epsilon();
        return filter * std::pow(delta, T(3));
    }

    template <bool is_vertex_face, typename T>
    static inline void numerical_error_bound(const T s0[3],
                                             const T s1[3],
                                             const T s2[3],
                                             const T s3[3],
                                             const T e0[3],
                                             const T e1[3],
                                             const T e2[3],
                                             const T e3[3],
                                             T error[3]) {
        for (int d = 0; d < 3; ++d) {
            error[d] = numerical_error_bound_component<is_vertex_face, T>(
                s0[d], s1[d], s2[d], s3[d], e0[d], e1[d], e2[d], e3[d]);
        }
    }

}  // namespace sccd

#endif
