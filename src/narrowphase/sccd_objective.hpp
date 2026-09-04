// VF
/* Auto-generated with SymPy (CSE + pow expansion) */
#ifndef SCCD_OBJECTIVE_HPP
#define SCCD_OBJECTIVE_HPP

#include "sccd_math.hpp"

namespace sccd {


template <typename T>
static inline void vf_objective(const T sv[3],
                                const T s1[3],
                                const T s2[3],
                                const T s3[3],
                                const T s4[3],
                                const T ev[3],
                                const T e1[3],
                                const T e2[3],
                                const T e3[3],
                                const T e4[3],
                                const T t,
                                const T u,
                                const T v,
                                T *out_f) {
    const T ssa0 = t - 1;
    const T ssa1 = u + v - 1;
    *out_f = (1.0 / 2.0) * sccd::pow2<T>(-ev[0] * t + ssa0 * sv[0] - ssa1 * (e1[0] * t - s1[0] * ssa0) +
                                         u * (e2[0] * t - s2[0] * ssa0) + v * (e3[0] * t - s3[0] * ssa0)) +
             (1.0 / 2.0) * sccd::pow2<T>(-ev[1] * t + ssa0 * sv[1] - ssa1 * (e1[1] * t - s1[1] * ssa0) +
                                         u * (e2[1] * t - s2[1] * ssa0) + v * (e3[1] * t - s3[1] * ssa0)) +
             (1.0 / 2.0) * sccd::pow2<T>(-ev[2] * t + ssa0 * sv[2] - ssa1 * (e1[2] * t - s1[2] * ssa0) +
                                         u * (e2[2] * t - s2[2] * ssa0) + v * (e3[2] * t - s3[2] * ssa0));
}


template <typename T>
static inline void vf_objective_dir(const T sv[3],
                                    const T s1[3],
                                    const T s2[3],
                                    const T s3[3],
                                    const T s4[3],
                                    const T ev[3],
                                    const T e1[3],
                                    const T e2[3],
                                    const T e3[3],
                                    const T e4[3],
                                    const T t,
                                    const T u,
                                    const T v,
                                    T *out_f,
                                    T out_p[3]) {
    const T ssa0 = t - 1;
    const T ssa1 = e2[0] * t;
    const T ssa2 = s2[0] * ssa0;
    const T ssa3 = e3[0] * t;
    const T ssa4 = s3[0] * ssa0;
    const T ssa5 = u + v - 1;
    const T ssa6 = e1[0] * t - s1[0] * ssa0;
    const T ssa7 = -ev[0] * t + ssa0 * sv[0] - ssa5 * ssa6 + u * (ssa1 - ssa2) + v * (ssa3 - ssa4);
    const T ssa8 = e2[1] * t;
    const T ssa9 = s2[1] * ssa0;
    const T ssa10 = e3[1] * t;
    const T ssa11 = s3[1] * ssa0;
    const T ssa12 = e1[1] * t - s1[1] * ssa0;
    const T ssa13 = -ev[1] * t + ssa0 * sv[1] - ssa12 * ssa5 + u * (ssa8 - ssa9) + v * (ssa10 - ssa11);
    const T ssa14 = e2[2] * t;
    const T ssa15 = s2[2] * ssa0;
    const T ssa16 = e3[2] * t;
    const T ssa17 = s3[2] * ssa0;
    const T ssa18 = e1[2] * t - s1[2] * ssa0;
    const T ssa19 = -ev[2] * t + ssa0 * sv[2] - ssa18 * ssa5 + u * (ssa14 - ssa15) + v * (ssa16 - ssa17);
    const T ssa20 = e1[0] - s1[0];
    const T ssa21 = -ev[0] - ssa20 * ssa5 + sv[0] + u * (e2[0] - s2[0]) + v * (e3[0] - s3[0]);
    const T ssa22 = e1[1] - s1[1];
    const T ssa23 = -ev[1] - ssa22 * ssa5 + sv[1] + u * (e2[1] - s2[1]) + v * (e3[1] - s3[1]);
    const T ssa24 = e1[2] - s1[2];
    const T ssa25 = -ev[2] - ssa24 * ssa5 + sv[2] + u * (e2[2] - s2[2]) + v * (e3[2] - s3[2]);
    const T ssa26 = 1.0 / (sccd::pow2<T>(ssa21) + sccd::pow2<T>(ssa23) + sccd::pow2<T>(ssa25));
    const T ssa27 = -ssa1 + ssa2 + ssa6;
    const T ssa28 = ssa12 - ssa8 + ssa9;
    const T ssa29 = -ssa14 + ssa15 + ssa18;
    const T ssa30 = -ssa13 * (-e2[1] + s2[1] + ssa22) - ssa19 * (-e2[2] + s2[2] + ssa24) - ssa21 * ssa27 -
                    ssa23 * ssa28 - ssa25 * ssa29 - ssa7 * (-e2[0] + s2[0] + ssa20);
    const T ssa31 =
        1.0 / (-ssa26 * sccd::pow2<T>(ssa30) + sccd::pow2<T>(ssa27) + sccd::pow2<T>(ssa28) + sccd::pow2<T>(ssa29));
    const T ssa32 = -ssa3 + ssa4 + ssa6;
    const T ssa33 = -ssa10 + ssa11 + ssa12;
    const T ssa34 = -ssa16 + ssa17 + ssa18;
    const T ssa35 = -ssa13 * (-e3[1] + s3[1] + ssa22) - ssa19 * (-e3[2] + s3[2] + ssa24) - ssa21 * ssa32 -
                    ssa23 * ssa33 - ssa25 * ssa34 - ssa7 * (-e3[0] + s3[0] + ssa20);
    const T ssa36 = ssa26 * ssa30;
    const T ssa37 = ssa27 * ssa32 + ssa28 * ssa33 + ssa29 * ssa34 - ssa35 * ssa36;
    const T ssa38 = ssa13 * ssa23 + ssa19 * ssa25 + ssa21 * ssa7;
    const T ssa39 = ssa13 * ssa28 + ssa19 * ssa29 + ssa27 * ssa7 + ssa36 * ssa38;
    const T ssa40 = (-ssa13 * ssa33 - ssa19 * ssa34 - ssa26 * ssa35 * ssa38 + ssa31 * ssa37 * ssa39 - ssa32 * ssa7) /
                    (-ssa26 * sccd::pow2<T>(ssa35) - ssa31 * sccd::pow2<T>(ssa37) + sccd::pow2<T>(ssa32) +
                     sccd::pow2<T>(ssa33) + sccd::pow2<T>(ssa34));
    const T ssa41 = ssa31 * (-ssa37 * ssa40 - ssa39);
    *out_f =
        (1.0 / 2.0) * sccd::pow2<T>(ssa13) + (1.0 / 2.0) * sccd::pow2<T>(ssa19) + (1.0 / 2.0) * sccd::pow2<T>(ssa7);
    out_p[0] = ssa26 * (-ssa30 * ssa41 - ssa35 * ssa40 + ssa38);
    out_p[1] = ssa41;
    out_p[2] = ssa40;
}



}  // namespace sccd

#endif