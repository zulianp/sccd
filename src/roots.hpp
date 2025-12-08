// VF
/* Auto-generated with SymPy (CSE + pow expansion) */
#ifndef ROOTS_HPP
#define ROOTS_HPP

#include "smath.hpp"

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
static inline void vf_gradient(const T sv[3],
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
                               T out_g[3]) {
    const T ssa0 = u + v - 1;
    const T ssa1 = t - 1;
    const T ssa2 = e2[0] * t;
    const T ssa3 = s2[0] * ssa1;
    const T ssa4 = e3[0] * t;
    const T ssa5 = s3[0] * ssa1;
    const T ssa6 = e1[0] * t - s1[0] * ssa1;
    const T ssa7 = -ev[0] * t - ssa0 * ssa6 + ssa1 * sv[0] + u * (ssa2 - ssa3) + v * (ssa4 - ssa5);
    const T ssa8 = e2[1] * t;
    const T ssa9 = s2[1] * ssa1;
    const T ssa10 = e3[1] * t;
    const T ssa11 = s3[1] * ssa1;
    const T ssa12 = e1[1] * t - s1[1] * ssa1;
    const T ssa13 = -ev[1] * t - ssa0 * ssa12 + ssa1 * sv[1] + u * (ssa8 - ssa9) + v * (ssa10 - ssa11);
    const T ssa14 = e2[2] * t;
    const T ssa15 = s2[2] * ssa1;
    const T ssa16 = e3[2] * t;
    const T ssa17 = s3[2] * ssa1;
    const T ssa18 = e1[2] * t - s1[2] * ssa1;
    const T ssa19 = -ev[2] * t - ssa0 * ssa18 + ssa1 * sv[2] + u * (ssa14 - ssa15) + v * (ssa16 - ssa17);
    out_g[0] = ssa13 * (-ev[1] - ssa0 * (e1[1] - s1[1]) + sv[1] + u * (e2[1] - s2[1]) + v * (e3[1] - s3[1])) +
               ssa19 * (-ev[2] - ssa0 * (e1[2] - s1[2]) + sv[2] + u * (e2[2] - s2[2]) + v * (e3[2] - s3[2])) +
               ssa7 * (-ev[0] - ssa0 * (e1[0] - s1[0]) + sv[0] + u * (e2[0] - s2[0]) + v * (e3[0] - s3[0]));
    out_g[1] = -ssa13 * (ssa12 - ssa8 + ssa9) - ssa19 * (-ssa14 + ssa15 + ssa18) - ssa7 * (-ssa2 + ssa3 + ssa6);
    out_g[2] = -ssa13 * (-ssa10 + ssa11 + ssa12) - ssa19 * (-ssa16 + ssa17 + ssa18) - ssa7 * (-ssa4 + ssa5 + ssa6);
}

template <typename T>
static inline void vf_hessian(const T sv[3],
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
                              T out_H[9]) {
    const T ssa0 = u + v - 1;
    const T ssa1 = e1[0] - s1[0];
    const T ssa2 = -ev[0] - ssa0 * ssa1 + sv[0] + u * (e2[0] - s2[0]) + v * (e3[0] - s3[0]);
    const T ssa3 = e1[1] - s1[1];
    const T ssa4 = -ev[1] - ssa0 * ssa3 + sv[1] + u * (e2[1] - s2[1]) + v * (e3[1] - s3[1]);
    const T ssa5 = e1[2] - s1[2];
    const T ssa6 = -ev[2] - ssa0 * ssa5 + sv[2] + u * (e2[2] - s2[2]) + v * (e3[2] - s3[2]);
    const T ssa7 = t - 1;
    const T ssa8 = s2[0] * ssa7;
    const T ssa9 = e2[0] * t;
    const T ssa10 = e1[0] * t - s1[0] * ssa7;
    const T ssa11 = ssa10 + ssa8 - ssa9;
    const T ssa12 = s2[1] * ssa7;
    const T ssa13 = e2[1] * t;
    const T ssa14 = e1[1] * t - s1[1] * ssa7;
    const T ssa15 = ssa12 - ssa13 + ssa14;
    const T ssa16 = s2[2] * ssa7;
    const T ssa17 = e2[2] * t;
    const T ssa18 = e1[2] * t - s1[2] * ssa7;
    const T ssa19 = ssa16 - ssa17 + ssa18;
    const T ssa20 = e3[0] * t;
    const T ssa21 = s3[0] * ssa7;
    const T ssa22 = -ev[0] * t - ssa0 * ssa10 + ssa7 * sv[0] + u * (-ssa8 + ssa9) + v * (ssa20 - ssa21);
    const T ssa23 = e3[1] * t;
    const T ssa24 = s3[1] * ssa7;
    const T ssa25 = -ev[1] * t - ssa0 * ssa14 + ssa7 * sv[1] + u * (-ssa12 + ssa13) + v * (ssa23 - ssa24);
    const T ssa26 = e3[2] * t;
    const T ssa27 = s3[2] * ssa7;
    const T ssa28 = -ev[2] * t - ssa0 * ssa18 + ssa7 * sv[2] + u * (-ssa16 + ssa17) + v * (ssa26 - ssa27);
    const T ssa29 = -ssa11 * ssa2 - ssa15 * ssa4 - ssa19 * ssa6 - ssa22 * (-e2[0] + s2[0] + ssa1) -
                    ssa25 * (-e2[1] + s2[1] + ssa3) - ssa28 * (-e2[2] + s2[2] + ssa5);
    const T ssa30 = ssa10 - ssa20 + ssa21;
    const T ssa31 = ssa14 - ssa23 + ssa24;
    const T ssa32 = ssa18 - ssa26 + ssa27;
    const T ssa33 = -ssa2 * ssa30 - ssa22 * (-e3[0] + s3[0] + ssa1) - ssa25 * (-e3[1] + s3[1] + ssa3) -
                    ssa28 * (-e3[2] + s3[2] + ssa5) - ssa31 * ssa4 - ssa32 * ssa6;
    const T ssa34 = ssa11 * ssa30 + ssa15 * ssa31 + ssa19 * ssa32;
    out_H[0] = sccd::pow2<T>(ssa2) + sccd::pow2<T>(ssa4) + sccd::pow2<T>(ssa6);
    out_H[1] = ssa29;
    out_H[2] = ssa33;
    out_H[3] = ssa29;
    out_H[4] = sccd::pow2<T>(ssa11) + sccd::pow2<T>(ssa15) + sccd::pow2<T>(ssa19);
    out_H[5] = ssa34;
    out_H[6] = ssa33;
    out_H[7] = ssa34;
    out_H[8] = sccd::pow2<T>(ssa30) + sccd::pow2<T>(ssa31) + sccd::pow2<T>(ssa32);
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

template <typename T>
static inline void vf_all(const T sv[3],
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
                          T out_g[3],
                          T out_H[9]) {
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
    const T ssa26 = -ssa1 + ssa2 + ssa6;
    const T ssa27 = ssa12 - ssa8 + ssa9;
    const T ssa28 = -ssa14 + ssa15 + ssa18;
    const T ssa29 = -ssa3 + ssa4 + ssa6;
    const T ssa30 = -ssa10 + ssa11 + ssa12;
    const T ssa31 = -ssa16 + ssa17 + ssa18;
    const T ssa32 = -ssa13 * (-e2[1] + s2[1] + ssa22) - ssa19 * (-e2[2] + s2[2] + ssa24) - ssa21 * ssa26 -
                    ssa23 * ssa27 - ssa25 * ssa28 - ssa7 * (-e2[0] + s2[0] + ssa20);
    const T ssa33 = -ssa13 * (-e3[1] + s3[1] + ssa22) - ssa19 * (-e3[2] + s3[2] + ssa24) - ssa21 * ssa29 -
                    ssa23 * ssa30 - ssa25 * ssa31 - ssa7 * (-e3[0] + s3[0] + ssa20);
    const T ssa34 = ssa26 * ssa29 + ssa27 * ssa30 + ssa28 * ssa31;
    *out_f =
        (1.0 / 2.0) * sccd::pow2<T>(ssa13) + (1.0 / 2.0) * sccd::pow2<T>(ssa19) + (1.0 / 2.0) * sccd::pow2<T>(ssa7);
    out_g[0] = ssa13 * ssa23 + ssa19 * ssa25 + ssa21 * ssa7;
    out_g[1] = -ssa13 * ssa27 - ssa19 * ssa28 - ssa26 * ssa7;
    out_g[2] = -ssa13 * ssa30 - ssa19 * ssa31 - ssa29 * ssa7;
    out_H[0] = sccd::pow2<T>(ssa21) + sccd::pow2<T>(ssa23) + sccd::pow2<T>(ssa25);
    out_H[1] = ssa32;
    out_H[2] = ssa33;
    out_H[3] = ssa32;
    out_H[4] = sccd::pow2<T>(ssa26) + sccd::pow2<T>(ssa27) + sccd::pow2<T>(ssa28);
    out_H[5] = ssa34;
    out_H[6] = ssa33;
    out_H[7] = ssa34;
    out_H[8] = sccd::pow2<T>(ssa29) + sccd::pow2<T>(ssa30) + sccd::pow2<T>(ssa31);
}
// EE
/* Auto-generated with SymPy (CSE + pow expansion) */

#include "smath.hpp"

template <typename T>
static inline void ee_objective(const T sv[3],
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
    const T ssa1 = -s3[0] * ssa0;
    const T ssa2 = v - 1;
    const T ssa3 = u - 1;
    const T ssa4 = -s3[1] * ssa0;
    const T ssa5 = -s3[2] * ssa0;
    *out_f = (1.0 / 2.0) * sccd::pow2<T>(ssa2 * (e3[0] * t + ssa1) - ssa3 * (e1[0] * t - s1[0] * ssa0) +
                                         u * (e2[0] * t - s2[0] * ssa0) - v * (e4[0] * t + ssa1)) +
             (1.0 / 2.0) * sccd::pow2<T>(ssa2 * (e3[1] * t + ssa4) - ssa3 * (e1[1] * t - s1[1] * ssa0) +
                                         u * (e2[1] * t - s2[1] * ssa0) - v * (e4[1] * t + ssa4)) +
             (1.0 / 2.0) * sccd::pow2<T>(ssa2 * (e3[2] * t + ssa5) - ssa3 * (e1[2] * t - s1[2] * ssa0) +
                                         u * (e2[2] * t - s2[2] * ssa0) - v * (e4[2] * t + ssa5));
}

template <typename T>
static inline void ee_gradient(const T sv[3],
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
                               T out_g[3]) {
    const T ssa0 = -s3[0];
    const T ssa1 = v - 1;
    const T ssa2 = u - 1;
    const T ssa3 = e2[0] * t;
    const T ssa4 = t - 1;
    const T ssa5 = s2[0] * ssa4;
    const T ssa6 = -s3[0] * ssa4;
    const T ssa7 = e1[0] * t - s1[0] * ssa4;
    const T ssa8 = ssa1 * (e3[0] * t + ssa6) - ssa2 * ssa7 + u * (ssa3 - ssa5) - v * (e4[0] * t + ssa6);
    const T ssa9 = -s3[1];
    const T ssa10 = e2[1] * t;
    const T ssa11 = s2[1] * ssa4;
    const T ssa12 = -s3[1] * ssa4;
    const T ssa13 = e1[1] * t - s1[1] * ssa4;
    const T ssa14 = ssa1 * (e3[1] * t + ssa12) - ssa13 * ssa2 + u * (ssa10 - ssa11) - v * (e4[1] * t + ssa12);
    const T ssa15 = -s3[2];
    const T ssa16 = e2[2] * t;
    const T ssa17 = s2[2] * ssa4;
    const T ssa18 = -s3[2] * ssa4;
    const T ssa19 = e1[2] * t - s1[2] * ssa4;
    const T ssa20 = ssa1 * (e3[2] * t + ssa18) - ssa19 * ssa2 + u * (ssa16 - ssa17) - v * (e4[2] * t + ssa18);
    out_g[0] = ssa14 * (ssa1 * (e3[1] + ssa9) - ssa2 * (e1[1] - s1[1]) + u * (e2[1] - s2[1]) - v * (e4[1] + ssa9)) +
               ssa20 * (ssa1 * (e3[2] + ssa15) - ssa2 * (e1[2] - s1[2]) + u * (e2[2] - s2[2]) - v * (e4[2] + ssa15)) +
               ssa8 * (ssa1 * (e3[0] + ssa0) - ssa2 * (e1[0] - s1[0]) + u * (e2[0] - s2[0]) - v * (e4[0] + ssa0));
    out_g[1] = -ssa14 * (-ssa10 + ssa11 + ssa13) - ssa20 * (-ssa16 + ssa17 + ssa19) - ssa8 * (-ssa3 + ssa5 + ssa7);
    out_g[2] = t * (ssa14 * (e3[1] - e4[1]) + ssa20 * (e3[2] - e4[2]) + ssa8 * (e3[0] - e4[0]));
}

template <typename T>
static inline void ee_hessian(const T sv[3],
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
                              T out_H[9]) {
    const T ssa0 = -s3[0];
    const T ssa1 = v - 1;
    const T ssa2 = u - 1;
    const T ssa3 = e1[0] - s1[0];
    const T ssa4 = ssa1 * (e3[0] + ssa0) - ssa2 * ssa3 + u * (e2[0] - s2[0]) - v * (e4[0] + ssa0);
    const T ssa5 = -s3[1];
    const T ssa6 = e1[1] - s1[1];
    const T ssa7 = ssa1 * (e3[1] + ssa5) - ssa2 * ssa6 + u * (e2[1] - s2[1]) - v * (e4[1] + ssa5);
    const T ssa8 = -s3[2];
    const T ssa9 = e1[2] - s1[2];
    const T ssa10 = ssa1 * (e3[2] + ssa8) - ssa2 * ssa9 + u * (e2[2] - s2[2]) - v * (e4[2] + ssa8);
    const T ssa11 = t - 1;
    const T ssa12 = s2[0] * ssa11;
    const T ssa13 = e2[0] * t;
    const T ssa14 = e1[0] * t - s1[0] * ssa11;
    const T ssa15 = ssa12 - ssa13 + ssa14;
    const T ssa16 = s2[1] * ssa11;
    const T ssa17 = e2[1] * t;
    const T ssa18 = e1[1] * t - s1[1] * ssa11;
    const T ssa19 = ssa16 - ssa17 + ssa18;
    const T ssa20 = s2[2] * ssa11;
    const T ssa21 = e2[2] * t;
    const T ssa22 = e1[2] * t - s1[2] * ssa11;
    const T ssa23 = ssa20 - ssa21 + ssa22;
    const T ssa24 = -s3[0] * ssa11;
    const T ssa25 = ssa1 * (e3[0] * t + ssa24) - ssa14 * ssa2 + u * (-ssa12 + ssa13) - v * (e4[0] * t + ssa24);
    const T ssa26 = -s3[1] * ssa11;
    const T ssa27 = ssa1 * (e3[1] * t + ssa26) - ssa18 * ssa2 + u * (-ssa16 + ssa17) - v * (e4[1] * t + ssa26);
    const T ssa28 = -s3[2] * ssa11;
    const T ssa29 = ssa1 * (e3[2] * t + ssa28) - ssa2 * ssa22 + u * (-ssa20 + ssa21) - v * (e4[2] * t + ssa28);
    const T ssa30 = -ssa10 * ssa23 - ssa15 * ssa4 - ssa19 * ssa7 - ssa25 * (-e2[0] + s2[0] + ssa3) -
                    ssa27 * (-e2[1] + s2[1] + ssa6) - ssa29 * (-e2[2] + s2[2] + ssa9);
    const T ssa31 = e3[0] - e4[0];
    const T ssa32 = e3[1] - e4[1];
    const T ssa33 = e3[2] - e4[2];
    const T ssa34 =
        ssa10 * ssa33 * t + ssa25 * ssa31 + ssa27 * ssa32 + ssa29 * ssa33 + ssa31 * ssa4 * t + ssa32 * ssa7 * t;
    const T ssa35 = t * (-ssa15 * ssa31 - ssa19 * ssa32 - ssa23 * ssa33);
    out_H[0] = sccd::pow2<T>(ssa10) + sccd::pow2<T>(ssa4) + sccd::pow2<T>(ssa7);
    out_H[1] = ssa30;
    out_H[2] = ssa34;
    out_H[3] = ssa30;
    out_H[4] = sccd::pow2<T>(ssa15) + sccd::pow2<T>(ssa19) + sccd::pow2<T>(ssa23);
    out_H[5] = ssa35;
    out_H[6] = ssa34;
    out_H[7] = ssa35;
    out_H[8] = (sccd::pow2<T>(ssa31) + sccd::pow2<T>(ssa32) + sccd::pow2<T>(ssa33)) * (t * t);
}

template <typename T>
static inline void ee_objective_dir(const T sv[3],
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
    const T ssa0 = e2[0] * t;
    const T ssa1 = t - 1;
    const T ssa2 = s2[0] * ssa1;
    const T ssa3 = -s3[0] * ssa1;
    const T ssa4 = v - 1;
    const T ssa5 = u - 1;
    const T ssa6 = e1[0] * t - s1[0] * ssa1;
    const T ssa7 = ssa4 * (e3[0] * t + ssa3) - ssa5 * ssa6 + u * (ssa0 - ssa2) - v * (e4[0] * t + ssa3);
    const T ssa8 = e2[1] * t;
    const T ssa9 = s2[1] * ssa1;
    const T ssa10 = -s3[1] * ssa1;
    const T ssa11 = e1[1] * t - s1[1] * ssa1;
    const T ssa12 = -ssa11 * ssa5 + ssa4 * (e3[1] * t + ssa10) + u * (ssa8 - ssa9) - v * (e4[1] * t + ssa10);
    const T ssa13 = e2[2] * t;
    const T ssa14 = s2[2] * ssa1;
    const T ssa15 = -s3[2] * ssa1;
    const T ssa16 = e1[2] * t - s1[2] * ssa1;
    const T ssa17 = -ssa16 * ssa5 + ssa4 * (e3[2] * t + ssa15) + u * (ssa13 - ssa14) - v * (e4[2] * t + ssa15);
    const T ssa18 = -s3[0];
    const T ssa19 = e1[0] - s1[0];
    const T ssa20 = -ssa19 * ssa5 + ssa4 * (e3[0] + ssa18) + u * (e2[0] - s2[0]) - v * (e4[0] + ssa18);
    const T ssa21 = -s3[1];
    const T ssa22 = e1[1] - s1[1];
    const T ssa23 = -ssa22 * ssa5 + ssa4 * (e3[1] + ssa21) + u * (e2[1] - s2[1]) - v * (e4[1] + ssa21);
    const T ssa24 = -s3[2];
    const T ssa25 = e1[2] - s1[2];
    const T ssa26 = -ssa25 * ssa5 + ssa4 * (e3[2] + ssa24) + u * (e2[2] - s2[2]) - v * (e4[2] + ssa24);
    const T ssa27 = 1.0 / (sccd::pow2<T>(ssa20) + sccd::pow2<T>(ssa23) + sccd::pow2<T>(ssa26));
    const T ssa28 = -ssa0 + ssa2 + ssa6;
    const T ssa29 = ssa11 - ssa8 + ssa9;
    const T ssa30 = -ssa13 + ssa14 + ssa16;
    const T ssa31 = -ssa12 * (-e2[1] + s2[1] + ssa22) - ssa17 * (-e2[2] + s2[2] + ssa25) - ssa20 * ssa28 -
                    ssa23 * ssa29 - ssa26 * ssa30 - ssa7 * (-e2[0] + s2[0] + ssa19);
    const T ssa32 =
        1.0 / (-ssa27 * sccd::pow2<T>(ssa31) + sccd::pow2<T>(ssa28) + sccd::pow2<T>(ssa29) + sccd::pow2<T>(ssa30));
    const T ssa33 = e3[0] - e4[0];
    const T ssa34 = e3[1] - e4[1];
    const T ssa35 = e3[2] - e4[2];
    const T ssa36 = ssa12 * ssa34 + ssa17 * ssa35 + ssa33 * ssa7;
    const T ssa37 = ssa20 * ssa33 * t + ssa23 * ssa34 * t + ssa26 * ssa35 * t + ssa36;
    const T ssa38 = ssa27 * ssa31;
    const T ssa39 = -ssa37 * ssa38 + t * (-ssa28 * ssa33 - ssa29 * ssa34 - ssa30 * ssa35);
    const T ssa40 = ssa12 * ssa23 + ssa17 * ssa26 + ssa20 * ssa7;
    const T ssa41 = ssa12 * ssa29 + ssa17 * ssa30 + ssa28 * ssa7 + ssa38 * ssa40;
    const T ssa42 = -ssa27 * ssa37 * ssa40 + ssa32 * ssa39 * ssa41 + ssa36 * t;
    const T ssa43 = -1.0 / (ssa27 * sccd::pow2<T>(ssa37) + ssa32 * sccd::pow2<T>(ssa39) -
                            (sccd::pow2<T>(ssa33) + sccd::pow2<T>(ssa34) + sccd::pow2<T>(ssa35)) * t * t);
    const T ssa44 = -ssa39 * ssa42 * ssa43 - ssa41;
    *out_f =
        (1.0 / 2.0) * sccd::pow2<T>(ssa12) + (1.0 / 2.0) * sccd::pow2<T>(ssa17) + (1.0 / 2.0) * sccd::pow2<T>(ssa7);
    out_p[0] = ssa27 * (-ssa31 * ssa32 * ssa44 - ssa37 * ssa42 * ssa43 + ssa40);
    out_p[1] = ssa32 * ssa44;
    out_p[2] = ssa42 * ssa43;
}

template <typename T>
static inline void ee_all(const T sv[3],
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
                          T out_g[3],
                          T out_H[9]) {
    const T ssa0 = e2[0] * t;
    const T ssa1 = t - 1;
    const T ssa2 = s2[0] * ssa1;
    const T ssa3 = -s3[0] * ssa1;
    const T ssa4 = v - 1;
    const T ssa5 = u - 1;
    const T ssa6 = e1[0] * t - s1[0] * ssa1;
    const T ssa7 = ssa4 * (e3[0] * t + ssa3) - ssa5 * ssa6 + u * (ssa0 - ssa2) - v * (e4[0] * t + ssa3);
    const T ssa8 = e2[1] * t;
    const T ssa9 = s2[1] * ssa1;
    const T ssa10 = -s3[1] * ssa1;
    const T ssa11 = e1[1] * t - s1[1] * ssa1;
    const T ssa12 = -ssa11 * ssa5 + ssa4 * (e3[1] * t + ssa10) + u * (ssa8 - ssa9) - v * (e4[1] * t + ssa10);
    const T ssa13 = e2[2] * t;
    const T ssa14 = s2[2] * ssa1;
    const T ssa15 = -s3[2] * ssa1;
    const T ssa16 = e1[2] * t - s1[2] * ssa1;
    const T ssa17 = -ssa16 * ssa5 + ssa4 * (e3[2] * t + ssa15) + u * (ssa13 - ssa14) - v * (e4[2] * t + ssa15);
    const T ssa18 = -s3[0];
    const T ssa19 = e1[0] - s1[0];
    const T ssa20 = -ssa19 * ssa5 + ssa4 * (e3[0] + ssa18) + u * (e2[0] - s2[0]) - v * (e4[0] + ssa18);
    const T ssa21 = -s3[1];
    const T ssa22 = e1[1] - s1[1];
    const T ssa23 = -ssa22 * ssa5 + ssa4 * (e3[1] + ssa21) + u * (e2[1] - s2[1]) - v * (e4[1] + ssa21);
    const T ssa24 = -s3[2];
    const T ssa25 = e1[2] - s1[2];
    const T ssa26 = -ssa25 * ssa5 + ssa4 * (e3[2] + ssa24) + u * (e2[2] - s2[2]) - v * (e4[2] + ssa24);
    const T ssa27 = -ssa0 + ssa2 + ssa6;
    const T ssa28 = ssa11 - ssa8 + ssa9;
    const T ssa29 = -ssa13 + ssa14 + ssa16;
    const T ssa30 = e3[0] - e4[0];
    const T ssa31 = e3[1] - e4[1];
    const T ssa32 = e3[2] - e4[2];
    const T ssa33 = ssa12 * ssa31 + ssa17 * ssa32 + ssa30 * ssa7;
    const T ssa34 = -ssa12 * (-e2[1] + s2[1] + ssa22) - ssa17 * (-e2[2] + s2[2] + ssa25) - ssa20 * ssa27 -
                    ssa23 * ssa28 - ssa26 * ssa29 - ssa7 * (-e2[0] + s2[0] + ssa19);
    const T ssa35 = ssa20 * ssa30 * t + ssa23 * ssa31 * t + ssa26 * ssa32 * t + ssa33;
    const T ssa36 = t * (-ssa27 * ssa30 - ssa28 * ssa31 - ssa29 * ssa32);
    *out_f =
        (1.0 / 2.0) * sccd::pow2<T>(ssa12) + (1.0 / 2.0) * sccd::pow2<T>(ssa17) + (1.0 / 2.0) * sccd::pow2<T>(ssa7);
    out_g[0] = ssa12 * ssa23 + ssa17 * ssa26 + ssa20 * ssa7;
    out_g[1] = -ssa12 * ssa28 - ssa17 * ssa29 - ssa27 * ssa7;
    out_g[2] = ssa33 * t;
    out_H[0] = sccd::pow2<T>(ssa20) + sccd::pow2<T>(ssa23) + sccd::pow2<T>(ssa26);
    out_H[1] = ssa34;
    out_H[2] = ssa35;
    out_H[3] = ssa34;
    out_H[4] = sccd::pow2<T>(ssa27) + sccd::pow2<T>(ssa28) + sccd::pow2<T>(ssa29);
    out_H[5] = ssa36;
    out_H[6] = ssa35;
    out_H[7] = ssa36;
    out_H[8] = (sccd::pow2<T>(ssa30) + sccd::pow2<T>(ssa31) + sccd::pow2<T>(ssa32)) * (t * t);
}


#endif