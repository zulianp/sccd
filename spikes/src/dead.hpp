#ifndef SCCD_SPIKES_DEAD_HPP
#define SCCD_SPIKES_DEAD_HPP

// Code that shipped in the library and was called by nothing.
//
// Not spikes in the usual sense -- none of these lost a measurement. They simply
// had no caller in src/, benchmark/, demo/ or the tests, and several sat in
// installed public headers: sccd_objective.hpp was 78% dead by line count, and
// sccd_narrowphase.hpp carried a declaration with no definition anywhere.
//
// Kept rather than deleted because some are complete, correct reference
// implementations worth checking a fast path against -- the scalar_*_range_self
// pair is the readable form of what the vectorised sweep does.
//
// spikes/README.md applies: not built by default, not installed, not covered by
// the correctness gate, deletable without notice. Nothing that ships may include
// this file. Recovering one means moving it back and giving it a caller.
//
// Device-side dead code is not here yet -- see wip/TODO.md.

#include "sccd_aabb.hpp"
#include "sccd_broadphase_sweep.hpp"
#include "sccd_base.hpp"
#include "sccd_math.hpp"

#include <cstddef>
#include <cstdint>

namespace sccd {
    namespace dead {

        // ---- from sccd_objective.hpp (was src/roots.hpp)
        // The objective's gradient, Hessian and fused variants, plus the whole
        // edge-edge family. The splitters use finite differences, so none of these
        // ever acquired a caller.

        //// vf_gradient + vf_hessian
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

        //// vf_all
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

        //// ee_objective, ee_gradient, ee_hessian, ee_objective_dir, ee_all
        // EE
        /* Auto-generated with SymPy (CSE + pow expansion) */


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

        // ---- from sccd_aabb.hpp (was src/vaabb.hpp)
        // vdisjoint is the AVX-512 / AVX2 / NEON one-to-many disjoint kernel,
        // superseded by vaabb_overlap_one_to_many_bits, which returns a bitmask
        // instead of a word per lane. vaabb_broadcast fed it.

        //// vdisjoint
            template <typename T>
            inline static void vdisjoint(const T* const SCCD_RESTRICT aminx,
                                         const T* const SCCD_RESTRICT aminy,
                                         const T* const SCCD_RESTRICT aminz,
                                         const T* const SCCD_RESTRICT amaxx,
                                         const T* const SCCD_RESTRICT amaxy,
                                         const T* const SCCD_RESTRICT amaxz,
                                         const T* const SCCD_RESTRICT bminx,
                                         const T* const SCCD_RESTRICT bminy,
                                         const T* const SCCD_RESTRICT bminz,
                                         const T* const SCCD_RESTRICT bmaxx,
                                         const T* const SCCD_RESTRICT bmaxy,
                                         const T* const SCCD_RESTRICT bmaxz,
                                         uint32_t* SCCD_RESTRICT mask) {
                if constexpr (std::is_same<T, double>::value)  //
                {
        #if defined(__AVX512F__)
                    for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 8) {
                        const __m512d a_minx = _mm512_loadu_pd(aminx + i);
                        const __m512d a_miny = _mm512_loadu_pd(aminy + i);
                        const __m512d a_minz = _mm512_loadu_pd(aminz + i);
                        const __m512d a_maxx = _mm512_loadu_pd(amaxx + i);
                        const __m512d a_maxy = _mm512_loadu_pd(amaxy + i);
                        const __m512d a_maxz = _mm512_loadu_pd(amaxz + i);

                        const __m512d b_minx = _mm512_loadu_pd(bminx + i);
                        const __m512d b_miny = _mm512_loadu_pd(bminy + i);
                        const __m512d b_minz = _mm512_loadu_pd(bminz + i);
                        const __m512d b_maxx = _mm512_loadu_pd(bmaxx + i);
                        const __m512d b_maxy = _mm512_loadu_pd(bmaxy + i);
                        const __m512d b_maxz = _mm512_loadu_pd(bmaxz + i);

                        __mmask8 k =
                            _mm512_cmp_pd_mask(a_minx, b_maxx, _CMP_GT_OQ) | _mm512_cmp_pd_mask(a_miny, b_maxy, _CMP_GT_OQ) |
                            _mm512_cmp_pd_mask(a_minz, b_maxz, _CMP_GT_OQ) | _mm512_cmp_pd_mask(b_minx, a_maxx, _CMP_GT_OQ) |
                            _mm512_cmp_pd_mask(b_miny, a_maxy, _CMP_GT_OQ) | _mm512_cmp_pd_mask(b_minz, a_maxz, _CMP_GT_OQ);

                        __m512i k_as_epi64 = _mm512_movm_epi64(k);
                        __m512i k_01 = _mm512_srli_epi64(k_as_epi64, 63);
                        alignas(64) uint64_t tmp[8];
                        _mm512_storeu_si512((__m512i*)tmp, k_01);
                        for (int lane = 0; lane < 8; ++lane) {
                            mask[i + lane] = static_cast<uint32_t>(tmp[lane]);
                        }
                    }
                    return;
        #elif defined(__AVX2__)
                    for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 4) {
                        const __m256d a_minx = _mm256_loadu_pd(aminx + i);
                        const __m256d a_miny = _mm256_loadu_pd(aminy + i);
                        const __m256d a_minz = _mm256_loadu_pd(aminz + i);
                        const __m256d a_maxx = _mm256_loadu_pd(amaxx + i);
                        const __m256d a_maxy = _mm256_loadu_pd(amaxy + i);
                        const __m256d a_maxz = _mm256_loadu_pd(amaxz + i);

                        const __m256d b_minx = _mm256_loadu_pd(bminx + i);
                        const __m256d b_miny = _mm256_loadu_pd(bminy + i);
                        const __m256d b_minz = _mm256_loadu_pd(bminz + i);
                        const __m256d b_maxx = _mm256_loadu_pd(bmaxx + i);
                        const __m256d b_maxy = _mm256_loadu_pd(bmaxy + i);
                        const __m256d b_maxz = _mm256_loadu_pd(bmaxz + i);

                        __m256d m = _mm256_or_pd(
                            _mm256_or_pd(_mm256_cmp_pd(a_minx, b_maxx, _CMP_GT_OQ), _mm256_cmp_pd(a_miny, b_maxy, _CMP_GT_OQ)),
                            _mm256_cmp_pd(a_minz, b_maxz, _CMP_GT_OQ));
                        m = _mm256_or_pd(
                            m,
                            _mm256_or_pd(_mm256_cmp_pd(b_minx, a_maxx, _CMP_GT_OQ), _mm256_cmp_pd(b_miny, a_maxy, _CMP_GT_OQ)));
                        m = _mm256_or_pd(m, _mm256_cmp_pd(b_minz, a_maxz, _CMP_GT_OQ));

                        const __m256i m_i = _mm256_castpd_si256(m);
                        const __m256i m_01 = _mm256_srli_epi64(m_i, 63);
                        alignas(32) uint64_t tmp[4];
                        _mm256_storeu_si256((__m256i*)tmp, m_01);
                        for (int lane = 0; lane < 4; ++lane) {
                            mask[i + lane] = static_cast<uint32_t>(tmp[lane]);
                        }
                    }
                    return;
        #elif defined(__ARM_NEON) || defined(__ARM_NEON__)
                    for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 2) {
                        const float64x2_t a_minx = vld1q_f64(aminx + i);
                        const float64x2_t a_miny = vld1q_f64(aminy + i);
                        const float64x2_t a_minz = vld1q_f64(aminz + i);
                        const float64x2_t a_maxx = vld1q_f64(amaxx + i);
                        const float64x2_t a_maxy = vld1q_f64(amaxy + i);
                        const float64x2_t a_maxz = vld1q_f64(amaxz + i);

                        const float64x2_t b_minx = vld1q_f64(bminx + i);
                        const float64x2_t b_miny = vld1q_f64(bminy + i);
                        const float64x2_t b_minz = vld1q_f64(bminz + i);
                        const float64x2_t b_maxx = vld1q_f64(bmaxx + i);
                        const float64x2_t b_maxy = vld1q_f64(bmaxy + i);
                        const float64x2_t b_maxz = vld1q_f64(bmaxz + i);

                        uint64x2_t m = vorrq_u64(vorrq_u64(vcgtq_f64(a_minx, b_maxx), vcgtq_f64(a_miny, b_maxy)),
                                                 vcgtq_f64(a_minz, b_maxz));
                        m = vorrq_u64(m, vorrq_u64(vcgtq_f64(b_minx, a_maxx), vcgtq_f64(b_miny, a_maxy)));
                        m = vorrq_u64(m, vcgtq_f64(b_minz, a_maxz));

                        const uint64x2_t m_01 = vshrq_n_u64(m, 63);
                        uint64_t tmp[2];
                        vst1q_u64(tmp, m_01);
                        mask[i] = static_cast<uint32_t>(tmp[0]);
                        mask[i + 1] = static_cast<uint32_t>(tmp[1]);
                    }
                    return;
        #endif
                }
        #pragma omp simd aligned(aminx, aminy, aminz, amaxx, amaxy, amaxz, bminx, bminy, bminz, bmaxx, bmaxy, bmaxz, mask : 64)
                for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i++) {
                    mask[i] = disjoint<T>(aminx[i],
                                          aminy[i],
                                          aminz[i],
                                          amaxx[i],
                                          amaxy[i],
                                          amaxz[i],
                                          bminx[i],
                                          bminy[i],
                                          bminz[i],
                                          bmaxx[i],
                                          bmaxy[i],
                                          bmaxz[i]);
                }
            }

        //// vaabb_broadcast
            /**
             * \brief Broadcast AABB at \p fi into SoA buffers sized for SIMD chunking.
             * \param aabbs SoA AABB arrays [6][...].
             * \param fi Index of the AABB to broadcast.
             * \param A_minx..A_maxz Output arrays of length SCCD_AABB_DISJOINT_CHUNK_SIZE.
             */
            template <typename T>
            inline static void vaabb_broadcast(T** const SCCD_RESTRICT aabbs,
                                               const size_t fi,
                                               T* const SCCD_RESTRICT A_minx,
                                               T* const SCCD_RESTRICT A_miny,
                                               T* const SCCD_RESTRICT A_minz,
                                               T* const SCCD_RESTRICT A_maxx,
                                               T* const SCCD_RESTRICT A_maxy,
                                               T* const SCCD_RESTRICT A_maxz) {
                const T aminx = aabbs[0][fi];
                const T aminy = aabbs[1][fi];
                const T aminz = aabbs[2][fi];
                const T amaxx = aabbs[3][fi];
                const T amaxy = aabbs[4][fi];
                const T amaxz = aabbs[5][fi];
                for (int k = 0; k < SCCD_AABB_DISJOINT_CHUNK_SIZE; ++k) {
                    A_minx[k] = aminx;
                    A_miny[k] = aminy;
                    A_minz[k] = aminz;
                    A_maxx[k] = amaxx;
                    A_maxy[k] = amaxy;
                    A_maxz[k] = amaxz;
                }
            }

        // ---- from sccd_broadphase_sweep.hpp (was src/broadphase.hpp)
        // Scalar reference forms of the count and collect passes, and two mask
        // helpers. Readable statements of what the vectorised sweep does.

        //// remap_idx
            /**
             * \brief Remap indices in-place through a permutation table.
             * \param n Number of entries.
             * \param idx Permutation table mapping old index -> new index.
             * \param remapped Array of indices to update; each entry is replaced by
             * idx[entry].
             */
            template <typename I>
            static void remap_idx(const ptrdiff_t n, const I *const SCCD_RESTRICT idx, I *const SCCD_RESTRICT remapped) {
                sccd::parallel_for_br(0, n, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                    for (ptrdiff_t i = rbegin; i < rend; i++) {
                        remapped[i] = idx[remapped[i]];
                    }
                });
            }

        //// mask_out_shared_self
                /**
                 * \brief Mark lanes where elements share a vertex in self-overlap path.
                 * \tparam N Number of vertices per element.
                 * \param dmask In/out lane mask; set to 1 when a vertex is shared.
                 * \param chunk_len Number of valid lanes.
                 * \param noffset Starting j index (j > i).
                 * \param ev Vertex indices of element i.
                 * \param idx Mapping from sorted position to element id.
                 * \param elements SoA vertex arrays.
                 * \param stride Stride between elements in the arrays.
                 */
                template <int N, typename I>
                static inline void mask_out_shared_self(uint32_t *const SCCD_RESTRICT dmask,
                                                        const ptrdiff_t chunk_len,
                                                        const ptrdiff_t noffset,
                                                        const I (&ev)[N],
                                                        const I *const SCCD_RESTRICT idx,
                                                        I **const SCCD_RESTRICT elements,
                                                        const ptrdiff_t stride) {
                    for (ptrdiff_t lane = 0; lane < chunk_len; ++lane) {
                        if (dmask[lane]) continue;

                        const ptrdiff_t j = noffset + lane;
                        const I jidx = idx[j];
                        I sev[N];
                        detail::load_ev<N>(elements, jidx, stride, sev);
                        if (detail::shares_vertex<N, N>(ev, sev)) {
                            dmask[lane] = 1;
                        }
                    }
                }

        //// scalar_count_range_self
                /**
                 * \brief Scalar reference: count self-overlaps in [begin,end) for element
                 * i.
                 * \return Number of non-disjoint, non-shared-vertex candidates with j>i.
                 */
                template <int N, typename T, typename I>
                static inline ptrdiff_t scalar_count_range_self(T **const SCCD_RESTRICT aabbs,
                                                                const ptrdiff_t fi,
                                                                I **const SCCD_RESTRICT elements,
                                                                const I *const SCCD_RESTRICT idx,
                                                                const ptrdiff_t stride,
                                                                const I (&ev)[N],
                                                                const ptrdiff_t begin,
                                                                const ptrdiff_t end) {
                    ptrdiff_t count = 0;
                    const T aminx = aabbs[0][fi];
                    const T aminy = aabbs[1][fi];
                    const T aminz = aabbs[2][fi];
                    const T amaxx = aabbs[3][fi];
                    const T amaxy = aabbs[4][fi];
                    const T amaxz = aabbs[5][fi];
                    for (ptrdiff_t j = begin; j < end; ++j) {
                        if (disjoint(aminx,
                                     aminy,
                                     aminz,
                                     amaxx,
                                     amaxy,
                                     amaxz,
                                     aabbs[0][j],
                                     aabbs[1][j],
                                     aabbs[2][j],
                                     aabbs[3][j],
                                     aabbs[4][j],
                                     aabbs[5][j])) {
                            continue;
                        }
                        const I jidx = idx[j];
                        I sev[N];
                        detail::load_ev<N>(elements, jidx, stride, sev);
                        const bool share = detail::shares_vertex<N, N>(ev, sev);
                        count += share ? 0 : 1;
                    }
                    return count;
                }

        //// scalar_collect_range_self
                /**
                 * \brief Scalar reference: collect self-overlaps in [begin,end) for element
                 * i.
                 * \return Number of pairs written to outputs, with (min(idxi,jidx),
                 * max(...)).
                 */
                template <int N, typename T, typename I>
                static inline ptrdiff_t scalar_collect_range_self(T **const SCCD_RESTRICT aabbs,
                                                                  const ptrdiff_t fi,
                                                                  const I idxi,
                                                                  I **const SCCD_RESTRICT elements,
                                                                  const I *const SCCD_RESTRICT idx,
                                                                  const ptrdiff_t stride,
                                                                  const I (&ev)[N],
                                                                  const ptrdiff_t begin,
                                                                  const ptrdiff_t end,
                                                                  I *const SCCD_RESTRICT first_out,
                                                                  I *const SCCD_RESTRICT second_out) {
                    ptrdiff_t count = 0;
                    const T aminx = aabbs[0][fi];
                    const T aminy = aabbs[1][fi];
                    const T aminz = aabbs[2][fi];
                    const T amaxx = aabbs[3][fi];
                    const T amaxy = aabbs[4][fi];
                    const T amaxz = aabbs[5][fi];
                    for (ptrdiff_t j = begin; j < end; ++j) {
                        if (disjoint(aminx,
                                     aminy,
                                     aminz,
                                     amaxx,
                                     amaxy,
                                     amaxz,
                                     aabbs[0][j],
                                     aabbs[1][j],
                                     aabbs[2][j],
                                     aabbs[3][j],
                                     aabbs[4][j],
                                     aabbs[5][j])) {
                            continue;
                        }
                        const I jidx = idx[j];
                        I sev[N];
                        detail::load_ev<N>(elements, jidx, stride, sev);
                        if (!detail::shares_vertex<N, N>(ev, sev)) {
                            first_out[count] = sccd::min(idxi, jidx);
                            second_out[count] = sccd::max(idxi, jidx);
                            count += 1;
                        }
                    }
                    return count;
                }

        // ---- from sccd_math.hpp (was src/smath.hpp)
        // Never called.

        //// array_min
        template <typename T> static inline T array_min(const int n, const T* arr) {
          T min = arr[0];
          for (int i = 1; i < n; i++) {
            if (arr[i] < min) {
              min = arr[i];
            }
          }
          return min;
        }

        //// array_max
        template <typename T> static inline T array_max(const int n, const T* arr) {
          T max = arr[0];
          for (int i = 1; i < n; i++) {
            if (arr[i] > max) {
              max = arr[i];
            }
          }
          return max;
        }

        // ---- from sccd_narrowphase.hpp (was src/narrowphase.hpp)
        // A declaration with no definition anywhere in the repository. It could not
        // have been called -- linking a call would have failed.

        //// narrow_phase_newton_pass_vf (declaration only)

            template <int nxe, typename T, typename I>
            int narrow_phase_newton_pass_vf(const size_t noverlaps,
                                            const I* const SCCD_RESTRICT voveralp,
                                            const I* const SCCD_RESTRICT foveralp,
                                            T** const SCCD_RESTRICT v0,
                                            T** const SCCD_RESTRICT v1,
                                            const size_t face_stride,
                                            I** const SCCD_RESTRICT faces,
                                            const T max_toi,
                                            T* const SCCD_RESTRICT toi,
                                            const T tol,
                                            const int toi_stride = 0);


    }  // namespace dead
}  // namespace sccd

#endif  // SCCD_SPIKES_DEAD_HPP
