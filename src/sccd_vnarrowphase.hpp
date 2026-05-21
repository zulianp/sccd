#ifndef SCCD_VNARROWPHASE_HPP
#define SCCD_VNARROWPHASE_HPP

#include "sccd_base.hpp"
#include "smath.hpp"
#include "snumtol.hpp"
#include "sparallel.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <vector>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#ifndef SCCD_VNARROWPHASE_ADAPTIVE_SPLIT
#define SCCD_VNARROWPHASE_ADAPTIVE_SPLIT 1
#endif

namespace sccd {

    template <typename T>
    static inline void atomic_min_relaxed(std::atomic<T>& target, const T value) {
        T current = target.load(std::memory_order_relaxed);
        while (value < current && !target.compare_exchange_weak(current, value, std::memory_order_relaxed)) {
        }
    }

    template <typename T, int VSIZE>
    struct VPoint {
        T x[VSIZE];
        T y[VSIZE];
        T z[VSIZE];
    };

    template <typename T, int VSIZE>
    struct VQuery {
        VPoint<T, VSIZE> s0;
        VPoint<T, VSIZE> s1;
        VPoint<T, VSIZE> s2;
        VPoint<T, VSIZE> s3;

        VPoint<T, VSIZE> e0;
        VPoint<T, VSIZE> e1;
        VPoint<T, VSIZE> e2;
        VPoint<T, VSIZE> e3;
    };

    template <typename T, int VSIZE>
    struct VInterval {
        T lower[VSIZE];
        T upper[VSIZE];
    };

    template <typename T, int VSIZE>
    struct VCodomain {
        VInterval<T, VSIZE> xyz[3];
    };

    template <typename T, int VSIZE>
    struct VTolerances {
        T axis[3][VSIZE];
    };

    template <typename T, typename I, int VSIZE>
    struct VDomain {
        I query_index[VSIZE];
        int depth[VSIZE];
        uint8_t active[VSIZE];
        VInterval<T, VSIZE> tuv[3];
    };

    template <typename T, typename I>
    struct VActiveBox {
        I query_index;
        int depth;
        T lower[3];
        T upper[3];
    };

    template <typename T, typename I, int VSIZE>
    static inline void deactivate_lane(VDomain<T, I, VSIZE>& domain, const ptrdiff_t i) {
        domain.query_index[i] = I(0);
        domain.depth[i] = 0;
        domain.active[i] = 0;
        for (int d = 0; d < 3; ++d) {
            domain.tuv[d].lower[i] = T(0);
            domain.tuv[d].upper[i] = T(0);
        }
    }

    template <typename Q, typename T, typename I, int VSIZE>
    static void load_query(const ptrdiff_t block_index,
                           const ptrdiff_t block_size,
                           const I* const SCCD_RESTRICT voveralp,
                           const I* const SCCD_RESTRICT foveralp,
                           T** const SCCD_RESTRICT v0,
                           T** const SCCD_RESTRICT v1,
                           const size_t face_stride,
                           I** const SCCD_RESTRICT faces,
                           VQuery<Q, VSIZE>& query) {
        const ptrdiff_t query_begin = block_index * VSIZE;

        const T* const SCCD_RESTRICT v0x = v0[0];
        const T* const SCCD_RESTRICT v0y = v0[1];
        const T* const SCCD_RESTRICT v0z = v0[2];
        const T* const SCCD_RESTRICT v1x = v1[0];
        const T* const SCCD_RESTRICT v1y = v1[1];
        const T* const SCCD_RESTRICT v1z = v1[2];

        const I* const SCCD_RESTRICT faces0 = faces[0];
        const I* const SCCD_RESTRICT faces1 = faces[1];
        const I* const SCCD_RESTRICT faces2 = faces[2];

        for (ptrdiff_t i = 0; i < block_size; ++i) {
            const ptrdiff_t qi = query_begin + i;
            const I vi = voveralp[qi];
            const I fi = foveralp[qi];
            const size_t face_offset = static_cast<size_t>(fi) * face_stride;

            const I n0 = faces0[face_offset];
            const I n1 = faces1[face_offset];
            const I n2 = faces2[face_offset];

            query.s0.x[i] = v0x[vi];
            query.s0.y[i] = v0y[vi];
            query.s0.z[i] = v0z[vi];
            query.e0.x[i] = v1x[vi];
            query.e0.y[i] = v1y[vi];
            query.e0.z[i] = v1z[vi];

            query.s1.x[i] = v0x[n0];
            query.s1.y[i] = v0y[n0];
            query.s1.z[i] = v0z[n0];
            query.e1.x[i] = v1x[n0];
            query.e1.y[i] = v1y[n0];
            query.e1.z[i] = v1z[n0];

            query.s2.x[i] = v0x[n1];
            query.s2.y[i] = v0y[n1];
            query.s2.z[i] = v0z[n1];
            query.e2.x[i] = v1x[n1];
            query.e2.y[i] = v1y[n1];
            query.e2.z[i] = v1z[n1];

            query.s3.x[i] = v0x[n2];
            query.s3.y[i] = v0y[n2];
            query.s3.z[i] = v0z[n2];
            query.e3.x[i] = v1x[n2];
            query.e3.y[i] = v1y[n2];
            query.e3.z[i] = v1z[n2];
        }

        for (ptrdiff_t i = block_size; i < VSIZE; ++i) {
            query.s0.x[i] = T(0);
            query.s0.y[i] = T(0);
            query.s0.z[i] = T(0);
            query.e0.x[i] = T(0);
            query.e0.y[i] = T(0);
            query.e0.z[i] = T(0);

            query.s1.x[i] = T(1);
            query.s1.y[i] = T(0);
            query.s1.z[i] = T(0);
            query.e1.x[i] = T(1);
            query.e1.y[i] = T(0);
            query.e1.z[i] = T(0);

            query.s2.x[i] = T(1);
            query.s2.y[i] = T(0);
            query.s2.z[i] = T(0);
            query.e2.x[i] = T(1);
            query.e2.y[i] = T(0);
            query.e2.z[i] = T(0);

            query.s3.x[i] = T(1);
            query.s3.y[i] = T(0);
            query.s3.z[i] = T(0);
            query.e3.x[i] = T(1);
            query.e3.y[i] = T(0);
            query.e3.z[i] = T(0);
        }
    }

    template <typename T, typename I, int VSIZE>
    static void init_domain(const ptrdiff_t block_size, const T max_toi, VDomain<T, I, VSIZE>& domain) {
        const T tupper = sccd::min<T>(T(1), max_toi);
        for (ptrdiff_t i = 0; i < block_size; ++i) {
            domain.query_index[i] = static_cast<I>(i);
            domain.depth[i] = 0;
            domain.active[i] = tupper > T(0);

            domain.tuv[0].lower[i] = T(0);
            domain.tuv[0].upper[i] = tupper;
            domain.tuv[1].lower[i] = T(0);
            domain.tuv[1].upper[i] = T(1);
            domain.tuv[2].lower[i] = T(0);
            domain.tuv[2].upper[i] = T(1);
        }

        for (ptrdiff_t i = block_size; i < VSIZE; ++i) {
            deactivate_lane(domain, i);
        }
    }

    template <typename T, int VSIZE>
    static void compute_tolerances(const T tol, const VQuery<T, VSIZE>& query, VTolerances<T, VSIZE>& tols) {
        for (int i = 0; i < VSIZE; ++i) {
            compute_face_vertex_tolerance_soa<T>(tol,
                                                 query.s0.x[i],
                                                 query.s0.y[i],
                                                 query.s0.z[i],
                                                 query.s1.x[i],
                                                 query.s1.y[i],
                                                 query.s1.z[i],
                                                 query.s2.x[i],
                                                 query.s2.y[i],
                                                 query.s2.z[i],
                                                 query.s3.x[i],
                                                 query.s3.y[i],
                                                 query.s3.z[i],
                                                 query.e0.x[i],
                                                 query.e0.y[i],
                                                 query.e0.z[i],
                                                 query.e1.x[i],
                                                 query.e1.y[i],
                                                 query.e1.z[i],
                                                 query.e2.x[i],
                                                 query.e2.y[i],
                                                 query.e2.z[i],
                                                 query.e3.x[i],
                                                 query.e3.y[i],
                                                 query.e3.z[i],
                                                 &tols.axis[0][i],
                                                 &tols.axis[1][i],
                                                 &tols.axis[2][i]);
        }
    }

    template <typename T>
    static inline void diff_vf_component(const T s0,
                                         const T s1,
                                         const T s2,
                                         const T s3,
                                         const T e0,
                                         const T e1,
                                         const T e2,
                                         const T e3,
                                         const T t,
                                         const T u,
                                         const T v,
                                         T& out) {
        const T omt = T(1) - t;
        const T o = T(1) - u - v;
        out = (omt * s0 + t * e0) - (omt * (o * s1 + u * s2 + v * s3) + t * (o * e1 + u * e2 + v * e3));
    }

    template <typename T, typename I, int VSIZE>
    static void compute_codomain_scalar(const VDomain<T, I, VSIZE>& domain,
                                        const VQuery<T, VSIZE>& query,
                                        VCodomain<T, VSIZE>& codomain) {
        for (int d = 0; d < 3; ++d) {
            for (int i = 0; i < VSIZE; ++i) {
                codomain.xyz[d].lower[i] = std::numeric_limits<T>::max();
                codomain.xyz[d].upper[i] = std::numeric_limits<T>::lowest();
            }
        }

        for (int mask = 0; mask < 8; ++mask) {
            const bool mt = (mask & 1) != 0;
            const bool mu = (mask & 2) != 0;
            const bool mv = (mask & 4) != 0;
#pragma omp simd
            for (int i = 0; i < VSIZE; ++i) {
                const int q = static_cast<int>(domain.query_index[i]);
                const T t = mt ? domain.tuv[0].upper[i] : domain.tuv[0].lower[i];
                const T u = mu ? domain.tuv[1].upper[i] : domain.tuv[1].lower[i];
                const T v = mv ? domain.tuv[2].upper[i] : domain.tuv[2].lower[i];
                T f;

                diff_vf_component<T>(query.s0.x[q],
                                     query.s1.x[q],
                                     query.s2.x[q],
                                     query.s3.x[q],
                                     query.e0.x[q],
                                     query.e1.x[q],
                                     query.e2.x[q],
                                     query.e3.x[q],
                                     t,
                                     u,
                                     v,
                                     f);
                codomain.xyz[0].lower[i] = sccd::min<T>(codomain.xyz[0].lower[i], f);
                codomain.xyz[0].upper[i] = sccd::max<T>(codomain.xyz[0].upper[i], f);

                diff_vf_component<T>(query.s0.y[q],
                                     query.s1.y[q],
                                     query.s2.y[q],
                                     query.s3.y[q],
                                     query.e0.y[q],
                                     query.e1.y[q],
                                     query.e2.y[q],
                                     query.e3.y[q],
                                     t,
                                     u,
                                     v,
                                     f);
                codomain.xyz[1].lower[i] = sccd::min<T>(codomain.xyz[1].lower[i], f);
                codomain.xyz[1].upper[i] = sccd::max<T>(codomain.xyz[1].upper[i], f);

                diff_vf_component<T>(query.s0.z[q],
                                     query.s1.z[q],
                                     query.s2.z[q],
                                     query.s3.z[q],
                                     query.e0.z[q],
                                     query.e1.z[q],
                                     query.e2.z[q],
                                     query.e3.z[q],
                                     t,
                                     u,
                                     v,
                                     f);
                codomain.xyz[2].lower[i] = sccd::min<T>(codomain.xyz[2].lower[i], f);
                codomain.xyz[2].upper[i] = sccd::max<T>(codomain.xyz[2].upper[i], f);
            }
        }
    }

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    static inline float64x2_t load2_indexed_f64(const double* const SCCD_RESTRICT a, const int q0, const int q1) {
        return vsetq_lane_f64(a[q1], vdupq_n_f64(a[q0]), 1);
    }

    static inline float64x2_t diff_vf_component_neon_indexed(const double* const SCCD_RESTRICT s0,
                                                             const double* const SCCD_RESTRICT s1,
                                                             const double* const SCCD_RESTRICT s2,
                                                             const double* const SCCD_RESTRICT s3,
                                                             const double* const SCCD_RESTRICT e0,
                                                             const double* const SCCD_RESTRICT e1,
                                                             const double* const SCCD_RESTRICT e2,
                                                             const double* const SCCD_RESTRICT e3,
                                                             const int q0,
                                                             const int q1,
                                                             const float64x2_t t,
                                                             const float64x2_t u,
                                                             const float64x2_t v) {
        const float64x2_t one = vdupq_n_f64(1.0);
        const float64x2_t omt = vsubq_f64(one, t);
        const float64x2_t o = vsubq_f64(vsubq_f64(one, u), v);

        const float64x2_t vs0 = load2_indexed_f64(s0, q0, q1);
        const float64x2_t vs1 = load2_indexed_f64(s1, q0, q1);
        const float64x2_t vs2 = load2_indexed_f64(s2, q0, q1);
        const float64x2_t vs3 = load2_indexed_f64(s3, q0, q1);
        const float64x2_t ve0 = load2_indexed_f64(e0, q0, q1);
        const float64x2_t ve1 = load2_indexed_f64(e1, q0, q1);
        const float64x2_t ve2 = load2_indexed_f64(e2, q0, q1);
        const float64x2_t ve3 = load2_indexed_f64(e3, q0, q1);

        const float64x2_t vertex = vfmaq_f64(vmulq_f64(omt, vs0), t, ve0);
        float64x2_t sface = vmulq_f64(o, vs1);
        sface = vfmaq_f64(sface, u, vs2);
        sface = vfmaq_f64(sface, v, vs3);
        float64x2_t eface = vmulq_f64(o, ve1);
        eface = vfmaq_f64(eface, u, ve2);
        eface = vfmaq_f64(eface, v, ve3);
        return vsubq_f64(vertex, vfmaq_f64(vmulq_f64(omt, sface), t, eface));
    }

    template <typename I, int VSIZE>
    static void compute_codomain_neon(const VDomain<double, I, VSIZE>& domain,
                                      const VQuery<double, VSIZE>& query,
                                      VCodomain<double, VSIZE>& codomain) {
        static_assert((VSIZE % 2) == 0, "double NEON codomain kernel expects an even vector size");
        for (int d = 0; d < 3; ++d) {
            for (int i = 0; i < VSIZE; i += 2) {
                vst1q_f64(codomain.xyz[d].lower + i, vdupq_n_f64(std::numeric_limits<double>::max()));
                vst1q_f64(codomain.xyz[d].upper + i, vdupq_n_f64(std::numeric_limits<double>::lowest()));
            }
        }

        for (int mask = 0; mask < 8; ++mask) {
            const double* const tptr = (mask & 1) ? domain.tuv[0].upper : domain.tuv[0].lower;
            const double* const uptr = (mask & 2) ? domain.tuv[1].upper : domain.tuv[1].lower;
            const double* const vptr = (mask & 4) ? domain.tuv[2].upper : domain.tuv[2].lower;
            for (int i = 0; i < VSIZE; i += 2) {
                const float64x2_t t = vld1q_f64(tptr + i);
                const float64x2_t u = vld1q_f64(uptr + i);
                const float64x2_t v = vld1q_f64(vptr + i);
                const int q0 = static_cast<int>(domain.query_index[i]);
                const int q1 = static_cast<int>(domain.query_index[i + 1]);

                float64x2_t f = diff_vf_component_neon_indexed(query.s0.x,
                                                               query.s1.x,
                                                               query.s2.x,
                                                               query.s3.x,
                                                               query.e0.x,
                                                               query.e1.x,
                                                               query.e2.x,
                                                               query.e3.x,
                                                               q0,
                                                               q1,
                                                               t,
                                                               u,
                                                               v);
                vst1q_f64(codomain.xyz[0].lower + i, vminq_f64(vld1q_f64(codomain.xyz[0].lower + i), f));
                vst1q_f64(codomain.xyz[0].upper + i, vmaxq_f64(vld1q_f64(codomain.xyz[0].upper + i), f));

                f = diff_vf_component_neon_indexed(query.s0.y,
                                                   query.s1.y,
                                                   query.s2.y,
                                                   query.s3.y,
                                                   query.e0.y,
                                                   query.e1.y,
                                                   query.e2.y,
                                                   query.e3.y,
                                                   q0,
                                                   q1,
                                                   t,
                                                   u,
                                                   v);
                vst1q_f64(codomain.xyz[1].lower + i, vminq_f64(vld1q_f64(codomain.xyz[1].lower + i), f));
                vst1q_f64(codomain.xyz[1].upper + i, vmaxq_f64(vld1q_f64(codomain.xyz[1].upper + i), f));

                f = diff_vf_component_neon_indexed(query.s0.z,
                                                   query.s1.z,
                                                   query.s2.z,
                                                   query.s3.z,
                                                   query.e0.z,
                                                   query.e1.z,
                                                   query.e2.z,
                                                   query.e3.z,
                                                   q0,
                                                   q1,
                                                   t,
                                                   u,
                                                   v);
                vst1q_f64(codomain.xyz[2].lower + i, vminq_f64(vld1q_f64(codomain.xyz[2].lower + i), f));
                vst1q_f64(codomain.xyz[2].upper + i, vmaxq_f64(vld1q_f64(codomain.xyz[2].upper + i), f));
            }
        }
    }
#endif

    template <typename T, typename I, int VSIZE>
    static void compute_codomain(const VDomain<T, I, VSIZE>& domain,
                                 const VQuery<T, VSIZE>& query,
                                 VCodomain<T, VSIZE>& codomain) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        if constexpr (std::is_same<T, double>::value) {
            compute_codomain_neon<I, VSIZE>(domain, query, codomain);
            return;
        }
#endif
        compute_codomain_scalar<T, I, VSIZE>(domain, query, codomain);
    }

    template <typename T, typename I, int VSIZE>
    static void compute_masks_scalar(const VDomain<T, I, VSIZE>& domain,
                                     const VCodomain<T, VSIZE>& codomain,
                                     const VTolerances<T, VSIZE>& tols,
                                     const T tol,
                                     uint8_t* const SCCD_RESTRICT contains_origin,
                                     uint8_t* const SCCD_RESTRICT accepted) {
        for (int i = 0; i < VSIZE; ++i) {
            const int q = static_cast<int>(domain.query_index[i]);
            bool contains = domain.active[i] != 0;
            bool box_in_eps = true;
            bool widths_within_tol = true;

            for (int d = 0; d < 3; ++d) {
                const T fmin = codomain.xyz[d].lower[i];
                const T fmax = codomain.xyz[d].upper[i];
                contains = contains && (fmin <= tol) && (fmax >= -tol);
                box_in_eps = box_in_eps && (fmin >= -tol) && (fmax <= tol);
                widths_within_tol =
                    widths_within_tol && ((domain.tuv[d].upper[i] - domain.tuv[d].lower[i]) <= tols.axis[d][q]);
            }

            const bool positive_time = domain.tuv[0].lower[i] > T(0);
            contains_origin[i] = contains;
            accepted[i] = contains && positive_time && (box_in_eps || widths_within_tol);
        }
    }

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    template <typename I, int VSIZE>
    static void compute_masks_neon(const VDomain<double, I, VSIZE>& domain,
                                   const VCodomain<double, VSIZE>& codomain,
                                   const VTolerances<double, VSIZE>& tols,
                                   const double tol,
                                   uint8_t* const SCCD_RESTRICT contains_origin,
                                   uint8_t* const SCCD_RESTRICT accepted) {
        static_assert((VSIZE % 2) == 0, "double NEON mask kernel expects an even vector size");
        const float64x2_t vtol = vdupq_n_f64(tol);
        const float64x2_t vntol = vdupq_n_f64(-tol);
        const float64x2_t vzero = vdupq_n_f64(0.0);

        for (int i = 0; i < VSIZE; i += 2) {
            uint64x2_t contains = vdupq_n_u64(~uint64_t(0));
            uint64x2_t box_in_eps = vdupq_n_u64(~uint64_t(0));
            uint64x2_t widths_within_tol = vdupq_n_u64(~uint64_t(0));
            const int q0 = static_cast<int>(domain.query_index[i]);
            const int q1 = static_cast<int>(domain.query_index[i + 1]);

            for (int d = 0; d < 3; ++d) {
                const float64x2_t fmin = vld1q_f64(codomain.xyz[d].lower + i);
                const float64x2_t fmax = vld1q_f64(codomain.xyz[d].upper + i);
                const float64x2_t width =
                    vsubq_f64(vld1q_f64(domain.tuv[d].upper + i), vld1q_f64(domain.tuv[d].lower + i));
                contains = vandq_u64(contains, vandq_u64(vcleq_f64(fmin, vtol), vcgeq_f64(fmax, vntol)));
                box_in_eps = vandq_u64(box_in_eps, vandq_u64(vcgeq_f64(fmin, vntol), vcleq_f64(fmax, vtol)));
                widths_within_tol =
                    vandq_u64(widths_within_tol, vcleq_f64(width, load2_indexed_f64(tols.axis[d], q0, q1)));
            }

            const uint64x2_t positive_time = vcgtq_f64(vld1q_f64(domain.tuv[0].lower + i), vzero);
            const uint64x2_t accept =
                vandq_u64(vandq_u64(contains, positive_time), vorrq_u64(box_in_eps, widths_within_tol));
            uint64_t co_tmp[2];
            uint64_t ac_tmp[2];
            vst1q_u64(co_tmp, contains);
            vst1q_u64(ac_tmp, accept);
            contains_origin[i] = domain.active[i] && co_tmp[0];
            contains_origin[i + 1] = domain.active[i + 1] && co_tmp[1];
            accepted[i] = domain.active[i] && ac_tmp[0];
            accepted[i + 1] = domain.active[i + 1] && ac_tmp[1];
        }
    }
#endif

    template <typename T, typename I, int VSIZE>
    static void compute_masks(const VDomain<T, I, VSIZE>& domain,
                              const VCodomain<T, VSIZE>& codomain,
                              const VTolerances<T, VSIZE>& tols,
                              const T tol,
                              uint8_t* const SCCD_RESTRICT contains_origin,
                              uint8_t* const SCCD_RESTRICT accepted) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        if constexpr (std::is_same<T, double>::value) {
            compute_masks_neon<I, VSIZE>(domain, codomain, tols, tol, contains_origin, accepted);
            return;
        }
#endif
        compute_masks_scalar<T, I, VSIZE>(domain, codomain, tols, tol, contains_origin, accepted);
    }

    template <typename T, typename I, int VSIZE>
    static void detect_longest_axis_to_split_scalar(const VDomain<T, I, VSIZE>& domain,
                                                    int* const SCCD_RESTRICT longest_axis_to_split) {
        for (int i = 0; i < VSIZE; ++i) {
            const T dt = domain.tuv[0].upper[i] - domain.tuv[0].lower[i];
            const T du = domain.tuv[1].upper[i] - domain.tuv[1].lower[i];
            const T dv = domain.tuv[2].upper[i] - domain.tuv[2].lower[i];
            if (du > dt && du >= dv) {
                longest_axis_to_split[i] = 1;
            } else if (dv > dt && dv > du) {
                longest_axis_to_split[i] = 2;
            } else {
                longest_axis_to_split[i] = 0;
            }
        }
    }

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    template <typename I, int VSIZE>
    static void detect_longest_axis_to_split_neon(const VDomain<double, I, VSIZE>& domain,
                                                  int* const SCCD_RESTRICT longest_axis_to_split) {
        static_assert((VSIZE % 2) == 0, "double NEON axis kernel expects an even vector size");
        for (int i = 0; i < VSIZE; i += 2) {
            const float64x2_t dt = vsubq_f64(vld1q_f64(domain.tuv[0].upper + i), vld1q_f64(domain.tuv[0].lower + i));
            const float64x2_t du = vsubq_f64(vld1q_f64(domain.tuv[1].upper + i), vld1q_f64(domain.tuv[1].lower + i));
            const float64x2_t dv = vsubq_f64(vld1q_f64(domain.tuv[2].upper + i), vld1q_f64(domain.tuv[2].lower + i));

            const uint64x2_t u_axis = vandq_u64(vcgtq_f64(du, dt), vcgeq_f64(du, dv));
            const uint64x2_t v_axis = vandq_u64(vcgtq_f64(dv, dt), vcgtq_f64(dv, du));
            uint64_t u_tmp[2];
            uint64_t v_tmp[2];
            vst1q_u64(u_tmp, u_axis);
            vst1q_u64(v_tmp, v_axis);
            longest_axis_to_split[i] = u_tmp[0] ? 1 : (v_tmp[0] ? 2 : 0);
            longest_axis_to_split[i + 1] = u_tmp[1] ? 1 : (v_tmp[1] ? 2 : 0);
        }
    }
#endif

    template <typename T, typename I, int VSIZE>
    static void detect_longest_axis_to_split(const VDomain<T, I, VSIZE>& domain,
                                             int* const SCCD_RESTRICT longest_axis_to_split) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        if constexpr (std::is_same<T, double>::value) {
            detect_longest_axis_to_split_neon<I, VSIZE>(domain, longest_axis_to_split);
            return;
        }
#endif
        detect_longest_axis_to_split_scalar<T, I, VSIZE>(domain, longest_axis_to_split);
    }

    template <typename T, typename I, int VSIZE>
    static inline bool lane_valid(const VDomain<T, I, VSIZE>& domain,
                                  const int i,
                                  const T toi,
                                  const VTolerances<T, VSIZE>&) {
        return domain.active[i] && domain.tuv[0].lower[i] < toi &&
               domain.tuv[1].lower[i] + domain.tuv[2].lower[i] <= T(1);
    }

    template <int SplitDim, int N, typename T, int VSIZE>
    static inline void normal_equation_axis_splitters_vf_lane(const VQuery<T, VSIZE>& query,
                                                              const int q,
                                                              const T tlower,
                                                              const T tupper,
                                                              const T ulower,
                                                              const T uupper,
                                                              const T vlower,
                                                              const T vupper,
                                                              T* const SCCD_RESTRICT splitters) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);

        const T lower[3] = {tlower, ulower, vlower};
        const T upper[3] = {tupper, uupper, vupper};
        const T sv[3] = {query.s0.x[q], query.s0.y[q], query.s0.z[q]};
        const T s1[3] = {query.s1.x[q], query.s1.y[q], query.s1.z[q]};
        const T s2[3] = {query.s2.x[q], query.s2.y[q], query.s2.z[q]};
        const T s3[3] = {query.s3.x[q], query.s3.y[q], query.s3.z[q]};
        const T ev[3] = {query.e0.x[q], query.e0.y[q], query.e0.z[q]};
        const T e1[3] = {query.e1.x[q], query.e1.y[q], query.e1.z[q]};
        const T e2[3] = {query.e2.x[q], query.e2.y[q], query.e2.z[q]};
        const T e3[3] = {query.e3.x[q], query.e3.y[q], query.e3.z[q]};

        const T lo = lower[SplitDim];
        const T hi = upper[SplitDim];
        const T h = (hi - lo) / T(N + 1);
        const T radius = h * T(0.45);
        const T mid_t = (lower[0] + upper[0]) * T(0.5);
        const T mid_u = (lower[1] + upper[1]) * T(0.5);
        const T mid_v = (lower[2] + upper[2]) * T(0.5);
        const T eps = std::numeric_limits<T>::epsilon();

        T F_base[3];
        T J_axis[3];
        T H_axis = T(0);

        const T base_t = SplitDim == 0 ? T(0) : mid_t;
        const T base_u = SplitDim == 1 ? T(0) : mid_u;
        const T base_v = SplitDim == 2 ? T(0) : mid_v;
        const T base_omt = T(1) - base_t;
        const T base_o = T(1) - base_u - base_v;

        for (int d = 0; d < 3; ++d) {
            const T vertex = base_omt * sv[d] + base_t * ev[d];
            const T face = base_omt * (base_o * s1[d] + base_u * s2[d] + base_v * s3[d]) +
                           base_t * (base_o * e1[d] + base_u * e2[d] + base_v * e3[d]);
            F_base[d] = vertex - face;

            if constexpr (SplitDim == 0) {
                const T o = T(1) - mid_u - mid_v;
                J_axis[d] = (ev[d] - sv[d]) - (o * (e1[d] - s1[d]) + mid_u * (e2[d] - s2[d]) + mid_v * (e3[d] - s3[d]));
            } else if constexpr (SplitDim == 1) {
                J_axis[d] = -((T(1) - mid_t) * (s2[d] - s1[d]) + mid_t * (e2[d] - e1[d]));
            } else {
                J_axis[d] = -((T(1) - mid_t) * (s3[d] - s1[d]) + mid_t * (e3[d] - e1[d]));
            }
            H_axis += J_axis[d] * J_axis[d];
        }

        const T step_scale = H_axis > eps ? T(1) / H_axis : T(0.00001);
#pragma omp simd aligned(splitters)
        for (int i = 0; i < N; ++i) {
            const T x0 = lo + h * T(i + 1);
            auto xmin = sccd::max<T>(lo, x0 - radius);
            auto xmax = sccd::min<T>(hi, x0 + radius);

            T g = T(0);
            for (int d = 0; d < 3; ++d) {
                const T J = J_axis[d];
                g += (F_base[d] + x0 * J) * J;
            }

            const T step = g * step_scale;
            splitters[i] = sccd::min<T>(xmax, sccd::max<T>(xmin, x0 - step));
        }
    }

    template <int N, typename T, int VSIZE>
    static inline void adaptive_splitters_vf(const VQuery<T, VSIZE>& query,
                                             const int q,
                                             const int split_axis,
                                             const T tlower,
                                             const T tupper,
                                             const T ulower,
                                             const T uupper,
                                             const T vlower,
                                             const T vupper,
                                             T* const SCCD_RESTRICT splitters) {
        if (split_axis == 0) {
            normal_equation_axis_splitters_vf_lane<0, N, T, VSIZE>(
                query, q, tlower, tupper, ulower, uupper, vlower, vupper, splitters);
        } else if (split_axis == 1) {
            normal_equation_axis_splitters_vf_lane<1, N, T, VSIZE>(
                query, q, tlower, tupper, ulower, uupper, vlower, vupper, splitters);
        } else {
            normal_equation_axis_splitters_vf_lane<2, N, T, VSIZE>(
                query, q, tlower, tupper, ulower, uupper, vlower, vupper, splitters);
        }
    }

    template <typename T>
    static inline bool codomain_acceptance_grid_vf(const T fmin[3],
                                                   const T fmax[3],
                                                   const T tol,
                                                   const T tols[3],
                                                   bool& accept) {
        bool contains_zero = true;
        bool smaller_than_axis_tol = true;
        bool inside_epsilon_box = true;
        bool last_axis_smaller_than_scalar_tol = false;
        bool degenerate_interval = true;

        for (int d = 0; d < 3; ++d) {
            const T interval_width = fmax[d] - fmin[d];
            contains_zero = contains_zero && (fmin[d] <= tol) && (fmax[d] >= -tol);
            smaller_than_axis_tol = smaller_than_axis_tol && (interval_width <= tols[d]);
            inside_epsilon_box = inside_epsilon_box && !((fmin[d] < tol) || (fmax[d] > -tol));
            last_axis_smaller_than_scalar_tol = interval_width < tol;
            degenerate_interval = degenerate_interval && (fmin[d] >= fmax[d]);
        }

        accept = contains_zero && (smaller_than_axis_tol || inside_epsilon_box || last_axis_smaller_than_scalar_tol ||
                                   degenerate_interval);
        return contains_zero;
    }

    template <typename T, typename I, int VSIZE>
    static inline bool active_box_codomain_acceptance_vf(const VActiveBox<T, I>& box,
                                                         const VQuery<T, VSIZE>& query,
                                                         const int q,
                                                         const T tol,
                                                         const VTolerances<T, VSIZE>& tols,
                                                         bool& accepted) {
        T fmin[3] = {std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
        T fmax[3] = {
            std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};

        for (int mask = 0; mask < 8; ++mask) {
            const T t = (mask & 1) ? box.upper[0] : box.lower[0];
            const T u = (mask & 2) ? box.upper[1] : box.lower[1];
            const T v = (mask & 4) ? box.upper[2] : box.lower[2];
            T f;

            diff_vf_component<T>(query.s0.x[q],
                                 query.s1.x[q],
                                 query.s2.x[q],
                                 query.s3.x[q],
                                 query.e0.x[q],
                                 query.e1.x[q],
                                 query.e2.x[q],
                                 query.e3.x[q],
                                 t,
                                 u,
                                 v,
                                 f);
            fmin[0] = sccd::min<T>(fmin[0], f);
            fmax[0] = sccd::max<T>(fmax[0], f);

            diff_vf_component<T>(query.s0.y[q],
                                 query.s1.y[q],
                                 query.s2.y[q],
                                 query.s3.y[q],
                                 query.e0.y[q],
                                 query.e1.y[q],
                                 query.e2.y[q],
                                 query.e3.y[q],
                                 t,
                                 u,
                                 v,
                                 f);
            fmin[1] = sccd::min<T>(fmin[1], f);
            fmax[1] = sccd::max<T>(fmax[1], f);

            diff_vf_component<T>(query.s0.z[q],
                                 query.s1.z[q],
                                 query.s2.z[q],
                                 query.s3.z[q],
                                 query.e0.z[q],
                                 query.e1.z[q],
                                 query.e2.z[q],
                                 query.e3.z[q],
                                 t,
                                 u,
                                 v,
                                 f);
            fmin[2] = sccd::min<T>(fmin[2], f);
            fmax[2] = sccd::max<T>(fmax[2], f);
        }

        const T query_tols[3] = {tols.axis[0][q], tols.axis[1][q], tols.axis[2][q]};
        return codomain_acceptance_grid_vf<T>(fmin, fmax, tol, query_tols, accepted);
    }

    template <typename T, typename I, int VSIZE>
    static inline void report_contained_domain(const VDomain<T, I, VSIZE>& domain,
                                               const int i,
                                               T* const SCCD_RESTRICT toi_packed,
                                               const VTolerances<T, VSIZE>& tols) {
        const int q = static_cast<int>(domain.query_index[i]);
        if (!lane_valid<T, I, VSIZE>(domain, i, toi_packed[q], tols)) {
            return;
        }

        const T candidate = domain.tuv[0].lower[i];
        if (candidate < toi_packed[q]) {
            toi_packed[q] = candidate;
        }
    }

    template <typename T, typename I, int VSIZE>
    static inline void push_child_lane(const VDomain<T, I, VSIZE>& src,
                                       const int src_lane,
                                       const int split_axis,
                                       const T lower,
                                       const T upper,
                                       const T toi,
                                       const VTolerances<T, VSIZE>& tols,
                                       VDomain<T, I, VSIZE>& pack,
                                       int& pack_size,
                                       std::vector<VDomain<T, I, VSIZE>>& stack) {
        if (!(lower < upper)) {
            return;
        }

        const int dst = pack_size++;
        pack.query_index[dst] = src.query_index[src_lane];
        pack.depth[dst] = src.depth[src_lane] + 1;
        pack.active[dst] = 1;
        for (int d = 0; d < 3; ++d) {
            pack.tuv[d].lower[dst] = src.tuv[d].lower[src_lane];
            pack.tuv[d].upper[dst] = src.tuv[d].upper[src_lane];
        }
        pack.tuv[split_axis].lower[dst] = lower;
        pack.tuv[split_axis].upper[dst] = upper;
        pack.tuv[0].upper[dst] = sccd::min<T>(pack.tuv[0].upper[dst], toi);

        if (!lane_valid<T, I, VSIZE>(pack, dst, toi, tols)) {
            --pack_size;
            return;
        }

        if (pack_size == VSIZE) {
            stack.push_back(pack);
            for (int i = 0; i < VSIZE; ++i) {
                deactivate_lane(pack, i);
            }
            pack_size = 0;
        }
    }

    template <typename T, typename I, int VSIZE>
    static void split_domain_along_longest_axis(const VDomain<T, I, VSIZE>& domain,
                                                const VQuery<T, VSIZE>& query,
                                                const uint8_t* const SCCD_RESTRICT contains_origin,
                                                const uint8_t* const SCCD_RESTRICT accepted,
                                                const int* const SCCD_RESTRICT longest_axis_to_split,
                                                T* const SCCD_RESTRICT toi_packed,
                                                const VTolerances<T, VSIZE>& tols,
                                                const int max_depth,
                                                std::vector<VDomain<T, I, VSIZE>>& stack) {
        VDomain<T, I, VSIZE> pack;
        for (int i = 0; i < VSIZE; ++i) {
            deactivate_lane(pack, i);
        }
        int pack_size = 0;

        for (int i = 0; i < VSIZE; ++i) {
            if (!domain.active[i] || !contains_origin[i] || accepted[i] || domain.depth[i] >= max_depth) {
                continue;
            }

            const int q = static_cast<int>(domain.query_index[i]);
            const T toi = toi_packed[q];
            if (!lane_valid<T, I, VSIZE>(domain, i, toi, tols)) {
                continue;
            }

            const int split_axis = longest_axis_to_split[i];
            const T lo = domain.tuv[split_axis].lower[i];
            const T hi = domain.tuv[split_axis].upper[i];

            if constexpr (SCCD_VNARROWPHASE_ADAPTIVE_SPLIT) {
                constexpr int NSPLIT = 4;
                alignas(64) T splitters[NSPLIT];
                adaptive_splitters_vf<NSPLIT, T, VSIZE>(query,
                                                        q,
                                                        split_axis,
                                                        domain.tuv[0].lower[i],
                                                        domain.tuv[0].upper[i],
                                                        domain.tuv[1].lower[i],
                                                        domain.tuv[1].upper[i],
                                                        domain.tuv[2].lower[i],
                                                        domain.tuv[2].upper[i],
                                                        splitters);

                T sample_min = lo;
                bool split = false;
                for (int s = 0; s <= NSPLIT; ++s) {
                    const T sample_max = s == NSPLIT ? hi : splitters[s];
                    if (sample_min < sample_max) {
                        push_child_lane<T, I, VSIZE>(
                            domain, i, split_axis, sample_min, sample_max, toi, tols, pack, pack_size, stack);
                        split = true;
                    }
                    sample_min = sample_max;
                }
                if (!split) {
                    report_contained_domain<T, I, VSIZE>(domain, i, toi_packed, tols);
                }
                continue;
            }

            const T mid = (lo + hi) * T(0.5);

            if (!(lo < mid && mid < hi)) {
                report_contained_domain<T, I, VSIZE>(domain, i, toi_packed, tols);
                continue;
            }

            push_child_lane<T, I, VSIZE>(domain, i, split_axis, lo, mid, toi, tols, pack, pack_size, stack);

            bool push_upper = true;
            if (split_axis == 0) {
                push_upper = mid < toi;
            } else if (split_axis == 1) {
                push_upper = mid + domain.tuv[2].lower[i] <= T(1);
            } else {
                push_upper = domain.tuv[1].lower[i] + mid <= T(1);
            }

            if (push_upper) {
                push_child_lane<T, I, VSIZE>(domain, i, split_axis, mid, hi, toi, tols, pack, pack_size, stack);
            }
        }

        if (pack_size != 0) {
            for (int i = pack_size; i < VSIZE; ++i) {
                deactivate_lane(pack, i);
            }
            stack.push_back(pack);
        }
    }

    template <typename T, typename I, int VSIZE>
    static void process_accepted_domains(const VDomain<T, I, VSIZE>& domain,
                                         const uint8_t* const SCCD_RESTRICT contains_origin,
                                         const uint8_t* const SCCD_RESTRICT accepted,
                                         T* const SCCD_RESTRICT toi_packed,
                                         const VTolerances<T, VSIZE>& tols,
                                         const int max_depth) {
        for (int i = 0; i < VSIZE; ++i) {
            if (!domain.active[i]) {
                continue;
            }

            const bool depth_limit = contains_origin[i] && domain.depth[i] >= max_depth;
            if (accepted[i] || depth_limit) {
                report_contained_domain<T, I, VSIZE>(domain, i, toi_packed, tols);
            }
        }
    }

    template <typename T, typename I, int VSIZE>
    static inline bool has_active_lanes(const VDomain<T, I, VSIZE>& domain) {
        for (int i = 0; i < VSIZE; ++i) {
            if (domain.active[i]) {
                return true;
            }
        }
        return false;
    }

    template <typename T, typename I, int VSIZE>
    static inline bool has_pending_split_lanes(const VDomain<T, I, VSIZE>& domain,
                                               const uint8_t* const SCCD_RESTRICT contains_origin,
                                               const uint8_t* const SCCD_RESTRICT accepted,
                                               const int max_depth) {
        for (int i = 0; i < VSIZE; ++i) {
            if (domain.active[i] && contains_origin[i] && !accepted[i] && domain.depth[i] < max_depth) {
                return true;
            }
        }
        return false;
    }

    template <typename T, typename I, int VSIZE>
    static inline bool active_box_valid(const VActiveBox<T, I>& box, const T toi, const VTolerances<T, VSIZE>&) {
        return box.lower[0] < toi && box.lower[1] + box.lower[2] <= T(1);
    }

    template <typename T, typename I, int VSIZE>
    static inline int active_box_next_split_axis(const VActiveBox<T, I>& box, const VTolerances<T, VSIZE>& tols) {
        const int q = static_cast<int>(box.query_index);
        int axis = 0;
        T best = -std::numeric_limits<T>::infinity();
        for (int d = 0; d < 3; ++d) {
            const T width = box.upper[d] - box.lower[d];
            if (width <= tols.axis[d][q]) {
                continue;
            }
            const T score = width / tols.axis[d][q];
            if (score > best) {
                best = score;
                axis = d;
            }
        }
        return axis;
    }

    template <typename T, typename I>
    static inline int active_box_widest_axis(const VActiveBox<T, I>& box) {
        const T dt = box.upper[0] - box.lower[0];
        const T du = box.upper[1] - box.lower[1];
        const T dv = box.upper[2] - box.lower[2];
        if (du > dt && du >= dv) {
            return 1;
        }
        if (dv > dt && dv > du) {
            return 2;
        }
        return 0;
    }

    template <typename T, typename I, int VSIZE>
    static inline void set_lane_from_active_box(VDomain<T, I, VSIZE>& domain,
                                                const int lane,
                                                const VActiveBox<T, I>& box,
                                                const T toi,
                                                const VTolerances<T, VSIZE>& tols) {
        domain.query_index[lane] = box.query_index;
        domain.depth[lane] = box.depth;
        domain.active[lane] = 1;
        for (int d = 0; d < 3; ++d) {
            domain.tuv[d].lower[lane] = box.lower[d];
            domain.tuv[d].upper[lane] = d == 0 ? sccd::min<T>(box.upper[d], toi) : box.upper[d];
        }
        if (!lane_valid<T, I, VSIZE>(domain, lane, toi, tols)) {
            deactivate_lane<T, I, VSIZE>(domain, lane);
        }
    }

    template <typename T, typename I, int VSIZE>
    static inline VActiveBox<T, I> active_box_from_lane(const VDomain<T, I, VSIZE>& domain, const int lane) {
        VActiveBox<T, I> box;
        box.query_index = domain.query_index[lane];
        box.depth = domain.depth[lane];
        for (int d = 0; d < 3; ++d) {
            box.lower[d] = domain.tuv[d].lower[lane];
            box.upper[d] = domain.tuv[d].upper[lane];
        }
        return box;
    }

    template <typename T, typename I, int VSIZE>
    static inline void refill_inactive_lanes(VDomain<T, I, VSIZE>& domain,
                                             std::vector<VActiveBox<T, I>>& stack,
                                             const T* const SCCD_RESTRICT toi_packed,
                                             const VTolerances<T, VSIZE>& tols) {
        for (int i = 0; i < VSIZE && !stack.empty(); ++i) {
            if (domain.active[i]) {
                continue;
            }

            VActiveBox<T, I> box = stack.back();
            stack.pop_back();
            set_lane_from_active_box<T, I, VSIZE>(domain, i, box, toi_packed[static_cast<int>(box.query_index)], tols);
        }
    }

    template <typename T, typename I, int VSIZE>
    static inline void split_or_deactivate_lane(VDomain<T, I, VSIZE>& domain,
                                                const int lane,
                                                const VQuery<T, VSIZE>& query,
                                                T* const SCCD_RESTRICT toi_packed,
                                                const VTolerances<T, VSIZE>& tols,
                                                const T tol,
                                                const int max_depth,
                                                std::vector<VActiveBox<T, I>>& stack) {
        VActiveBox<T, I> parent = active_box_from_lane<T, I, VSIZE>(domain, lane);
        const int q = static_cast<int>(parent.query_index);
        T toi = toi_packed[q];
        const int split_axis = [](const VActiveBox<T, I>& box, const VTolerances<T, VSIZE>& tolerances) {
            if constexpr (SCCD_VNARROWPHASE_ADAPTIVE_SPLIT) {
                return active_box_widest_axis<T, I>(box);
            }
            return active_box_next_split_axis<T, I, VSIZE>(box, tolerances);
        }(parent, tols);

        if constexpr (SCCD_VNARROWPHASE_ADAPTIVE_SPLIT) {
            constexpr int NSPLIT = 4;
            alignas(64) T splitters[NSPLIT];
            adaptive_splitters_vf<NSPLIT, T, VSIZE>(query,
                                                    q,
                                                    split_axis,
                                                    parent.lower[0],
                                                    parent.upper[0],
                                                    parent.lower[1],
                                                    parent.upper[1],
                                                    parent.lower[2],
                                                    parent.upper[2],
                                                    splitters);

            VActiveBox<T, I> keep;
            bool has_keep = false;
            T sample_min = parent.lower[split_axis];
            for (int s = 0; s <= NSPLIT; ++s) {
                const T sample_max = s == NSPLIT ? parent.upper[split_axis] : splitters[s];
                if (!(sample_min < sample_max)) {
                    sample_min = sample_max;
                    continue;
                }

                VActiveBox<T, I> child = parent;
                child.depth = parent.depth + 1;
                child.lower[split_axis] = sample_min;
                child.upper[split_axis] = sample_max;
                child.upper[0] = sccd::min<T>(child.upper[0], toi);

                if (child.lower[0] >= toi || child.lower[1] + child.lower[2] >= T(1) + tol) {
                    sample_min = sample_max;
                    continue;
                }

                bool accepted = false;
                if (active_box_codomain_acceptance_vf<T, I, VSIZE>(child, query, q, tol, tols, accepted)) {
                    accepted = accepted && (child.lower[0] > T(0));
                    if (accepted || child.depth > max_depth) {
                        if (child.lower[0] < toi &&
                            child.lower[1] + child.lower[2] < T(1) + tols.axis[1][q] + tols.axis[2][q]) {
                            toi = sccd::min<T>(toi, child.lower[0]);
                            toi_packed[q] = toi;
                        }
                    } else {
                        if (!has_keep || child.lower[0] < keep.lower[0]) {
                            if (has_keep) {
                                stack.push_back(keep);
                            }
                            keep = child;
                            has_keep = true;
                        } else {
                            stack.push_back(child);
                        }
                    }
                }
                sample_min = sample_max;
            }

            if (has_keep) {
                set_lane_from_active_box<T, I, VSIZE>(domain, lane, keep, toi, tols);
            } else {
                deactivate_lane<T, I, VSIZE>(domain, lane);
            }
            return;
        }

        const T mid = (parent.lower[split_axis] + parent.upper[split_axis]) * T(0.5);

        if (!(parent.lower[split_axis] < mid && mid < parent.upper[split_axis])) {
            report_contained_domain<T, I, VSIZE>(domain, lane, toi_packed, tols);
            deactivate_lane<T, I, VSIZE>(domain, lane);
            return;
        }

        VActiveBox<T, I> left = parent;
        VActiveBox<T, I> right = parent;
        left.depth = parent.depth + 1;
        right.depth = parent.depth + 1;
        left.upper[split_axis] = mid;
        right.lower[split_axis] = mid;
        left.upper[0] = sccd::min<T>(left.upper[0], toi);
        right.upper[0] = sccd::min<T>(right.upper[0], toi);

        const bool left_valid = active_box_valid<T, I, VSIZE>(left, toi, tols);
        const bool right_valid = active_box_valid<T, I, VSIZE>(right, toi, tols);

        if (left_valid && right_valid) {
            const bool keep_left = left.lower[0] <= right.lower[0];
            const VActiveBox<T, I>& keep = keep_left ? left : right;
            const VActiveBox<T, I>& push = keep_left ? right : left;
            set_lane_from_active_box<T, I, VSIZE>(domain, lane, keep, toi, tols);
            stack.push_back(push);
        } else if (left_valid) {
            set_lane_from_active_box<T, I, VSIZE>(domain, lane, left, toi, tols);
        } else if (right_valid) {
            set_lane_from_active_box<T, I, VSIZE>(domain, lane, right, toi, tols);
        } else {
            deactivate_lane<T, I, VSIZE>(domain, lane);
        }
    }

    template <typename T, typename I, int VSIZE>
    static inline void advance_lanes_after_evaluation(VDomain<T, I, VSIZE>& domain,
                                                      const uint8_t* const SCCD_RESTRICT contains_origin,
                                                      const uint8_t* const SCCD_RESTRICT accepted,
                                                      const VQuery<T, VSIZE>& query,
                                                      T* const SCCD_RESTRICT toi_packed,
                                                      const VTolerances<T, VSIZE>& tols,
                                                      const T tol,
                                                      const int max_depth,
                                                      std::vector<VActiveBox<T, I>>& stack) {
        for (int i = 0; i < VSIZE; ++i) {
            if (!domain.active[i]) {
                continue;
            }

            const int q = static_cast<int>(domain.query_index[i]);
            if (!contains_origin[i] || accepted[i] || domain.depth[i] >= max_depth ||
                !lane_valid<T, I, VSIZE>(domain, i, toi_packed[q], tols)) {
                deactivate_lane<T, I, VSIZE>(domain, i);
            } else {
                split_or_deactivate_lane<T, I, VSIZE>(domain, i, query, toi_packed, tols, tol, max_depth, stack);
            }
        }
    }

    template <int nxe, typename T, typename I>
    int v_narrow_phase_vf(const size_t noverlaps,
                          const I* const SCCD_RESTRICT voveralp,
                          const I* const SCCD_RESTRICT foveralp,
                          T** const SCCD_RESTRICT v0,
                          T** const SCCD_RESTRICT v1,
                          const size_t face_stride,
                          I** const SCCD_RESTRICT faces,
                          const T max_toi,
                          T* const SCCD_RESTRICT toi,
                          const int max_depth,
                          const T tol,
                          const int toi_stride) {
        using T_HP = double;

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) {
                toi[0] = sccd::min<T>(T(1), max_toi);
            }
            return 0;
        }
        assert(toi != nullptr);

        const T_HP max_domain_toi = sccd::min<T_HP>(T_HP(1), T_HP(max_toi));

        int constexpr VSIZE = 8;
        using VQueryT = VQuery<T_HP, VSIZE>;
        using VDomainT = VDomain<T_HP, I, VSIZE>;
        using VCodomainT = VCodomain<T_HP, VSIZE>;
        using VTolerancesT = VTolerances<T_HP, VSIZE>;
        using VActiveBoxT = VActiveBox<T_HP, I>;

        const ptrdiff_t nblocks = static_cast<ptrdiff_t>((noverlaps + VSIZE - 1) / VSIZE);
        std::atomic<T_HP> global_min{max_domain_toi};

        sccd::parallel_for_br(0, nblocks, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            std::vector<VActiveBoxT> stack;
            stack.reserve(1024);

            for (ptrdiff_t ib = rbegin; ib < rend; ++ib) {
                const ptrdiff_t block_begin = ib * VSIZE;
                const ptrdiff_t block_size =
                    std::min<ptrdiff_t>(VSIZE, static_cast<ptrdiff_t>(noverlaps) - block_begin);

                VQueryT query;
                load_query<T_HP, T, I, VSIZE>(ib, block_size, voveralp, foveralp, v0, v1, face_stride, faces, query);

                VTolerancesT tols;
                compute_tolerances<T_HP, VSIZE>(tol, query, tols);

                T_HP toi_packed[VSIZE];
                const T_HP initial_toi = toi_stride == 0 ? global_min.load(std::memory_order_relaxed) : max_domain_toi;
                for (int i = 0; i < VSIZE; ++i) {
                    toi_packed[i] = initial_toi;
                }

                stack.clear();

                for (ptrdiff_t i = block_size - 1; i >= 0; --i) {
                    VActiveBoxT root;
                    root.query_index = static_cast<I>(i);
                    root.depth = 0;
                    root.lower[0] = T_HP(0);
                    root.upper[0] = toi_packed[i];
                    root.lower[1] = T_HP(0);
                    root.upper[1] = T_HP(1);
                    root.lower[2] = T_HP(0);
                    root.upper[2] = T_HP(1);
                    stack.push_back(root);
                }

                VDomainT domain;
                for (int i = 0; i < VSIZE; ++i) {
                    deactivate_lane<T_HP, I, VSIZE>(domain, i);
                }
                refill_inactive_lanes<T_HP, I, VSIZE>(domain, stack, toi_packed, tols);

                while (has_active_lanes<T_HP, I, VSIZE>(domain) || !stack.empty()) {
                    refill_inactive_lanes<T_HP, I, VSIZE>(domain, stack, toi_packed, tols);

                    const T_HP current_min_toi =
                        toi_stride == 0 ? global_min.load(std::memory_order_relaxed) : max_domain_toi;

                    for (int i = 0; i < VSIZE; ++i) {
                        if (!domain.active[i]) {
                            continue;
                        }
                        const int q = static_cast<int>(domain.query_index[i]);
                        if (toi_stride == 0) {
                            toi_packed[q] = sccd::min<T_HP>(toi_packed[q], current_min_toi);
                        }
                        domain.tuv[0].upper[i] = sccd::min<T_HP>(domain.tuv[0].upper[i], toi_packed[q]);
                        if (!lane_valid<T_HP, I, VSIZE>(domain, i, toi_packed[q], tols)) {
                            domain.active[i] = 0;
                        }
                    }

                    if (!has_active_lanes<T_HP, I, VSIZE>(domain)) {
                        continue;
                    }

                    VCodomainT codomain;
                    compute_codomain<T_HP, I, VSIZE>(domain, query, codomain);

                    uint8_t contains_origin_mask[VSIZE];
                    uint8_t acceptance_mask[VSIZE];
                    compute_masks<T_HP, I, VSIZE>(
                        domain, codomain, tols, tol, contains_origin_mask, acceptance_mask);

                    process_accepted_domains<T_HP, I, VSIZE>(
                        domain, contains_origin_mask, acceptance_mask, toi_packed, tols, max_depth);

                    if (toi_stride == 0) {
                        T_HP local_min = max_domain_toi;
                        for (int i = 0; i < VSIZE; ++i) {
                            local_min = sccd::min<T_HP>(local_min, toi_packed[i]);
                        }
                        if (local_min < max_domain_toi) {
                            atomic_min_relaxed<T_HP>(global_min, local_min);
                        }
                    }

                    if (!has_pending_split_lanes<T_HP, I, VSIZE>(
                            domain, contains_origin_mask, acceptance_mask, max_depth)) {
                        for (int i = 0; i < VSIZE; ++i) {
                            if (domain.active[i]) {
                                deactivate_lane<T_HP, I, VSIZE>(domain, i);
                            }
                        }
                        continue;
                    }

                    advance_lanes_after_evaluation<T_HP, I, VSIZE>(domain,
                                                                   contains_origin_mask,
                                                                   acceptance_mask,
                                                                   query,
                                                                   toi_packed,
                                                                   tols,
                                                                   tol,
                                                                   max_depth,
                                                                   stack);
                }

                for (ptrdiff_t i = 0; i < block_size; ++i) {
                    if (toi_stride == 1) {
                        toi[block_begin + i] = T(toi_packed[i]);
                    }
                }
            }
        });

        if (toi_stride == 0) {
            toi[0] = T(global_min.load(std::memory_order_relaxed));
        }

        return 0;
    }
}  // namespace sccd

#endif
