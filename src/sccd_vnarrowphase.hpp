#ifndef SCCD_VNARROWPHASE_HPP
#define SCCD_VNARROWPHASE_HPP

#include "sccd_base.hpp"
#include "smath.hpp"
#include "snumtol.hpp"
#include "sparallel.hpp"

#include <algorithm>
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

namespace sccd {

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

    template <typename T, int VSIZE>
    static inline void copy_query_lane(const VQuery<T, VSIZE>& src,
                                       const int src_lane,
                                       VQuery<T, VSIZE>& dst,
                                       const int dst_lane) {
        dst.s0.x[dst_lane] = src.s0.x[src_lane];
        dst.s0.y[dst_lane] = src.s0.y[src_lane];
        dst.s0.z[dst_lane] = src.s0.z[src_lane];
        dst.s1.x[dst_lane] = src.s1.x[src_lane];
        dst.s1.y[dst_lane] = src.s1.y[src_lane];
        dst.s1.z[dst_lane] = src.s1.z[src_lane];
        dst.s2.x[dst_lane] = src.s2.x[src_lane];
        dst.s2.y[dst_lane] = src.s2.y[src_lane];
        dst.s2.z[dst_lane] = src.s2.z[src_lane];
        dst.s3.x[dst_lane] = src.s3.x[src_lane];
        dst.s3.y[dst_lane] = src.s3.y[src_lane];
        dst.s3.z[dst_lane] = src.s3.z[src_lane];
        dst.e0.x[dst_lane] = src.e0.x[src_lane];
        dst.e0.y[dst_lane] = src.e0.y[src_lane];
        dst.e0.z[dst_lane] = src.e0.z[src_lane];
        dst.e1.x[dst_lane] = src.e1.x[src_lane];
        dst.e1.y[dst_lane] = src.e1.y[src_lane];
        dst.e1.z[dst_lane] = src.e1.z[src_lane];
        dst.e2.x[dst_lane] = src.e2.x[src_lane];
        dst.e2.y[dst_lane] = src.e2.y[src_lane];
        dst.e2.z[dst_lane] = src.e2.z[src_lane];
        dst.e3.x[dst_lane] = src.e3.x[src_lane];
        dst.e3.y[dst_lane] = src.e3.y[src_lane];
        dst.e3.z[dst_lane] = src.e3.z[src_lane];
    }

    template <typename T, typename I, int VSIZE>
    static void gather_query(const VDomain<T, I, VSIZE>& domain,
                             const VQuery<T, VSIZE>& query,
                             VQuery<T, VSIZE>& packed_query) {
        for (int i = 0; i < VSIZE; ++i) {
            copy_query_lane<T, VSIZE>(query, static_cast<int>(domain.query_index[i]), packed_query, i);
        }
    }

    template <typename T, typename I, int VSIZE>
    static void gather_tolerances(const VDomain<T, I, VSIZE>& domain,
                                  const VTolerances<T, VSIZE>& tols,
                                  VTolerances<T, VSIZE>& packed_tols) {
        for (int i = 0; i < VSIZE; ++i) {
            const int q = static_cast<int>(domain.query_index[i]);
            packed_tols.axis[0][i] = tols.axis[0][q];
            packed_tols.axis[1][i] = tols.axis[1][q];
            packed_tols.axis[2][i] = tols.axis[2][q];
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
    static inline float64x2_t diff_vf_component_neon(const double* const SCCD_RESTRICT s0,
                                                     const double* const SCCD_RESTRICT s1,
                                                     const double* const SCCD_RESTRICT s2,
                                                     const double* const SCCD_RESTRICT s3,
                                                     const double* const SCCD_RESTRICT e0,
                                                     const double* const SCCD_RESTRICT e1,
                                                     const double* const SCCD_RESTRICT e2,
                                                     const double* const SCCD_RESTRICT e3,
                                                     const float64x2_t t,
                                                     const float64x2_t u,
                                                     const float64x2_t v,
                                                     const int i) {
        const float64x2_t one = vdupq_n_f64(1.0);
        const float64x2_t omt = vsubq_f64(one, t);
        const float64x2_t o = vsubq_f64(vsubq_f64(one, u), v);

        const float64x2_t vs0 = vld1q_f64(s0 + i);
        const float64x2_t vs1 = vld1q_f64(s1 + i);
        const float64x2_t vs2 = vld1q_f64(s2 + i);
        const float64x2_t vs3 = vld1q_f64(s3 + i);
        const float64x2_t ve0 = vld1q_f64(e0 + i);
        const float64x2_t ve1 = vld1q_f64(e1 + i);
        const float64x2_t ve2 = vld1q_f64(e2 + i);
        const float64x2_t ve3 = vld1q_f64(e3 + i);

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
        VQuery<double, VSIZE> packed_query;
        gather_query<double, I, VSIZE>(domain, query, packed_query);

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

                float64x2_t f = diff_vf_component_neon(packed_query.s0.x,
                                                       packed_query.s1.x,
                                                       packed_query.s2.x,
                                                       packed_query.s3.x,
                                                       packed_query.e0.x,
                                                       packed_query.e1.x,
                                                       packed_query.e2.x,
                                                       packed_query.e3.x,
                                                       t,
                                                       u,
                                                       v,
                                                       i);
                vst1q_f64(codomain.xyz[0].lower + i, vminq_f64(vld1q_f64(codomain.xyz[0].lower + i), f));
                vst1q_f64(codomain.xyz[0].upper + i, vmaxq_f64(vld1q_f64(codomain.xyz[0].upper + i), f));

                f = diff_vf_component_neon(packed_query.s0.y,
                                           packed_query.s1.y,
                                           packed_query.s2.y,
                                           packed_query.s3.y,
                                           packed_query.e0.y,
                                           packed_query.e1.y,
                                           packed_query.e2.y,
                                           packed_query.e3.y,
                                           t,
                                           u,
                                           v,
                                           i);
                vst1q_f64(codomain.xyz[1].lower + i, vminq_f64(vld1q_f64(codomain.xyz[1].lower + i), f));
                vst1q_f64(codomain.xyz[1].upper + i, vmaxq_f64(vld1q_f64(codomain.xyz[1].upper + i), f));

                f = diff_vf_component_neon(packed_query.s0.z,
                                           packed_query.s1.z,
                                           packed_query.s2.z,
                                           packed_query.s3.z,
                                           packed_query.e0.z,
                                           packed_query.e1.z,
                                           packed_query.e2.z,
                                           packed_query.e3.z,
                                           t,
                                           u,
                                           v,
                                           i);
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
            bool smaller_than_axis_tol = true;
            bool inside_epsilon_box = true;
            T true_tol = T(0);

            for (int d = 0; d < 3; ++d) {
                const T fmin = codomain.xyz[d].lower[i];
                const T fmax = codomain.xyz[d].upper[i];
                const T width = fmax - fmin;
                contains = contains && (fmin <= tol) && (fmax >= -tol);
                smaller_than_axis_tol = smaller_than_axis_tol && (width <= tols.axis[d][q]);
                inside_epsilon_box = inside_epsilon_box && (fmin >= -tol) && (fmax <= tol);
                true_tol = sccd::max<T>(true_tol, width);
            }

            const bool real_tol_smaller = (domain.tuv[0].lower[i] > T(0)) && (true_tol < tol);
            const bool terminal = (domain.tuv[0].lower[i] >= domain.tuv[0].upper[i]) ||
                                  (domain.tuv[1].lower[i] >= domain.tuv[1].upper[i]) ||
                                  (domain.tuv[2].lower[i] >= domain.tuv[2].upper[i]);
            const bool positive_time = domain.tuv[0].lower[i] > T(0);

            contains_origin[i] = contains;
            accepted[i] = contains && positive_time &&
                          (smaller_than_axis_tol || inside_epsilon_box || real_tol_smaller || terminal);
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
        VTolerances<double, VSIZE> packed_tols;
        gather_tolerances<double, I, VSIZE>(domain, tols, packed_tols);

        const float64x2_t vtol = vdupq_n_f64(tol);
        const float64x2_t vntol = vdupq_n_f64(-tol);
        const float64x2_t vzero = vdupq_n_f64(0.0);

        for (int i = 0; i < VSIZE; i += 2) {
            uint64x2_t contains = vdupq_n_u64(~uint64_t(0));
            uint64x2_t smaller = vdupq_n_u64(~uint64_t(0));
            uint64x2_t inside = vdupq_n_u64(~uint64_t(0));
            float64x2_t true_tol = vdupq_n_f64(0.0);

            for (int d = 0; d < 3; ++d) {
                const float64x2_t fmin = vld1q_f64(codomain.xyz[d].lower + i);
                const float64x2_t fmax = vld1q_f64(codomain.xyz[d].upper + i);
                const float64x2_t width = vsubq_f64(fmax, fmin);
                contains = vandq_u64(contains, vandq_u64(vcleq_f64(fmin, vtol), vcgeq_f64(fmax, vntol)));
                smaller = vandq_u64(smaller, vcleq_f64(width, vld1q_f64(packed_tols.axis[d] + i)));
                inside = vandq_u64(inside, vandq_u64(vcgeq_f64(fmin, vntol), vcleq_f64(fmax, vtol)));
                true_tol = vmaxq_f64(true_tol, width);
            }

            const float64x2_t tl = vld1q_f64(domain.tuv[0].lower + i);
            const uint64x2_t positive_time = vcgtq_f64(tl, vzero);
            const uint64x2_t real_small = vandq_u64(vcgtq_f64(tl, vzero), vcltq_f64(true_tol, vtol));
            uint64x2_t terminal = vcgeq_f64(vld1q_f64(domain.tuv[0].lower + i), vld1q_f64(domain.tuv[0].upper + i));
            terminal =
                vorrq_u64(terminal, vcgeq_f64(vld1q_f64(domain.tuv[1].lower + i), vld1q_f64(domain.tuv[1].upper + i)));
            terminal =
                vorrq_u64(terminal, vcgeq_f64(vld1q_f64(domain.tuv[2].lower + i), vld1q_f64(domain.tuv[2].upper + i)));

            const uint64x2_t accept =
                vandq_u64(vandq_u64(contains, positive_time),
                          vorrq_u64(vorrq_u64(smaller, inside), vorrq_u64(real_small, terminal)));
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
                                  const VTolerances<T, VSIZE>& tols) {
        const int q = static_cast<int>(domain.query_index[i]);
        return domain.active[i] && domain.tuv[0].lower[i] < toi &&
               domain.tuv[1].lower[i] + domain.tuv[2].lower[i] < T(1) + tols.axis[1][q] + tols.axis[2][q];
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
            const T mid = (lo + hi) * T(0.5);

            if (!(lo < mid && mid < hi)) {
                if (domain.tuv[0].lower[i] > T(0)) {
                    toi_packed[q] = sccd::min<T>(toi_packed[q], domain.tuv[0].lower[i]);
                }
                continue;
            }

            push_child_lane<T, I, VSIZE>(domain, i, split_axis, lo, mid, toi, tols, pack, pack_size, stack);

            bool push_upper = true;
            if (split_axis == 0) {
                push_upper = mid < toi;
            } else if (split_axis == 1) {
                push_upper = mid + domain.tuv[2].lower[i] < T(1) + tols.axis[1][q] + tols.axis[2][q];
            } else {
                push_upper = domain.tuv[1].lower[i] + mid < T(1) + tols.axis[1][q] + tols.axis[2][q];
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

            const int q = static_cast<int>(domain.query_index[i]);
            const bool depth_limit = contains_origin[i] && domain.depth[i] >= max_depth;
            const bool positive_time = domain.tuv[0].lower[i] > T(0);
            if ((accepted[i] || (depth_limit && positive_time)) &&
                lane_valid<T, I, VSIZE>(domain, i, toi_packed[q], tols)) {
                toi_packed[q] = sccd::min<T>(toi_packed[q], domain.tuv[0].lower[i]);
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
                          const int toi_stride) {
        using T_HP = double;

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) {
                toi[0] = max_toi;
            }
            return 0;
        }
        assert(toi != nullptr);

        int SCCD_MAX_DEPTH = 32;
        SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);

        T_HP SCCD_TOL = std::is_same<float, T_HP>::value ? T_HP(1e-8) : T_HP(1e-14);
        SCCD_READ_ENV(SCCD_TOL, atof);

        int constexpr VSIZE = 16;
        using VQueryT = VQuery<T_HP, VSIZE>;
        using VDomainT = VDomain<T_HP, I, VSIZE>;
        using VCodomainT = VCodomain<T_HP, VSIZE>;
        using VTolerancesT = VTolerances<T_HP, VSIZE>;

        const ptrdiff_t nblocks = static_cast<ptrdiff_t>((noverlaps + VSIZE - 1) / VSIZE);
        std::vector<T> block_min(toi_stride == 0 ? nblocks : 0, max_toi);

        sccd::parallel_for_br(0, nblocks, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            std::vector<VDomainT> stack;
            stack.reserve(1024);

            for (ptrdiff_t ib = rbegin; ib < rend; ++ib) {
                const ptrdiff_t block_begin = ib * VSIZE;
                const ptrdiff_t block_size =
                    std::min<ptrdiff_t>(VSIZE, static_cast<ptrdiff_t>(noverlaps) - block_begin);

                VQueryT query;
                load_query<T_HP, T, I, VSIZE>(ib, block_size, voveralp, foveralp, v0, v1, face_stride, faces, query);

                VTolerancesT tols;
                compute_tolerances<T_HP, VSIZE>(SCCD_TOL, query, tols);

                VDomainT domain;
                init_domain<T_HP, I, VSIZE>(block_size, T_HP(max_toi), domain);

                T_HP toi_packed[VSIZE];
                for (int i = 0; i < VSIZE; ++i) {
                    toi_packed[i] = T_HP(max_toi);
                }

                stack.clear();
                stack.push_back(domain);

                while (!stack.empty()) {
                    domain = stack.back();
                    stack.pop_back();

                    for (int i = 0; i < VSIZE; ++i) {
                        if (!domain.active[i]) {
                            continue;
                        }
                        const int q = static_cast<int>(domain.query_index[i]);
                        domain.tuv[0].upper[i] = sccd::min<T_HP>(domain.tuv[0].upper[i], toi_packed[q]);
                        if (!lane_valid<T_HP, I, VSIZE>(domain, i, toi_packed[q], tols)) {
                            domain.active[i] = 0;
                        }
                    }

                    VCodomainT codomain;
                    compute_codomain<T_HP, I, VSIZE>(domain, query, codomain);

                    uint8_t contains_origin_mask[VSIZE];
                    uint8_t acceptance_mask[VSIZE];
                    compute_masks<T_HP, I, VSIZE>(
                        domain, codomain, tols, SCCD_TOL, contains_origin_mask, acceptance_mask);

                    process_accepted_domains<T_HP, I, VSIZE>(
                        domain, contains_origin_mask, acceptance_mask, toi_packed, tols, SCCD_MAX_DEPTH);

                    int longest_axis_to_split[VSIZE];
                    detect_longest_axis_to_split<T_HP, I, VSIZE>(domain, longest_axis_to_split);

                    split_domain_along_longest_axis<T_HP, I, VSIZE>(domain,
                                                                    contains_origin_mask,
                                                                    acceptance_mask,
                                                                    longest_axis_to_split,
                                                                    toi_packed,
                                                                    tols,
                                                                    SCCD_MAX_DEPTH,
                                                                    stack);
                }

                T_HP local_min = T_HP(max_toi);
                for (ptrdiff_t i = 0; i < block_size; ++i) {
                    local_min = sccd::min<T_HP>(local_min, toi_packed[i]);
                    if (toi_stride == 1) {
                        toi[block_begin + i] = T(toi_packed[i]);
                    }
                }

                if (toi_stride == 0) {
                    block_min[ib] = T(local_min);
                }
            }
        });

        if (toi_stride == 0) {
            T min_t = max_toi;
            for (ptrdiff_t ib = 0; ib < nblocks; ++ib) {
                min_t = sccd::min<T>(min_t, block_min[ib]);
            }
            toi[0] = min_t;
        }

        return 0;
    }
}  // namespace sccd

#endif
