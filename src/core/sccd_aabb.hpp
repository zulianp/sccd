#ifndef SCCD_AABB_HPP
#define SCCD_AABB_HPP

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "sccd_base.hpp"
#include "sccd_math.hpp"
#include "sccd_parallel.hpp"

#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#define SCCD_AABB_DISJOINT_CHUNK_SIZE 32
#define SCCD_AABB_DISJOINT_NOVECTORIZE_THRESHOLD 16

namespace sccd {
    /**
     * \brief Return the next representable value toward +infinity.
     * \param x Input value.
     * \return Next representable float greater than or equal to x.
     */

    template <typename T>
    inline static T nextafter_up(const T x) noexcept {
        if constexpr (std::is_same<T, float>::value) {
            return nextafterf(x, FLT_MAX);
        } else {
            return nextafter(x, DBL_MAX);
        }
    }

    /**
     * \brief Return the next representable value toward -infinity.
     * \param x Input value.
     * \return Next representable float less than or equal to x.
     */
    template <typename T>
    inline static T nextafter_down(const T x) noexcept {
        if constexpr (std::is_same<T, float>::value) {
            return nextafterf(x, -FLT_MAX);
        } else {
            return nextafter(x, -DBL_MAX);
        }
    }

    template <typename T>
    inline static uint32_t disjoint(const T aminx,
                                    const T aminy,
                                    const T aminz,
                                    const T amaxx,
                                    const T amaxy,
                                    const T amaxz,
                                    const T bminx,
                                    const T bminy,
                                    const T bminz,
                                    const T bmaxx,
                                    const T bmaxy,
                                    const T bmaxz) {
        return aminx > bmaxx | aminy > bmaxy | aminz > bmaxz | bminx > amaxx | bminy > amaxy | bminz > amaxz;
    }



    /**
     * \brief Overlap bitmask of one AABB against up to 32 consecutive AABBs.
     *
     * Returns a packed mask whose bit \p i is set when B[start+i] is NOT
     * axis-disjoint from A, i.e. the pair survives to the next filter. Bits at
     * or above \p count are always zero.
     *
     * This is the bit-packed counterpart of vaabb_disjoint_one_to_many. Keeping
     * the result in a register rather than expanding it to one uint32_t per lane
     * lets callers count survivors with a popcount and visit them with a
     * bit-scan, instead of storing 32 words and re-reading them in a scalar
     * loop. It also takes an explicit \p count and handles the remainder with
     * scalar code, so no padded scratch buffers are needed for partial chunks.
     *
     * \param count Number of valid lanes, 0..SCCD_AABB_DISJOINT_CHUNK_SIZE.
     */
    template <typename T>
    inline static uint32_t vaabb_overlap_one_to_many_bits(const T aminx,
                                                          const T aminy,
                                                          const T aminz,
                                                          const T amaxx,
                                                          const T amaxy,
                                                          const T amaxz,
                                                          const T* const SCCD_RESTRICT bminx,
                                                          const T* const SCCD_RESTRICT bminy,
                                                          const T* const SCCD_RESTRICT bminz,
                                                          const T* const SCCD_RESTRICT bmaxx,
                                                          const T* const SCCD_RESTRICT bmaxy,
                                                          const T* const SCCD_RESTRICT bmaxz,
                                                          const int count) {
        assert(count >= 0 && count <= SCCD_AABB_DISJOINT_CHUNK_SIZE);

        uint32_t disjoint_bits = 0;
        int i = 0;

        if constexpr (std::is_same<T, double>::value) {
#if defined(__AVX512F__)
            const __m512d a_minx = _mm512_set1_pd(aminx);
            const __m512d a_miny = _mm512_set1_pd(aminy);
            const __m512d a_minz = _mm512_set1_pd(aminz);
            const __m512d a_maxx = _mm512_set1_pd(amaxx);
            const __m512d a_maxy = _mm512_set1_pd(amaxy);
            const __m512d a_maxz = _mm512_set1_pd(amaxz);

            for (; i + 8 <= count; i += 8) {
                const __mmask8 k = _mm512_cmp_pd_mask(a_minx, _mm512_loadu_pd(bmaxx + i), _CMP_GT_OQ) |
                                   _mm512_cmp_pd_mask(a_miny, _mm512_loadu_pd(bmaxy + i), _CMP_GT_OQ) |
                                   _mm512_cmp_pd_mask(a_minz, _mm512_loadu_pd(bmaxz + i), _CMP_GT_OQ) |
                                   _mm512_cmp_pd_mask(_mm512_loadu_pd(bminx + i), a_maxx, _CMP_GT_OQ) |
                                   _mm512_cmp_pd_mask(_mm512_loadu_pd(bminy + i), a_maxy, _CMP_GT_OQ) |
                                   _mm512_cmp_pd_mask(_mm512_loadu_pd(bminz + i), a_maxz, _CMP_GT_OQ);
                disjoint_bits |= static_cast<uint32_t>(k) << i;
            }
#elif defined(__AVX2__)
            const __m256d a_minx = _mm256_set1_pd(aminx);
            const __m256d a_miny = _mm256_set1_pd(aminy);
            const __m256d a_minz = _mm256_set1_pd(aminz);
            const __m256d a_maxx = _mm256_set1_pd(amaxx);
            const __m256d a_maxy = _mm256_set1_pd(amaxy);
            const __m256d a_maxz = _mm256_set1_pd(amaxz);

            for (; i + 4 <= count; i += 4) {
                __m256d m = _mm256_or_pd(_mm256_or_pd(_mm256_cmp_pd(a_minx, _mm256_loadu_pd(bmaxx + i), _CMP_GT_OQ),
                                                      _mm256_cmp_pd(a_miny, _mm256_loadu_pd(bmaxy + i), _CMP_GT_OQ)),
                                         _mm256_cmp_pd(a_minz, _mm256_loadu_pd(bmaxz + i), _CMP_GT_OQ));
                m = _mm256_or_pd(m,
                                 _mm256_or_pd(_mm256_cmp_pd(_mm256_loadu_pd(bminx + i), a_maxx, _CMP_GT_OQ),
                                              _mm256_cmp_pd(_mm256_loadu_pd(bminy + i), a_maxy, _CMP_GT_OQ)));
                m = _mm256_or_pd(m, _mm256_cmp_pd(_mm256_loadu_pd(bminz + i), a_maxz, _CMP_GT_OQ));
                disjoint_bits |= static_cast<uint32_t>(_mm256_movemask_pd(m)) << i;
            }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
            const float64x2_t a_minx = vdupq_n_f64(aminx);
            const float64x2_t a_miny = vdupq_n_f64(aminy);
            const float64x2_t a_minz = vdupq_n_f64(aminz);
            const float64x2_t a_maxx = vdupq_n_f64(amaxx);
            const float64x2_t a_maxy = vdupq_n_f64(amaxy);
            const float64x2_t a_maxz = vdupq_n_f64(amaxz);

            for (; i + 2 <= count; i += 2) {
                uint64x2_t m = vorrq_u64(vorrq_u64(vcgtq_f64(a_minx, vld1q_f64(bmaxx + i)),
                                                   vcgtq_f64(a_miny, vld1q_f64(bmaxy + i))),
                                         vcgtq_f64(a_minz, vld1q_f64(bmaxz + i)));
                m = vorrq_u64(m,
                              vorrq_u64(vcgtq_f64(vld1q_f64(bminx + i), a_maxx),
                                        vcgtq_f64(vld1q_f64(bminy + i), a_maxy)));
                m = vorrq_u64(m, vcgtq_f64(vld1q_f64(bminz + i), a_maxz));

                const uint32_t lo = static_cast<uint32_t>(vgetq_lane_u64(m, 0) & 1u);
                const uint32_t hi = static_cast<uint32_t>(vgetq_lane_u64(m, 1) & 1u);
                disjoint_bits |= (lo | (hi << 1)) << i;
            }
#endif
        } else if constexpr (std::is_same<T, float>::value) {
#if defined(__AVX512F__)
            const __m512 a_minx = _mm512_set1_ps(aminx);
            const __m512 a_miny = _mm512_set1_ps(aminy);
            const __m512 a_minz = _mm512_set1_ps(aminz);
            const __m512 a_maxx = _mm512_set1_ps(amaxx);
            const __m512 a_maxy = _mm512_set1_ps(amaxy);
            const __m512 a_maxz = _mm512_set1_ps(amaxz);

            for (; i + 16 <= count; i += 16) {
                const __mmask16 k = _mm512_cmp_ps_mask(a_minx, _mm512_loadu_ps(bmaxx + i), _CMP_GT_OQ) |
                                    _mm512_cmp_ps_mask(a_miny, _mm512_loadu_ps(bmaxy + i), _CMP_GT_OQ) |
                                    _mm512_cmp_ps_mask(a_minz, _mm512_loadu_ps(bmaxz + i), _CMP_GT_OQ) |
                                    _mm512_cmp_ps_mask(_mm512_loadu_ps(bminx + i), a_maxx, _CMP_GT_OQ) |
                                    _mm512_cmp_ps_mask(_mm512_loadu_ps(bminy + i), a_maxy, _CMP_GT_OQ) |
                                    _mm512_cmp_ps_mask(_mm512_loadu_ps(bminz + i), a_maxz, _CMP_GT_OQ);
                disjoint_bits |= static_cast<uint32_t>(k) << i;
            }
#elif defined(__AVX2__)
            const __m256 a_minx = _mm256_set1_ps(aminx);
            const __m256 a_miny = _mm256_set1_ps(aminy);
            const __m256 a_minz = _mm256_set1_ps(aminz);
            const __m256 a_maxx = _mm256_set1_ps(amaxx);
            const __m256 a_maxy = _mm256_set1_ps(amaxy);
            const __m256 a_maxz = _mm256_set1_ps(amaxz);

            for (; i + 8 <= count; i += 8) {
                __m256 m = _mm256_or_ps(_mm256_or_ps(_mm256_cmp_ps(a_minx, _mm256_loadu_ps(bmaxx + i), _CMP_GT_OQ),
                                                     _mm256_cmp_ps(a_miny, _mm256_loadu_ps(bmaxy + i), _CMP_GT_OQ)),
                                        _mm256_cmp_ps(a_minz, _mm256_loadu_ps(bmaxz + i), _CMP_GT_OQ));
                m = _mm256_or_ps(m,
                                 _mm256_or_ps(_mm256_cmp_ps(_mm256_loadu_ps(bminx + i), a_maxx, _CMP_GT_OQ),
                                              _mm256_cmp_ps(_mm256_loadu_ps(bminy + i), a_maxy, _CMP_GT_OQ)));
                m = _mm256_or_ps(m, _mm256_cmp_ps(_mm256_loadu_ps(bminz + i), a_maxz, _CMP_GT_OQ));
                disjoint_bits |= static_cast<uint32_t>(_mm256_movemask_ps(m)) << i;
            }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
            const float32x4_t a_minx = vdupq_n_f32(aminx);
            const float32x4_t a_miny = vdupq_n_f32(aminy);
            const float32x4_t a_minz = vdupq_n_f32(aminz);
            const float32x4_t a_maxx = vdupq_n_f32(amaxx);
            const float32x4_t a_maxy = vdupq_n_f32(amaxy);
            const float32x4_t a_maxz = vdupq_n_f32(amaxz);

            for (; i + 4 <= count; i += 4) {
                uint32x4_t m = vorrq_u32(vorrq_u32(vcgtq_f32(a_minx, vld1q_f32(bmaxx + i)),
                                                   vcgtq_f32(a_miny, vld1q_f32(bmaxy + i))),
                                         vcgtq_f32(a_minz, vld1q_f32(bmaxz + i)));
                m = vorrq_u32(m,
                              vorrq_u32(vcgtq_f32(vld1q_f32(bminx + i), a_maxx),
                                        vcgtq_f32(vld1q_f32(bminy + i), a_maxy)));
                m = vorrq_u32(m, vcgtq_f32(vld1q_f32(bminz + i), a_maxz));

                const uint32x4_t m_01 = vshrq_n_u32(m, 31);
                disjoint_bits |= (vgetq_lane_u32(m_01, 0) | (vgetq_lane_u32(m_01, 1) << 1) |
                                  (vgetq_lane_u32(m_01, 2) << 2) | (vgetq_lane_u32(m_01, 3) << 3))
                                 << i;
            }
#endif
        }

        for (; i < count; ++i) {
            const uint32_t d = disjoint<T>(
                aminx, aminy, aminz, amaxx, amaxy, amaxz, bminx[i], bminy[i], bminz[i], bmaxx[i], bmaxy[i], bmaxz[i]);
            disjoint_bits |= (d != 0 ? 1u : 0u) << i;
        }

        const uint32_t valid = count >= 32 ? ~uint32_t(0) : ((uint32_t(1) << count) - 1);
        return ~disjoint_bits & valid;
    }

    template <typename T>
    inline static void vaabb_disjoint_one_to_many(const T aminx,
                                                  const T aminy,
                                                  const T aminz,
                                                  const T amaxx,
                                                  const T amaxy,
                                                  const T amaxz,
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
            const __m512d a_minx = _mm512_set1_pd(aminx);
            const __m512d a_miny = _mm512_set1_pd(aminy);
            const __m512d a_minz = _mm512_set1_pd(aminz);
            const __m512d a_maxx = _mm512_set1_pd(amaxx);
            const __m512d a_maxy = _mm512_set1_pd(amaxy);
            const __m512d a_maxz = _mm512_set1_pd(amaxz);

            for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 8) {
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

                const __m512i k_as_epi64 = _mm512_movm_epi64(k);
                const __m512i k_01 = _mm512_srli_epi64(k_as_epi64, 63);
                alignas(64) uint64_t tmp[8];
                _mm512_storeu_si512((__m512i*)tmp, k_01);
                for (int lane = 0; lane < 8; ++lane) {
                    mask[i + lane] = static_cast<uint32_t>(tmp[lane]);
                }
            }
            return;
#elif defined(__AVX2__)
            const __m256d a_minx = _mm256_set1_pd(aminx);
            const __m256d a_miny = _mm256_set1_pd(aminy);
            const __m256d a_minz = _mm256_set1_pd(aminz);
            const __m256d a_maxx = _mm256_set1_pd(amaxx);
            const __m256d a_maxy = _mm256_set1_pd(amaxy);
            const __m256d a_maxz = _mm256_set1_pd(amaxz);

            for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 4) {
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
            const float64x2_t a_minx = vdupq_n_f64(aminx);
            const float64x2_t a_miny = vdupq_n_f64(aminy);
            const float64x2_t a_minz = vdupq_n_f64(aminz);
            const float64x2_t a_maxx = vdupq_n_f64(amaxx);
            const float64x2_t a_maxy = vdupq_n_f64(amaxy);
            const float64x2_t a_maxz = vdupq_n_f64(amaxz);

            for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 2) {
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
        } else if constexpr (std::is_same<T, float>::value) {
#if defined(__AVX512F__)
            const __m512 a_minx = _mm512_set1_ps(aminx);
            const __m512 a_miny = _mm512_set1_ps(aminy);
            const __m512 a_minz = _mm512_set1_ps(aminz);
            const __m512 a_maxx = _mm512_set1_ps(amaxx);
            const __m512 a_maxy = _mm512_set1_ps(amaxy);
            const __m512 a_maxz = _mm512_set1_ps(amaxz);

            for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 16) {
                const __m512 b_minx = _mm512_loadu_ps(bminx + i);
                const __m512 b_miny = _mm512_loadu_ps(bminy + i);
                const __m512 b_minz = _mm512_loadu_ps(bminz + i);
                const __m512 b_maxx = _mm512_loadu_ps(bmaxx + i);
                const __m512 b_maxy = _mm512_loadu_ps(bmaxy + i);
                const __m512 b_maxz = _mm512_loadu_ps(bmaxz + i);

                __mmask16 k =
                    _mm512_cmp_ps_mask(a_minx, b_maxx, _CMP_GT_OQ) | _mm512_cmp_ps_mask(a_miny, b_maxy, _CMP_GT_OQ) |
                    _mm512_cmp_ps_mask(a_minz, b_maxz, _CMP_GT_OQ) | _mm512_cmp_ps_mask(b_minx, a_maxx, _CMP_GT_OQ) |
                    _mm512_cmp_ps_mask(b_miny, a_maxy, _CMP_GT_OQ) | _mm512_cmp_ps_mask(b_minz, a_maxz, _CMP_GT_OQ);

                const __m512i k_as_epi32 = _mm512_movm_epi32(k);
                const __m512i k_01 = _mm512_srli_epi32(k_as_epi32, 31);
                _mm512_storeu_si512((__m512i*)(mask + i), k_01);
            }
            return;
#elif defined(__AVX2__)
            const __m256 a_minx = _mm256_set1_ps(aminx);
            const __m256 a_miny = _mm256_set1_ps(aminy);
            const __m256 a_minz = _mm256_set1_ps(aminz);
            const __m256 a_maxx = _mm256_set1_ps(amaxx);
            const __m256 a_maxy = _mm256_set1_ps(amaxy);
            const __m256 a_maxz = _mm256_set1_ps(amaxz);

            for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 8) {
                const __m256 b_minx = _mm256_loadu_ps(bminx + i);
                const __m256 b_miny = _mm256_loadu_ps(bminy + i);
                const __m256 b_minz = _mm256_loadu_ps(bminz + i);
                const __m256 b_maxx = _mm256_loadu_ps(bmaxx + i);
                const __m256 b_maxy = _mm256_loadu_ps(bmaxy + i);
                const __m256 b_maxz = _mm256_loadu_ps(bmaxz + i);

                __m256 m = _mm256_or_ps(
                    _mm256_or_ps(_mm256_cmp_ps(a_minx, b_maxx, _CMP_GT_OQ), _mm256_cmp_ps(a_miny, b_maxy, _CMP_GT_OQ)),
                    _mm256_cmp_ps(a_minz, b_maxz, _CMP_GT_OQ));
                m = _mm256_or_ps(
                    m,
                    _mm256_or_ps(_mm256_cmp_ps(b_minx, a_maxx, _CMP_GT_OQ), _mm256_cmp_ps(b_miny, a_maxy, _CMP_GT_OQ)));
                m = _mm256_or_ps(m, _mm256_cmp_ps(b_minz, a_maxz, _CMP_GT_OQ));

                const __m256i m_i = _mm256_castps_si256(m);
                const __m256i m_01 = _mm256_srli_epi32(m_i, 31);
                _mm256_storeu_si256((__m256i*)(mask + i), m_01);
            }
            return;
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
            const float32x4_t a_minx = vdupq_n_f32(aminx);
            const float32x4_t a_miny = vdupq_n_f32(aminy);
            const float32x4_t a_minz = vdupq_n_f32(aminz);
            const float32x4_t a_maxx = vdupq_n_f32(amaxx);
            const float32x4_t a_maxy = vdupq_n_f32(amaxy);
            const float32x4_t a_maxz = vdupq_n_f32(amaxz);

            for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i += 4) {
                const float32x4_t b_minx = vld1q_f32(bminx + i);
                const float32x4_t b_miny = vld1q_f32(bminy + i);
                const float32x4_t b_minz = vld1q_f32(bminz + i);
                const float32x4_t b_maxx = vld1q_f32(bmaxx + i);
                const float32x4_t b_maxy = vld1q_f32(bmaxy + i);
                const float32x4_t b_maxz = vld1q_f32(bmaxz + i);

                uint32x4_t m = vorrq_u32(vorrq_u32(vcgtq_f32(a_minx, b_maxx), vcgtq_f32(a_miny, b_maxy)),
                                         vcgtq_f32(a_minz, b_maxz));
                m = vorrq_u32(m, vorrq_u32(vcgtq_f32(b_minx, a_maxx), vcgtq_f32(b_miny, a_maxy)));
                m = vorrq_u32(m, vcgtq_f32(b_minz, a_maxz));

                const uint32x4_t m_01 = vshrq_n_u32(m, 31);
                vst1q_u32(mask + i, m_01);
            }
            return;
#endif
        }
#pragma omp simd
        for (int i = 0; i < SCCD_AABB_DISJOINT_CHUNK_SIZE; i++) {
            mask[i] = disjoint<T>(
                aminx, aminy, aminz, amaxx, amaxy, amaxz, bminx[i], bminy[i], bminz[i], bmaxx[i], bmaxy[i], bmaxz[i]);
        }
    }

    template <typename idx_t, typename geom_t, typename aabb_t>
    void compute_aabbs(const int nxe,
                       const ptrdiff_t n_elements,
                       const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                       const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                       const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                       aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                       const BoxRounding rounding = BoxRounding::Exact) {
        // One sweep per dimension, reducing over the element's vertices in
        // registers. The previous shape ran nxe+1 separate parallel passes per
        // dimension, each one a full read-modify-write of the min/max rows.
        for (int d = 0; d < SCCD_DIM; d++) {
            aabb_t *const SCCD_RESTRICT amin = aabbs[d];
            aabb_t *const SCCD_RESTRICT amax = aabbs[SCCD_DIM + d];
            const geom_t *const SCCD_RESTRICT p0d = points0[d];
            const geom_t *const SCCD_RESTRICT p1d = points1[d];

            sccd::parallel_for_br(0, n_elements, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                for (ptrdiff_t i = rbegin; i < rend; i++) {
                    aabb_t e_min = std::numeric_limits<aabb_t>::max();
                    aabb_t e_max = std::numeric_limits<aabb_t>::lowest();

                    for (int v = 0; v < nxe; v++) {
                        const idx_t ii = elements[v][i];
                        const geom_t p0 = p0d[ii];
                        const geom_t p1 = p1d[ii];
                        e_min = std::min(e_min, aabb_t(std::min(p0, p1)));
                        e_max = std::max(e_max, aabb_t(std::max(p0, p1)));
                    }

                    if (rounding == BoxRounding::OutwardUlp) {
                        e_min = nextafter_down<aabb_t>(e_min);
                        e_max = nextafter_up<aabb_t>(e_max);
                    }

                    amin[i] = e_min;
                    amax[i] = e_max;
                }
            });
        }
    }

    template <typename geom_t, typename aabb_t>
    void compute_aabbs(const ptrdiff_t n_nodes,
                       const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                       const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                       aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                       const BoxRounding rounding = BoxRounding::Exact) {
        for (int d = 0; d < SCCD_DIM; d++) {
            aabb_t *const SCCD_RESTRICT amin = aabbs[d];
            aabb_t *const SCCD_RESTRICT amax = aabbs[SCCD_DIM + d];
            const geom_t *const SCCD_RESTRICT p0d = points0[d];
            const geom_t *const SCCD_RESTRICT p1d = points1[d];

            sccd::parallel_for_br(0, n_nodes, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                for (ptrdiff_t i = rbegin; i < rend; i++) {
                    aabb_t p_min = std::min(p0d[i], p1d[i]);
                    aabb_t p_max = std::max(p0d[i], p1d[i]);
                    if (rounding == BoxRounding::OutwardUlp) {
                        p_min = nextafter_down<aabb_t>(p_min);
                        p_max = nextafter_up<aabb_t>(p_max);
                    }
                    amin[i] = p_min;
                    amax[i] = p_max;
                }
            });
        }
    }

    template <typename idx_t, typename geom_t, typename disp_t, typename aabb_t>
    void compute_aabbs(const int nxe,
                       const ptrdiff_t n_elements,
                       const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                       const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points,
                       const ptrdiff_t stride_disp,
                       const disp_t* const SCCD_RESTRICT* const SCCD_RESTRICT disp0,
                       const disp_t* const SCCD_RESTRICT* const SCCD_RESTRICT disp1,
                       aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                       const BoxRounding rounding = BoxRounding::Exact) {
        for (int d = 0; d < SCCD_DIM; d++) {
            aabb_t *const SCCD_RESTRICT amin = aabbs[d];
            aabb_t *const SCCD_RESTRICT amax = aabbs[SCCD_DIM + d];
            const geom_t *const SCCD_RESTRICT pd = points[d];
            const disp_t *const SCCD_RESTRICT d0 = disp0[d];
            const disp_t *const SCCD_RESTRICT d1 = disp1[d];

            sccd::parallel_for_br(0, n_elements, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                for (ptrdiff_t i = rbegin; i < rend; i++) {
                    aabb_t e_min = std::numeric_limits<aabb_t>::max();
                    aabb_t e_max = std::numeric_limits<aabb_t>::lowest();

                    for (int v = 0; v < nxe; v++) {
                        const idx_t ii = elements[v][i];
                        const geom_t p = pd[ii];
                        const geom_t disp0_i = p + d0[ii * stride_disp];
                        const geom_t disp1_i = p + d1[ii * stride_disp];
                        e_min = std::min(e_min, aabb_t(std::min(disp0_i, disp1_i)));
                        e_max = std::max(e_max, aabb_t(std::max(disp0_i, disp1_i)));
                    }

                    if (rounding == BoxRounding::OutwardUlp) {
                        e_min = nextafter_down<aabb_t>(e_min);
                        e_max = nextafter_up<aabb_t>(e_max);
                    }

                    amin[i] = e_min;
                    amax[i] = e_max;
                }
            });
        }
    }

    template <typename geom_t, typename disp_t, typename aabb_t>
    void compute_aabbs(const ptrdiff_t n_nodes,
                       const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points,
                       const ptrdiff_t stride_disp,
                       const disp_t* const SCCD_RESTRICT* const SCCD_RESTRICT disp0,
                       const disp_t* const SCCD_RESTRICT* const SCCD_RESTRICT disp1,
                       aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                       const BoxRounding rounding = BoxRounding::Exact) {
        for (int d = 0; d < SCCD_DIM; d++) {
            aabb_t *const SCCD_RESTRICT amin = aabbs[d];
            aabb_t *const SCCD_RESTRICT amax = aabbs[SCCD_DIM + d];
            const geom_t *const SCCD_RESTRICT pd = points[d];
            const disp_t *const SCCD_RESTRICT d0 = disp0[d];
            const disp_t *const SCCD_RESTRICT d1 = disp1[d];

            sccd::parallel_for_br(0, n_nodes, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                for (ptrdiff_t i = rbegin; i < rend; i++) {
                    const geom_t p = pd[i];
                    const geom_t disp0_i = p + d0[i * stride_disp];
                    const geom_t disp1_i = p + d1[i * stride_disp];
                    aabb_t p_min = std::min(disp0_i, disp1_i);
                    aabb_t p_max = std::max(disp0_i, disp1_i);
                    if (rounding == BoxRounding::OutwardUlp) {
                        p_min = nextafter_down<aabb_t>(p_min);
                        p_max = nextafter_up<aabb_t>(p_max);
                    }
                    amin[i] = p_min;
                    amax[i] = p_max;
                }
            });
        }
    }

    template <typename idx_t, typename geom_t, typename aabb_t>
    void compute_aabbs(const int nxe,
                       const ptrdiff_t n_elements,
                       const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                       const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points,
                       aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                       const BoxRounding rounding = BoxRounding::Exact) {
        for (int d = 0; d < SCCD_DIM; d++) {
            aabb_t *const SCCD_RESTRICT amin = aabbs[d];
            aabb_t *const SCCD_RESTRICT amax = aabbs[SCCD_DIM + d];
            const geom_t *const SCCD_RESTRICT pd = points[d];

            sccd::parallel_for_br(0, n_elements, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                for (ptrdiff_t i = rbegin; i < rend; i++) {
                    aabb_t e_min = std::numeric_limits<aabb_t>::max();
                    aabb_t e_max = std::numeric_limits<aabb_t>::lowest();

                    for (int v = 0; v < nxe; v++) {
                        const idx_t ii = elements[v][i];
                        const aabb_t p = aabb_t(pd[ii]);
                        e_min = std::min(e_min, p);
                        e_max = std::max(e_max, p);
                    }

                    if (rounding == BoxRounding::OutwardUlp) {
                        e_min = nextafter_down<aabb_t>(e_min);
                        e_max = nextafter_up<aabb_t>(e_max);
                    }

                    amin[i] = e_min;
                    amax[i] = e_max;
                }
            });
        }
    }
}  // namespace sccd

#endif  // SCCD_AABB_HPP
