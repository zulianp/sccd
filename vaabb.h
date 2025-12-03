#ifndef VAABB_H
#define VAABB_H

#include <stdlib.h>
#include <float.h>

#define SFEM_READ_ENV(name, conversion)                                        \
    do {                                                                       \
        char* var = getenv(#name);                                             \
        if (var) {                                                             \
            name = conversion(var);                                            \
        }                                                                      \
    } while (0)

#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#define AABB_DISJOINT_CHUNK_SIZE 32
#define AABB_DISJOINT_NOVECTORIZE_THRESHOLD 16

#ifndef _WIN32
#define SFEM_RESTRICT __restrict__
#else
#define SFEM_RESTRICT __restrict
#endif

namespace sccd {

/// Geometry scalar type used for coordinates and AABB values.
using geom_t = float;
/// Integer type used for element/vertex indices.
using idx_t = int;
/// Integer type used for counts and prefix sums.
using count_t = int;

/**
 * \brief Return the next representable value toward +infinity.
 * \param x Input value.
 * \return Next representable float greater than or equal to x.
 */
inline static geom_t nextafter_up(const geom_t x)
{
    if constexpr (std::is_same<geom_t, float>::value) {
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
inline static geom_t nextafter_down(const geom_t x)
{
    if constexpr (std::is_same<geom_t, float>::value) {
        return nextafterf(x, -FLT_MAX);
    } else {
        return nextafter(x, -DBL_MAX);
    }
}

inline static uint32_t disjoint(
    const geom_t aminx,
    const geom_t aminy,
    const geom_t aminz,
    const geom_t amaxx,
    const geom_t amaxy,
    const geom_t amaxz,
    const geom_t bminx,
    const geom_t bminy,
    const geom_t bminz,
    const geom_t bmaxx,
    const geom_t bmaxy,
    const geom_t bmaxz)
{
    return aminx > bmaxx | aminy > bmaxy | aminz > bmaxz | bminx > amaxx
        | bminy > amaxy | bminz > amaxz;
}

inline static void vdisjoint(
    const geom_t* const SFEM_RESTRICT aminx,
    const geom_t* const SFEM_RESTRICT aminy,
    const geom_t* const SFEM_RESTRICT aminz,
    const geom_t* const SFEM_RESTRICT amaxx,
    const geom_t* const SFEM_RESTRICT amaxy,
    const geom_t* const SFEM_RESTRICT amaxz,
    const geom_t* const SFEM_RESTRICT bminx,
    const geom_t* const SFEM_RESTRICT bminy,
    const geom_t* const SFEM_RESTRICT bminz,
    const geom_t* const SFEM_RESTRICT bmaxx,
    const geom_t* const SFEM_RESTRICT bmaxy,
    const geom_t* const SFEM_RESTRICT bmaxz,
    uint32_t* SFEM_RESTRICT mask)
{
#if defined(__AVX512F__)
    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i += 16) {
        const __m512 a_minx = _mm512_loadu_ps(aminx + i);
        const __m512 a_miny = _mm512_loadu_ps(aminy + i);
        const __m512 a_minz = _mm512_loadu_ps(aminz + i);
        const __m512 a_maxx = _mm512_loadu_ps(amaxx + i);
        const __m512 a_maxy = _mm512_loadu_ps(amaxy + i);
        const __m512 a_maxz = _mm512_loadu_ps(amaxz + i);

        const __m512 b_minx = _mm512_loadu_ps(bminx + i);
        const __m512 b_miny = _mm512_loadu_ps(bminy + i);
        const __m512 b_minz = _mm512_loadu_ps(bminz + i);
        const __m512 b_maxx = _mm512_loadu_ps(bmaxx + i);
        const __m512 b_maxy = _mm512_loadu_ps(bmaxy + i);
        const __m512 b_maxz = _mm512_loadu_ps(bmaxz + i);

        __mmask16 k = _mm512_cmp_ps_mask(a_minx, b_maxx, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(a_miny, b_maxy, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(a_minz, b_maxz, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(b_minx, a_maxx, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(b_miny, a_maxy, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(b_minz, a_maxz, _CMP_GT_OQ);

        __m512i k_as_epi32 = _mm512_movm_epi32(k);
        __m512i k_01 = _mm512_srli_epi32(k_as_epi32, 31);
        _mm512_storeu_si512((__m512i*)(mask + i), k_01);
    }
#elif defined(__AVX2__)
    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i += 8) {
        const __m256 a_minx = _mm256_loadu_ps(aminx + i);
        const __m256 a_miny = _mm256_loadu_ps(aminy + i);
        const __m256 a_minz = _mm256_loadu_ps(aminz + i);
        const __m256 a_maxx = _mm256_loadu_ps(amaxx + i);
        const __m256 a_maxy = _mm256_loadu_ps(amaxy + i);
        const __m256 a_maxz = _mm256_loadu_ps(amaxz + i);

        const __m256 b_minx = _mm256_loadu_ps(bminx + i);
        const __m256 b_miny = _mm256_loadu_ps(bminy + i);
        const __m256 b_minz = _mm256_loadu_ps(bminz + i);
        const __m256 b_maxx = _mm256_loadu_ps(bmaxx + i);
        const __m256 b_maxy = _mm256_loadu_ps(bmaxy + i);
        const __m256 b_maxz = _mm256_loadu_ps(bmaxz + i);

        __m256 m = _mm256_or_ps(
            _mm256_or_ps(
                _mm256_cmp_ps(a_minx, b_maxx, _CMP_GT_OQ),
                _mm256_cmp_ps(a_miny, b_maxy, _CMP_GT_OQ)),
            _mm256_cmp_ps(a_minz, b_maxz, _CMP_GT_OQ));
        m = _mm256_or_ps(
            m,
            _mm256_or_ps(
                _mm256_cmp_ps(b_minx, a_maxx, _CMP_GT_OQ),
                _mm256_cmp_ps(b_miny, a_maxy, _CMP_GT_OQ)));
        m = _mm256_or_ps(m, _mm256_cmp_ps(b_minz, a_maxz, _CMP_GT_OQ));

        const __m256i m_i = _mm256_castps_si256(m);
        const __m256i m_01 = _mm256_srli_epi32(m_i, 31);
        _mm256_storeu_si256((__m256i*)(mask + i), m_01);
    }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i += 4) {
        const float32x4_t a_minx = vld1q_f32(aminx + i);
        const float32x4_t a_miny = vld1q_f32(aminy + i);
        const float32x4_t a_minz = vld1q_f32(aminz + i);
        const float32x4_t a_maxx = vld1q_f32(amaxx + i);
        const float32x4_t a_maxy = vld1q_f32(amaxy + i);
        const float32x4_t a_maxz = vld1q_f32(amaxz + i);

        const float32x4_t b_minx = vld1q_f32(bminx + i);
        const float32x4_t b_miny = vld1q_f32(bminy + i);
        const float32x4_t b_minz = vld1q_f32(bminz + i);
        const float32x4_t b_maxx = vld1q_f32(bmaxx + i);
        const float32x4_t b_maxy = vld1q_f32(bmaxy + i);
        const float32x4_t b_maxz = vld1q_f32(bmaxz + i);

        uint32x4_t m = vorrq_u32(
            vorrq_u32(vcgtq_f32(a_minx, b_maxx), vcgtq_f32(a_miny, b_maxy)),
            vcgtq_f32(a_minz, b_maxz));
        m = vorrq_u32(
            m, vorrq_u32(vcgtq_f32(b_minx, a_maxx), vcgtq_f32(b_miny, a_maxy)));
        m = vorrq_u32(m, vcgtq_f32(b_minz, a_maxz));

        const uint32x4_t m_01 = vshrq_n_u32(m, 31);
        vst1q_u32(mask + i, m_01);
    }
#else
#pragma omp simd aligned(                                                      \
    aminx, aminy, aminz, amaxx, amaxy, amaxz, bminx, bminy, bminz, bmaxx,      \
    bmaxy, bmaxz, mask : 64)
    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i++) {
        mask[i] = disjoint(
            aminx[i], aminy[i], aminz[i], amaxx[i], amaxy[i], amaxz[i],
            bminx[i], bminy[i], bminz[i], bmaxx[i], bmaxy[i], bmaxz[i]);
    }
#endif
}

/**
 * \brief Broadcast AABB at \p fi into SoA buffers sized for SIMD chunking.
 * \param aabbs SoA AABB arrays [6][...].
 * \param fi Index of the AABB to broadcast.
 * \param A_minx..A_maxz Output arrays of length AABB_DISJOINT_CHUNK_SIZE.
 */
inline static void vaabb_broadcast(
    geom_t** const SFEM_RESTRICT aabbs,
    const size_t fi,
    geom_t* const SFEM_RESTRICT A_minx,
    geom_t* const SFEM_RESTRICT A_miny,
    geom_t* const SFEM_RESTRICT A_minz,
    geom_t* const SFEM_RESTRICT A_maxx,
    geom_t* const SFEM_RESTRICT A_maxy,
    geom_t* const SFEM_RESTRICT A_maxz)
{
    const geom_t aminx = aabbs[0][fi];
    const geom_t aminy = aabbs[1][fi];
    const geom_t aminz = aabbs[2][fi];
    const geom_t amaxx = aabbs[3][fi];
    const geom_t amaxy = aabbs[4][fi];
    const geom_t amaxz = aabbs[5][fi];
    for (int k = 0; k < AABB_DISJOINT_CHUNK_SIZE; ++k) {
        A_minx[k] = aminx;
        A_miny[k] = aminy;
        A_minz[k] = aminz;
        A_maxx[k] = amaxx;
        A_maxy[k] = amaxy;
        A_maxz[k] = amaxz;
    }
}

inline static void vdisjoint_one_to_many(
    const geom_t aminx,
    const geom_t aminy,
    const geom_t aminz,
    const geom_t amaxx,
    const geom_t amaxy,
    const geom_t amaxz,
    const geom_t* const SFEM_RESTRICT bminx,
    const geom_t* const SFEM_RESTRICT bminy,
    const geom_t* const SFEM_RESTRICT bminz,
    const geom_t* const SFEM_RESTRICT bmaxx,
    const geom_t* const SFEM_RESTRICT bmaxy,
    const geom_t* const SFEM_RESTRICT bmaxz,
    uint32_t* SFEM_RESTRICT mask)
{

#if defined(__AVX512F__)
    const __m512 a_minx = _mm512_set1_ps(aminx);
    const __m512 a_miny = _mm512_set1_ps(aminy);
    const __m512 a_minz = _mm512_set1_ps(aminz);
    const __m512 a_maxx = _mm512_set1_ps(amaxx);
    const __m512 a_maxy = _mm512_set1_ps(amaxy);
    const __m512 a_maxz = _mm512_set1_ps(amaxz);

    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i += 16) {
        const __m512 b_minx = _mm512_loadu_ps(bminx + i);
        const __m512 b_miny = _mm512_loadu_ps(bminy + i);
        const __m512 b_minz = _mm512_loadu_ps(bminz + i);
        const __m512 b_maxx = _mm512_loadu_ps(bmaxx + i);
        const __m512 b_maxy = _mm512_loadu_ps(bmaxy + i);
        const __m512 b_maxz = _mm512_loadu_ps(bmaxz + i);

        __mmask16 k = _mm512_cmp_ps_mask(a_minx, b_maxx, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(a_miny, b_maxy, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(a_minz, b_maxz, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(b_minx, a_maxx, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(b_miny, a_maxy, _CMP_GT_OQ)
            | _mm512_cmp_ps_mask(b_minz, a_maxz, _CMP_GT_OQ);

        const __m512i k_as_epi32 = _mm512_movm_epi32(k);
        const __m512i k_01 = _mm512_srli_epi32(k_as_epi32, 31);
        _mm512_storeu_si512((__m512i*)(mask + i), k_01);
    }
#elif defined(__AVX2__)
    const __m256 a_minx = _mm256_set1_ps(aminx);
    const __m256 a_miny = _mm256_set1_ps(aminy);
    const __m256 a_minz = _mm256_set1_ps(aminz);
    const __m256 a_maxx = _mm256_set1_ps(amaxx);
    const __m256 a_maxy = _mm256_set1_ps(amaxy);
    const __m256 a_maxz = _mm256_set1_ps(amaxz);

    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i += 8) {
        const __m256 b_minx = _mm256_loadu_ps(bminx + i);
        const __m256 b_miny = _mm256_loadu_ps(bminy + i);
        const __m256 b_minz = _mm256_loadu_ps(bminz + i);
        const __m256 b_maxx = _mm256_loadu_ps(bmaxx + i);
        const __m256 b_maxy = _mm256_loadu_ps(bmaxy + i);
        const __m256 b_maxz = _mm256_loadu_ps(bmaxz + i);

        __m256 m = _mm256_or_ps(
            _mm256_or_ps(
                _mm256_cmp_ps(a_minx, b_maxx, _CMP_GT_OQ),
                _mm256_cmp_ps(a_miny, b_maxy, _CMP_GT_OQ)),
            _mm256_cmp_ps(a_minz, b_maxz, _CMP_GT_OQ));
        m = _mm256_or_ps(
            m,
            _mm256_or_ps(
                _mm256_cmp_ps(b_minx, a_maxx, _CMP_GT_OQ),
                _mm256_cmp_ps(b_miny, a_maxy, _CMP_GT_OQ)));
        m = _mm256_or_ps(m, _mm256_cmp_ps(b_minz, a_maxz, _CMP_GT_OQ));

        const __m256i m_i = _mm256_castps_si256(m);
        const __m256i m_01 = _mm256_srli_epi32(m_i, 31);
        _mm256_storeu_si256((__m256i*)(mask + i), m_01);
    }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    const float32x4_t a_minx = vdupq_n_f32(aminx);
    const float32x4_t a_miny = vdupq_n_f32(aminy);
    const float32x4_t a_minz = vdupq_n_f32(aminz);
    const float32x4_t a_maxx = vdupq_n_f32(amaxx);
    const float32x4_t a_maxy = vdupq_n_f32(amaxy);
    const float32x4_t a_maxz = vdupq_n_f32(amaxz);

    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i += 4) {
        const float32x4_t b_minx = vld1q_f32(bminx + i);
        const float32x4_t b_miny = vld1q_f32(bminy + i);
        const float32x4_t b_minz = vld1q_f32(bminz + i);
        const float32x4_t b_maxx = vld1q_f32(bmaxx + i);
        const float32x4_t b_maxy = vld1q_f32(bmaxy + i);
        const float32x4_t b_maxz = vld1q_f32(bmaxz + i);

        uint32x4_t m = vorrq_u32(
            vorrq_u32(vcgtq_f32(a_minx, b_maxx), vcgtq_f32(a_miny, b_maxy)),
            vcgtq_f32(a_minz, b_maxz));
        m = vorrq_u32(
            m, vorrq_u32(vcgtq_f32(b_minx, a_maxx), vcgtq_f32(b_miny, a_maxy)));
        m = vorrq_u32(m, vcgtq_f32(b_minz, a_maxz));

        const uint32x4_t m_01 = vshrq_n_u32(m, 31);
        vst1q_u32(mask + i, m_01);
    }
#else
#pragma omp simd
    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i++) {
        mask[i] = disjoint(
            aminx, aminy, aminz, amaxx, amaxy, amaxz, bminx[i], bminy[i],
            bminz[i], bmaxx[i], bmaxy[i], bmaxz[i]);
    }
#endif
}

} // namespace sccd

#endif // VAABB_H
