#ifndef VAABB_H
#define VAABB_H

#include <stdint.h>

#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif

#define AABB_DISJOINT_CHUNK_SIZE 32

#ifndef SFEM_RESTRICT // defined in sccd.hpp
#define SFEM_RESTRICT __restrict
#endif

typedef float geom_t;

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

void vdisjoint(
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
        aminx, aminy, aminz, amaxx, amaxy, amaxz, bminx, bminy, bminz, bmaxx,  \
            bmaxy, bmaxz, mask : 64)
    for (int i = 0; i < AABB_DISJOINT_CHUNK_SIZE; i++) {
        mask[i] = disjoint(
            aminx[i], aminy[i], aminz[i], amaxx[i], amaxy[i], amaxz[i],
            bminx[i], bminy[i], bminz[i], bmaxx[i], bmaxy[i], bmaxz[i]);
    }
#endif
}

// Return bitmask of disjoint lanes (1 = disjoint), only low `lanes` bits used
// static inline uint32_t vdisjoint_mask(
//     const geom_t* SFEM_RESTRICT aminx,
//     const geom_t* SFEM_RESTRICT aminy,
//     const geom_t* SFEM_RESTRICT aminz,
//     const geom_t* SFEM_RESTRICT amaxx,
//     const geom_t* SFEM_RESTRICT amaxy,
//     const geom_t* SFEM_RESTRICT amaxz,
//     const geom_t* SFEM_RESTRICT bminx,
//     const geom_t* SFEM_RESTRICT bminy,
//     const geom_t* SFEM_RESTRICT bminz,
//     const geom_t* SFEM_RESTRICT bmaxx,
//     const geom_t* SFEM_RESTRICT bmaxy,
//     const geom_t* SFEM_RESTRICT bmaxz,
//     size_t lanes)
// {
// #if defined(__AVX512F__)
//   const __m512 a_minx = _mm512_loadu_ps(aminx);
//   const __m512 a_miny = _mm512_loadu_ps(aminy);
//   const __m512 a_minz = _mm512_loadu_ps(aminz);
//   const __m512 a_maxx = _mm512_loadu_ps(amaxx);
//   const __m512 a_maxy = _mm512_loadu_ps(amaxy);
//   const __m512 a_maxz = _mm512_loadu_ps(amaxz);
//   const __m512 b_minx = _mm512_loadu_ps(bminx);
//   const __m512 b_miny = _mm512_loadu_ps(bminy);
//   const __m512 b_minz = _mm512_loadu_ps(bminz);
//   const __m512 b_maxx = _mm512_loadu_ps(bmaxx);
//   const __m512 b_maxy = _mm512_loadu_ps(bmaxy);
//   const __m512 b_maxz = _mm512_loadu_ps(bmaxz);
//   const __mmask16 k =
//       _mm512_cmp_ps_mask(a_minx, b_maxx, _CMP_GT_OQ)
//     | _mm512_cmp_ps_mask(a_miny, b_maxy, _CMP_GT_OQ)
//     | _mm512_cmp_ps_mask(a_minz, b_maxz, _CMP_GT_OQ)
//     | _mm512_cmp_ps_mask(b_minx, a_maxx, _CMP_GT_OQ)
//     | _mm512_cmp_ps_mask(b_miny, a_maxy, _CMP_GT_OQ)
//     | _mm512_cmp_ps_mask(b_minz, a_maxz, _CMP_GT_OQ);
//   return (uint32_t)k & ((lanes < 32) ? ((1u << lanes) - 1u) : 0xFFFFFFFFu);
// #elif defined(__AVX2__)
//   uint32_t bits = 0;
//   size_t processed = 0;
//   while (processed < lanes) {
//     const size_t n = (lanes - processed) >= 8 ? 8 : (lanes - processed);
//     const __m256 a_minx = _mm256_loadu_ps(aminx + processed);
//     const __m256 a_miny = _mm256_loadu_ps(aminy + processed);
//     const __m256 a_minz = _mm256_loadu_ps(aminz + processed);
//     const __m256 a_maxx = _mm256_loadu_ps(amaxx + processed);
//     const __m256 a_maxy = _mm256_loadu_ps(amaxy + processed);
//     const __m256 a_maxz = _mm256_loadu_ps(amaxz + processed);
//     const __m256 b_minx = _mm256_loadu_ps(bminx + processed);
//     const __m256 b_miny = _mm256_loadu_ps(bminy + processed);
//     const __m256 b_minz = _mm256_loadu_ps(bminz + processed);
//     const __m256 b_maxx = _mm256_loadu_ps(bmaxx + processed);
//     const __m256 b_maxy = _mm256_loadu_ps(bmaxy + processed);
//     const __m256 b_maxz = _mm256_loadu_ps(bmaxz + processed);
//     __m256 m =
//         _mm256_or_ps(
//             _mm256_or_ps(_mm256_cmp_ps(a_minx, b_maxx, _CMP_GT_OQ),
//                          _mm256_cmp_ps(a_miny, b_maxy, _CMP_GT_OQ)),
//             _mm256_cmp_ps(a_minz, b_maxz, _CMP_GT_OQ));
//     m = _mm256_or_ps(
//         m, _mm256_or_ps(_mm256_cmp_ps(b_minx, a_maxx, _CMP_GT_OQ),
//                         _mm256_cmp_ps(b_miny, a_maxy, _CMP_GT_OQ)));
//     m = _mm256_or_ps(m, _mm256_cmp_ps(b_minz, a_maxz, _CMP_GT_OQ));
//     uint32_t mbits = (uint32_t)_mm256_movemask_ps(m) & ((1u << n) - 1u);
//     bits |= (mbits << processed);
//     processed += n;
//   }
//   return bits;
// #else
//   uint32_t bits = 0;
//   for (size_t i = 0; i < lanes; ++i) {
//     const uint32_t d = disjoint(
//         aminx[i], aminy[i], aminz[i], amaxx[i], amaxy[i], amaxz[i],
//         bminx[i], bminy[i], bminz[i], bmaxx[i], bmaxy[i], bmaxz[i]);
//     bits |= (d ? (1u << i) : 0u);
//   }
//   return bits;
// #endif
// }

// Return a bitmask (LSB -> lane 0) of lanes that SHARE at least one vertex with
// ev0/ev1 Lanes >= chunk_len are ignored and masked to 0.
// static inline uint32_t vshare_edge_shares_vertex_mask(
//     const int* SFEM_RESTRICT idx,
//     const size_t base_index,
//     const size_t chunk_len,
//     const int stride,
//     const int* SFEM_RESTRICT elements0,
//     const int* SFEM_RESTRICT elements1,
//     const int ev0,
//     const int ev1)
// {
// #if defined(__AVX512F__)
//     const size_t lanes = chunk_len > 16 ? 16 : chunk_len;
//     alignas(64) int offsets[16];
//     for (size_t l = 0; l < 16; ++l) {
//         offsets[l] = (l < lanes) ? idx[base_index + l] * stride : 0;
//     }
//     const __m512i off = _mm512_loadu_si512((const void*)offsets);
//     const __m512i sec0 = _mm512_i32gather_epi32(off, (const int*)elements0, 4);
//     const __m512i sec1 = _mm512_i32gather_epi32(off, (const int*)elements1, 4);
//     const __m512i ev0v = _mm512_set1_epi32(ev0);
//     const __m512i ev1v = _mm512_set1_epi32(ev1);
//     const __mmask16 m0 = _mm512_cmpeq_epi32_mask(sec0, ev0v);
//     const __mmask16 m1 = _mm512_cmpeq_epi32_mask(sec0, ev1v);
//     const __mmask16 m2 = _mm512_cmpeq_epi32_mask(sec1, ev0v);
//     const __mmask16 m3 = _mm512_cmpeq_epi32_mask(sec1, ev1v);
//     const uint32_t share_bits = (uint32_t)(m0 | m1 | m2 | m3)
//         & ((lanes < 32) ? ((1u << lanes) - 1u) : 0xFFFFFFFFu);
//     return share_bits;
// #elif defined(__AVX2__)
//     uint32_t share_bits = 0;
//     size_t processed = 0;
//     while (processed < chunk_len) {
//         const size_t lanes =
//             ((chunk_len - processed) >= 8) ? 8 : (chunk_len - processed);
//         alignas(32) int offsets8[8];
//         for (size_t l = 0; l < lanes; ++l) {
//             const size_t lane = processed + l;
//             offsets8[l] = idx[base_index + lane] * stride;
//         }
//         for (size_t l = lanes; l < 8; ++l)
//             offsets8[l] = 0;

//         const __m256i off = _mm256_loadu_si256((const __m256i*)offsets8);
//         const __m256i sec0 =
//             _mm256_i32gather_epi32((const int*)elements0, off, 4);
//         const __m256i sec1 =
//             _mm256_i32gather_epi32((const int*)elements1, off, 4);
//         const __m256i ev0v = _mm256_set1_epi32(ev0);
//         const __m256i ev1v = _mm256_set1_epi32(ev1);

//         const __m256i e0eq0 = _mm256_cmpeq_epi32(sec0, ev0v);
//         const __m256i e0eq1 = _mm256_cmpeq_epi32(sec0, ev1v);
//         const __m256i e1eq0 = _mm256_cmpeq_epi32(sec1, ev0v);
//         const __m256i e1eq1 = _mm256_cmpeq_epi32(sec1, ev1v);
//         const __m256i anyeq = _mm256_or_si256(
//             _mm256_or_si256(e0eq0, e0eq1), _mm256_or_si256(e1eq0, e1eq1));
//         const uint32_t half_bits =
//             (uint32_t)_mm256_movemask_ps(_mm256_castsi256_ps(anyeq))
//             & ((1u << lanes) - 1u);
//         share_bits |= (half_bits << processed);
//         processed += lanes;
//     }
//     return share_bits;
// #else
//     uint32_t share_bits = 0;
//     for (size_t lane = 0; lane < chunk_len; ++lane) {
//         const int jidx = idx[base_index + lane];
//         const int s0 = elements0[jidx * stride];
//         const int s1 = elements1[jidx * stride];
//         const bool share =
//             (s0 == ev0) | (s0 == ev1) | (s1 == ev0) | (s1 == ev1);
//         if (share)
//             share_bits |= (1u << lane);
//     }
//     return share_bits;
// #endif
// }

// Fill mask array (1 = shares a vertex with (ev0,ev1), 0 = no share)
// static inline void vshare_edge_shares_vertex(
//     const int* SFEM_RESTRICT idx,
//     const size_t base_index,
//     const size_t chunk_len,
//     const int stride,
//     const int* SFEM_RESTRICT elements0,
//     const int* SFEM_RESTRICT elements1,
//     const int ev0,
//     const int ev1,
//     uint32_t* SFEM_RESTRICT mask_out)
// {
// #if defined(__AVX512F__)
//   alignas(64) int offsets[16];
//   const size_t lanes = chunk_len > 16 ? 16 : chunk_len;
//   for (size_t l = 0; l < 16; ++l) {
//     offsets[l] = (l < lanes) ? idx[base_index + l] * stride : 0;
//   }
//   const __m512i off = _mm512_loadu_si512((const void*)offsets);
//   const __m512i sec0 = _mm512_i32gather_epi32(off, (const int*)elements0, 4);
//   const __m512i sec1 = _mm512_i32gather_epi32(off, (const int*)elements1, 4);
//   const __m512i ev0v = _mm512_set1_epi32(ev0);
//   const __m512i ev1v = _mm512_set1_epi32(ev1);
//   __mmask16 k =
//       _mm512_cmpeq_epi32_mask(sec0, ev0v)
//     | _mm512_cmpeq_epi32_mask(sec0, ev1v)
//     | _mm512_cmpeq_epi32_mask(sec1, ev0v)
//     | _mm512_cmpeq_epi32_mask(sec1, ev1v);
//   // Zero out lanes >= chunk_len
//   if (lanes < 16) {
//     const __mmask16 lane_mask = (lanes < 32) ? (__mmask16)((1u << lanes) - 1u) : (__mmask16)0xFFFF;
//     k &= lane_mask;
//   }
//   const __m512i k_i = _mm512_movm_epi32(k);
//   const __m512i k_01 = _mm512_srli_epi32(k_i, 31);
//   _mm512_storeu_si512((__m512i*)mask_out, k_01);
// #elif defined(__AVX2__)
//   size_t processed = 0;
//   while (processed < chunk_len) {
//     const size_t lanes = ((chunk_len - processed) >= 8) ? 8 : (chunk_len - processed);
//     alignas(32) int offsets8[8];
//     for (size_t l = 0; l < lanes; ++l) {
//       const size_t lane = processed + l;
//       offsets8[l] = idx[base_index + lane] * stride;
//     }
//     for (size_t l = lanes; l < 8; ++l) offsets8[l] = 0;

//     const __m256i off = _mm256_loadu_si256((const __m256i*)offsets8);
//     const __m256i sec0 = _mm256_i32gather_epi32((const int*)elements0, off, 4);
//     const __m256i sec1 = _mm256_i32gather_epi32((const int*)elements1, off, 4);
//     const __m256i ev0v = _mm256_set1_epi32(ev0);
//     const __m256i ev1v = _mm256_set1_epi32(ev1);
//     const __m256i e0eq0 = _mm256_cmpeq_epi32(sec0, ev0v);
//     const __m256i e0eq1 = _mm256_cmpeq_epi32(sec0, ev1v);
//     const __m256i e1eq0 = _mm256_cmpeq_epi32(sec1, ev0v);
//     const __m256i e1eq1 = _mm256_cmpeq_epi32(sec1, ev1v);
//     const __m256i anyeq = _mm256_or_si256(_mm256_or_si256(e0eq0, e0eq1),
//                                           _mm256_or_si256(e1eq0, e1eq1));
//     const __m256i m_01 = _mm256_srli_epi32(_mm256_castps_si256(_mm256_castsi256_ps(anyeq)), 31);
//     _mm256_storeu_si256((__m256i*)(mask_out + processed), m_01);
//     processed += lanes;
//   }
//   // Zero-pad remaining entries in mask_out up to 16 (optional)
// #else
//   for (size_t lane = 0; lane < chunk_len; ++lane) {
//     const int jidx = idx[base_index + lane];
//     const int s0 = elements0[jidx * stride];
//     const int s1 = elements1[jidx * stride];
//     mask_out[lane] = (s0 == ev0) | (s0 == ev1) | (s1 == ev0) | (s1 == ev1);
//   }
// #endif
// }

#endif // VAABB_H
