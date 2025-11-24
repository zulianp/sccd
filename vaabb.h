#ifndef VAABB_H
#define VAABB_H

#include <stdint.h>

#define CHUNK_SIZE 16
#define SFEM_RESTRICT __restrict

typedef float geom_t;

inline static uint32_t disjoint(const geom_t aminx, const geom_t aminy,
                                const geom_t aminz, const geom_t amaxx,
                                const geom_t amaxy, const geom_t amaxz,
                                const geom_t bminx, const geom_t bminy,
                                const geom_t bminz, const geom_t bmaxx,
                                const geom_t bmaxy, const geom_t bmaxz) {
  return aminx > bmaxx | aminy > bmaxy | aminz > bmaxz | bminx > amaxx |
         bminy > amaxy | bminz > amaxz;
}

void vdisjoint(const geom_t *const SFEM_RESTRICT aminx,
               const geom_t *const SFEM_RESTRICT aminy,
               const geom_t *const SFEM_RESTRICT aminz,
               const geom_t *const SFEM_RESTRICT amaxx,
               const geom_t *const SFEM_RESTRICT amaxy,
               const geom_t *const SFEM_RESTRICT amaxz,
               const geom_t *const SFEM_RESTRICT bminx,
               const geom_t *const SFEM_RESTRICT bminy,
               const geom_t *const SFEM_RESTRICT bminz,
               const geom_t *const SFEM_RESTRICT bmaxx,
               const geom_t *const SFEM_RESTRICT bmaxy,
               const geom_t *const SFEM_RESTRICT bmaxz,
               uint32_t *SFEM_RESTRICT mask) {
#pragma omp simd aligned(aminx, aminy, aminz, amaxx, amaxy, amaxz, bminx,      \
                         bminy, bminz, bmaxx, bmaxy, bmaxz, mask : 64)
  for (int i = 0; i < CHUNK_SIZE; i++) {
    mask[i] =
        disjoint(aminx[i], aminy[i], aminz[i], amaxx[i], amaxy[i], amaxz[i],
                 bminx[i], bminy[i], bminz[i], bmaxx[i], bmaxy[i], bmaxz[i]);
  }
}

#endif // VAABB_H
