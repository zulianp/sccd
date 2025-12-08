#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "assert.h"

#include "roots.hpp"
#include "srootfinder.hpp"
#include "vaabb.hpp"

namespace sccd {

    template <int nxe, typename T, typename I>
    T narrow_phase_vf(const size_t noverlaps,
                      const I* const SCCD_RESTRICT voveralp,
                      const I* const SCCD_RESTRICT foveralp,
                      // Geometric data
                      T** const SCCD_RESTRICT v0,
                      T** const SCCD_RESTRICT v1,
                      const size_t face_stride,
                      I** const SCCD_RESTRICT faces,
                      // Output
                      T* const SCCD_RESTRICT toi) {
        using T_HP = double;
        const T infty = 100000;

        int USE_TI = 0;
        SFEM_READ_ENV(USE_TI, atoi);

        int SCCD_MAX_ITER = 12;
        SFEM_READ_ENV(SCCD_MAX_ITER, atoi);

        T_HP tol = 1e-12;

        sccd::parallel_for_br(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            std::vector<Box<T_HP>> stack;
            for (ptrdiff_t i = rbegin; i < rend; i++) {
                const I vi = voveralp[i];
                const I fi = foveralp[i];

                I nodes[3] = {faces[0][fi * face_stride], faces[1][fi * face_stride], faces[2][fi * face_stride]};

                const T_HP sv[3] = {v0[0][vi], v0[1][vi], v0[2][vi]};
                const T_HP ev[3] = {v1[0][vi], v1[1][vi], v1[2][vi]};

                const T_HP s1[3] = {v0[0][nodes[0]], v0[1][nodes[0]], v0[2][nodes[0]]};
                const T_HP s2[3] = {v0[0][nodes[1]], v0[1][nodes[1]], v0[2][nodes[1]]};
                const T_HP s3[3] = {v0[0][nodes[2]], v0[1][nodes[2]], v0[2][nodes[2]]};

                const T_HP e1[3] = {v1[0][nodes[0]], v1[1][nodes[0]], v1[2][nodes[0]]};
                const T_HP e2[3] = {v1[0][nodes[1]], v1[1][nodes[1]], v1[2][nodes[1]]};
                const T_HP e3[3] = {v1[0][nodes[2]], v1[1][nodes[2]], v1[2][nodes[2]]};

                // Iteration variables
                T_HP t = infty;
                T_HP u = 0;
                T_HP v = 0;

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (USE_TI) {
                    if (find_root_tight_inclusion_vf<T_HP>(
                            SCCD_MAX_ITER * 1000, tol, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
                        toi[i] = t;
                    } else {
                        toi[i] = infty;
                    }
                    continue;
                }
#endif
                if (find_root_grid_rotate_vf<T_HP>(SCCD_MAX_ITER, tol, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack))
                {
                    toi[i] = t;
                } else {
                    toi[i] = infty;
                }
            }
        });

        T min_t = infty;
        for (size_t i = 0; i < noverlaps; i++) {
            min_t = sccd::min<T>(toi[i], min_t);
        }


        return min_t;
    }

    template <typename T, typename I>
    T narrow_phase_ee(const size_t noverlaps,
                      const I* const SCCD_RESTRICT e0overalp,
                      const I* const SCCD_RESTRICT e1overalp,
                      // Geometric data
                      T** const SCCD_RESTRICT v0,
                      T** const SCCD_RESTRICT v1,
                      const size_t edge_stride,
                      I** const SCCD_RESTRICT edges,
                      // Output
                      T* const SCCD_RESTRICT toi) {
        using T_HP = double;
        const T infty = 100000;
        

        int USE_TI = 0;
        SFEM_READ_ENV(USE_TI, atoi);

        int SCCD_MAX_ITER = 12;
        SFEM_READ_ENV(SCCD_MAX_ITER, atoi);

        T_HP tol = 1e-12;

        sccd::parallel_for_br(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            std::vector<Box<T_HP>> stack;
            for (ptrdiff_t i = rbegin; i < rend; i++) {
                const I i0 = e0overalp[i];
                const I i1 = e1overalp[i];

                I nodes0[2] = {edges[0][i0 * edge_stride], edges[1][i0 * edge_stride]};
                I nodes1[2] = {edges[0][i1 * edge_stride], edges[1][i1 * edge_stride]};

                const T_HP s1[3] = {v0[0][nodes0[0]], v0[1][nodes0[0]], v0[2][nodes0[0]]};
                const T_HP s2[3] = {v0[0][nodes0[1]], v0[1][nodes0[1]], v0[2][nodes0[1]]};

                const T_HP s3[3] = {v0[0][nodes1[2]], v0[1][nodes1[2]], v0[2][nodes1[2]]};
                const T_HP s4[3] = {v0[0][nodes1[2]], v0[1][nodes1[2]], v0[2][nodes1[2]]};

                const T_HP e1[3] = {v1[0][nodes0[0]], v1[1][nodes0[0]], v1[2][nodes0[0]]};
                const T_HP e2[3] = {v1[0][nodes0[1]], v1[1][nodes0[1]], v1[2][nodes0[1]]};

                const T_HP e3[3] = {v1[0][nodes1[2]], v1[1][nodes1[2]], v1[2][nodes1[2]]};
                const T_HP e4[3] = {v1[0][nodes1[2]], v1[1][nodes1[2]], v1[2][nodes1[2]]};

                // Iteration variables
                T_HP t = infty;
                T_HP u = 0;
                T_HP v = 0;

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (USE_TI) {
                    if (find_root_tight_inclusion_ee<T_HP>(
                            SCCD_MAX_ITER * 1000, tol, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v)) {
                        toi[i] = t;
                    } else {
                        toi[i] = infty;
                    }
                    continue;
                }
#endif
                if (find_root_grid_ee<T_HP>(SCCD_MAX_ITER, tol, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack)) {
                    toi[i] = t;
                } else {
                    toi[i] = infty;
                }
            }
        });

        T min_t = infty;
        for (size_t i = 0; i < noverlaps; i++) {
            min_t = sccd::min<T>(toi[i], min_t);
        }

        return min_t;
    }

}  // namespace sccd

#endif