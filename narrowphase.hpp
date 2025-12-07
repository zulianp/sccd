#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "assert.h"
#include "roots.hpp"
#include "srootfinder.hpp"
#include "vaabb.h"

namespace sccd {

    template <int nxe, typename T>
    T narrow_phase_vf(const size_t noverlaps,
                      const idx_t* const SFEM_RESTRICT voveralp,
                      const idx_t* const SFEM_RESTRICT foveralp,
                      // Geometric data
                      T** const SFEM_RESTRICT v0,
                      T** const SFEM_RESTRICT v1,
                      const size_t face_stride,
                      idx_t** const SFEM_RESTRICT faces,
                      // Output
                      T* const SFEM_RESTRICT toi) {
        const T infty = 100000;
        T min_t = infty;

        int USE_TI = 0;
        SFEM_READ_ENV(USE_TI, atoi);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, noverlaps), [&](const tbb::blocked_range<size_t>& r) {
            std::vector<Box<double>> stack;
            for (size_t i = r.begin(); i < r.end(); i++) {
                // for (size_t i = 0; i < noverlaps; i++) {
                const idx_t vi = voveralp[i];
                const idx_t fi = foveralp[i];

                idx_t nodes[3] = {faces[0][fi * face_stride], faces[1][fi * face_stride], faces[2][fi * face_stride]};

                const double sv[3] = {v0[0][vi], v0[1][vi], v0[2][vi]};
                const double ev[3] = {v1[0][vi], v1[1][vi], v1[2][vi]};

                const double s1[3] = {v0[0][nodes[0]], v0[1][nodes[0]], v0[2][nodes[0]]};
                const double s2[3] = {v0[0][nodes[1]], v0[1][nodes[1]], v0[2][nodes[1]]};
                const double s3[3] = {v0[0][nodes[2]], v0[1][nodes[2]], v0[2][nodes[2]]};

                const double e1[3] = {v1[0][nodes[0]], v1[1][nodes[0]], v1[2][nodes[0]]};
                const double e2[3] = {v1[0][nodes[1]], v1[1][nodes[1]], v1[2][nodes[1]]};
                const double e3[3] = {v1[0][nodes[2]], v1[1][nodes[2]], v1[2][nodes[2]]};

                // Iteration variables
                double t = 0;
                double u = 0;
                double v = 0;
                double tol = 1e-12;
                int max_iter = 16;

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (USE_TI) {
                    if (find_root_tight_inclusion_vf<double>(
                            max_iter * 1000, tol, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
                        toi[i] = t;
                        min_t = sccd::min<T>(t, min_t);
                    } else {
                        toi[i] = infty;
                    }
                    continue;
                }
#endif
                if (find_root_grid_vf<double>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack))
                // if (find_root_bisection<double>(max_iter * 3*3*3, tol, sv,
                // s1, s2, s3, ev, e1, e2, e3, t, u, v))
                {
                    toi[i] = t;
                    min_t = sccd::min<T>(t, min_t);
                } else {
                    toi[i] = infty;
                }
            }
        });

        return min_t;
    }


        template <typename T>
        T narrow_phase_ee(const size_t noverlaps,
                          const idx_t* const SFEM_RESTRICT e0overalp,
                          const idx_t* const SFEM_RESTRICT e1overalp,
                          // Geometric data
                          T** const SFEM_RESTRICT v0,
                          T** const SFEM_RESTRICT v1,
                          const size_t edge_stride,
                          idx_t** const SFEM_RESTRICT edges,
                          // Output
                          T* const SFEM_RESTRICT toi) {
            const T infty = 100000;
            T min_t = infty;

            int USE_TI = 0;
            SFEM_READ_ENV(USE_TI, atoi);

            tbb::parallel_for(tbb::blocked_range<size_t>(0, noverlaps), [&](const tbb::blocked_range<size_t>& r) {
                std::vector<Box<double>> stack;
                for (size_t i = r.begin(); i < r.end(); i++) {
                    // for (size_t i = 0; i < noverlaps; i++) {
                    const idx_t i0 = e0overalp[i];
                    const idx_t i1 = e1overalp[i];

                    idx_t nodes0[2] = {edges[0][i0 * edge_stride], edges[1][i0 * edge_stride]};
                    idx_t nodes1[2] = {edges[0][i1 * edge_stride], edges[1][i1 * edge_stride]};

                    const double s1[3] = {v0[0][nodes0[0]], v0[1][nodes0[0]], v0[2][nodes0[0]]};
                    const double s2[3] = {v0[0][nodes0[1]], v0[1][nodes0[1]], v0[2][nodes0[1]]};

                    const double s3[3] = {v0[0][nodes1[2]], v0[1][nodes1[2]], v0[2][nodes1[2]]};
                    const double s4[3] = {v0[0][nodes1[2]], v0[1][nodes1[2]], v0[2][nodes1[2]]};

                    const double e1[3] = {v1[0][nodes0[0]], v1[1][nodes0[0]], v1[2][nodes0[0]]};
                    const double e2[3] = {v1[0][nodes0[1]], v1[1][nodes0[1]], v1[2][nodes0[1]]};

                    const double e3[3] = {v1[0][nodes1[2]], v1[1][nodes1[2]], v1[2][nodes1[2]]};
                    const double e4[3] = {v1[0][nodes1[2]], v1[1][nodes1[2]], v1[2][nodes1[2]]};


                    // Iteration variables
                    double t = 0;
                    double u = 0;
                    double v = 0;
                    double tol = 1e-12;
                    int max_iter = 16;

    #ifdef SCCD_ENABLE_TIGHT_INCLUSION
    #warning "SCCD_ENABLE_TIGHT_INCLUSION"
                    if (USE_TI) {
                        if (find_root_tight_inclusion_ee<double>(
                                max_iter * 1000, tol, 
                                s1, s2, s3, s4, 
                                e1, e2, e3, e4, 
                                t, u, v)) 
                        {
                            toi[i] = t;
                            min_t = sccd::min<T>(t, min_t);
                        } else {
                            toi[i] = infty;
                        }
                        continue;
                    }
    #endif
                    if (find_root_grid_ee<double>(max_iter, tol, 
                        s1, s2, s3, s4, 
                        e1, e2, e3, e4, 
                        t, u, v, stack))
                    {
                        toi[i] = t;
                        min_t = sccd::min<T>(t, min_t);
                    } else {
                        toi[i] = infty;
                    }
                }
            });

            return min_t;
        }

}  // namespace sccd

#endif