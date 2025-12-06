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
                double tol = 1e-6;
                int max_iter = 4000;

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (find_root_tight_inclusion<double>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
                    toi[i] = t;
                    min_t = sccd::min<T>(t, min_t);
                } else {
                    toi[i] = infty;
                }
#else
                if (find_root_grid<double>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack)) 
                // if (find_root_bisection<double>(4000, 1e-10, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) 
                {
                    toi[i] = t;
                    min_t = sccd::min<T>(t, min_t);
                } else {
                    toi[i] = infty;
                }
#endif
            }
        });

        return min_t;
    }

}  // namespace sccd

#endif