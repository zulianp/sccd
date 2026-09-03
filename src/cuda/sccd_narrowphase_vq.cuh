#ifndef SCCD_NARROWPHASE_VQ_CUH
#define SCCD_NARROWPHASE_VQ_CUH

#include "sccd_base.hpp"

#include <cstddef>

/**
 * \file
 * \brief Device vertex-quad narrow phase.
 *
 * Quads had no device narrow phase at all: `sccd_smesh_CCD.hpp` raised
 * `SMESH_ERROR("CUDA QUADSHELL4 narrow phase is not implemented")`, so a quad
 * mesh could be broad-phased on the GPU and then not finished there.
 *
 * ## Why this is a separate kernel
 *
 * The triangle kernel in `sccd_narrowphase.cu` is built end to end around a
 * four-point query packed into a `Vec4` -- vertex plus three corners for
 * vertex-face, four endpoints for edge-edge -- and on a barycentric domain, which
 * is why its `is_domain_valid` rejects `u + v >= 1`. A vertex-quad query is five
 * points on a domain where `u` and `v` range independently. Extending that kernel
 * would mean changing the load, the tolerances, the error bound, the corner
 * evaluation, the domain test and the splitting all at once. This is the same
 * branch and bound over the quad's own parameter domain, written directly.
 *
 * ## What it guarantees
 *
 * The same three properties the host path rests on, and for the same reasons:
 *
 *  - the rejection test pads by the certified numerical error bound, never by
 *    the user's tolerance;
 *  - an accepted box reports its `t` lower bound, which is at or before any root
 *    inside it;
 *  - exhaustion accepts. The depth cap and a full stack both accept rather than
 *    drop, because dropping a box that may contain a root is the one way this
 *    algorithm loses a collision.
 *
 * The per-thread stack is sized so overflow cannot occur, and the search depth
 * is clamped to what it holds. Those two numbers are tied together deliberately.
 * Neither limit can cost a collision -- exhaustion accepts, and accepting is
 * always safe -- so what an undersized stack buys is false positives and a time
 * of impact further before the true one than it needs to be. That is an accuracy
 * cost, not a safety one, and it is silent, which is why the two numbers are
 * derived from each other rather than chosen independently.
 *
 * Arithmetic is double throughout regardless of the storage type: the error
 * bound and the tolerances that terminate the search are too close together in
 * single precision for the guarantee to survive. A float result is narrowed
 * toward negative infinity, since reporting an earlier time of impact is safe
 * and reporting a later one is exactly the failure this exists to prevent.
 */

namespace sccd {
    namespace device {

        /**
         * \brief Vertex-quad times of impact for \p noverlaps candidate pairs.
         *
         * \p toi_stride 0 writes a single shared minimum to `toi[0]`; 1 writes one
         * time of impact per pair. Mirrors the host `narrow_phase_vq`.
         */
        template <typename T, typename I>
        int narrow_phase_vq(const size_t noverlaps,
                            const I* const SCCD_RESTRICT voverlap,
                            const I* const SCCD_RESTRICT qoverlap,
                            T** const SCCD_RESTRICT v0,
                            T** const SCCD_RESTRICT v1,
                            const size_t quad_stride,
                            I** const SCCD_RESTRICT quads,
                            const T max_toi,
                            T* const SCCD_RESTRICT toi,
                            const int max_depth,
                            const T tol,
                            const int toi_stride = 0);

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_NARROWPHASE_VQ_CUH
