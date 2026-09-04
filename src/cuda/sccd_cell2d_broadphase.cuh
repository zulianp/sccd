#ifndef SCCD_CELL2D_BROADPHASE_CUH
#define SCCD_CELL2D_BROADPHASE_CUH

#include "sccd_base.hpp"

#include <cstddef>

/**
 * \file
 * \brief Device port of the two-dimensional cell-list broad phase.
 *
 * The host version is `src/cell2d_broadphase.hpp`; the reasoning for two axes
 * rather than three, and for binning rather than sorting, is written up there and
 * in `wip/BROADPHASE.md`. On the host it took the broad phase from 4,467 ms
 * to 1,063 ms at 1.5M triangles.
 *
 * The shape suits a GPU better than the sweep does. Binning is count, prefix sum,
 * scatter -- the same idiom the rest of this file already uses for the pair CRS --
 * and the query has no sequential window walk, so every element is an independent
 * thread with no dependence on a sorted order it has to march along.
 */

namespace sccd {
    namespace device {

        /** \brief A uniform 2D cell list over two of the three coordinate axes. */
        template <typename T>
        struct Cell2DGridD {
            int axis0;
            int axis1;
            T min0, min1;
            T inv0, inv1;
            int n0, n1;

            // The grid is passed to kernels by value but is also built and
            // queried from host code, including from tests compiled as plain C++,
            // so the qualifiers only exist when nvcc is the compiler.
#ifdef __CUDACC__
            __host__ __device__
#endif
                ptrdiff_t
                ncells() const {
                return (ptrdiff_t)n0 * (ptrdiff_t)n1;
            }
        };

        /**
         * \brief Size the grid and bin \p n boxes into it.
         *
         * \p cellptr must have ncells + 1 entries and \p cellidx must have room
         * for the total span count; the caller gets that count back through
         * \p span_count so it can size \p cellidx between the two calls, exactly
         * as it does for the pair arrays.
         */
        template <typename T, typename I>
        void cell2d_setup_and_count(const int dim,
                                    const ptrdiff_t n,
                                    T** const SCCD_RESTRICT aabbs,
                                    Cell2DGridD<T>& grid,
                                    ptrdiff_t* const SCCD_RESTRICT cellptr,
                                    ptrdiff_t* const SCCD_RESTRICT span_count);

        template <typename T, typename I>
        void cell2d_fill(const ptrdiff_t n,
                         T** const SCCD_RESTRICT aabbs,
                         const Cell2DGridD<T>& grid,
                         const ptrdiff_t* const SCCD_RESTRICT cellptr,
                         I* const SCCD_RESTRICT cellidx,
                         ptrdiff_t* const SCCD_RESTRICT cursor);

        template <int first_nxe, int second_nxe, typename T, typename I>
        void cell2d_count_overlaps(const ptrdiff_t first_count,
                                   T** const SCCD_RESTRICT first_aabbs,
                                   I* const SCCD_RESTRICT first_idx,
                                   const ptrdiff_t first_stride,
                                   I** const SCCD_RESTRICT first_elements,
                                   T** const SCCD_RESTRICT second_aabbs,
                                   I* const SCCD_RESTRICT second_idx,
                                   const ptrdiff_t second_stride,
                                   I** const SCCD_RESTRICT second_elements,
                                   const Cell2DGridD<T>& grid,
                                   const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                   const I* const SCCD_RESTRICT cellidx,
                                   ptrdiff_t* const SCCD_RESTRICT ccdptr);

        template <int first_nxe, int second_nxe, typename T, typename I>
        void cell2d_collect_overlaps(const ptrdiff_t first_count,
                                     T** const SCCD_RESTRICT first_aabbs,
                                     I* const SCCD_RESTRICT first_idx,
                                     const ptrdiff_t first_stride,
                                     I** const SCCD_RESTRICT first_elements,
                                     T** const SCCD_RESTRICT second_aabbs,
                                     I* const SCCD_RESTRICT second_idx,
                                     const ptrdiff_t second_stride,
                                     I** const SCCD_RESTRICT second_elements,
                                     const Cell2DGridD<T>& grid,
                                     const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                     const I* const SCCD_RESTRICT cellidx,
                                     const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                     I* const SCCD_RESTRICT first_out,
                                     I* const SCCD_RESTRICT second_out);

        template <int nxe, typename T, typename I>
        void cell2d_count_self_overlaps(const ptrdiff_t element_count,
                                        T** const SCCD_RESTRICT aabbs,
                                        I* const SCCD_RESTRICT idx,
                                        const ptrdiff_t stride,
                                        I** const SCCD_RESTRICT elements,
                                        const Cell2DGridD<T>& grid,
                                        const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                        const I* const SCCD_RESTRICT cellidx,
                                        ptrdiff_t* const SCCD_RESTRICT ccdptr);

        template <int nxe, typename T, typename I>
        void cell2d_collect_self_overlaps(const ptrdiff_t element_count,
                                          T** const SCCD_RESTRICT aabbs,
                                          I* const SCCD_RESTRICT idx,
                                          const ptrdiff_t stride,
                                          I** const SCCD_RESTRICT elements,
                                          const Cell2DGridD<T>& grid,
                                          const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                          const I* const SCCD_RESTRICT cellidx,
                                          const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                          I* const SCCD_RESTRICT out0,
                                          I* const SCCD_RESTRICT out1);

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_CELL2D_BROADPHASE_CUH
