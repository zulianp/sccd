#ifndef SCCD_CELL2D_BROADPHASE_CUH
#define SCCD_CELL2D_BROADPHASE_CUH

#include "sccd_base.hpp"

#include <cstddef>

/**
 * \file
 * \brief Device port of the two-dimensional cell-list broad phase.
 *
 * The host version is `src/broadphase/sccd_broadphase_cell2d.hpp`; the reasoning for two axes
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
         * \brief Write `idx[i] = i` over \p n entries.
         *
         * The cell list does not sort, so it needs the identity permutation
         * where the sweep leaves the order its sort produced. A step that
         * switched strategies would otherwise inherit the other one's
         * permutation, so this is written every step rather than once.
         */
        template <typename I>
        void fill_identity(const ptrdiff_t n, I* const SCCD_RESTRICT idx);

        /**
         * \brief Size the grid and bin \p n boxes into it.
         *
         * \p ranges must have `4 * n + 4` entries: the cell range each box
         * covers, which the scatter pass reads back so the two passes cannot
         * disagree about it, and the grid it was computed for.
         *
         * \p cellptr must have ncells + 1 entries and \p cellidx must have room
         * for the total span count; the caller gets that count back through
         * \p span_count so it can size \p cellidx between the two calls, exactly
         * as it does for the pair arrays.
         */
        template <typename T, typename I>
        void cell2d_setup_and_count(const ptrdiff_t n,
                                    T** const SCCD_RESTRICT aabbs,
                                    Cell2DGridD<T>& grid,
                                    ptrdiff_t* const SCCD_RESTRICT cellptr,
                                    int* const SCCD_RESTRICT ranges,
                                    ptrdiff_t* const SCCD_RESTRICT span_count);

        /**
         * \brief Scatter each box into its cells.
         *
         * \p capacity is how many entries \p cellidx holds, and \p cursor must
         * have `ncells + 1` of them: the last counts writes the kernel refused
         * because they fell outside the array. The counting pass reserved
         * exactly one slot per span, so a non-zero count means the two passes
         * disagreed and the result is incomplete -- it is reported through
         * \p out_rejected (which may be null) rather than written past the end.
         */
        template <typename T, typename I>
        void cell2d_fill(const ptrdiff_t n,
                         T** const SCCD_RESTRICT aabbs,
                         const Cell2DGridD<T>& grid,
                         const ptrdiff_t* const SCCD_RESTRICT cellptr,
                         const int* const SCCD_RESTRICT ranges,
                         I* const SCCD_RESTRICT cellidx,
                         const ptrdiff_t capacity,
                         ptrdiff_t* const SCCD_RESTRICT cursor,
                         ptrdiff_t* const SCCD_RESTRICT out_rejected);

        template <int first_nxe, int second_nxe, typename T, typename I>
        void cell2d_count_overlaps(const ptrdiff_t first_count,
                                   T** const SCCD_RESTRICT first_aabbs,
                                   I* const SCCD_RESTRICT first_idx,
                                   const ptrdiff_t first_element_stride,
                                   I** const SCCD_RESTRICT first_elements,
                                   T** const SCCD_RESTRICT second_aabbs,
                                   I* const SCCD_RESTRICT second_idx,
                                   const ptrdiff_t second_element_stride,
                                   I** const SCCD_RESTRICT second_elements,
                                   const Cell2DGridD<T>& grid,
                                   const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                   const I* const SCCD_RESTRICT cellidx,
                                   ptrdiff_t* const SCCD_RESTRICT ccdptr);

        template <int first_nxe, int second_nxe, typename T, typename I>
        void cell2d_collect_overlaps(const ptrdiff_t first_count,
                                     T** const SCCD_RESTRICT first_aabbs,
                                     I* const SCCD_RESTRICT first_idx,
                                     const ptrdiff_t first_element_stride,
                                     I** const SCCD_RESTRICT first_elements,
                                     T** const SCCD_RESTRICT second_aabbs,
                                     I* const SCCD_RESTRICT second_idx,
                                     const ptrdiff_t second_element_stride,
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
                                        const ptrdiff_t element_stride,
                                        I** const SCCD_RESTRICT elements,
                                        const Cell2DGridD<T>& grid,
                                        const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                        const I* const SCCD_RESTRICT cellidx,
                                        ptrdiff_t* const SCCD_RESTRICT ccdptr);

        template <int nxe, typename T, typename I>
        void cell2d_collect_self_overlaps(const ptrdiff_t element_count,
                                          T** const SCCD_RESTRICT aabbs,
                                          I* const SCCD_RESTRICT idx,
                                          const ptrdiff_t element_stride,
                                          I** const SCCD_RESTRICT elements,
                                          const Cell2DGridD<T>& grid,
                                          const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                          const I* const SCCD_RESTRICT cellidx,
                                          const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                          I* const SCCD_RESTRICT first_out,
                                          I* const SCCD_RESTRICT second_out);


        /**
         * \brief Size the grid and bin \p n boxes by their minimum corner alone.
         *
         * The edge-edge form of the cell list: a box enters the one cell holding
         * its minimum corner, so the cell array holds exactly \p n entries and the
         * caller needs no span count to size it. No \p ranges array either --
         * the cell a box lands in is one clamp of one coordinate pair, a pure
         * function of the box and the grid, so the counting and the scatter
         * cannot disagree about it the way two nested loops could.
         */
        template <typename T, typename I>
        void cell2dmin_setup_and_count(const ptrdiff_t n,
                                       T** const SCCD_RESTRICT aabbs,
                                       Cell2DGridD<T>& grid,
                                       ptrdiff_t* const SCCD_RESTRICT cellptr);

        /** \brief Scatter each box index into the cell its minimum corner is in. */
        template <typename T, typename I>
        void cell2dmin_fill(const ptrdiff_t n,
                            T** const SCCD_RESTRICT aabbs,
                            const Cell2DGridD<T>& grid,
                            const ptrdiff_t* const SCCD_RESTRICT cellptr,
                            I* const SCCD_RESTRICT cellidx,
                            ptrdiff_t* const SCCD_RESTRICT cursor);

        /**
         * \brief The two bounds the forward walk prunes with.
         *
         * \p cell_hi1 is the largest upper bound on the second grid axis of the
         * boxes in a cell, and \p row_prefix the running maximum along each row of
         * the largest upper bound on the first. Both have `ncells` entries. The
         * prefix is what lets a query skip a run of columns in one binary search
         * rather than a test per cell, and it is why the binning needs no bound
         * on how many cells a box may span.
         */
        template <typename T, typename I>
        void cell2dmin_bounds(const Cell2DGridD<T>& grid,
                              T** const SCCD_RESTRICT aabbs,
                              const ptrdiff_t* const SCCD_RESTRICT cellptr,
                              const I* const SCCD_RESTRICT cellidx,
                              T* const SCCD_RESTRICT row_prefix,
                              T* const SCCD_RESTRICT cell_hi1);

        template <int nxe, typename T, typename I>
        void cell2dmin_count_self_overlaps(const ptrdiff_t element_count,
                                           T** const SCCD_RESTRICT aabbs,
                                           I* const SCCD_RESTRICT idx,
                                           const ptrdiff_t element_stride,
                                           I** const SCCD_RESTRICT elements,
                                           const Cell2DGridD<T>& grid,
                                           const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                           const I* const SCCD_RESTRICT cellidx,
                                           const T* const SCCD_RESTRICT row_prefix,
                                           const T* const SCCD_RESTRICT cell_hi1,
                                           ptrdiff_t* const SCCD_RESTRICT ccdptr);

        template <int nxe, typename T, typename I>
        void cell2dmin_collect_self_overlaps(const ptrdiff_t element_count,
                                             T** const SCCD_RESTRICT aabbs,
                                             I* const SCCD_RESTRICT idx,
                                             const ptrdiff_t element_stride,
                                             I** const SCCD_RESTRICT elements,
                                             const Cell2DGridD<T>& grid,
                                             const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                             const I* const SCCD_RESTRICT cellidx,
                                             const T* const SCCD_RESTRICT row_prefix,
                                             const T* const SCCD_RESTRICT cell_hi1,
                                             const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                             I* const SCCD_RESTRICT first_out,
                                             I* const SCCD_RESTRICT second_out);

        /**
         * \brief The most rows on the second grid axis any of \p n boxes spans.
         *
         * The cross-list walk bounds its row range with this. It is written to
         * \p krow, a single `int` of device memory that the query kernels read
         * directly, so the value stays on the device and no step synchronises
         * for it.
         */
        template <typename T>
        void cell2dmin_max_row_span(const ptrdiff_t n,
                                    T** const SCCD_RESTRICT aabbs,
                                    const Cell2DGridD<T>& grid,
                                    int* const SCCD_RESTRICT krow);

        /**
         * \brief Vertex-face count over a minimum-corner binning of the faces.
         *
         * The mirror image of \ref cell2d_count_overlaps for this query: the
         * faces are binned and the vertices walk, so the cells are face-sized and
         * a vertex reads a fraction of the cells a face does. The pair set is the
         * same -- this changes where a pair is found, not which pairs exist.
         *
         * \p krow points at the single device `int` \ref cell2dmin_max_row_span
         * filled over the same faces, and
         * the grid, \p cellptr, \p cellidx, \p row_prefix and \p cell_hi1 from
         * binning the faces with \ref cell2dmin_setup_and_count,
         * \ref cell2dmin_fill and \ref cell2dmin_bounds.
         *
         * \p ccdptr holds `vertex_count + 1` entries, one per vertex, since the
         * vertices are what the walk is indexed by here.
         */
        template <int S, typename T, typename I>
        void cell2dmin_count_vf_overlaps(const ptrdiff_t vertex_count,
                                         T** const SCCD_RESTRICT vaabbs,
                                         const I* const SCCD_RESTRICT vertex_idx,
                                         T** const SCCD_RESTRICT faabbs,
                                         const I* const SCCD_RESTRICT face_idx,
                                         const ptrdiff_t face_element_stride,
                                         I** const SCCD_RESTRICT face_elements,
                                         const Cell2DGridD<T>& grid,
                                         const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                         const I* const SCCD_RESTRICT cellidx,
                                         const T* const SCCD_RESTRICT row_prefix,
                                         const T* const SCCD_RESTRICT cell_hi1,
                                         const int* const SCCD_RESTRICT krow,
                                         ptrdiff_t* const SCCD_RESTRICT ccdptr);

        /** \brief Write those pairs, face first, to match the shipped query. */
        template <int S, typename T, typename I>
        void cell2dmin_collect_vf_overlaps(const ptrdiff_t vertex_count,
                                           T** const SCCD_RESTRICT vaabbs,
                                           const I* const SCCD_RESTRICT vertex_idx,
                                           T** const SCCD_RESTRICT faabbs,
                                           const I* const SCCD_RESTRICT face_idx,
                                           const ptrdiff_t face_element_stride,
                                           I** const SCCD_RESTRICT face_elements,
                                           const Cell2DGridD<T>& grid,
                                           const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                           const I* const SCCD_RESTRICT cellidx,
                                           const T* const SCCD_RESTRICT row_prefix,
                                           const T* const SCCD_RESTRICT cell_hi1,
                                           const int* const SCCD_RESTRICT krow,
                                           const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                           I* const SCCD_RESTRICT face_out,
                                           I* const SCCD_RESTRICT vertex_out);

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_CELL2D_BROADPHASE_CUH
