#ifndef SCCD_BROADPHASE_CUH
#define SCCD_BROADPHASE_CUH

#include "sccd_base.hpp"

#include <cstddef>

namespace sccd {
    namespace device {

        template <typename T>
        int choose_axis(const int dim, const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);

        template <typename T>
        void enumerate(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx);

        template <typename T, typename I>
        void sort_along_axis(const int dim,
                             const ptrdiff_t n,
                             const int sort_axis,
                             T** const SCCD_RESTRICT arrays,
                             I* const SCCD_RESTRICT idx,
                             T* const SCCD_RESTRICT scratch);

        template <int nxe, typename T, typename I>
        void count_self_overlaps(const int sort_axis,
                                 const ptrdiff_t element_count,
                                 T** const SCCD_RESTRICT aabbs,
                                 I* const SCCD_RESTRICT idx,
                                 const ptrdiff_t stride,
                                 I** const SCCD_RESTRICT elements,
                                 ptrdiff_t* const SCCD_RESTRICT ccdptr);

        template <int nxe, typename T, typename I>
        void collect_self_overlaps(const int sort_axis,
                                   const ptrdiff_t element_count,
                                   T** const SCCD_RESTRICT aabbs,
                                   I* const SCCD_RESTRICT idx,
                                   const ptrdiff_t stride,
                                   I** const elements,
                                   const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                   I* SCCD_RESTRICT foverlap,
                                   I* SCCD_RESTRICT noverlap);

        template <typename T>
        void cummax(const ptrdiff_t n, const T* const SCCD_RESTRICT in, T* const SCCD_RESTRICT out);

        template <int first_nxe, int second_nxe, typename T, typename I>
        bool count_overlaps(const int sort_axis,
                            const ptrdiff_t first_count,
                            T** const SCCD_RESTRICT first_aabbs,
                            I* const SCCD_RESTRICT first_idx,
                            const ptrdiff_t first_stride,
                            I** const SCCD_RESTRICT first_elements,
                            const ptrdiff_t second_count,
                            T** const SCCD_RESTRICT second_aabbs,
                            I* const SCCD_RESTRICT second_idx,
                            const ptrdiff_t second_stride,
                            I** const SCCD_RESTRICT second_elements,
                            ptrdiff_t* const SCCD_RESTRICT ccdptr,
                            const T* const SCCD_RESTRICT cummax);

        template <int first_nxe, int second_nxe, typename T, typename I>
        void collect_overlaps(const int sort_axis,
                              const ptrdiff_t first_count,
                              T** const SCCD_RESTRICT first_aabbs,
                              I* const SCCD_RESTRICT first_idx,
                              const ptrdiff_t first_stride,
                              I** SCCD_RESTRICT const first_elements,
                              const ptrdiff_t second_count,
                              T** const SCCD_RESTRICT second_aabbs,
                              I* const SCCD_RESTRICT second_idx,
                              const ptrdiff_t second_stride,
                              I** SCCD_RESTRICT const second_elements,
                              const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                              const T* const SCCD_RESTRICT cummax,
                              I* SCCD_RESTRICT foverlap,
                              I* SCCD_RESTRICT noverlap);

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_BROADPHASE_CUH