#ifndef SCCD_NARROWPHASE_CUH
#define SCCD_NARROWPHASE_CUH

#include "sccd_base.hpp"

#include <cstddef>

namespace sccd {
    namespace device {

        // toi_stride: 1 -> per-candidate toi (toi has length noverlaps);
        //             0 -> shared toi across all candidates (toi has length 1).
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
                          const T max_toi,
                          const int toi_stride = 0);

        template <int nxe, typename T, typename I>
        T narrow_phase_vf(const size_t noverlaps,
                          const I* const SCCD_RESTRICT voveralp,
                          const I* const SCCD_RESTRICT foveralp,
                          // Geometric data
                          T** const SCCD_RESTRICT v0,
                          T** const SCCD_RESTRICT v1,
                          const size_t face_stride,
                          I** const SCCD_RESTRICT faces,
                          const T max_toi);

    }  // namespace device
}  // namespace sccd

#endif
