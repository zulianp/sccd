#include "sccd_narrowphase.cuh"

namespace sccd {
    namespace device {

        __device__ void sample_f_ee(const T tl,
                                    const T tu,
                                    const T ul,
                                    const T uu,
                                    const T vl,
                                    const T vu,
                                    const Vec4 sx,
                                    const Vec4 ex,
                                    T* const f) {
            // Compute spatial displacements range
            Vec4 dx = ex - sx;

            // Compute temporal displacements for lower bound
            {
                Vec4 xt = tl * dx + sx;
                f[0] = ((xt[1] - xt[0]) * ul + xt[0] - (xt[3] - xt[2]) * vl + xt[2]);
                f[1] = ((xt[1] - xt[0]) * ul + xt[0] - (xt[3] - xt[2]) * vu + xt[2]);
                f[2] = ((xt[1] - xt[0]) * uu + xt[0] - (xt[3] - xt[2]) * vl + xt[2]);
                f[3] = ((xt[1] - xt[0]) * uu + xt[0] - (xt[3] - xt[2]) * vu + xt[2]);
            }

            // Compute temporal displacements for upper bound
            {
                Vec4 xt = tu * dy + sy;
                f[4] = ((xt[1] - xt[0]) * ul + xt[0] - (xt[3] - xt[2]) * vl + xt[2]);
                f[5] = ((xt[1] - xt[0]) * ul + xt[0] - (xt[3] - xt[2]) * vu + xt[2]);
                f[6] = ((xt[1] - xt[0]) * uu + xt[0] - (xt[3] - xt[2]) * vl + xt[2]);
                f[7] = ((xt[1] - xt[0]) * uu + xt[0] - (xt[3] - xt[2]) * vu + xt[2]);
            }
        }

        __device__ void fminmax(const T* const f, T& fmin, T& fmax) {
            fmin = f[0];
            fmax = f[0];
            for (int i = 1; i < 8; i++) {
                fmin = MIN(fmin, f[i]);
                fmax = MAX(fmax, f[i]);
            }
        }

#define MASK_DOMAIN_SMALLER_THAN_TOL 1
#define MASK_BOX_INSIDE_EPSILON_BOX 2
#define MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE 4
#define MASK_INTERVAL_TERMINAL 8

        __device__ uint8_t cond_mask(const T fmin, const T fmax, const T tol, const T adaptive_tol) {
            uint8_t cond1 = (fmax - fmin <= adaptive_tol);
            uint8_t cond2 = !((fmin < tol) | (fmax > -tol));
            uint8_t cond3 = (fmax - fmin < tol);
            uint8_t cond4 = (fmin >= fmax);

            uint8_t cond_mask = cond1 & MASK_DOMAIN_SMALLER_THAN_TOL;
            cond_mask |= cond2 & MASK_BOX_INSIDE_EPSILON_BOX;
            cond_mask |= cond3 & MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE;
            cond_mask |= cond4 & MASK_INTERVAL_TERMINAL;
            cond_mask &= ((fmin <= tol) & ((fmax >= -tol) ? 0xf : 0));
            return cond_mask;
        }

        template <typename T, typename Vec4>
        __device__ void contains_origin_ee(const int nx,
                                           const int ny,
                                           const int x,
                                           const int y,
                                           const T tlower,
                                           const T tupper,
                                           const T ulower,
                                           const T uupper,
                                           const T vlower,
                                           const T vupper,
                                           const Vec4 sx,
                                           const Vec4 sy,
                                           const Vec4 sz,
                                           const Vec4 ex,
                                           const Vec4 ey,
                                           const Vec4 ez,
                                           const T tol,
                                           const T adaptive_tol,
                                           bool* const SCCD_RESTRICT contains_origin,
                                           bool* const SCCD_RESTRICT accept) {
            const T u_h = (uupper - ulower) / nx;
            const T v_h = (vupper - vlower) / ny;

            const T ul = ulower + x * u_h;
            const T vl = vlower + y * v_h;
            const T uu = ulower + (x + 1) * u_h;
            const T vu = vlower + (y + 1) * v_h;

            T fmin, fmax;
            T f[8];

            sample_f_ee(tlower, tupper, ul, uu, vl, vu, sx, ex, f);
            fminmax(f, fmin, fmax);
            const uint8_t x_mask = cond_mask(fmin, fmax, tol, adaptive_tol);
            *contains_origin = (fmin <= tol) & (fmax >= -tol);

            sample_f_ee(tlower, tupper, ul, uu, vl, vu, sy, ey, f);
            fminmax(f, fmin, fmax);
            const uint8_t y_mask = cond_mask(fmin, fmax, tol, adaptive_tol);
            *contains_origin &= (fmin <= tol) & (fmax >= -tol);

            sample_f_ee(tlower, tupper, ul, uu, vl, vu, sz, ez, f);
            fminmax(f, fmin, fmax);
            const uint8_t z_mask = cond_mask(fmin, fmax, tol, adaptive_tol);
            *contains_origin &= (fmin <= tol) & (fmax >= -tol);

            const uint8_t and_mask = (x_mask & y_mask & z_mask);
            const uint8_t or_mask = (x_mask | y_mask | z_mask);

            const bool cond1 = and_mask & MASK_DOMAIN_SMALLER_THAN_TOL;
            const bool cond2 = and_mask & MASK_BOX_INSIDE_EPSILON_BOX;
            const bool cond3 = or_mask & MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE;
            const bool cond4 = and_mask & MASK_INTERVAL_TERMINAL;

            *accept = *contains_origin && (cond1 || cond2 || cond3 || cond4);
            // If contains origin and does not accept, it is a nutcase (should we handle it separately?)
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
                          const T max_toi) {
            return max_toi;
        }

        template <int nxe, typename T, typename I>
        T narrow_phase_vf(const size_t noverlaps,
                          const I* const SCCD_RESTRICT voveralp,
                          const I* const SCCD_RESTRICT foveralp,
                          // Geometric data
                          T** const SCCD_RESTRICT v0,
                          T** const SCCD_RESTRICT v1,
                          const size_t face_stride,
                          I** const SCCD_RESTRICT faces,
                          const T max_toi) {
            return max_toi;
        }

    }  // namespace device
}  // namespace sccd

#define INSTANTIATE_NARROW_PHASE_EE(T, I)                                                  \
    template T sccd::device::narrow_phase_ee<T, I>(const size_t noverlaps,                 \
                                                   const I* const SCCD_RESTRICT e0overalp, \
                                                   const I* const SCCD_RESTRICT e1overalp, \
                                                   T** const SCCD_RESTRICT v0,             \
                                                   T** const SCCD_RESTRICT v1,             \
                                                   const size_t edge_stride,               \
                                                   I** const SCCD_RESTRICT edges,          \
                                                   const T max_toi);

#define INSTANTIATE_NARROW_PHASE_VF(NXE, T, I)                                                 \
    template T sccd::device::narrow_phase_vf<NXE, T, I>(const size_t noverlaps,                \
                                                        const I* const SCCD_RESTRICT voveralp, \
                                                        const I* const SCCD_RESTRICT foveralp, \
                                                        T** const SCCD_RESTRICT v0,            \
                                                        T** const SCCD_RESTRICT v1,            \
                                                        const size_t face_stride,              \
                                                        I** const SCCD_RESTRICT faces,         \
                                                        const T max_toi);

INSTANTIATE_NARROW_PHASE_EE(float, int32_t);
INSTANTIATE_NARROW_PHASE_EE(float, int64_t);
INSTANTIATE_NARROW_PHASE_EE(double, int32_t);
INSTANTIATE_NARROW_PHASE_EE(double, int64_t);

INSTANTIATE_NARROW_PHASE_VF(3, float, int32_t);
INSTANTIATE_NARROW_PHASE_VF(3, float, int64_t);
INSTANTIATE_NARROW_PHASE_VF(3, double, int32_t);
INSTANTIATE_NARROW_PHASE_VF(3, double, int64_t);

#undef INSTANTIATE_NARROW_PHASE_EE
#undef INSTANTIATE_NARROW_PHASE_VF