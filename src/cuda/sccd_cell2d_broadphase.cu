#include "sccd_cell2d_broadphase.cuh"

// soa_device_row: the SoA row accessor shared with the sweep implementation.
#include "sccd_broadphase.cuh"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include <cub/device/device_reduce.cuh>
#include <cub/device/device_scan.cuh>

#include "sccd_base.hpp"
#include "sccd_cuda_base.cuh"
#include "sccd_device_workspace.cuh"

#define SCCD_C2D_N_WARPS_PER_BLOCK 8

namespace sccd {
    namespace device {

        namespace detail {

            template <typename T>
            static __device__ __forceinline__ int clamp_cell(const T v, const T lo, const T inv, const int n) {
                const T f = (v - lo) * inv;
                // The negated form also rejects NaN, which a plain f < 0 would let
                // through and turn into an out-of-range cell index.
                if (!(f > T(0))) return 0;
                const int c = (int)f;
                return c >= n ? n - 1 : c;
            }

            template <typename T>
            static __device__ __forceinline__ int cell0(const Cell2DGridD<T>& g, const T v) {
                return clamp_cell<T>(v, g.min0, g.inv0, g.n0);
            }

            template <typename T>
            static __device__ __forceinline__ int cell1(const Cell2DGridD<T>& g, const T v) {
                return clamp_cell<T>(v, g.min1, g.inv1, g.n1);
            }

            template <typename T>
            static __device__ __forceinline__ ptrdiff_t cell_of(const Cell2DGridD<T>& g, const int c0, const int c1) {
                return (ptrdiff_t)c1 * g.n0 + c0;
            }

            template <typename T>
            static __device__ __forceinline__ bool disjoint(const T aminx,
                                                            const T aminy,
                                                            const T aminz,
                                                            const T amaxx,
                                                            const T amaxy,
                                                            const T amaxz,
                                                            const T bminx,
                                                            const T bminy,
                                                            const T bminz,
                                                            const T bmaxx,
                                                            const T bmaxy,
                                                            const T bmaxz) {
                return aminx > bmaxx || aminy > bmaxy || aminz > bmaxz || bminx > amaxx || bminy > amaxy ||
                       bminz > amaxz;
            }

            template <int nxe, typename I>
            static __device__ __forceinline__ void load_ev(I** const SCCD_RESTRICT elements,
                                                           const I elem_idx,
                                                           const ptrdiff_t element_stride,
                                                           I (&out)[nxe]) {
                for (int v = 0; v < nxe; ++v) {
                    out[v] = elements[v][elem_idx * element_stride];
                }
            }

            template <int n1, int n2, typename I>
            static __device__ __forceinline__ bool shares_vertex(const I (&a)[n1], const I (&b)[n2]) {
                for (int i = 0; i < n1; ++i) {
                    for (int j = 0; j < n2; ++j) {
                        if (a[i] == b[j]) return true;
                    }
                }
                return false;
            }

            /**
             * \brief Per-axis min, max and mean extent, for sizing the grid.
             *
             * One block with a shared-array tree reduction rather than atomics:
             * there is no atomic_max for floating point here, and this runs once
             * per axis over arrays the later passes traverse many times, so its
             * cost is not worth optimising.
             */
            template <typename T>
            __global__ void grid_stats_kernel(const ptrdiff_t n,
                                              const T* const SCCD_RESTRICT lo,
                                              const T* const SCCD_RESTRICT hi,
                                              T* const SCCD_RESTRICT out_min,
                                              T* const SCCD_RESTRICT out_max,
                                              T* const SCCD_RESTRICT out_sum) {
                extern __shared__ char s_raw[];
                T* const s_min = (T*)s_raw;
                T* const s_max = s_min + blockDim.x;
                T* const s_sum = s_max + blockDim.x;

                T tmin = lo[0], tmax = hi[0], tsum = T(0);
                for (ptrdiff_t i = threadIdx.x; i < n; i += blockDim.x) {
                    tmin = lo[i] < tmin ? lo[i] : tmin;
                    tmax = hi[i] > tmax ? hi[i] : tmax;
                    tsum += hi[i] - lo[i];
                }
                s_min[threadIdx.x] = tmin;
                s_max[threadIdx.x] = tmax;
                s_sum[threadIdx.x] = tsum;
                __syncthreads();

                for (unsigned element_stride = blockDim.x / 2; element_stride > 0; element_stride >>= 1) {
                    if (threadIdx.x < element_stride) {
                        const unsigned o = threadIdx.x + element_stride;
                        s_min[threadIdx.x] = s_min[o] < s_min[threadIdx.x] ? s_min[o] : s_min[threadIdx.x];
                        s_max[threadIdx.x] = s_max[o] > s_max[threadIdx.x] ? s_max[o] : s_max[threadIdx.x];
                        s_sum[threadIdx.x] += s_sum[o];
                    }
                    __syncthreads();
                }

                if (threadIdx.x == 0) {
                    *out_min = s_min[0];
                    *out_max = s_max[0];
                    *out_sum = s_sum[0];
                }
            }

            template <typename T>
            __global__ void bin_count_kernel(const ptrdiff_t n,
                                             const T* const SCCD_RESTRICT lo0,
                                             const T* const SCCD_RESTRICT hi0,
                                             const T* const SCCD_RESTRICT lo1,
                                             const T* const SCCD_RESTRICT hi1,
                                             const Cell2DGridD<T> grid,
                                             ptrdiff_t* const SCCD_RESTRICT cellptr,
                                             int* const SCCD_RESTRICT ranges) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (i >= n) return;

                const int a = cell0<T>(grid, lo0[i]);
                const int b = cell0<T>(grid, hi0[i]);
                const int c = cell1<T>(grid, lo1[i]);
                const int d = cell1<T>(grid, hi1[i]);

                // The scatter pass reads these back rather than recomputing
                // them, so the two passes cannot disagree about which cells a
                // box covers however the arithmetic is scheduled.
                if (i == 0) {
                    ranges[4 * n + 0] = grid.n0;
                    ranges[4 * n + 1] = grid.n1;
                    ranges[4 * n + 2] = grid.axis0;
                    ranges[4 * n + 3] = grid.axis1;
                }
                ranges[4 * i + 0] = a;
                ranges[4 * i + 1] = b;
                ranges[4 * i + 2] = c;
                ranges[4 * i + 3] = d;

                for (int j = c; j <= d; ++j) {
                    for (int k = a; k <= b; ++k) {
                        // +1 so the array is already the CRS row pointer after an
                        // inclusive scan, matching how ccdptr is built.
                        atomicAdd((unsigned long long*)&cellptr[cell_of<T>(grid, k, j) + 1], 1ull);
                    }
                }
            }

            template <typename I>
            __global__ void fill_identity_kernel(const ptrdiff_t n, I* const SCCD_RESTRICT idx) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (i < n) idx[i] = (I)i;
            }

            template <typename T, typename I>
            __global__ void bin_fill_kernel(const ptrdiff_t n,
                                            const int* const SCCD_RESTRICT ranges,
                                            const Cell2DGridD<T> grid,
                                            const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                            I* const SCCD_RESTRICT cellidx,
                                            const ptrdiff_t capacity,
                                            ptrdiff_t* const SCCD_RESTRICT cursor) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (i >= n) return;

                const int a = ranges[4 * i + 0];
                const int b = ranges[4 * i + 1];
                const int c = ranges[4 * i + 2];
                const int d = ranges[4 * i + 3];
                const ptrdiff_t ncells = grid.ncells();
                if (i == 0 && (ranges[4 * n + 0] != grid.n0 || ranges[4 * n + 1] != grid.n1 ||
                               ranges[4 * n + 2] != grid.axis0 || ranges[4 * n + 3] != grid.axis1)) {
                    printf("GRID MISMATCH count=%dx%d axes %d,%d  fill=%dx%d axes %d,%d\n",
                           ranges[4 * n + 0], ranges[4 * n + 1], ranges[4 * n + 2], ranges[4 * n + 3],
                           grid.n0, grid.n1, grid.axis0, grid.axis1);
                }

                for (int j = c; j <= d; ++j) {
                    for (int k = a; k <= b; ++k) {
                        const ptrdiff_t cell = cell_of<T>(grid, k, j);
                        // A write outside the array is counted rather than made.
                        // The counting pass reserved exactly one slot per span,
                        // so a bad index here means the two passes disagreed,
                        // which the caller reports instead of corrupting memory.
                        if (cell < 0 || cell >= ncells) {
                            atomicAdd((unsigned long long*)&cursor[ncells], 1ull);
                            continue;
                        }
                        const unsigned long long at =
                            atomicAdd((unsigned long long*)&cursor[cell], 1ull);
                        const ptrdiff_t pos = cellptr[cell] + (ptrdiff_t)at;
                        if (pos < 0 || pos >= capacity) {
                            atomicAdd((unsigned long long*)&cursor[ncells], 1ull);
                            continue;
                        }
                        cellidx[pos] = (I)i;
                    }
                }
            }

            /**
             * \brief Visit each partner of box \p fi exactly once.
             *
             * A box spans several cells, so a pair can be met in several of them.
             * It is attributed to the cell holding the minimum corner of the two
             * boxes' overlap: that corner lies inside both boxes, so both are
             * binned there, and the cell is unique. Two clamps, and no per-pair
             * state -- which is what makes it usable on a GPU, where a hash set or
             * a mark array per thread would not be.
             */
            template <int F, int S, typename T, typename I, typename Visit>
            static __device__ __forceinline__ void for_each_partner(T** const SCCD_RESTRICT first_aabbs,
                                                                    const ptrdiff_t fi,
                                                                    T** const SCCD_RESTRICT second_aabbs,
                                                                    const I* const SCCD_RESTRICT second_idx,
                                                                    I** const SCCD_RESTRICT second_elements,
                                                                    const ptrdiff_t second_element_stride,
                                                                    const I (&ev)[F],
                                                                    const Cell2DGridD<T>& grid,
                                                                    const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                                    const I* const SCCD_RESTRICT cellidx,
                                                                    const bool self_mode,
                                                                    Visit visit) {
                const T aminx = first_aabbs[0][fi], aminy = first_aabbs[1][fi], aminz = first_aabbs[2][fi];
                const T amaxx = first_aabbs[3][fi], amaxy = first_aabbs[4][fi], amaxz = first_aabbs[5][fi];

                const T amin0 = first_aabbs[grid.axis0][fi];
                const T amax0 = first_aabbs[3 + grid.axis0][fi];
                const T amin1 = first_aabbs[grid.axis1][fi];
                const T amax1 = first_aabbs[3 + grid.axis1][fi];

                const int c0b = cell0<T>(grid, amin0), c0e = cell0<T>(grid, amax0);
                const int c1b = cell1<T>(grid, amin1), c1e = cell1<T>(grid, amax1);

                for (int c1 = c1b; c1 <= c1e; ++c1) {
                    for (int c0 = c0b; c0 <= c0e; ++c0) {
                        const ptrdiff_t cell = cell_of<T>(grid, c0, c1);
                        const ptrdiff_t begin = cellptr[cell];
                        const ptrdiff_t end = cellptr[cell + 1];

                        for (ptrdiff_t k = begin; k < end; ++k) {
                            const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                            if (self_mode && j <= fi) continue;

                            if (disjoint<T>(aminx,
                                            aminy,
                                            aminz,
                                            amaxx,
                                            amaxy,
                                            amaxz,
                                            second_aabbs[0][j],
                                            second_aabbs[1][j],
                                            second_aabbs[2][j],
                                            second_aabbs[3][j],
                                            second_aabbs[4][j],
                                            second_aabbs[5][j])) {
                                continue;
                            }

                            const T o0 = amin0 > second_aabbs[grid.axis0][j] ? amin0 : second_aabbs[grid.axis0][j];
                            const T o1 = amin1 > second_aabbs[grid.axis1][j] ? amin1 : second_aabbs[grid.axis1][j];
                            if (cell0<T>(grid, o0) != c0 || cell1<T>(grid, o1) != c1) continue;

                            const I jidx = second_idx[j];
                            bool share = false;
                            if (S > 1) {
                                I sev[S > 1 ? S : 1];
                                load_ev<(S > 1 ? S : 1), I>(second_elements, jidx, second_element_stride, sev);
                                share = shares_vertex<F, (S > 1 ? S : 1), I>(ev, sev);
                            } else {
                                for (int a = 0; a < F; ++a) {
                                    if (ev[a] == jidx) {
                                        share = true;
                                        break;
                                    }
                                }
                            }
                            if (share) continue;

                            visit(j, jidx);
                        }
                    }
                }
            }

            template <int first_nxe, int second_nxe, typename T, typename I>
            __global__ void count_overlaps_kernel(const ptrdiff_t first_count,
                                                  T** const SCCD_RESTRICT first_aabbs,
                                                  I* const SCCD_RESTRICT first_idx,
                                                  const ptrdiff_t first_element_stride,
                                                  I** const SCCD_RESTRICT first_elements,
                                                  T** const SCCD_RESTRICT second_aabbs,
                                                  I* const SCCD_RESTRICT second_idx,
                                                  const ptrdiff_t second_element_stride,
                                                  I** const SCCD_RESTRICT second_elements,
                                                  const Cell2DGridD<T> grid,
                                                  const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                  const I* const SCCD_RESTRICT cellidx,
                                                  ptrdiff_t* const SCCD_RESTRICT ccdptr) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi == 0) ccdptr[0] = 0;
                if (fi >= first_count) return;

                I ev[first_nxe];
                load_ev<first_nxe, I>(first_elements, first_idx[fi], first_element_stride, ev);

                ptrdiff_t count = 0;
                for_each_partner<first_nxe, second_nxe, T, I>(first_aabbs,
                                                              fi,
                                                              second_aabbs,
                                                              second_idx,
                                                              second_element_stride ? second_elements : second_elements,
                                                              second_element_stride,
                                                              ev,
                                                              grid,
                                                              cellptr,
                                                              cellidx,
                                                              false,
                                                              [&](const ptrdiff_t, const I) { ++count; });
                ccdptr[fi + 1] = count;
            }

            template <int first_nxe, int second_nxe, typename T, typename I>
            __global__ void collect_overlaps_kernel(const ptrdiff_t first_count,
                                                    T** const SCCD_RESTRICT first_aabbs,
                                                    I* const SCCD_RESTRICT first_idx,
                                                    const ptrdiff_t first_element_stride,
                                                    I** const SCCD_RESTRICT first_elements,
                                                    T** const SCCD_RESTRICT second_aabbs,
                                                    I* const SCCD_RESTRICT second_idx,
                                                    const ptrdiff_t second_element_stride,
                                                    I** const SCCD_RESTRICT second_elements,
                                                    const Cell2DGridD<T> grid,
                                                    const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                    const I* const SCCD_RESTRICT cellidx,
                                                    const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                    I* const SCCD_RESTRICT first_out,
                                                    I* const SCCD_RESTRICT second_out) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi >= first_count) return;

                const I first_idxi = first_idx[fi];
                I ev[first_nxe];
                load_ev<first_nxe, I>(first_elements, first_idxi, first_element_stride, ev);

                ptrdiff_t at = ccdptr[fi];
                for_each_partner<first_nxe, second_nxe, T, I>(first_aabbs,
                                                              fi,
                                                              second_aabbs,
                                                              second_idx,
                                                              second_elements,
                                                              second_element_stride,
                                                              ev,
                                                              grid,
                                                              cellptr,
                                                              cellidx,
                                                              false,
                                                              [&](const ptrdiff_t, const I jidx) {
                                                                  first_out[at] = first_idxi;
                                                                  second_out[at] = jidx;
                                                                  ++at;
                                                              });
            }

            template <int nxe, typename T, typename I>
            __global__ void count_self_kernel(const ptrdiff_t element_count,
                                              T** const SCCD_RESTRICT aabbs,
                                              I* const SCCD_RESTRICT idx,
                                              const ptrdiff_t element_stride,
                                              I** const SCCD_RESTRICT elements,
                                              const Cell2DGridD<T> grid,
                                              const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                              const I* const SCCD_RESTRICT cellidx,
                                              ptrdiff_t* const SCCD_RESTRICT ccdptr) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi == 0) ccdptr[0] = 0;
                if (fi >= element_count) return;

                I ev[nxe];
                load_ev<nxe, I>(elements, idx[fi], element_stride, ev);

                ptrdiff_t count = 0;
                for_each_partner<nxe, nxe, T, I>(aabbs,
                                                 fi,
                                                 aabbs,
                                                 idx,
                                                 elements,
                                                 element_stride,
                                                 ev,
                                                 grid,
                                                 cellptr,
                                                 cellidx,
                                                 true,
                                                 [&](const ptrdiff_t, const I) { ++count; });
                ccdptr[fi + 1] = count;
            }

            template <int nxe, typename T, typename I>
            __global__ void collect_self_kernel(const ptrdiff_t element_count,
                                                T** const SCCD_RESTRICT aabbs,
                                                I* const SCCD_RESTRICT idx,
                                                const ptrdiff_t element_stride,
                                                I** const SCCD_RESTRICT elements,
                                                const Cell2DGridD<T> grid,
                                                const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                const I* const SCCD_RESTRICT cellidx,
                                                const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                I* const SCCD_RESTRICT first_out,
                                                I* const SCCD_RESTRICT second_out) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi >= element_count) return;

                const I idxi = idx[fi];
                I ev[nxe];
                load_ev<nxe, I>(elements, idxi, element_stride, ev);

                ptrdiff_t at = ccdptr[fi];
                for_each_partner<nxe, nxe, T, I>(aabbs,
                                                 fi,
                                                 aabbs,
                                                 idx,
                                                 elements,
                                                 element_stride,
                                                 ev,
                                                 grid,
                                                 cellptr,
                                                 cellidx,
                                                 true,
                                                 [&](const ptrdiff_t, const I jidx) {
                                                     first_out[at] = idxi < jidx ? idxi : jidx;
                                                     second_out[at] = idxi < jidx ? jidx : idxi;
                                                     ++at;
                                                 });
            }


            // ------------------------------------------------------------------
            // The edge-edge form: one entry per box, at its minimum corner.
            // ------------------------------------------------------------------

            /**
             * \brief The value an empty cell's bound carries.
             *
             * `std::numeric_limits` is not available in device code, so the two
             * types the library instantiates say it themselves. An empty cell
             * holds this, and every comparison against it fails, which is how an
             * empty cell is skipped without a separate case.
             */
            template <typename T>
            static __device__ __forceinline__ T bound_floor();
            template <>
            __device__ __forceinline__ float bound_floor<float>() {
                return -3.402823466e+38f;
            }
            template <>
            __device__ __forceinline__ double bound_floor<double>() {
                return -1.7976931348623157e+308;
            }

            template <typename T>
            __global__ void min_bin_count_kernel(const ptrdiff_t n,
                                                 const T* const SCCD_RESTRICT lo0,
                                                 const T* const SCCD_RESTRICT lo1,
                                                 const Cell2DGridD<T> grid,
                                                 ptrdiff_t* const SCCD_RESTRICT cellptr) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (i >= n) return;
                // One cell, so no range to record for the scatter to read back:
                // both passes call the same clamp on the same two coordinates.
                const ptrdiff_t cell = cell_of<T>(grid, cell0<T>(grid, lo0[i]), cell1<T>(grid, lo1[i]));
                atomicAdd((unsigned long long*)&cellptr[cell + 1], 1ull);
            }

            template <typename T, typename I>
            __global__ void min_bin_fill_kernel(const ptrdiff_t n,
                                                const T* const SCCD_RESTRICT lo0,
                                                const T* const SCCD_RESTRICT lo1,
                                                const Cell2DGridD<T> grid,
                                                const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                I* const SCCD_RESTRICT cellidx,
                                                ptrdiff_t* const SCCD_RESTRICT cursor) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (i >= n) return;
                const ptrdiff_t cell = cell_of<T>(grid, cell0<T>(grid, lo0[i]), cell1<T>(grid, lo1[i]));
                const unsigned long long at = atomicAdd((unsigned long long*)&cursor[cell], 1ull);
                cellidx[cellptr[cell] + (ptrdiff_t)at] = (I)i;
            }

            /** \brief Per cell, the largest upper bound of the boxes it holds, on both grid axes. */
            template <typename T, typename I>
            __global__ void min_cell_bounds_kernel(const ptrdiff_t ncells,
                                                   const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                   const I* const SCCD_RESTRICT cellidx,
                                                   const T* const SCCD_RESTRICT hi0,
                                                   const T* const SCCD_RESTRICT hi1,
                                                   T* const SCCD_RESTRICT cell_hi0,
                                                   T* const SCCD_RESTRICT cell_hi1) {
                const ptrdiff_t c = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (c >= ncells) return;
                const T floor = bound_floor<T>();
                T m0 = floor, m1 = floor;
                for (ptrdiff_t k = cellptr[c]; k < cellptr[c + 1]; ++k) {
                    const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                    if (hi0[j] > m0) m0 = hi0[j];
                    if (hi1[j] > m1) m1 = hi1[j];
                }
                cell_hi0[c] = m0;
                cell_hi1[c] = m1;
            }

            /**
             * \brief Turn each row of per-cell bounds into its running maximum.
             *
             * One thread per row, walking the row. The parallelism is the number
             * of rows rather than of cells, which is the square root of the work,
             * but the work itself is one pass over an array the query then reads
             * many times.
             */
            template <typename T>
            __global__ void min_row_prefix_kernel(const int n0,
                                                  const int n1,
                                                  T* const SCCD_RESTRICT row_prefix) {
                const int r = blockIdx.x * blockDim.x + threadIdx.x;
                if (r >= n1) return;
                const ptrdiff_t base = (ptrdiff_t)r * n0;
                T run = row_prefix[base];
                for (int c = 1; c < n0; ++c) {
                    const T v = row_prefix[base + c];
                    run = v > run ? v : run;
                    row_prefix[base + c] = run;
                }
            }

            /** \brief First column of a row whose running maximum reaches \p v. */
            template <typename T>
            static __device__ __forceinline__ int lower_bound_row(const T* const SCCD_RESTRICT p,
                                                                  const int count,
                                                                  const T v) {
                int lo = 0, hi = count;
                while (lo < hi) {
                    const int mid = (lo + hi) >> 1;
                    if (p[mid] < v) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                return lo;
            }

            /**
             * \brief Visit each partner of box \p fi exactly once, walking forward.
             *
             * The host form is `for_each_forward_self_partner` in
             * `src/broadphase/sccd_broadphase_cell2d.hpp` and the reasoning is
             * written up there: one entry per box means a partner is met at most
             * once, row-major linear index is a total order so a forward walk is
             * complete over it, and the box's own cell is the only one where two
             * boxes share an index and the index has to decide.
             */
            template <int NXE, typename T, typename I, typename Visit>
            static __device__ __forceinline__ void for_each_forward_self_partner(
                T** const SCCD_RESTRICT aabbs,
                const ptrdiff_t fi,
                I* const SCCD_RESTRICT idx,
                I** const SCCD_RESTRICT elements,
                const ptrdiff_t element_stride,
                const I (&ev)[NXE],
                const Cell2DGridD<T>& grid,
                const ptrdiff_t* const SCCD_RESTRICT cellptr,
                const I* const SCCD_RESTRICT cellidx,
                const T* const SCCD_RESTRICT row_prefix,
                const T* const SCCD_RESTRICT cell_hi1,
                Visit&& visit) {
                const T aminx = aabbs[0][fi], aminy = aabbs[1][fi], aminz = aabbs[2][fi];
                const T amaxx = aabbs[3][fi], amaxy = aabbs[4][fi], amaxz = aabbs[5][fi];

                const T amin0 = aabbs[grid.axis0][fi];
                const T amin1 = aabbs[grid.axis1][fi];

                const int c0b = cell0<T>(grid, amin0);
                const int c0e = cell0<T>(grid, aabbs[SCCD_DIM + grid.axis0][fi]);
                const int c1b = cell1<T>(grid, amin1);
                const int c1e = cell1<T>(grid, aabbs[SCCD_DIM + grid.axis1][fi]);
                const ptrdiff_t own = cell_of<T>(grid, c0b, c1b);

                // A cell's entries, with OWN deciding whether the index still has
                // to break a tie. Everything the walk reads after its own cell is
                // ahead of it, so there the question does not arise.
                const auto scan = [&](const ptrdiff_t cell, const bool is_own) {
                    for (ptrdiff_t k = cellptr[cell]; k < cellptr[cell + 1]; ++k) {
                        const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                        if (is_own && j <= fi) continue;
                        if (disjoint<T>(aminx, aminy, aminz, amaxx, amaxy, amaxz,
                                        aabbs[0][j], aabbs[1][j], aabbs[2][j],
                                        aabbs[3][j], aabbs[4][j], aabbs[5][j])) {
                            continue;
                        }
                        const I jidx = idx[j];
                        I sev[NXE];
                        load_ev<NXE, I>(elements, jidx, element_stride, sev);
                        if (shares_vertex<NXE, NXE>(ev, sev)) continue;
                        visit(j, jidx);
                    }
                };

                // Its own cell first. Both bounds hold this box's own maximum, so
                // neither can rule the cell out and neither is worth reading.
                scan(own, true);

                for (int c1 = c1b; c1 <= c1e; ++c1) {
                    const ptrdiff_t row = (ptrdiff_t)c1 * grid.n0;
                    int c0;
                    if (c1 == c1b) {
                        // This box's own maximum is in its row's prefix at column
                        // c0b, so a search here can never return a later column
                        // than the cell just read.
                        c0 = c0b + 1;
                    } else {
                        c0 = lower_bound_row<T>(row_prefix + row, c0e + 1, amin0);
                    }
                    for (; c0 <= c0e; ++c0) {
                        const ptrdiff_t cell = row + c0;
                        if (cell_hi1[cell] < amin1) continue;
                        scan(cell, false);
                    }
                }
            }

            template <int nxe, typename T, typename I>
            __global__ void min_count_self_kernel(const ptrdiff_t element_count,
                                                  T** const SCCD_RESTRICT aabbs,
                                                  I* const SCCD_RESTRICT idx,
                                                  const ptrdiff_t element_stride,
                                                  I** const SCCD_RESTRICT elements,
                                                  const Cell2DGridD<T> grid,
                                                  const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                  const I* const SCCD_RESTRICT cellidx,
                                                  const T* const SCCD_RESTRICT row_prefix,
                                                  const T* const SCCD_RESTRICT cell_hi1,
                                                  ptrdiff_t* const SCCD_RESTRICT ccdptr) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi == 0) ccdptr[0] = 0;
                if (fi >= element_count) return;

                I ev[nxe];
                load_ev<nxe, I>(elements, idx[fi], element_stride, ev);

                ptrdiff_t count = 0;
                for_each_forward_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx,
                    row_prefix, cell_hi1, [&](const ptrdiff_t, const I) { ++count; });
                ccdptr[fi + 1] = count;
            }

            template <int nxe, typename T, typename I>
            __global__ void min_collect_self_kernel(const ptrdiff_t element_count,
                                                    T** const SCCD_RESTRICT aabbs,
                                                    I* const SCCD_RESTRICT idx,
                                                    const ptrdiff_t element_stride,
                                                    I** const SCCD_RESTRICT elements,
                                                    const Cell2DGridD<T> grid,
                                                    const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                    const I* const SCCD_RESTRICT cellidx,
                                                    const T* const SCCD_RESTRICT row_prefix,
                                                    const T* const SCCD_RESTRICT cell_hi1,
                                                    const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                    I* const SCCD_RESTRICT first_out,
                                                    I* const SCCD_RESTRICT second_out) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi >= element_count) return;

                const I idxi = idx[fi];
                I ev[nxe];
                load_ev<nxe, I>(elements, idxi, element_stride, ev);

                ptrdiff_t at = ccdptr[fi];
                for_each_forward_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx,
                    row_prefix, cell_hi1, [&](const ptrdiff_t, const I jidx) {
                        first_out[at] = idxi < jidx ? idxi : jidx;
                        second_out[at] = idxi < jidx ? jidx : idxi;
                        ++at;
                    });
            }

            template <typename T>
            static void inclusive_sum(ptrdiff_t* const data, const ptrdiff_t n) {
                size_t bytes = 0;
                SCCD_CHECK_CUDA(cub::DeviceScan::InclusiveSum(nullptr, bytes, data, data, n));
                void* const tmp = workspace(WorkspaceSlot::TempStorage).get(bytes);
                SCCD_CHECK_CUDA(cub::DeviceScan::InclusiveSum(tmp, bytes, data, data, n));
            }

        }  // namespace detail

        /**
         * \brief Shrink a grid until it holds at most `4 * n` cells.
         *
         * The caller sizes its cell array from this bound, so it has to hold
         * exactly, not approximately. The host cell list enforces it the same
         * way, in `cap_cells`.
         */
        static void cap_cells_device(const ptrdiff_t n, int& n0, int& n1) {
            ptrdiff_t cap = 4 * n;
            if (cap < 1) cap = 1;
            if (n0 < 1) n0 = 1;
            if (n1 < 1) n1 = 1;
            if ((ptrdiff_t)n0 > cap) n0 = (int)cap;
            if ((ptrdiff_t)n1 > cap) n1 = (int)cap;
            if ((ptrdiff_t)n0 * (ptrdiff_t)n1 <= cap) return;
            if (n0 >= n1) {
                const ptrdiff_t v = cap / (ptrdiff_t)n1;
                n0 = (int)(v < 1 ? 1 : v);
            } else {
                const ptrdiff_t v = cap / (ptrdiff_t)n0;
                n1 = (int)(v < 1 ? 1 : v);
            }
        }

        /**
         * \brief Choose the two axes and the cell size for \p n boxes.
         *
         * Shared by both binnings, which differ only in what they then put
         * into the grid: every cell a box touches, or the one its minimum
         * corner is in.
         */
        template <typename T>
        static void size_grid_device(const ptrdiff_t n,
                                     T** const SCCD_RESTRICT aabbs,
                                     Cell2DGridD<T>& grid) {
            // Axis choice is on the host: it needs three reductions and then a
            // decision, and the arrays are small enough that a kernel per axis is
            // cheaper than the launch overhead of doing it any other way.
            T* const stats = workspace(WorkspaceSlot::Scratch).get_as<T>(9);
            T h_min[3], h_max[3], h_sum[3];
            for (int d = 0; d < 3; ++d) {
                detail::grid_stats_kernel<T><<<1, 256, 3 * 256 * sizeof(T)>>>(n,
                                                      soa_device_row<T>(aabbs, d),
                                                      soa_device_row<T>(aabbs, SCCD_DIM + d),
                                                      stats + 3 * d,
                                                      stats + 3 * d + 1,
                                                      stats + 3 * d + 2);
            }
            SCCD_CUDA_LAST_ERROR();
            T host_stats[9];
            SCCD_CHECK_CUDA(cudaMemcpy(host_stats, stats, sizeof(T) * 9, cudaMemcpyDeviceToHost));
            for (int d = 0; d < 3; ++d) {
                h_min[d] = host_stats[3 * d];
                h_max[d] = host_stats[3 * d + 1];
                h_sum[d] = host_stats[3 * d + 2];
            }

            // Two widest axes by extent. The host version uses centre variance;
            // extent is the same decision on this geometry and needs no second
            // pass over the data.
            int order[3] = {0, 1, 2};
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 3; ++j) {
                    if ((h_max[order[j]] - h_min[order[j]]) > (h_max[order[i]] - h_min[order[i]])) {
                        const int t = order[i];
                        order[i] = order[j];
                        order[j] = t;
                    }
                }
            }
            grid.axis0 = order[0];
            grid.axis1 = order[1];

            const T span0 = h_max[grid.axis0] - h_min[grid.axis0];
            const T span1 = h_max[grid.axis1] - h_min[grid.axis1];
            const T s0 = span0 > T(0) ? span0 : T(1);
            const T s1 = span1 > T(0) ? span1 : T(1);
            const T mean0 = h_sum[grid.axis0] / (T)n > T(0) ? h_sum[grid.axis0] / (T)n : s0;
            const T mean1 = h_sum[grid.axis1] / (T)n > T(0) ? h_sum[grid.axis1] / (T)n : s1;

            double want0 = (double)(s0 / mean0);
            double want1 = (double)(s1 / mean1);
            // Negated so a NaN extent lands on one cell rather than on whatever
            // the cast to int makes of it.
            if (!(want0 >= 1)) want0 = 1;
            if (!(want1 >= 1)) want1 = 1;
            const double cap = 4.0 * (double)n;
            if (want0 * want1 > cap) {
                const double sc = sqrt(cap / (want0 * want1));
                want0 = want0 * sc < 1 ? 1 : want0 * sc;
                want1 = want1 * sc < 1 ? 1 : want1 * sc;
            }
            // An infinite or still-huge side would be undefined as an int, so
            // both are brought inside the cap before the cast.
            if (want0 > cap) want0 = cap;
            if (want1 > cap) want1 = cap;
            grid.n0 = (int)want0;
            grid.n1 = (int)want1;
            // Scaling both sides by the same factor does not enforce the cap on
            // a scene whose two axes differ wildly: a side scaled below one is
            // raised back to one, and the product grows past the cap again by
            // exactly that much. n-body-simulation, puffer-ball and rod-twist
            // are all such scenes, and the caller sizes cellptr from this bound,
            // so the bound is enforced on the integers it is sized by.
            cap_cells_device(n, grid.n0, grid.n1);
            grid.min0 = h_min[grid.axis0];
            grid.min1 = h_min[grid.axis1];
            grid.inv0 = (T)grid.n0 / (s0 * (T)1.0000001);
            grid.inv1 = (T)grid.n1 / (s1 * (T)1.0000001);

        }

        template <typename T, typename I>
        void cell2d_setup_and_count(const ptrdiff_t n,
                                    T** const SCCD_RESTRICT aabbs,
                                    Cell2DGridD<T>& grid,
                                    ptrdiff_t* const SCCD_RESTRICT cellptr,
                                    int* const SCCD_RESTRICT ranges,
                                    ptrdiff_t* const SCCD_RESTRICT span_count) {
            SCCD_CUDA_LAST_ERROR();
            if (n <= 0) {
                *span_count = 0;
                return;
            }

            size_grid_device<T>(n, aabbs, grid);

            const ptrdiff_t ncells = grid.ncells();
            SCCD_CHECK_CUDA(cudaMemset(cellptr, 0, sizeof(ptrdiff_t) * (size_t)(ncells + 1)));

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::bin_count_kernel<T><<<gridsz, block>>>(n,
                                                        soa_device_row<T>(aabbs, grid.axis0),
                                                        soa_device_row<T>(aabbs, SCCD_DIM + grid.axis0),
                                                        soa_device_row<T>(aabbs, grid.axis1),
                                                        soa_device_row<T>(aabbs, SCCD_DIM + grid.axis1),
                                                        grid,
                                                        cellptr,
                                                        ranges);
            SCCD_CUDA_LAST_ERROR();

            detail::inclusive_sum<T>(cellptr, ncells + 1);
            SCCD_CHECK_CUDA(cudaMemcpy(span_count, cellptr + ncells, sizeof(ptrdiff_t), cudaMemcpyDeviceToHost));
        }

        template <typename I>
        void fill_identity(const ptrdiff_t n, I* const SCCD_RESTRICT idx) {
            if (n <= 0) return;
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::fill_identity_kernel<I><<<gridsz, block>>>(n, idx);
            SCCD_CUDA_LAST_ERROR();
        }

        template <typename T, typename I>
        void cell2d_fill(const ptrdiff_t n,
                         T** const SCCD_RESTRICT aabbs,
                         const Cell2DGridD<T>& grid,
                         const ptrdiff_t* const SCCD_RESTRICT cellptr,
                         const int* const SCCD_RESTRICT ranges,
                         I* const SCCD_RESTRICT cellidx,
                         const ptrdiff_t capacity,
                         ptrdiff_t* const SCCD_RESTRICT cursor,
                         ptrdiff_t* const SCCD_RESTRICT out_rejected) {
            if (n <= 0) {
                if (out_rejected) *out_rejected = 0;
                return;
            }
            // One past the cells: the slot the kernel counts rejected writes in.
            SCCD_CHECK_CUDA(cudaMemset(cursor, 0, sizeof(ptrdiff_t) * (size_t)(grid.ncells() + 1)));

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::bin_fill_kernel<T, I><<<gridsz, block>>>(n,
                                                          ranges,
                                                          grid,
                                                          cellptr,
                                                          cellidx,
                                                          capacity,
                                                          cursor);
            SCCD_CUDA_LAST_ERROR();
            if (out_rejected) {
                SCCD_CHECK_CUDA(cudaMemcpy(out_rejected, cursor + grid.ncells(),
                                           sizeof(ptrdiff_t), cudaMemcpyDeviceToHost));
            }
        }

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
                                   ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            SCCD_CUDA_LAST_ERROR();
            if (first_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((first_count + block.x - 1) / block.x);
            detail::count_overlaps_kernel<first_nxe, second_nxe, T, I><<<gridsz, block>>>(first_count,
                                                                                       first_aabbs,
                                                                                       first_idx,
                                                                                       first_element_stride,
                                                                                       first_elements,
                                                                                       second_aabbs,
                                                                                       second_idx,
                                                                                       second_element_stride,
                                                                                       second_elements,
                                                                                       grid,
                                                                                       cellptr,
                                                                                       cellidx,
                                                                                       ccdptr);
            SCCD_CUDA_LAST_ERROR();
            detail::inclusive_sum<T>(ccdptr, first_count + 1);
        }

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
                                     I* const SCCD_RESTRICT second_out) {
            if (first_count <= 0) return;
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((first_count + block.x - 1) / block.x);
            detail::collect_overlaps_kernel<first_nxe, second_nxe, T, I><<<gridsz, block>>>(first_count,
                                                                                         first_aabbs,
                                                                                         first_idx,
                                                                                         first_element_stride,
                                                                                         first_elements,
                                                                                         second_aabbs,
                                                                                         second_idx,
                                                                                         second_element_stride,
                                                                                         second_elements,
                                                                                         grid,
                                                                                         cellptr,
                                                                                         cellidx,
                                                                                         ccdptr,
                                                                                         first_out,
                                                                                         second_out);
            SCCD_CUDA_LAST_ERROR();
        }

        template <int nxe, typename T, typename I>
        void cell2d_count_self_overlaps(const ptrdiff_t element_count,
                                        T** const SCCD_RESTRICT aabbs,
                                        I* const SCCD_RESTRICT idx,
                                        const ptrdiff_t element_stride,
                                        I** const SCCD_RESTRICT elements,
                                        const Cell2DGridD<T>& grid,
                                        const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                        const I* const SCCD_RESTRICT cellidx,
                                        ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            SCCD_CUDA_LAST_ERROR();
            if (element_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((element_count + block.x - 1) / block.x);
            detail::count_self_kernel<nxe, T, I><<<gridsz, block>>>(
                element_count, aabbs, idx, element_stride, elements, grid, cellptr, cellidx, ccdptr);
            SCCD_CUDA_LAST_ERROR();
            detail::inclusive_sum<T>(ccdptr, element_count + 1);
        }

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
                                          I* const SCCD_RESTRICT second_out) {
            if (element_count <= 0) return;
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((element_count + block.x - 1) / block.x);
            detail::collect_self_kernel<nxe, T, I><<<gridsz, block>>>(
                element_count, aabbs, idx, element_stride, elements, grid, cellptr, cellidx, ccdptr, first_out, second_out);
            SCCD_CUDA_LAST_ERROR();
        }

        template <typename T, typename I>
        void cell2dmin_setup_and_count(const ptrdiff_t n,
                                       T** const SCCD_RESTRICT aabbs,
                                       Cell2DGridD<T>& grid,
                                       ptrdiff_t* const SCCD_RESTRICT cellptr) {
            SCCD_CUDA_LAST_ERROR();
            if (n <= 0) return;

            size_grid_device<T>(n, aabbs, grid);

            const ptrdiff_t ncells = grid.ncells();
            SCCD_CHECK_CUDA(cudaMemset(cellptr, 0, sizeof(ptrdiff_t) * (size_t)(ncells + 1)));

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::min_bin_count_kernel<T><<<gridsz, block>>>(
                n, soa_device_row<T>(aabbs, grid.axis0), soa_device_row<T>(aabbs, grid.axis1),
                grid, cellptr);
            SCCD_CUDA_LAST_ERROR();

            detail::inclusive_sum<T>(cellptr, ncells + 1);
        }

        template <typename T, typename I>
        void cell2dmin_fill(const ptrdiff_t n,
                            T** const SCCD_RESTRICT aabbs,
                            const Cell2DGridD<T>& grid,
                            const ptrdiff_t* const SCCD_RESTRICT cellptr,
                            I* const SCCD_RESTRICT cellidx,
                            ptrdiff_t* const SCCD_RESTRICT cursor) {
            if (n <= 0) return;
            SCCD_CHECK_CUDA(cudaMemset(cursor, 0, sizeof(ptrdiff_t) * (size_t)grid.ncells()));

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::min_bin_fill_kernel<T, I><<<gridsz, block>>>(
                n, soa_device_row<T>(aabbs, grid.axis0), soa_device_row<T>(aabbs, grid.axis1),
                grid, cellptr, cellidx, cursor);
            SCCD_CUDA_LAST_ERROR();
        }

        template <typename T, typename I>
        void cell2dmin_bounds(const Cell2DGridD<T>& grid,
                              T** const SCCD_RESTRICT aabbs,
                              const ptrdiff_t* const SCCD_RESTRICT cellptr,
                              const I* const SCCD_RESTRICT cellidx,
                              T* const SCCD_RESTRICT row_prefix,
                              T* const SCCD_RESTRICT cell_hi1) {
            const ptrdiff_t ncells = grid.ncells();
            if (ncells <= 0) return;

            // The per-cell maximum on the first axis is written straight into
            // row_prefix, which the next kernel turns into the running maximum in
            // place: the two are the same array at two stages, as on the host.
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((ncells + block.x - 1) / block.x);
            detail::min_cell_bounds_kernel<T, I><<<gridsz, block>>>(
                ncells, cellptr, cellidx,
                soa_device_row<T>(aabbs, SCCD_DIM + grid.axis0),
                soa_device_row<T>(aabbs, SCCD_DIM + grid.axis1),
                row_prefix, cell_hi1);
            SCCD_CUDA_LAST_ERROR();

            dim3 rowsz((grid.n1 + block.x - 1) / block.x);
            detail::min_row_prefix_kernel<T><<<rowsz, block>>>(grid.n0, grid.n1, row_prefix);
            SCCD_CUDA_LAST_ERROR();
        }

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
                                           ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            if (element_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((element_count + block.x - 1) / block.x);
            detail::min_count_self_kernel<nxe, T, I><<<gridsz, block>>>(
                element_count, aabbs, idx, element_stride, elements, grid, cellptr, cellidx,
                row_prefix, cell_hi1, ccdptr);
            SCCD_CUDA_LAST_ERROR();
            detail::inclusive_sum<T>(ccdptr, element_count + 1);
        }

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
                                             I* const SCCD_RESTRICT second_out) {
            if (element_count <= 0) return;
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((element_count + block.x - 1) / block.x);
            detail::min_collect_self_kernel<nxe, T, I><<<gridsz, block>>>(
                element_count, aabbs, idx, element_stride, elements, grid, cellptr, cellidx,
                row_prefix, cell_hi1, ccdptr, first_out, second_out);
            SCCD_CUDA_LAST_ERROR();
        }

    }  // namespace device
}  // namespace sccd

// The face-vertex pair collectors are instantiated at both supported face
// widths. Quads are a supported topology, so the device path gets <4,1> next to
// <3,1> -- the host cell list was shipped with only the triangle instantiation
// and that made quad meshes fail on the default broad phase.
#define SCCD_C2D_INSTANTIATE_FV(NXE, T, I)                                                     \
    template void sccd::device::cell2d_count_overlaps<NXE, 1, T, I>(                           \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        ptrdiff_t*);                                                                           \
    template void sccd::device::cell2d_collect_overlaps<NXE, 1, T, I>(                         \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        const ptrdiff_t*,                                                                      \
        I*,                                                                                    \
        I*);

// Indexed by the index type alone, so it is instantiated once rather than once
// per scalar type like everything below it.
#define SCCD_C2D_INSTANTIATE_IDX(I) \
    template void sccd::device::fill_identity<I>(const ptrdiff_t, I*);

SCCD_C2D_INSTANTIATE_IDX(int32_t)

#define SCCD_C2D_INSTANTIATE(T, I)                                                             \
    template void sccd::device::cell2d_setup_and_count<T, I>(                                  \
        const ptrdiff_t, T**, sccd::device::Cell2DGridD<T>&, ptrdiff_t*, int*, ptrdiff_t*);          \
    template void sccd::device::cell2d_fill<T, I>(const ptrdiff_t,                             \
                                                  T**,                                         \
                                                  const sccd::device::Cell2DGridD<T>&,         \
                                                  const ptrdiff_t*,                            \
                                                  const int*,                                  \
                                                  I*,                                          \
                                                  const ptrdiff_t,                             \
                                                  ptrdiff_t*,                                  \
                                                  ptrdiff_t*);                                 \
    template void sccd::device::cell2dmin_setup_and_count<T, I>(                                \
        const ptrdiff_t, T**, sccd::device::Cell2DGridD<T>&, ptrdiff_t*);                      \
    template void sccd::device::cell2dmin_fill<T, I>(const ptrdiff_t,                          \
                                                     T**,                                      \
                                                     const sccd::device::Cell2DGridD<T>&,      \
                                                     const ptrdiff_t*,                         \
                                                     I*,                                       \
                                                     ptrdiff_t*);                              \
    template void sccd::device::cell2dmin_bounds<T, I>(const sccd::device::Cell2DGridD<T>&,    \
                                                       T**,                                    \
                                                       const ptrdiff_t*,                       \
                                                       const I*,                               \
                                                       T*,                                     \
                                                       T*);                                    \
    template void sccd::device::cell2dmin_count_self_overlaps<2, T, I>(                        \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        const T*,                                                                              \
        const T*,                                                                              \
        ptrdiff_t*);                                                                           \
    template void sccd::device::cell2dmin_collect_self_overlaps<2, T, I>(                      \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        const T*,                                                                              \
        const T*,                                                                              \
        const ptrdiff_t*,                                                                      \
        I*,                                                                                    \
        I*);                                                                                   \
    SCCD_C2D_INSTANTIATE_FV(3, T, I)                                                           \
    SCCD_C2D_INSTANTIATE_FV(4, T, I)                                                           \
    template void sccd::device::cell2d_count_self_overlaps<2, T, I>(                           \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        ptrdiff_t*);                                                                           \
    template void sccd::device::cell2d_collect_self_overlaps<2, T, I>(                         \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        I*,                                                                                    \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        const ptrdiff_t*,                                                                      \
        I*,                                                                                    \
        I*);

SCCD_C2D_INSTANTIATE(float, int32_t)
SCCD_C2D_INSTANTIATE(double, int32_t)

#undef SCCD_C2D_INSTANTIATE
#undef SCCD_C2D_INSTANTIATE_FV
#undef SCCD_C2D_INSTANTIATE_IDX
