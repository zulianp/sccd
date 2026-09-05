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
                                                           const ptrdiff_t stride,
                                                           I (&out)[nxe]) {
                for (int v = 0; v < nxe; ++v) {
                    out[v] = elements[v][elem_idx * stride];
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

                for (unsigned stride = blockDim.x / 2; stride > 0; stride >>= 1) {
                    if (threadIdx.x < stride) {
                        const unsigned o = threadIdx.x + stride;
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
                                             ptrdiff_t* const SCCD_RESTRICT cellptr) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (i >= n) return;

                const int a = cell0<T>(grid, lo0[i]);
                const int b = cell0<T>(grid, hi0[i]);
                const int c = cell1<T>(grid, lo1[i]);
                const int d = cell1<T>(grid, hi1[i]);

                for (int j = c; j <= d; ++j) {
                    for (int k = a; k <= b; ++k) {
                        // +1 so the array is already the CRS row pointer after an
                        // inclusive scan, matching how ccdptr is built.
                        atomicAdd((unsigned long long*)&cellptr[cell_of<T>(grid, k, j) + 1], 1ull);
                    }
                }
            }

            template <typename T, typename I>
            __global__ void bin_fill_kernel(const ptrdiff_t n,
                                            const T* const SCCD_RESTRICT lo0,
                                            const T* const SCCD_RESTRICT hi0,
                                            const T* const SCCD_RESTRICT lo1,
                                            const T* const SCCD_RESTRICT hi1,
                                            const Cell2DGridD<T> grid,
                                            const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                            I* const SCCD_RESTRICT cellidx,
                                            ptrdiff_t* const SCCD_RESTRICT cursor) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (i >= n) return;

                const int a = cell0<T>(grid, lo0[i]);
                const int b = cell0<T>(grid, hi0[i]);
                const int c = cell1<T>(grid, lo1[i]);
                const int d = cell1<T>(grid, hi1[i]);

                for (int j = c; j <= d; ++j) {
                    for (int k = a; k <= b; ++k) {
                        const ptrdiff_t cell = cell_of<T>(grid, k, j);
                        const unsigned long long at =
                            atomicAdd((unsigned long long*)&cursor[cell], 1ull);
                        cellidx[cellptr[cell] + (ptrdiff_t)at] = (I)i;
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
                                                                    const ptrdiff_t second_stride,
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
                                load_ev<(S > 1 ? S : 1), I>(second_elements, jidx, second_stride, sev);
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
                                                  const ptrdiff_t first_stride,
                                                  I** const SCCD_RESTRICT first_elements,
                                                  T** const SCCD_RESTRICT second_aabbs,
                                                  I* const SCCD_RESTRICT second_idx,
                                                  const ptrdiff_t second_stride,
                                                  I** const SCCD_RESTRICT second_elements,
                                                  const Cell2DGridD<T> grid,
                                                  const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                  const I* const SCCD_RESTRICT cellidx,
                                                  ptrdiff_t* const SCCD_RESTRICT ccdptr) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi == 0) ccdptr[0] = 0;
                if (fi >= first_count) return;

                I ev[first_nxe];
                load_ev<first_nxe, I>(first_elements, first_idx[fi], first_stride, ev);

                ptrdiff_t count = 0;
                for_each_partner<first_nxe, second_nxe, T, I>(first_aabbs,
                                                              fi,
                                                              second_aabbs,
                                                              second_idx,
                                                              second_stride ? second_elements : second_elements,
                                                              second_stride,
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
                                                    const ptrdiff_t first_stride,
                                                    I** const SCCD_RESTRICT first_elements,
                                                    T** const SCCD_RESTRICT second_aabbs,
                                                    I* const SCCD_RESTRICT second_idx,
                                                    const ptrdiff_t second_stride,
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
                load_ev<first_nxe, I>(first_elements, first_idxi, first_stride, ev);

                ptrdiff_t at = ccdptr[fi];
                for_each_partner<first_nxe, second_nxe, T, I>(first_aabbs,
                                                              fi,
                                                              second_aabbs,
                                                              second_idx,
                                                              second_elements,
                                                              second_stride,
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
                                              const ptrdiff_t stride,
                                              I** const SCCD_RESTRICT elements,
                                              const Cell2DGridD<T> grid,
                                              const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                              const I* const SCCD_RESTRICT cellidx,
                                              ptrdiff_t* const SCCD_RESTRICT ccdptr) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi == 0) ccdptr[0] = 0;
                if (fi >= element_count) return;

                I ev[nxe];
                load_ev<nxe, I>(elements, idx[fi], stride, ev);

                ptrdiff_t count = 0;
                for_each_partner<nxe, nxe, T, I>(aabbs,
                                                 fi,
                                                 aabbs,
                                                 idx,
                                                 elements,
                                                 stride,
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
                                                const ptrdiff_t stride,
                                                I** const SCCD_RESTRICT elements,
                                                const Cell2DGridD<T> grid,
                                                const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                const I* const SCCD_RESTRICT cellidx,
                                                const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                I* const SCCD_RESTRICT out0,
                                                I* const SCCD_RESTRICT out1) {
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                if (fi >= element_count) return;

                const I idxi = idx[fi];
                I ev[nxe];
                load_ev<nxe, I>(elements, idxi, stride, ev);

                ptrdiff_t at = ccdptr[fi];
                for_each_partner<nxe, nxe, T, I>(aabbs,
                                                 fi,
                                                 aabbs,
                                                 idx,
                                                 elements,
                                                 stride,
                                                 ev,
                                                 grid,
                                                 cellptr,
                                                 cellidx,
                                                 true,
                                                 [&](const ptrdiff_t, const I jidx) {
                                                     out0[at] = idxi < jidx ? idxi : jidx;
                                                     out1[at] = idxi < jidx ? jidx : idxi;
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

        template <typename T, typename I>
        void cell2d_setup_and_count(const ptrdiff_t n,
                                    T** const SCCD_RESTRICT aabbs,
                                    Cell2DGridD<T>& grid,
                                    ptrdiff_t* const SCCD_RESTRICT cellptr,
                                    ptrdiff_t* const SCCD_RESTRICT span_count) {
            SCCD_CUDA_LAST_ERROR();
            if (n <= 0) {
                *span_count = 0;
                return;
            }

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
            if (want0 < 1) want0 = 1;
            if (want1 < 1) want1 = 1;
            const double cap = 4.0 * (double)n;
            if (want0 * want1 > cap) {
                const double sc = sqrt(cap / (want0 * want1));
                want0 = want0 * sc < 1 ? 1 : want0 * sc;
                want1 = want1 * sc < 1 ? 1 : want1 * sc;
            }
            grid.n0 = (int)want0;
            grid.n1 = (int)want1;
            grid.min0 = h_min[grid.axis0];
            grid.min1 = h_min[grid.axis1];
            grid.inv0 = (T)grid.n0 / (s0 * (T)1.0000001);
            grid.inv1 = (T)grid.n1 / (s1 * (T)1.0000001);

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
                                                        cellptr);
            SCCD_CUDA_LAST_ERROR();

            detail::inclusive_sum<T>(cellptr, ncells + 1);
            SCCD_CHECK_CUDA(cudaMemcpy(span_count, cellptr + ncells, sizeof(ptrdiff_t), cudaMemcpyDeviceToHost));
        }

        template <typename T, typename I>
        void cell2d_fill(const ptrdiff_t n,
                         T** const SCCD_RESTRICT aabbs,
                         const Cell2DGridD<T>& grid,
                         const ptrdiff_t* const SCCD_RESTRICT cellptr,
                         I* const SCCD_RESTRICT cellidx,
                         ptrdiff_t* const SCCD_RESTRICT cursor) {
            if (n <= 0) return;
            SCCD_CHECK_CUDA(cudaMemset(cursor, 0, sizeof(ptrdiff_t) * (size_t)grid.ncells()));

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::bin_fill_kernel<T, I><<<gridsz, block>>>(n,
                                                          soa_device_row<T>(aabbs, grid.axis0),
                                                          soa_device_row<T>(aabbs, SCCD_DIM + grid.axis0),
                                                          soa_device_row<T>(aabbs, grid.axis1),
                                                          soa_device_row<T>(aabbs, SCCD_DIM + grid.axis1),
                                                          grid,
                                                          cellptr,
                                                          cellidx,
                                                          cursor);
            SCCD_CUDA_LAST_ERROR();
        }

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
                                                                                       first_stride,
                                                                                       first_elements,
                                                                                       second_aabbs,
                                                                                       second_idx,
                                                                                       second_stride,
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
                                     I* const SCCD_RESTRICT second_out) {
            if (first_count <= 0) return;
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((first_count + block.x - 1) / block.x);
            detail::collect_overlaps_kernel<first_nxe, second_nxe, T, I><<<gridsz, block>>>(first_count,
                                                                                         first_aabbs,
                                                                                         first_idx,
                                                                                         first_stride,
                                                                                         first_elements,
                                                                                         second_aabbs,
                                                                                         second_idx,
                                                                                         second_stride,
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
                                        const ptrdiff_t stride,
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
                element_count, aabbs, idx, stride, elements, grid, cellptr, cellidx, ccdptr);
            SCCD_CUDA_LAST_ERROR();
            detail::inclusive_sum<T>(ccdptr, element_count + 1);
        }

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
                                          I* const SCCD_RESTRICT out1) {
            if (element_count <= 0) return;
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((element_count + block.x - 1) / block.x);
            detail::collect_self_kernel<nxe, T, I><<<gridsz, block>>>(
                element_count, aabbs, idx, stride, elements, grid, cellptr, cellidx, ccdptr, out0, out1);
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

#define SCCD_C2D_INSTANTIATE(T, I)                                                             \
    template void sccd::device::cell2d_setup_and_count<T, I>(                                  \
        const ptrdiff_t, T**, sccd::device::Cell2DGridD<T>&, ptrdiff_t*, ptrdiff_t*);          \
    template void sccd::device::cell2d_fill<T, I>(const ptrdiff_t,                             \
                                                  T**,                                         \
                                                  const sccd::device::Cell2DGridD<T>&,         \
                                                  const ptrdiff_t*,                            \
                                                  I*,                                          \
                                                  ptrdiff_t*);                                 \
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
