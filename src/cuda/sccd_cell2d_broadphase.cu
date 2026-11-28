#include "sccd_cell2d_broadphase.cuh"

// soa_device_row: the SoA row accessor shared with the sweep implementation.
#include "sccd_broadphase.cuh"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include "sccd_gpu_cub.cuh"

#include "sccd_base.hpp"
#include "sccd_cuda_base.cuh"
#include "sccd_device_workspace.cuh"

#define SCCD_C2D_N_WARPS_PER_BLOCK 8

// The four warp-per-box launches (self count and collect, vertex-face count and
// collect) give every box a whole warp, so their grid holds 32 or 64 times as
// many threads as there are boxes. HIP bounds a grid dimension by its
// work-items, gridDim.x * blockDim.x < 2^32, and a 64-wide wavefront crosses it
// past 2^32 / 64 = 67M boxes: cloth-ball refined to level 5 (94M faces, 141M
// edges) failed its self-overlap count with "invalid configuration argument" on an
// MI300A. Set, such a launch spills into gridDim.y and those kernels number
// their warps over both dimensions. CUDA bounds blocks instead (2^31 - 1 in x),
// which a warp-per-box launch reaches only past 17 billion boxes, so the CUDA
// build keeps its one-dimensional launch and the H100's kernels are unchanged.
#ifndef SCCD_C2D_WARP_GRID_2D
#if defined(SCCD_ENABLE_HIP)
#define SCCD_C2D_WARP_GRID_2D 1
#else
#define SCCD_C2D_WARP_GRID_2D 0
#endif
#endif
#if SCCD_C2D_WARP_GRID_2D
// The most work-items a launch puts in x. Only a test lowers it, so that small
// inputs take the two-dimensional path.
#ifndef SCCD_C2D_WARP_GRID_MAX_ITEMS
#define SCCD_C2D_WARP_GRID_MAX_ITEMS ((((ptrdiff_t)1) << 32) - 1)
#endif
// A warp-per-box kernel's thread index over both grid dimensions.
#define SCCD_C2D_WARP_TID() \
    (((ptrdiff_t)blockIdx.y * gridDim.x + blockIdx.x) * (ptrdiff_t)blockDim.x + threadIdx.x)
// The grid of a warp-per-box launch over n boxes: one dimension while it fits,
// else as many rows of the widest x as it takes. The kernels return for a warp
// past n, so the last row's surplus does nothing.
static inline dim3 sccd_c2d_warp_grid(const ptrdiff_t n, const unsigned block) {
    const ptrdiff_t blocks = (n * SCCD_WARP_SIZE + block - 1) / block;
    const ptrdiff_t max_x = SCCD_C2D_WARP_GRID_MAX_ITEMS / block;
    if (blocks <= max_x) return dim3((unsigned)blocks);
    return dim3((unsigned)max_x, (unsigned)((blocks + max_x - 1) / max_x));
}
#else
#define SCCD_C2D_WARP_TID() ((ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x)
#endif

// The box arrays reach this file as a device table of row pointers (T** in
// device memory), and soa_device_row copies that table back to the host for
// every row it is asked for: a blocking 48-byte device-to-host copy, two to six
// per kernel launch, 26 per device broad-phase prep. On MI300A (ROCm 6.4) each
// costs tens of microseconds, and together they were most of the prep's gap to
// GH200 (cloth-funnel median 1.05 ms against 0.69 ms per case). With this set,
// each host function copies the table once and reads every row it needs from
// that copy: the same pointers, 9 copies per prep instead of 26 (MI300A,
// cloth-funnel median prep 1.20 -> 0.88 ms). A table its owner registered
// (register_soa_rows) is not copied at all, which takes the remaining 9 away
// for the smesh CCD's three box tables.
//
// The CUDA build keeps one copy per row unless this is set explicitly, so the
// H100 numbers stay those of the code that produced them. The macro itself is
// in sccd_cell2d_broadphase.cuh.

#if SCCD_DEVICE_SOA_ROWS_ONCE
#include <cstring>
#include <mutex>
#endif

// The forward walks test SCCD_WARP_SIZE cells of a row at once, one per lane, and then scan the
// survivors in ascending order -- the order the serial test visits them in, so the pairs and their
// order are unchanged. Tested one at a time, a box spanning many cells kept its warp running long
// after the rest of the launch had finished: Nsight Compute had cloth-funnel's edge-edge kernels
// active on 16-18% of the elapsed SM cycles. On a GH200 the broad-phase queries fall 26% on
// cloth-funnel, 12% on armadillo-rollers, 14% on n-body, 13% on puffer-ball and 7% on rod-twist, with
// the same candidate pairs on every case.
#ifndef SCCD_C2D_WARP_CELL_BALLOT
#define SCCD_C2D_WARP_CELL_BALLOT 1
#endif

// Grid sizing reduces all three axes in two launches instead of six; each block does the same
// arithmetic in the same order, so the grid it chooses is bit-identical. Four fewer launches per
// binning take 8% off the device prep on cloth-funnel and 12% on armadillo-rollers, the scenes made
// of many small steps (GH200), and nothing measurable on the large ones.
#ifndef SCCD_DEVICE_FUSED_GRID_STATS
#define SCCD_DEVICE_FUSED_GRID_STATS 1
#endif

namespace sccd {
    namespace device {

#if SCCD_DEVICE_SOA_ROWS_ONCE
        // Defined in sccd_broadphase.cu next to soa_device_row.
        bool is_ptr_device(const void* ptr);
#endif

#if SCCD_DEVICE_SOA_ROWS_ONCE && SCCD_DEVICE_SOA_ROWS_REGISTRY
        namespace {
            // A handful of tables at most -- three per live CCD -- so a short
            // array searched linearly under a lock.
            struct SoaRowsRegistry {
                static constexpr int kMax = 32;
                std::mutex mutex;
                int n = 0;
                const void* table[kMax];
                void* rows[kMax][2 * SCCD_DIM];
            };

            SoaRowsRegistry& soa_rows_registry() {
                static SoaRowsRegistry* r = new SoaRowsRegistry();  // never destroyed, like the runtime it outlives
                return *r;
            }

            bool lookup_soa_rows(const void* device_table, void** out) {
                SoaRowsRegistry& r = soa_rows_registry();
                std::lock_guard<std::mutex> lock(r.mutex);
                for (int i = 0; i < r.n; ++i) {
                    if (r.table[i] == device_table) {
                        std::memcpy(out, r.rows[i], sizeof(r.rows[i]));
                        return true;
                    }
                }
                return false;
            }
        }  // namespace

        void register_soa_rows(const void* device_table, void* const* host_rows) {
            SoaRowsRegistry& r = soa_rows_registry();
            std::lock_guard<std::mutex> lock(r.mutex);
            int i = 0;
            while (i < r.n && r.table[i] != device_table) ++i;
            if (i == r.n) {
                if (r.n == SoaRowsRegistry::kMax) return;  // full: the table is simply copied back
                ++r.n;
            }
            r.table[i] = device_table;
            std::memcpy(r.rows[i], host_rows, sizeof(r.rows[i]));
        }

        void unregister_soa_rows(const void* device_table) {
            SoaRowsRegistry& r = soa_rows_registry();
            std::lock_guard<std::mutex> lock(r.mutex);
            for (int i = 0; i < r.n; ++i) {
                if (r.table[i] == device_table) {
                    --r.n;
                    r.table[i] = r.table[r.n];
                    std::memcpy(r.rows[i], r.rows[r.n], sizeof(r.rows[i]));
                    return;
                }
            }
        }
#endif

        namespace detail {

            // The rows of one SoA box table, read the way SCCD_DEVICE_SOA_ROWS_ONCE
            // says: the whole table on the first access, or one row per access
            // through soa_device_row as before. It lives for one host function,
            // so it cannot outlast the table it copied.
            template <typename T>
            class SoaRows {
            public:
                explicit SoaRows(T** const arrays) : arrays_(arrays) {}

                T* operator[](const int row) {
#if SCCD_DEVICE_SOA_ROWS_ONCE
                    if (!fetched_) {
#if SCCD_DEVICE_SOA_ROWS_REGISTRY
                        if (lookup_soa_rows(arrays_, reinterpret_cast<void**>(rows_))) {
                            // registered by its owner: nothing to copy
                        } else
#endif
                        {
#if defined(__HIP_PLATFORM_SPIRV__)
                            // chipStar cannot say where a table it did not allocate lives, so it is
                            // read through the runtime's copy, as soa_rows_to_host does (sccd_broadphase.cu).
                            SCCD_CHECK_CUDA(hipMemcpy(rows_, arrays_, sizeof(rows_), hipMemcpyDefault));
#else
                            if (is_ptr_device(arrays_)) {
                                SCCD_CHECK_CUDA(cudaMemcpy(rows_, arrays_, sizeof(rows_), cudaMemcpyDeviceToHost));
                            } else {
                                for (int d = 0; d < 2 * SCCD_DIM; ++d) rows_[d] = arrays_[d];
                            }
#endif
                        }
                        fetched_ = true;
                    }
                    return rows_[row];
#else
                    return soa_device_row<T>(arrays_, row);
#endif
                }

            private:
                T** arrays_;
#if SCCD_DEVICE_SOA_ROWS_ONCE
                T* rows_[2 * SCCD_DIM] = {};
                bool fetched_ = false;
#endif
            };

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

            /** \brief Threads a block brings to the grid-sizing reduction. */
            static constexpr int kGridStatsBlock = 256;
            /** \brief Blocks it splits the array over, and so partials to combine. */
            static constexpr int kGridStatsBlocks = 512;

            /**
             * \brief Reduce a block's stride of the array to one (min, max, sum).
             *
             * Per-axis min, max and mean extent size the grid. A shared-array tree
             * rather than atomics, because there is no floating-point atomic
             * minimum or maximum here; over many blocks rather than one, because
             * one block is 256 threads of a device that has hundreds of times
             * that, and the reduction is then linear in the element count with
             * none of the machine working. Profiled on puffer-ball it was 10.1%
             * of all GPU time, and the largest single term in the difference
             * between binning the faces and binning the vertices -- not because
             * either is wrong, but because the face list is the longer one and
             * this kernel charges by length.
             */
            template <typename T>
            static __device__ __forceinline__ void grid_stats_partial_block(const ptrdiff_t n,
                                                                            const T* const SCCD_RESTRICT lo,
                                                                            const T* const SCCD_RESTRICT hi,
                                                                            T* const SCCD_RESTRICT part_min,
                                                                            T* const SCCD_RESTRICT part_max,
                                                                            T* const SCCD_RESTRICT part_sum) {
                extern __shared__ char s_raw[];
                T* const s_min = (T*)s_raw;
                T* const s_max = s_min + blockDim.x;
                T* const s_sum = s_max + blockDim.x;

                // lo[0] and hi[0] are the identity for a lane with no element:
                // they are a real box's bounds, so folding them in cannot widen
                // the result past what the array holds.
                T tmin = lo[0], tmax = hi[0], tsum = T(0);
                const ptrdiff_t stride = (ptrdiff_t)blockDim.x * (ptrdiff_t)gridDim.x;
                for (ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x; i < n; i += stride) {
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
                    part_min[blockIdx.x] = s_min[0];
                    part_max[blockIdx.x] = s_max[0];
                    part_sum[blockIdx.x] = s_sum[0];
                }
            }

            template <typename T>
            __global__ void grid_stats_partial_kernel(const ptrdiff_t n,
                                                      const T* const SCCD_RESTRICT lo,
                                                      const T* const SCCD_RESTRICT hi,
                                                      T* const SCCD_RESTRICT part_min,
                                                      T* const SCCD_RESTRICT part_max,
                                                      T* const SCCD_RESTRICT part_sum) {
                grid_stats_partial_block<T>(n, lo, hi, part_min, part_max, part_sum);
            }

            /** \brief The three axes' rows and partial arrays, for the fused launches. */
            template <typename T>
            struct GridStatsAxes {
                const T* lo[3];
                const T* hi[3];
                T* parts;  // 9 arrays of nslots: (min, max, sum) per axis
                T* out;    // 9 values: (min, max, sum) per axis
                int nslots;
            };

            /** \brief All three axes' partial reductions in one launch, one axis per
             *         grid row (blockIdx.y); each block does what the per-axis kernel does. */
            template <typename T>
            __global__ void grid_stats_partial3_kernel(const ptrdiff_t n, const GridStatsAxes<T> ax) {
                const int d = blockIdx.y;
                T* const p = ax.parts + (ptrdiff_t)(3 * d) * ax.nslots;
                grid_stats_partial_block<T>(n, ax.lo[d], ax.hi[d], p, p + ax.nslots, p + 2 * (ptrdiff_t)ax.nslots);
            }

            /**
             * \brief Combine the per-block partials into the axis's three values.
             *
             * One block, which is right here: there are at most
             * \ref kGridStatsBlocks of them however long the array was, so this
             * stage does not grow with the input.
             *
             * The sums arrive already summed, so this one adds where the first
             * pass subtracted -- the two stages are not the same reduction and
             * cannot share a kernel.
             */
            template <typename T>
            static __device__ __forceinline__ void grid_stats_final_block(const int nparts,
                                                                          const T* const SCCD_RESTRICT part_min,
                                                                          const T* const SCCD_RESTRICT part_max,
                                                                          const T* const SCCD_RESTRICT part_sum,
                                                                          T* const SCCD_RESTRICT out_min,
                                                                          T* const SCCD_RESTRICT out_max,
                                                                          T* const SCCD_RESTRICT out_sum) {
                extern __shared__ char s_raw[];
                T* const s_min = (T*)s_raw;
                T* const s_max = s_min + blockDim.x;
                T* const s_sum = s_max + blockDim.x;

                T tmin = part_min[0], tmax = part_max[0], tsum = T(0);
                for (int i = threadIdx.x; i < nparts; i += blockDim.x) {
                    tmin = part_min[i] < tmin ? part_min[i] : tmin;
                    tmax = part_max[i] > tmax ? part_max[i] : tmax;
                    tsum += part_sum[i];
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
            __global__ void grid_stats_final_kernel(const int nparts,
                                                    const T* const SCCD_RESTRICT part_min,
                                                    const T* const SCCD_RESTRICT part_max,
                                                    const T* const SCCD_RESTRICT part_sum,
                                                    T* const SCCD_RESTRICT out_min,
                                                    T* const SCCD_RESTRICT out_max,
                                                    T* const SCCD_RESTRICT out_sum) {
                grid_stats_final_block<T>(nparts, part_min, part_max, part_sum, out_min, out_max, out_sum);
            }

            /** \brief All three axes' final combines in one launch, one block per axis. */
            template <typename T>
            __global__ void grid_stats_final3_kernel(const int nparts, const GridStatsAxes<T> ax) {
                const int d = blockIdx.y;
                const T* const p = ax.parts + (ptrdiff_t)(3 * d) * ax.nslots;
                grid_stats_final_block<T>(nparts, p, p + ax.nslots, p + 2 * (ptrdiff_t)ax.nslots,
                                          ax.out + 3 * d, ax.out + 3 * d + 1, ax.out + 3 * d + 2);
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
             * and that looks like the reason this walk gains less on the device
             * than on the host. It is not: replacing this with a segmented scan
             * by key over every cell, which makes the parallelism the cell count,
             * moved the four measured scenes by 1.531x to 1.523x, 1.047x to
             * 0.963x, 0.665x to 0.673x and 1.638x to 1.710x -- noise in both
             * directions. The simple kernel is kept and the scan is not worth its
             * two cub iterators.
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

            /** \brief Index of the lowest set lane of \p m (which must be non-zero). */
            static __device__ __forceinline__ int lowest_lane(const device::lane_mask_t m) {
                if constexpr (sizeof(device::lane_mask_t) == 8) {
                    return __ffsll((long long)m) - 1;
                } else {
                    return __ffs((int)m) - 1;
                }
            }

            /**
             * \brief Visit the cells \p c0 .. \p c0e of row \p row whose bound reaches
             *        \p amin1, in ascending order, testing a warp's width of them at once.
             */
            template <typename T, typename Visit>
            static __device__ __forceinline__ void warp_row_cells(const ptrdiff_t row,
                                                                  int c0,
                                                                  const int c0e,
                                                                  const T amin1,
                                                                  const T* const SCCD_RESTRICT cell_hi1,
                                                                  Visit&& visit) {
                const int lane = (int)(threadIdx.x & (SCCD_WARP_SIZE - 1));
                for (; c0 <= c0e; c0 += SCCD_WARP_SIZE) {
                    const int mine = c0 + lane;
                    const bool ok = mine <= c0e && !(cell_hi1[row + mine] < amin1);
                    device::lane_mask_t m = device::warp_ballot(ok);
                    while (m) {
                        const int b = lowest_lane(m);
                        m &= m - 1;
                        visit(row + c0 + b);
                    }
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
             * \brief Walk forward over the cells, one warp to a box.
             *
             * The thread-per-box form of this walk diverges on everything it
             * does: the number of rows a box spans, the columns its binary
             * search leaves, and the number of entries in each cell all differ
             * between neighbouring boxes, so the lanes of a warp spend most of
             * their time masked off, and each lane reads a different part of the
             * cell array. A warp to a box inverts that. Every lane computes the
             * same walk -- the same rows, the same searches, redundantly and
             * branch-free -- and the lanes split each cell's entries between
             * them, which is the only loop left that can diverge and is also the
             * one whose loads are then consecutive.
             *
             * \p combine is handed one candidate per lane per step, with the
             * lanes that found nothing marked. It is a warp-wide operation
             * because the count and the write offset are: every lane iterates
             * the same number of times so that the ballot inside it is taken
             * with the whole warp converged.
             */
            template <int NXE, typename T, typename I, typename Combine>
            static __device__ __forceinline__ void warp_forward_self_partner(
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
                Combine&& combine) {
                const unsigned lane = threadIdx.x & (SCCD_WARP_SIZE - 1);

                const T aminx = aabbs[0][fi], aminy = aabbs[1][fi], aminz = aabbs[2][fi];
                const T amaxx = aabbs[3][fi], amaxy = aabbs[4][fi], amaxz = aabbs[5][fi];

                const T amin0 = aabbs[grid.axis0][fi];
                const T amin1 = aabbs[grid.axis1][fi];

                const int c0b = cell0<T>(grid, amin0);
                const int c0e = cell0<T>(grid, aabbs[SCCD_DIM + grid.axis0][fi]);
                const int c1b = cell1<T>(grid, amin1);
                const int c1e = cell1<T>(grid, aabbs[SCCD_DIM + grid.axis1][fi]);
                const ptrdiff_t own = cell_of<T>(grid, c0b, c1b);

                const auto scan = [&](const ptrdiff_t cell, const bool is_own) {
                    const ptrdiff_t begin = cellptr[cell];
                    const ptrdiff_t end = cellptr[cell + 1];
                    // Every lane runs the same trip count, so the warp is
                    // converged where combine takes its ballot.
                    for (ptrdiff_t base = begin; base < end; base += SCCD_WARP_SIZE) {
                        const ptrdiff_t k = base + (ptrdiff_t)lane;
                        bool keep = false;
                        ptrdiff_t j = 0;
                        I jidx = 0;
                        if (k < end) {
                            j = (ptrdiff_t)cellidx[k];
                            if (!(is_own && j <= fi) &&
                                !disjoint<T>(aminx, aminy, aminz, amaxx, amaxy, amaxz,
                                             aabbs[0][j], aabbs[1][j], aabbs[2][j],
                                             aabbs[3][j], aabbs[4][j], aabbs[5][j])) {
                                jidx = idx[j];
                                I sev[NXE];
                                load_ev<NXE, I>(elements, jidx, element_stride, sev);
                                keep = !shares_vertex<NXE, NXE>(ev, sev);
                            }
                        }
                        combine(keep, jidx);
                    }
                };

                // Its own cell first, where the index decides which of two boxes
                // sharing it reports the pair. Both bounds hold this box's own
                // maximum, so neither can rule the cell out.
                scan(own, true);

                for (int c1 = c1b; c1 <= c1e; ++c1) {
                    const ptrdiff_t row = (ptrdiff_t)c1 * grid.n0;
                    int c0;
                    if (c1 == c1b) {
                        c0 = c0b + 1;
                    } else {
                        c0 = lower_bound_row<T>(row_prefix + row, c0e + 1, amin0);
                    }
#if SCCD_C2D_WARP_CELL_BALLOT
                    warp_row_cells<T>(row, c0, c0e, amin1, cell_hi1,
                                      [&](const ptrdiff_t cell) { scan(cell, false); });
#else
                    for (; c0 <= c0e; ++c0) {
                        const ptrdiff_t cell = row + c0;
                        if (cell_hi1[cell] < amin1) continue;
                        scan(cell, false);
                    }
#endif
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
                const ptrdiff_t tid = SCCD_C2D_WARP_TID();
                const ptrdiff_t fi = tid / SCCD_WARP_SIZE;
                if (tid == 0) ccdptr[0] = 0;
                // Uniform over the warp, so the whole warp leaves together and
                // the ballots below are taken converged.
                if (fi >= element_count) return;

                I ev[nxe];
                load_ev<nxe, I>(elements, idx[fi], element_stride, ev);

                ptrdiff_t count = 0;
                warp_forward_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx,
                    row_prefix, cell_hi1, [&](const bool keep, const I) {
                        count += device::warp_popc(device::warp_ballot(keep));
                    });
                // Every lane accumulated the same total, so one writes it.
                if ((threadIdx.x & (SCCD_WARP_SIZE - 1)) == 0) ccdptr[fi + 1] = count;
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
                const ptrdiff_t tid = SCCD_C2D_WARP_TID();
                const ptrdiff_t fi = tid / SCCD_WARP_SIZE;
                if (fi >= element_count) return;

                const unsigned lane = threadIdx.x & (SCCD_WARP_SIZE - 1);
                const device::lane_mask_t below = device::lanes_below(lane);

                const I idxi = idx[fi];
                I ev[nxe];
                load_ev<nxe, I>(elements, idxi, element_stride, ev);

                // `at` stays the same on every lane: each step advances it by the
                // whole warp's count, and a lane writes at its rank within that.
                ptrdiff_t at = ccdptr[fi];
                warp_forward_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx,
                    row_prefix, cell_hi1, [&](const bool keep, const I jidx) {
                        const device::lane_mask_t mask = device::warp_ballot(keep);
                        if (keep) {
                            const ptrdiff_t w = at + (ptrdiff_t)device::warp_popc(mask & below);
                            first_out[w] = idxi < jidx ? idxi : jidx;
                            second_out[w] = idxi < jidx ? jidx : idxi;
                        }
                        at += (ptrdiff_t)device::warp_popc(mask);
                    });
            }

            /**
             * \brief The most rows any box of the binned list spans.
             *
             * The second axis has no prefix array to search, so the cross-list
             * walk bounds its row range with this instead: a box reaching the
             * querying box on that axis has its maximum at or above the querying
             * box's minimum row, and spans at most this many rows, so it is
             * binned no lower than that row less this.
             */
            template <typename T>
            __global__ void min_row_span_kernel(const ptrdiff_t n,
                                                const T* const SCCD_RESTRICT lo,
                                                const T* const SCCD_RESTRICT hi,
                                                const Cell2DGridD<T> grid,
                                                int* const SCCD_RESTRICT out) {
                const ptrdiff_t i = (ptrdiff_t)blockIdx.x * blockDim.x + threadIdx.x;
                int span = 0;
                if (i < n) {
                    span = cell1<T>(grid, hi[i]) - cell1<T>(grid, lo[i]);
                }
                // One atomic per warp rather than per thread: the maximum is a
                // single number and the whole grid contends for it.
                for (int off = SCCD_WARP_SIZE / 2; off > 0; off >>= 1) {
                    const int other = device::warp_shfl_down(span, off);
                    span = other > span ? other : span;
                }
                if ((threadIdx.x & (SCCD_WARP_SIZE - 1)) == 0 && span > 0) {
                    atomicMax(out, span);
                }
            }

            /**
             * \brief Walk a minimum-corner binning of the other list, one warp
             *        to a querying box.
             *
             * The device counterpart of \ref sccd::detail::for_each_min_corner_partner,
             * and it inherits the warp shape of \ref warp_forward_self_partner for
             * the same reason: the rows a box spans, the columns its search
             * leaves and the entries in each cell all differ between neighbouring
             * boxes. Every lane computes the same walk and the lanes split each
             * cell's entries between them.
             *
             * Two things the self walk needs are absent here. There is no own
             * cell to visit first and no index test, because the two lists are
             * distinct and a partner has exactly one cell, so a pair is met
             * exactly once. In exchange the row range is bounded from below by
             * \p krow instead of starting at this box's own row, and every row
             * pays a search -- the self walk skips the search on its own row.
             */
            template <int S, typename T, typename I, typename Combine>
            static __device__ __forceinline__ void warp_min_corner_partner(
                T** const SCCD_RESTRICT first_aabbs,
                const ptrdiff_t fi,
                const I vidx,
                T** const SCCD_RESTRICT second_aabbs,
                const I* const SCCD_RESTRICT second_idx,
                I** const SCCD_RESTRICT second_elements,
                const ptrdiff_t second_element_stride,
                const Cell2DGridD<T>& grid,
                const ptrdiff_t* const SCCD_RESTRICT cellptr,
                const I* const SCCD_RESTRICT cellidx,
                const T* const SCCD_RESTRICT row_prefix,
                const T* const SCCD_RESTRICT cell_hi1,
                const int* const SCCD_RESTRICT krow,
                Combine&& combine) {
                const unsigned lane = threadIdx.x & (SCCD_WARP_SIZE - 1);

                const T aminx = first_aabbs[0][fi];
                const T aminy = first_aabbs[1][fi];
                const T aminz = first_aabbs[2][fi];
                const T amaxx = first_aabbs[SCCD_DIM + 0][fi];
                const T amaxy = first_aabbs[SCCD_DIM + 1][fi];
                const T amaxz = first_aabbs[SCCD_DIM + 2][fi];

                const T amin0 = first_aabbs[grid.axis0][fi];
                const T amin1 = first_aabbs[grid.axis1][fi];

                // Above and to the right needs nothing stored: a partner
                // overlapping this box has its minimum corner at or below this
                // box's maximum on both axes.
                const int c0e = cell0<T>(grid, first_aabbs[SCCD_DIM + grid.axis0][fi]);
                const int c1e = cell1<T>(grid, first_aabbs[SCCD_DIM + grid.axis1][fi]);
                // Read from the device rather than taken as a launch argument:
                // it is computed by a kernel, and fetching it to pass it back
                // would put a full synchronisation in every step.
                int c1b = cell1<T>(grid, amin1) - *krow;
                if (c1b < 0) c1b = 0;

                for (int c1 = c1b; c1 <= c1e; ++c1) {
                    const ptrdiff_t row = (ptrdiff_t)c1 * grid.n0;
                    int c0 = lower_bound_row<T>(row_prefix + row, c0e + 1, amin0);

                    const auto scan = [&](const ptrdiff_t cell) {
                        const ptrdiff_t begin = cellptr[cell];
                        const ptrdiff_t end = cellptr[cell + 1];
                        // Every lane runs the same trip count, so the warp is
                        // converged where combine takes its ballot.
                        for (ptrdiff_t base = begin; base < end; base += SCCD_WARP_SIZE) {
                            const ptrdiff_t k = base + (ptrdiff_t)lane;
                            bool keep = false;
                            I jidx = 0;
                            if (k < end) {
                                const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                                if (!disjoint<T>(aminx, aminy, aminz, amaxx, amaxy, amaxz,
                                                 second_aabbs[0][j],
                                                 second_aabbs[1][j],
                                                 second_aabbs[2][j],
                                                 second_aabbs[SCCD_DIM + 0][j],
                                                 second_aabbs[SCCD_DIM + 1][j],
                                                 second_aabbs[SCCD_DIM + 2][j])) {
                                    jidx = second_idx[j];
                                    keep = true;
                                    for (int v = 0; v < S; ++v) {
                                        if (second_elements[v][(ptrdiff_t)jidx * second_element_stride] ==
                                            vidx) {
                                            keep = false;
                                            break;
                                        }
                                    }
                                }
                            }
                            combine(keep, jidx);
                        }
                    };
#if SCCD_C2D_WARP_CELL_BALLOT
                    warp_row_cells<T>(row, c0, c0e, amin1, cell_hi1, scan);
#else
                    for (; c0 <= c0e; ++c0) {
                        const ptrdiff_t cell = row + c0;
                        if (cell_hi1[cell] < amin1) continue;
                        scan(cell);
                    }
#endif
                }
            }

            template <int S, typename T, typename I>
            __global__ void min_count_vf_kernel(const ptrdiff_t vertex_count,
                                                T** const SCCD_RESTRICT vaabbs,
                                                const I* const SCCD_RESTRICT vertex_idx,
                                                T** const SCCD_RESTRICT faabbs,
                                                const I* const SCCD_RESTRICT face_idx,
                                                const ptrdiff_t face_element_stride,
                                                I** const SCCD_RESTRICT face_elements,
                                                const Cell2DGridD<T> grid,
                                                const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                const I* const SCCD_RESTRICT cellidx,
                                                const T* const SCCD_RESTRICT row_prefix,
                                                const T* const SCCD_RESTRICT cell_hi1,
                                                const int* const SCCD_RESTRICT krow,
                                                ptrdiff_t* const SCCD_RESTRICT ccdptr) {
                const ptrdiff_t tid = SCCD_C2D_WARP_TID();
                const ptrdiff_t vi = tid / SCCD_WARP_SIZE;
                if (tid == 0) ccdptr[0] = 0;
                // Uniform over the warp, so the whole warp leaves together and
                // the ballots below are taken converged.
                if (vi >= vertex_count) return;

                ptrdiff_t count = 0;
                warp_min_corner_partner<S, T, I>(
                    vaabbs, vi, vertex_idx[vi], faabbs, face_idx, face_elements,
                    face_element_stride, grid, cellptr, cellidx, row_prefix, cell_hi1, krow,
                    [&](const bool keep, const I) {
                        count += device::warp_popc(device::warp_ballot(keep));
                    });
                // Every lane accumulated the same total, so one writes it.
                if ((threadIdx.x & (SCCD_WARP_SIZE - 1)) == 0) ccdptr[vi + 1] = count;
            }

            template <int S, typename T, typename I>
            __global__ void min_collect_vf_kernel(const ptrdiff_t vertex_count,
                                                  T** const SCCD_RESTRICT vaabbs,
                                                  const I* const SCCD_RESTRICT vertex_idx,
                                                  T** const SCCD_RESTRICT faabbs,
                                                  const I* const SCCD_RESTRICT face_idx,
                                                  const ptrdiff_t face_element_stride,
                                                  I** const SCCD_RESTRICT face_elements,
                                                  const Cell2DGridD<T> grid,
                                                  const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                  const I* const SCCD_RESTRICT cellidx,
                                                  const T* const SCCD_RESTRICT row_prefix,
                                                  const T* const SCCD_RESTRICT cell_hi1,
                                                  const int* const SCCD_RESTRICT krow,
                                                  const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                  I* const SCCD_RESTRICT face_out,
                                                  I* const SCCD_RESTRICT vertex_out) {
                const ptrdiff_t tid = SCCD_C2D_WARP_TID();
                const ptrdiff_t vi = tid / SCCD_WARP_SIZE;
                if (vi >= vertex_count) return;

                const unsigned lane = threadIdx.x & (SCCD_WARP_SIZE - 1);
                const device::lane_mask_t below = device::lanes_below(lane);

                const I vidx = vertex_idx[vi];

                // `at` stays the same on every lane: each step advances it by the
                // whole warp's count, and a lane writes at its rank within that.
                ptrdiff_t at = ccdptr[vi];
                warp_min_corner_partner<S, T, I>(
                    vaabbs, vi, vidx, faabbs, face_idx, face_elements, face_element_stride,
                    grid, cellptr, cellidx, row_prefix, cell_hi1, krow,
                    [&](const bool keep, const I jidx) {
                        const device::lane_mask_t mask = device::warp_ballot(keep);
                        if (keep) {
                            const ptrdiff_t w = at + (ptrdiff_t)device::warp_popc(mask & below);
                            face_out[w] = jidx;
                            vertex_out[w] = vidx;
                        }
                        at += (ptrdiff_t)device::warp_popc(mask);
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
            detail::SoaRows<T> rows(aabbs);
            // Axis choice is on the host: it needs three reductions and then a
            // decision, and the arrays are small enough that a kernel per axis is
            // cheaper than the launch overhead of doing it any other way.
            // Nine results and three partial arrays, in one allocation so the
            // three axes can be in flight together as they were before.
            constexpr int kB = detail::kGridStatsBlock;
            constexpr int kG = detail::kGridStatsBlocks;
            // At least one, since both callers return before this on an empty
            // list, and at most kG so the partial arrays are a fixed size.
            const ptrdiff_t want = (n + kB - 1) / kB;
            const int nparts = (int)(want < (ptrdiff_t)kG ? want : (ptrdiff_t)kG);
            T* const stats = workspace(WorkspaceSlot::Scratch).get_as<T>(9 + 9 * (ptrdiff_t)kG);
            T* const parts = stats + 9;
            T h_min[3], h_max[3], h_sum[3];
            const size_t shmem = 3 * kB * sizeof(T);
#if SCCD_DEVICE_FUSED_GRID_STATS
            {
                detail::GridStatsAxes<T> ax;
                for (int d = 0; d < 3; ++d) {
                    ax.lo[d] = rows[d];
                    ax.hi[d] = rows[SCCD_DIM + d];
                }
                ax.parts = parts;
                ax.out = stats;
                ax.nslots = kG;
                detail::grid_stats_partial3_kernel<T><<<dim3(nparts, 3), kB, shmem>>>(n, ax);
                detail::grid_stats_final3_kernel<T><<<dim3(1, 3), kB, shmem>>>(nparts, ax);
            }
#else
            for (int d = 0; d < 3; ++d) {
                // Each axis gets its own third of the partial space, so the
                // three reductions stay independent.
                T* const pmin = parts + (ptrdiff_t)(3 * d + 0) * kG;
                T* const pmax = parts + (ptrdiff_t)(3 * d + 1) * kG;
                T* const psum = parts + (ptrdiff_t)(3 * d + 2) * kG;
                detail::grid_stats_partial_kernel<T><<<nparts, kB, shmem>>>(
                    n,
                    rows[d],
                    rows[SCCD_DIM + d],
                    pmin, pmax, psum);
                detail::grid_stats_final_kernel<T><<<1, kB, shmem>>>(
                    nparts, pmin, pmax, psum, stats + 3 * d, stats + 3 * d + 1, stats + 3 * d + 2);
            }
#endif
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
            detail::SoaRows<T> rows(aabbs);
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
                                                        rows[grid.axis0],
                                                        rows[SCCD_DIM + grid.axis0],
                                                        rows[grid.axis1],
                                                        rows[SCCD_DIM + grid.axis1],
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
            detail::SoaRows<T> rows(aabbs);
            SCCD_CUDA_LAST_ERROR();
            if (n <= 0) return;

            size_grid_device<T>(n, aabbs, grid);

            const ptrdiff_t ncells = grid.ncells();
            SCCD_CHECK_CUDA(cudaMemset(cellptr, 0, sizeof(ptrdiff_t) * (size_t)(ncells + 1)));

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::min_bin_count_kernel<T><<<gridsz, block>>>(
                n, rows[grid.axis0], rows[grid.axis1],
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
            detail::SoaRows<T> rows(aabbs);
            if (n <= 0) return;
            SCCD_CHECK_CUDA(cudaMemset(cursor, 0, sizeof(ptrdiff_t) * (size_t)grid.ncells()));

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::min_bin_fill_kernel<T, I><<<gridsz, block>>>(
                n, rows[grid.axis0], rows[grid.axis1],
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
            detail::SoaRows<T> rows(aabbs);
            const ptrdiff_t ncells = grid.ncells();
            if (ncells <= 0) return;

            // The per-cell maximum on the first axis is written straight into
            // row_prefix, which the next kernel turns into the running maximum in
            // place: the two are the same array at two stages, as on the host.
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((ncells + block.x - 1) / block.x);
            detail::min_cell_bounds_kernel<T, I><<<gridsz, block>>>(
                ncells, cellptr, cellidx,
                rows[SCCD_DIM + grid.axis0],
                rows[SCCD_DIM + grid.axis1],
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
            // One warp to a box, so the launch is sized in warps.
#if SCCD_C2D_WARP_GRID_2D
            const dim3 gridsz = sccd_c2d_warp_grid(element_count, block.x);
#else
            dim3 gridsz((element_count * SCCD_WARP_SIZE + block.x - 1) / block.x);
#endif
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
#if SCCD_C2D_WARP_GRID_2D
            const dim3 gridsz = sccd_c2d_warp_grid(element_count, block.x);
#else
            dim3 gridsz((element_count * SCCD_WARP_SIZE + block.x - 1) / block.x);
#endif
            detail::min_collect_self_kernel<nxe, T, I><<<gridsz, block>>>(
                element_count, aabbs, idx, element_stride, elements, grid, cellptr, cellidx,
                row_prefix, cell_hi1, ccdptr, first_out, second_out);
            SCCD_CUDA_LAST_ERROR();
        }

        template <typename T>
        void cell2dmin_max_row_span(const ptrdiff_t n,
                                    T** const SCCD_RESTRICT aabbs,
                                    const Cell2DGridD<T>& grid,
                                    int* const SCCD_RESTRICT krow) {
            detail::SoaRows<T> rows(aabbs);
            SCCD_CHECK_CUDA(cudaMemset(krow, 0, sizeof(int)));
            if (n <= 0) return;

            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 gridsz((n + block.x - 1) / block.x);
            detail::min_row_span_kernel<T><<<gridsz, block>>>(
                n,
                rows[grid.axis1],
                rows[SCCD_DIM + grid.axis1],
                grid,
                krow);
            SCCD_CUDA_LAST_ERROR();
        }

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
                                         ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            if (vertex_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            // One warp to a vertex, so the launch is sized in warps.
#if SCCD_C2D_WARP_GRID_2D
            const dim3 gridsz = sccd_c2d_warp_grid(vertex_count, block.x);
#else
            dim3 gridsz((vertex_count * SCCD_WARP_SIZE + block.x - 1) / block.x);
#endif
            detail::min_count_vf_kernel<S, T, I><<<gridsz, block>>>(
                vertex_count, vaabbs, vertex_idx, faabbs, face_idx, face_element_stride,
                face_elements, grid, cellptr, cellidx, row_prefix, cell_hi1, krow, ccdptr);
            SCCD_CUDA_LAST_ERROR();
            detail::inclusive_sum<T>(ccdptr, vertex_count + 1);
        }

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
                                           I* const SCCD_RESTRICT vertex_out) {
            if (vertex_count <= 0) return;
            dim3 block(SCCD_C2D_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
#if SCCD_C2D_WARP_GRID_2D
            const dim3 gridsz = sccd_c2d_warp_grid(vertex_count, block.x);
#else
            dim3 gridsz((vertex_count * SCCD_WARP_SIZE + block.x - 1) / block.x);
#endif
            detail::min_collect_vf_kernel<S, T, I><<<gridsz, block>>>(
                vertex_count, vaabbs, vertex_idx, faabbs, face_idx, face_element_stride,
                face_elements, grid, cellptr, cellidx, row_prefix, cell_hi1, krow, ccdptr,
                face_out, vertex_out);
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
        I*);                                                                                   \
    template void sccd::device::cell2dmin_count_vf_overlaps<NXE, T, I>(                        \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        const I*,                                                                              \
        T**,                                                                                   \
        const I*,                                                                              \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        const T*,                                                                              \
        const T*,                                                                              \
        const int*,                                                                            \
        ptrdiff_t*);                                                                           \
    template void sccd::device::cell2dmin_collect_vf_overlaps<NXE, T, I>(                      \
        const ptrdiff_t,                                                                       \
        T**,                                                                                   \
        const I*,                                                                              \
        T**,                                                                                   \
        const I*,                                                                              \
        const ptrdiff_t,                                                                       \
        I**,                                                                                   \
        const sccd::device::Cell2DGridD<T>&,                                                   \
        const ptrdiff_t*,                                                                      \
        const I*,                                                                              \
        const T*,                                                                              \
        const T*,                                                                              \
        const int*,                                                                            \
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
    template void sccd::device::cell2dmin_max_row_span<T>(                                     \
        const ptrdiff_t, T**, const sccd::device::Cell2DGridD<T>&, int*);                      \
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
