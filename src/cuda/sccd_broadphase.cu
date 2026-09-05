#include "sccd_broadphase.cuh"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_scan.cuh>

#include "sccd_base.hpp"
#include "sccd_cuda_base.cuh"

#include "sccd_device_workspace.cuh"
#include "sccd_reduce.cuh"

#define SCCD_BP_N_WARPS_PER_BLOCK 8

namespace sccd {
    namespace device {

        template <typename T>
        __device__ __forceinline__ const T* lower_bound(const T* const SCCD_RESTRICT begin,
                                                        const T* const SCCD_RESTRICT end,
                                                        const T val) {
            ptrdiff_t count = end - begin;
            const T* it = begin;

            while (count > 0) {
                const ptrdiff_t step = count >> 1;
                const T* const mid = it + step;

                if (*mid < val) {
                    it = mid + 1;
                    count -= step + 1;
                } else {
                    count = step;
                }
            }

            return it;
        }

        template <typename T>
        __global__ void choose_axis_mean_kernel(                                                const ptrdiff_t n,
                                                const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                                T* const SCCD_RESTRICT mean) {
            ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= n) return;

            T local_mean[3] = {0};

            for (int d = 0; d < SCCD_DIM; d++) {
                const T p0 = aabbs[d][i];
                const T p1 = aabbs[SCCD_DIM + d][i];
                const T p = (p0 + p1) / 2;
                local_mean[d] += p;
            }

            __shared__ T block_accumulator[SCCD_BP_N_WARPS_PER_BLOCK];
            for (int d = 0; d < SCCD_DIM; d++) {
                block_reduce_to_gmem<T>(local_mean[d], block_accumulator, &mean[d]);
            }
        }

        template <typename T>
        __global__ void choose_axis_var_kernel(                                               const ptrdiff_t n,
                                               const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                               T* const SCCD_RESTRICT mean,
                                               T* const SCCD_RESTRICT var) {
            ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= n) return;

            T local_var[3] = {0};
            for (int d = 0; d < SCCD_DIM; d++) {
                const T m = mean[d] / n;
                const T p = (aabbs[d][i] + aabbs[SCCD_DIM + d][i]) / 2;
                local_var[d] += (p - m) * (p - m);
            }

            __shared__ T block_accumulator[SCCD_BP_N_WARPS_PER_BLOCK];
            for (int d = 0; d < SCCD_DIM; d++) {
                block_reduce_to_gmem<T>(local_var[d], block_accumulator, &var[d]);
            }
        }

        template <typename T>
        int choose_axis(const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs) {
            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((n + block.x - 1) / block.x);

            // One cached allocation holding [mean | var]; these used to be two
            // cudaMalloc/cudaFree pairs on every call.
            T* const mean = workspace(WorkspaceSlot::TempStorage).get_as<T>(2 * SCCD_DIM);
            T* const var = mean + SCCD_DIM;

            cudaMemset(mean, 0, 2 * SCCD_DIM * sizeof(T));

            choose_axis_mean_kernel<T><<<grid, block>>>(n, aabbs, mean);
            choose_axis_var_kernel<T><<<grid, block>>>(n, aabbs, mean, var);

            cudaError_t error = cudaGetLastError();

            T* hvar = (T*)malloc(SCCD_DIM * sizeof(T));
            cudaMemcpy(hvar, var, SCCD_DIM * sizeof(T), cudaMemcpyDeviceToHost);

            int fargmax = 0;
            T fmax = hvar[0];

            for (int d = 1; d < SCCD_DIM; d++) {
                if (fmax < hvar[d]) {
                    fmax = hvar[d];
                    fargmax = d;
                }
            }

            free(hvar);

            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                exit(1);
            }

            return fargmax;
        }

        template <typename T>
        __global__ void enumerate_kernel(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx) {
            ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= end - begin) return;
            idx[i] = begin + i;
        }

        template <typename T>
        void enumerate(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx) {
            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((end - begin + block.x - 1) / block.x);
            enumerate_kernel<T><<<grid, block>>>(begin, end, idx);
        }

        template <typename T, typename I>
        __global__ void permute_kernel(const ptrdiff_t n,
                                       const I* const SCCD_RESTRICT idx,
                                       const T* const SCCD_RESTRICT src,
                                       T* const SCCD_RESTRICT dst) {
            const ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= n) return;
            dst[i] = src[idx[i]];
        }

        bool is_ptr_device(const void* ptr) {
            cudaPointerAttributes attributes;
            cudaError_t err = cudaPointerGetAttributes(&attributes, ptr);

            if (err != cudaSuccess) {
                fprintf(stderr, "cudaPointerGetAttributes failed: %s\n", cudaGetErrorString(err));
            }

#if CUDART_VERSION >= 10000
            // CUDA 10.0 and newer
            return attributes.type == cudaMemoryTypeDevice;
#else
            return attributes.memoryType == cudaMemoryTypeDevice;
#endif
        }

        template <typename T>
        T* soa_device_row(T** const SCCD_RESTRICT arrays, const int row) {
            T* host_rows[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
            if (is_ptr_device(arrays)) {
                SCCD_CHECK_CUDA(cudaMemcpy(host_rows, arrays, 2 * SCCD_DIM * sizeof(T*), cudaMemcpyDeviceToHost));
            } else {
                for (int d = 0; d < 2 * SCCD_DIM; d++) {
                    host_rows[d] = arrays[d];
                }
            }
            return host_rows[row];
        }

        template <typename T, typename I>
        void sort_along_axis(const ptrdiff_t n,
                             const int sort_axis,
                             T** const SCCD_RESTRICT arrays,
                             I* const SCCD_RESTRICT idx,
                             T* const SCCD_RESTRICT scratch) {
            if (n <= 0) return;

            SCCD_CUDA_LAST_ERROR();

            if (n == 1) {
                enumerate<I>(0, n, idx);
                SCCD_CUDA_LAST_ERROR();
                return;
            }

            T* host_arrays[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

            if (is_ptr_device(arrays)) {
                SCCD_CHECK_CUDA(cudaMemcpy(host_arrays, arrays, 2 * SCCD_DIM * sizeof(T*), cudaMemcpyDeviceToHost));
            } else {
                for (int d = 0; d < 2 * SCCD_DIM; d++) {
                    host_arrays[d] = arrays[d];
                }
            }

            const T* const SCCD_RESTRICT x = host_arrays[sort_axis];

            I* const tmp_idx = workspace(WorkspaceSlot::SortIndex).get_as<I>(n);
            enumerate<I>(0, n, tmp_idx);
            SCCD_CUDA_LAST_ERROR();

            size_t tmp_storage_bytes = 0;
            SCCD_CHECK_CUDA(cub::DeviceRadixSort::SortPairs(nullptr, tmp_storage_bytes, x, scratch, tmp_idx, idx, n));
            void* const tmp_storage = workspace(WorkspaceSlot::TempStorage).get(tmp_storage_bytes);
            SCCD_CHECK_CUDA(
                cub::DeviceRadixSort::SortPairs(tmp_storage, tmp_storage_bytes, x, scratch, tmp_idx, idx, n));

            SCCD_CHECK_CUDA(cudaMemcpy(host_arrays[sort_axis], scratch, n * sizeof(T), cudaMemcpyDeviceToDevice));

            // NOTE: the permuted rows have to end up in the caller's own
            // buffers, so the copy back out of `scratch` cannot be avoided by
            // swapping pointers here -- `host_arrays` is only a local copy of
            // them. Removing these copies would require the permutation to be
            // done against caller-owned double buffers.
            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((n + block.x - 1) / block.x);
            for (int d = 0; d < 2 * SCCD_DIM; d++) {
                if (d == sort_axis) continue;
                permute_kernel<T, I><<<grid, block>>>(n, idx, host_arrays[d], scratch);
                SCCD_CUDA_LAST_ERROR();
                SCCD_CHECK_CUDA(cudaMemcpy(host_arrays[d], scratch, n * sizeof(T), cudaMemcpyDeviceToDevice));
            }

            SCCD_CUDA_LAST_ERROR();
        }

        template <typename T>
        static inline __device__ void compute_candidate_window_progressive(const T fimin,
                                                                           const T fimax,
                                                                           const T* const SCCD_RESTRICT second_xmax,
                                                                           const T* const SCCD_RESTRICT second_xmin,
                                                                           const ptrdiff_t second_count,
                                                                           ptrdiff_t& begin,
                                                                           ptrdiff_t& end) {
            for (; begin < second_count; ++begin) {
                // Inclusive: the overlap predicate counts touching boxes as
                // overlapping (it rejects only on a strict `amin > bmax`), so the
                // window that feeds it has to keep a box whose xmax equals fimin.
                // With a strict `<` here the two disagreed, and the sweep dropped
                // real pairs -- a false negative, which the conservativeness
                // invariant does not allow. Boxes are sorted by xmin, so this can
                // only bite when xmin == xmax == fimin, i.e. a zero-extent box
                // sitting exactly at another box's lower bound on the sort axis.
                // That is what an axis-aligned face produces, and it cost 20 of
                // 2220 edge-edge pairs on a refined cube.
                if (fimin <= second_xmax[begin]) {
                    break;
                }
            }
            end = begin;
            for (; end < second_count; ++end) {
                if (fimax < second_xmin[end]) {
                    break;
                }
            }
        }

        template <int nxe, typename I>
        static inline __device__ void load_ev(I** const SCCD_RESTRICT elements,
                                              const I elem_idx,
                                              const ptrdiff_t stride,
                                              I (&out)[nxe]) {
            for (int v = 0; v < nxe; ++v) {
                out[v] = elements[v][elem_idx * stride];
            }
        }

        template <int n1, int n2, typename I>
        static inline __device__ bool shares_vertex(const I (&a)[n1], const I (&b)[n2]) {
            for (int i = 0; i < n1; ++i) {
                for (int j = 0; j < n2; ++j) {
                    if (a[i] == b[j]) {
                        return true;
                    }
                }
            }
            return false;
        }

        template <typename T>
        static __device__ __forceinline__ uint32_t disjoint(const T aminx,
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
            return aminx > bmaxx | aminy > bmaxy | aminz > bmaxz | bminx > amaxx | bminy > amaxy | bminz > amaxz;
        }

        template <int nxe, typename T, typename I>
        __global__ void count_self_overlaps_kernel(const int sort_axis,
                                                   const ptrdiff_t element_count,
                                                   T** const SCCD_RESTRICT aabbs,
                                                   I* const SCCD_RESTRICT idx,
                                                   const ptrdiff_t stride,
                                                   I** const SCCD_RESTRICT elements,
                                                   ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= element_count) return;
            if (fi == 0) {
                ccdptr[0] = 0;
            }

            const T* const SCCD_RESTRICT xmin = aabbs[sort_axis];
            const T* const SCCD_RESTRICT xmax = aabbs[3 + sort_axis];

            const T fimin = xmin[fi];
            const T fimax = xmax[fi];
            const I idxi = idx[fi];

            assert(idxi >= 0);
            assert(idxi < element_count);

            I ev[nxe];
            for (int v = 0; v < nxe; v++) {
                ev[v] = elements[v][idxi * stride];
            }

            ptrdiff_t begin = fi + 1;
            ptrdiff_t end = begin;
            compute_candidate_window_progressive(fimin, fimax, xmax, xmin, element_count, begin, end);

            if (begin >= end) {
                ccdptr[fi + 1] = 0;
                return;
            }

            ptrdiff_t count = 0;
            const T aminx = aabbs[0][fi];
            const T aminy = aabbs[1][fi];
            const T aminz = aabbs[2][fi];
            const T amaxx = aabbs[3][fi];
            const T amaxy = aabbs[4][fi];
            const T amaxz = aabbs[5][fi];

            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
                             aminy,
                             aminz,
                             amaxx,
                             amaxy,
                             amaxz,
                             aabbs[0][j],
                             aabbs[1][j],
                             aabbs[2][j],
                             aabbs[3][j],
                             aabbs[4][j],
                             aabbs[5][j])) {
                    continue;
                }
                const I jidx = idx[j];
                I sev[nxe];
                load_ev<nxe>(elements, jidx, stride, sev);
                const bool share = shares_vertex<nxe, nxe>(ev, sev);
                count += share ? 0 : 1;
            }

            ccdptr[fi + 1] = count;
        }

        template <int nxe, typename T, typename I>
        void count_self_overlaps(const int sort_axis,
                                 const ptrdiff_t element_count,
                                 T** const SCCD_RESTRICT aabbs,
                                 I* const SCCD_RESTRICT idx,
                                 const ptrdiff_t stride,
                                 I** const SCCD_RESTRICT elements,
                                 ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            SCCD_CUDA_LAST_ERROR();

            if (element_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }

            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((element_count + block.x - 1) / block.x);
            count_self_overlaps_kernel<nxe, T, I>
                <<<grid, block>>>(sort_axis, element_count, aabbs, idx, stride, elements, ccdptr);
            SCCD_CUDA_LAST_ERROR();

            size_t tmp_storage_bytes = 0;
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(nullptr, tmp_storage_bytes, ccdptr, ccdptr, element_count + 1));
            void* const tmp_storage = workspace(WorkspaceSlot::TempStorage).get(tmp_storage_bytes);
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(tmp_storage, tmp_storage_bytes, ccdptr, ccdptr, element_count + 1));

            SCCD_CUDA_LAST_ERROR();
        }

        template <int nxe, typename T, typename I>
        __global__ void collect_self_overlaps_kernel(const int sort_axis,
                                                     const ptrdiff_t element_count,
                                                     T** const SCCD_RESTRICT aabbs,
                                                     I* const SCCD_RESTRICT idx,
                                                     const ptrdiff_t stride,
                                                     I** const elements,
                                                     const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                     I* SCCD_RESTRICT foverlap,
                                                     I* SCCD_RESTRICT noverlap) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= element_count) return;

            const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
            if (!expected_count) return;

            const T fimin = aabbs[sort_axis][fi];
            const T fimax = aabbs[3 + sort_axis][fi];
            const I idxi = idx[fi];

            I ev[nxe];
            for (int v = 0; v < nxe; v++) {
                ev[v] = elements[v][idxi * stride];
            }

            ptrdiff_t begin = fi + 1;
            ptrdiff_t end = begin;
            compute_candidate_window_progressive(
                fimin, fimax, aabbs[3 + sort_axis], aabbs[sort_axis], element_count, begin, end);

            if (begin >= end) {
                return;
            }

            I* SCCD_RESTRICT const first_local_elements = &foverlap[ccdptr[fi]];
            I* SCCD_RESTRICT const second_local_elements = &noverlap[ccdptr[fi]];

            ptrdiff_t count = 0;
            const T aminx = aabbs[0][fi];
            const T aminy = aabbs[1][fi];
            const T aminz = aabbs[2][fi];
            const T amaxx = aabbs[3][fi];
            const T amaxy = aabbs[4][fi];
            const T amaxz = aabbs[5][fi];

            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
                             aminy,
                             aminz,
                             amaxx,
                             amaxy,
                             amaxz,
                             aabbs[0][j],
                             aabbs[1][j],
                             aabbs[2][j],
                             aabbs[3][j],
                             aabbs[4][j],
                             aabbs[5][j])) {
                    continue;
                }
                const I jidx = idx[j];
                I sev[nxe];
                load_ev<nxe>(elements, jidx, stride, sev);
                if (shares_vertex<nxe, nxe>(ev, sev)) {
                    continue;
                }

                first_local_elements[count] = SCCD_MIN(idxi, jidx);
                second_local_elements[count] = SCCD_MAX(idxi, jidx);

                count++;
            }

            assert(expected_count == count);
        }

        template <int nxe, typename T, typename I>
        void collect_self_overlaps(const int sort_axis,
                                   const ptrdiff_t element_count,
                                   T** const SCCD_RESTRICT aabbs,
                                   I* const SCCD_RESTRICT idx,
                                   const ptrdiff_t stride,
                                   I** const elements,
                                   const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                   I* SCCD_RESTRICT foverlap,
                                   I* SCCD_RESTRICT noverlap) {
            SCCD_CUDA_LAST_ERROR();

            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((element_count + block.x - 1) / block.x);
            collect_self_overlaps_kernel<nxe, T, I>
                <<<grid, block>>>(sort_axis, element_count, aabbs, idx, stride, elements, ccdptr, foverlap, noverlap);

            SCCD_CUDA_LAST_ERROR();
        }

        template <typename T>
        void cummax(const ptrdiff_t n, const T* const SCCD_RESTRICT in, T* const SCCD_RESTRICT out) {
            if (n <= 0) return;
            SCCD_CUDA_LAST_ERROR();

            size_t tmp_storage_bytes = 0;

            SCCD_CHECK_CUDA(cub::DeviceScan::InclusiveScan(nullptr, tmp_storage_bytes, in, out, device_max_op<T>(), n));
            void* const tmp_storage = workspace(WorkspaceSlot::TempStorage).get(tmp_storage_bytes);
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveScan(tmp_storage, tmp_storage_bytes, in, out, device_max_op<T>(), n));

            SCCD_CUDA_LAST_ERROR();
        }

        template <int first_nxe, int second_nxe, typename T, typename I>
        static inline __device__ ptrdiff_t count_overlaps_kernel_aux(const ptrdiff_t fi,
                                                                     ptrdiff_t begin,
                                                                     const int sort_axis,
                                                                     const ptrdiff_t first_count,
                                                                     T** const SCCD_RESTRICT first_aabbs,
                                                                     I* const SCCD_RESTRICT first_idx,
                                                                     const ptrdiff_t first_stride,
                                                                     I** const SCCD_RESTRICT first_elements,
                                                                     const ptrdiff_t second_count,
                                                                     T** const SCCD_RESTRICT second_aabbs,
                                                                     I* const SCCD_RESTRICT second_idx,
                                                                     const ptrdiff_t second_stride,
                                                                     I** const SCCD_RESTRICT second_elements)

        {
            const T* const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
            const T* const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
            const T* const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
            const T* const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

            const T fimin = first_aabbs[sort_axis][fi];
            const T fimax = first_aabbs[3 + sort_axis][fi];
            const I first_idxi = first_idx[fi];

            assert(first_idxi >= 0);
            assert(first_idxi < first_count);

            I ev[first_nxe];
            for (int v = 0; v < first_nxe; v++) {
                ev[v] = first_elements[v][first_idxi * first_stride];
            }

            ptrdiff_t end = begin;
            compute_candidate_window_progressive(fimin, fimax, second_xmax, second_xmin, second_count, begin, end);

            if (begin >= end) {
                return 0;
            }

            ptrdiff_t count = 0;
            const T aminx = first_aabbs[0][fi];
            const T aminy = first_aabbs[1][fi];
            const T aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi];
            const T amaxy = first_aabbs[4][fi];
            const T amaxz = first_aabbs[5][fi];
            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
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
                bool share = false;
                const I jidx = second_idx[j];
                assert(jidx >= 0);
                assert(jidx < second_count);
                if constexpr (second_nxe > 1) {
                    I sev[second_nxe];
                    for (int v = 0; v < second_nxe; ++v) {
                        sev[v] = second_elements[v][jidx * second_stride];
                    }
                    share = shares_vertex<first_nxe, second_nxe>(ev, sev);
                } else {
                    for (int a = 0; a < first_nxe; ++a) {
                        if (ev[a] == jidx) {
                            share = true;
                            break;
                        }
                    }
                }
                count += share ? 0 : 1;
            }

            return count;
        }

        template <int first_nxe, int second_nxe, typename T, typename I>
        __global__ void count_overlaps_kernel(const int sort_axis,
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
                                              const T* const SCCD_RESTRICT second_xmax_running,
                                              ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= first_count) return;
            if (fi == 0) {
                ccdptr[0] = 0;
            }

            const T fimin = first_aabbs[sort_axis][fi];
            ptrdiff_t begin = lower_bound(second_xmax_running, second_xmax_running + second_count, fimin) - second_xmax_running;
            ccdptr[fi + 1] = count_overlaps_kernel_aux<first_nxe, second_nxe, T, I>(fi,
                                                                                    begin,
                                                                                    sort_axis,
                                                                                    first_count,
                                                                                    first_aabbs,
                                                                                    first_idx,
                                                                                    first_stride,
                                                                                    first_elements,
                                                                                    second_count,
                                                                                    second_aabbs,
                                                                                    second_idx,
                                                                                    second_stride,
                                                                                    second_elements);
        }

        template <int first_nxe, int second_nxe, typename T, typename I>
        void count_overlaps(const int sort_axis,
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
                            const T* const SCCD_RESTRICT second_xmax_running) {
            SCCD_CUDA_LAST_ERROR();

            if (first_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }

            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((first_count + block.x - 1) / block.x);
            count_overlaps_kernel<first_nxe, second_nxe, T, I><<<grid, block>>>(sort_axis,
                                                                                first_count,
                                                                                first_aabbs,
                                                                                first_idx,
                                                                                first_stride,
                                                                                first_elements,
                                                                                second_count,
                                                                                second_aabbs,
                                                                                second_idx,
                                                                                second_stride,
                                                                                second_elements,
                                                                                second_xmax_running,
                                                                                ccdptr);

            SCCD_CUDA_LAST_ERROR();

            size_t tmp_storage_bytes = 0;
            SCCD_CHECK_CUDA(cub::DeviceScan::InclusiveSum(nullptr, tmp_storage_bytes, ccdptr, ccdptr, first_count + 1));
            void* const tmp_storage = workspace(WorkspaceSlot::TempStorage).get(tmp_storage_bytes);
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(tmp_storage, tmp_storage_bytes, ccdptr, ccdptr, first_count + 1));

            SCCD_CUDA_LAST_ERROR();
        }

        template <int first_nxe, int second_nxe, typename T, typename I>
        static inline __device__ void collect_overlaps_kernel_aux(const ptrdiff_t fi,
                                                                  ptrdiff_t begin,
                                                                  const int sort_axis,
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
                                                                  I* SCCD_RESTRICT foverlap,
                                                                  I* SCCD_RESTRICT noverlap) {
            const T* const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
            const T* const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
            const T* const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
            const T* const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

#ifndef NDEBUG
            const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
#endif

            const T fimin = first_aabbs[sort_axis][fi];
            const T fimax = first_aabbs[3 + sort_axis][fi];
            const I first_idxi = first_idx[fi];

            assert(first_idxi >= 0);
            assert(first_idxi < first_count);

            I ev[first_nxe];
            for (int v = 0; v < first_nxe; v++) {
                ev[v] = first_elements[v][first_idxi * first_stride];
            }

            ptrdiff_t end = begin;
            compute_candidate_window_progressive(fimin, fimax, second_xmax, second_xmin, second_count, begin, end);

            if (begin >= end) {
                return;
            }

            I* SCCD_RESTRICT const first_local_elements = &foverlap[ccdptr[fi]];
            I* SCCD_RESTRICT const second_local_elements = &noverlap[ccdptr[fi]];

            ptrdiff_t count = 0;
            const T aminx = first_aabbs[0][fi];
            const T aminy = first_aabbs[1][fi];
            const T aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi];
            const T amaxy = first_aabbs[4][fi];
            const T amaxz = first_aabbs[5][fi];
            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
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
                bool share = false;
                const I jidx = second_idx[j];

                assert(jidx >= 0);
                assert(jidx < second_count);

                if constexpr (second_nxe > 1) {
                    I sev[second_nxe];
                    for (int v = 0; v < second_nxe; ++v) {
                        sev[v] = second_elements[v][jidx * second_stride];
                    }
                    share = shares_vertex<first_nxe, second_nxe>(ev, sev);
                } else {
                    for (int a = 0; a < first_nxe; ++a) {
                        if (ev[a] == jidx) {
                            share = true;
                            break;
                        }
                    }
                }

                if (share) continue;
                first_local_elements[count] = first_idxi;
                second_local_elements[count] = jidx;
                count += 1;
            }

            assert(expected_count == count);
        }

        template <int first_nxe, int second_nxe, typename T, typename I>
        __global__ void collect_overlaps_kernel(const int sort_axis,
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
                                                const T* const SCCD_RESTRICT second_xmax_running,
                                                I* SCCD_RESTRICT foverlap,
                                                I* SCCD_RESTRICT noverlap) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= first_count) return;
            const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
            if (!expected_count) return;

            const T fimin = first_aabbs[sort_axis][fi];

            ptrdiff_t begin = lower_bound(second_xmax_running, second_xmax_running + second_count, fimin) - second_xmax_running;

            collect_overlaps_kernel_aux<first_nxe, second_nxe, T, I>(fi,
                                                                     begin,
                                                                     sort_axis,
                                                                     first_count,
                                                                     first_aabbs,
                                                                     first_idx,
                                                                     first_stride,
                                                                     first_elements,
                                                                     second_count,
                                                                     second_aabbs,
                                                                     second_idx,
                                                                     second_stride,
                                                                     second_elements,
                                                                     ccdptr,
                                                                     foverlap,
                                                                     noverlap);
        }

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
                              const T* const SCCD_RESTRICT second_xmax_running,
                              I* SCCD_RESTRICT foverlap,
                              I* SCCD_RESTRICT noverlap) {
            SCCD_CUDA_LAST_ERROR();

            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((first_count + block.x - 1) / block.x);
            collect_overlaps_kernel<first_nxe, second_nxe, T, I><<<grid, block>>>(sort_axis,
                                                                                  first_count,
                                                                                  first_aabbs,
                                                                                  first_idx,
                                                                                  first_stride,
                                                                                  first_elements,
                                                                                  second_count,
                                                                                  second_aabbs,
                                                                                  second_idx,
                                                                                  second_stride,
                                                                                  second_elements,
                                                                                  ccdptr,
                                                                                  second_xmax_running,
                                                                                  foverlap,
                                                                                  noverlap);

            SCCD_CUDA_LAST_ERROR();
        }



    }  // namespace device
}  // namespace sccd

#define SCCD_BP_INSTANTIATE_CHOOSE_AXIS(T)             \
    template int sccd::device::choose_axis<T>( \
        const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs)

#define SCCD_BP_INSTANTIATE_ENUMERATE(T) \
    template void sccd::device::enumerate<T>(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx)

#define SCCD_BP_INSTANTIATE_SORT_ALONG_AXIS(T, I)                                             \
    template void sccd::device::sort_along_axis<T, I>(                                \
                                                      const ptrdiff_t n,              \
                                                      const int sort_axis,            \
                                                      T** const SCCD_RESTRICT arrays, \
                                                      I* const SCCD_RESTRICT idx,     \
                                                      T* const SCCD_RESTRICT scratch)

#define SCCD_BP_INSTANTIATE_COUNT_SELF_OVERLAPS(NXE, T, I)                                               \
    template void sccd::device::count_self_overlaps<NXE, T, I>(const int sort_axis,              \
                                                               const ptrdiff_t element_count,    \
                                                               T** const SCCD_RESTRICT aabbs,    \
                                                               I* const SCCD_RESTRICT idx,       \
                                                               const ptrdiff_t stride,           \
                                                               I** const SCCD_RESTRICT elements, \
                                                               ptrdiff_t* const SCCD_RESTRICT ccdptr)

#define SCCD_BP_INSTANTIATE_COLLECT_SELF_OVERLAPS(NXE, T, I)                                                          \
    template void sccd::device::collect_self_overlaps<NXE, T, I>(const int sort_axis,                         \
                                                                 const ptrdiff_t element_count,               \
                                                                 T** const SCCD_RESTRICT aabbs,               \
                                                                 I* const SCCD_RESTRICT idx,                  \
                                                                 const ptrdiff_t stride,                      \
                                                                 I** const SCCD_RESTRICT elements,            \
                                                                 const ptrdiff_t* const SCCD_RESTRICT ccdptr, \
                                                                 I* SCCD_RESTRICT foverlap,                   \
                                                                 I* SCCD_RESTRICT noverlap)

#define SCCD_BP_INSTANTIATE_CUMMAX(T)              \
    template void sccd::device::cummax<T>( \
        const ptrdiff_t n, const T* const SCCD_RESTRICT in, T* const SCCD_RESTRICT out)

#define SCCD_BP_INSTANTIATE_SOA_DEVICE_ROW(T) \
    template T* sccd::device::soa_device_row<T>(T * * const SCCD_RESTRICT arrays, const int row)

#define SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(FIRST_NXE, SECOND_NXE, T, I)                                                      \
    template void sccd::device::count_overlaps<FIRST_NXE, SECOND_NXE, T, I>(const int sort_axis,                     \
                                                                            const ptrdiff_t first_count,             \
                                                                            T** const SCCD_RESTRICT first_aabbs,     \
                                                                            I* const SCCD_RESTRICT first_idx,        \
                                                                            const ptrdiff_t first_stride,            \
                                                                            I** const SCCD_RESTRICT first_elements,  \
                                                                            const ptrdiff_t second_count,            \
                                                                            T** const SCCD_RESTRICT second_aabbs,    \
                                                                            I* const SCCD_RESTRICT second_idx,       \
                                                                            const ptrdiff_t second_stride,           \
                                                                            I** const SCCD_RESTRICT second_elements, \
                                                                            ptrdiff_t* const SCCD_RESTRICT ccdptr,   \
                                                                            const T* const SCCD_RESTRICT second_xmax_running)

#define SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(FIRST_NXE, SECOND_NXE, T, I)              \
    template void sccd::device::collect_overlaps<FIRST_NXE, SECOND_NXE, T, I>( \
        const int sort_axis,                                                   \
        const ptrdiff_t first_count,                                           \
        T** const SCCD_RESTRICT first_aabbs,                                   \
        I* const SCCD_RESTRICT first_idx,                                      \
        const ptrdiff_t first_stride,                                          \
        I** const SCCD_RESTRICT first_elements,                                \
        const ptrdiff_t second_count,                                          \
        T** const SCCD_RESTRICT second_aabbs,                                  \
        I* const SCCD_RESTRICT second_idx,                                     \
        const ptrdiff_t second_stride,                                         \
        I** const SCCD_RESTRICT second_elements,                               \
        const ptrdiff_t* const SCCD_RESTRICT ccdptr,                           \
        const T* const SCCD_RESTRICT second_xmax_running,                                   \
        I* SCCD_RESTRICT foverlap,                                             \
        I* SCCD_RESTRICT noverlap)



SCCD_BP_INSTANTIATE_CHOOSE_AXIS(float);
SCCD_BP_INSTANTIATE_CHOOSE_AXIS(double);
SCCD_BP_INSTANTIATE_ENUMERATE(int32_t);
SCCD_BP_INSTANTIATE_ENUMERATE(int64_t);
SCCD_BP_INSTANTIATE_SORT_ALONG_AXIS(float, int32_t);
SCCD_BP_INSTANTIATE_SORT_ALONG_AXIS(float, int64_t);
SCCD_BP_INSTANTIATE_SORT_ALONG_AXIS(double, int32_t);
SCCD_BP_INSTANTIATE_SORT_ALONG_AXIS(double, int64_t);

SCCD_BP_INSTANTIATE_COUNT_SELF_OVERLAPS(2, float, int32_t);
SCCD_BP_INSTANTIATE_COUNT_SELF_OVERLAPS(2, float, int64_t);
SCCD_BP_INSTANTIATE_COUNT_SELF_OVERLAPS(2, double, int32_t);
SCCD_BP_INSTANTIATE_COUNT_SELF_OVERLAPS(2, double, int64_t);

SCCD_BP_INSTANTIATE_COLLECT_SELF_OVERLAPS(2, float, int32_t);
SCCD_BP_INSTANTIATE_COLLECT_SELF_OVERLAPS(2, float, int64_t);
SCCD_BP_INSTANTIATE_COLLECT_SELF_OVERLAPS(2, double, int32_t);
SCCD_BP_INSTANTIATE_COLLECT_SELF_OVERLAPS(2, double, int64_t);

SCCD_BP_INSTANTIATE_CUMMAX(float);
SCCD_BP_INSTANTIATE_CUMMAX(double);

SCCD_BP_INSTANTIATE_SOA_DEVICE_ROW(float);
SCCD_BP_INSTANTIATE_SOA_DEVICE_ROW(double);

SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(3, 1, float, int32_t);
SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(3, 1, float, int64_t);
SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(3, 1, double, int32_t);
SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(3, 1, double, int64_t);

SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(3, 1, float, int32_t);
SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(3, 1, float, int64_t);
SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(3, 1, double, int32_t);
SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(3, 1, double, int64_t);

// Quads. broad_phase_fv_step_device_ dispatches QUADSHELL4 to <4, 1>, so
// without these the CUDA build does not link at all once a quad mesh is
// reachable -- which it always is, the dispatch being on a runtime element type.
SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(4, 1, float, int32_t);
SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(4, 1, float, int64_t);
SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(4, 1, double, int32_t);
SCCD_BP_INSTANTIATE_COUNT_OVERLAPS(4, 1, double, int64_t);

SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(4, 1, float, int32_t);
SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(4, 1, float, int64_t);
SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(4, 1, double, int32_t);
SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS(4, 1, double, int64_t);

#undef SCCD_BP_INSTANTIATE_CHOOSE_AXIS
#undef SCCD_BP_INSTANTIATE_ENUMERATE
#undef SCCD_BP_INSTANTIATE_SORT_ALONG_AXIS
#undef SCCD_BP_INSTANTIATE_COUNT_SELF_OVERLAPS
#undef SCCD_BP_INSTANTIATE_COLLECT_SELF_OVERLAPS
#undef SCCD_BP_INSTANTIATE_CUMMAX
#undef SCCD_BP_INSTANTIATE_SOA_DEVICE_ROW
#undef SCCD_BP_INSTANTIATE_COUNT_OVERLAPS
#undef SCCD_BP_INSTANTIATE_COLLECT_OVERLAPS

// Clean-up kernel macros
#undef SCCD_BP_N_WARPS_PER_BLOCK
