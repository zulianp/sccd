#ifndef SCCD_DEVICE_DFS_STACK_CUH
#define SCCD_DEVICE_DFS_STACK_CUH

#include "sccd_base.hpp"

#include <cuda_runtime.h>

/**
 * \file
 * \brief The branch-and-bound work stack shared by the device narrow phases.
 *
 * Every device narrow phase searches boxes in (t, u, v) by the same scheme: a
 * small per-block stack in shared memory, and a global queue for the boxes that
 * do not fit in it, double-buffered so that a launch reads one buffer and writes
 * the other. What differs between the vertex-face, edge-edge and vertex-quad
 * kernels is the inclusion function, not the bookkeeping.
 *
 * This lived inside sccd_narrowphase.cu, which meant the vertex-quad kernel in
 * its own translation unit could not reach it and kept its stack in thread-local
 * memory instead -- an 8 KB stack frame and 255 registers, the only kernel in
 * the build that spilled.
 */

// Entries in the per-block shared stack. An entry is six scalars plus two ints,
// so 1024 entries is 32,800 bytes of static __shared__ for float and 57,376 for
// double. See SharedStackCap below for why that is legal and what it costs.
#ifndef SCCD_NP_SHARED_STACK_CAP
#define SCCD_NP_SHARED_STACK_CAP 64
#endif

namespace sccd {
    namespace device {


        // Shared-stack capacity for scalar type T.
        //
        // Small on purpose. The block-local stack is where a query's subtree
        // lives, so a large one keeps a heavy query inside the block that seeded
        // it: measured per query on cloth-funnel, the conservative search needs
        // one box for most queries and over a million for 208 of them, and at
        // 1024 entries per block 99% of pushes never reach the global queue at
        // all. One block then grinds through a 19.7-million-box query while its
        // neighbours finish and idle.
        //
        // Spilling early hands those boxes to the global queue, which the next
        // drain round redistributes across every block. Measured on GH200, mode 2,
        // 16 cases (narrow-phase milliseconds):
        //
        //                        cap 1024   cap 256   cap 64
        //   cloth-funnel           4059.2     198.2     37.2
        //   armadillo-rollers      3133.7     277.2     75.8
        //   cloth-ball              129.8     110.2    105.9
        //
        // False positives and negatives are identical at every capacity, and mode
        // 0 -- which averages under four boxes per query and so never fills even a
        // small stack -- is unaffected.
        //
        // This only became available once the global queue stopped deadlocking
        // under sustained use; see the Stack comment above. Before that, shrinking
        // this number hung the kernel, which is why it was 1024.
        //
        // 64 entries is 3,584 bytes for double against 57,376 at 1024. That does
        // not buy occupancy -- at 238 registers per thread the kernel is capped at
        // two blocks per SM either way -- it buys balance.
        template <typename T>
        struct SharedStackCap {
            static constexpr int value = SCCD_NP_SHARED_STACK_CAP;
        };

        // Axis-aligned subdomain in (t, u, v) parameter space.
        template <typename T>
        struct Domain {
            T tlower;
            T tupper;
            T ulower;
            T uupper;
            T vlower;
            T vupper;
        };

        /**
         * \brief One half of the double-buffered global work queue.
         *
         * A launch reads from one buffer and writes to the other, never both, so
         * there is no producer/consumer handshake here and no spinning. Reading
         * uses `top` as a claim cursor over the `count` entries a previous launch
         * left; writing uses `top` as a bump allocator and `count` is unused.
         *
         * This replaces a single shared buffer in which producers and consumers
         * ran concurrently and coordinated through the `qid` array: a writer
         * spun on `atomicCAS(&qid[slot], EMPTY, WRITING)` until a reader released
         * the slot, and a reader spun until the writer committed. With the queue
         * lightly used that handshake was invisible; driven hard -- which is
         * exactly what using it to balance the load requires -- it hung the
         * kernel. Separating the buffers removes the interaction rather than
         * tuning around it.
         */
        template <typename T>
        struct Stack {
            T* tlower;
            T* tupper;
            T* ulower;
            T* uupper;
            T* vlower;
            T* vupper;
            int* level;
            int* qid;
            int* top;      // write cursor when writing, claim cursor when reading
            int* request;  // boxes that did not fit, for the host to size the retry
            int capacity;
            int count;     // entries available to readers; 0 for a write buffer
        };

        /**
         * \brief Claim one slot and write a box to it. No waiting, ever.
         *
         * The index comes from a bump allocator, so it belongs to this thread
         * alone and nobody reads the buffer until the launch ends. A box that does
         * not fit is counted in `request` and dropped, which is safe -- dropping
         * can only leave the time of impact earlier than the truth, never later --
         * and the host grows the queue and retries.
         */
        template <typename T>
        static inline __device__ void push_global(const Stack<T>& out,
                                                  const Domain<T>& box,
                                                  const int level,
                                                  const int qid) {
            const int slot = atomicAdd(out.top, 1);
            if (slot >= out.capacity) {
                atomicAdd(out.request, 1);
                return;
            }
            out.tlower[slot] = box.tlower;
            out.tupper[slot] = box.tupper;
            out.ulower[slot] = box.ulower;
            out.uupper[slot] = box.uupper;
            out.vlower[slot] = box.vlower;
            out.vupper[slot] = box.vupper;
            out.level[slot] = level;
            out.qid[slot] = qid;
        }

        /**
         * \brief Take the next box from the read buffer, or return 0 when empty.
         *
         * A plain claim on a cursor over entries a previous launch committed, so
         * the loads need no ordering against a concurrent writer -- there is none.
         */
        template <typename T>
        static inline __device__ int pop_global(const Stack<T>& in,
                                                Domain<T>& d,
                                                int& level,
                                                int& qid) {
            if (in.count <= 0) return 0;
            const int slot = atomicAdd(in.top, 1);
            if (slot >= in.count) return 0;

            d.tlower = in.tlower[slot];
            d.tupper = in.tupper[slot];
            d.ulower = in.ulower[slot];
            d.uupper = in.uupper[slot];
            d.vlower = in.vlower[slot];
            d.vupper = in.vupper[slot];
            level = in.level[slot];
            qid = in.qid[slot];
            return 1;
        }


        /**
         * \brief DFS stack buffers that persist between narrow-phase calls.
         *
         * cudaMalloc and cudaFree synchronize the whole device. Allocating and
         * releasing the eight stack arrays plus the two counters on every
         * invocation put twenty such stalls into each CCD step, and because the
         * capacity was a local it always restarted from zero and re-grew.
         * Keeping them alive across calls also preserves the capacity the
         * previous step converged on, so steady state does no allocation at all.
         *
         * One instance per scalar type; like the rest of the device pipeline it
         * assumes a single host thread drives it.
         */
        template <typename T>
        struct PersistentDfsStack {
            T* tlower{nullptr};
            T* tupper{nullptr};
            T* ulower{nullptr};
            T* uupper{nullptr};
            T* vlower{nullptr};
            T* vupper{nullptr};
            int* level{nullptr};
            int* qid{nullptr};
            int* counters{nullptr};       // device, [0]/[1] = write cursors for the
                                          // two buffers, [2] = g_request
            int* host_counters{nullptr};  // pinned staging for reading them back
            int cap{0};

            ~PersistentDfsStack() { release(); }

            void release() {
                cudaFree(tlower);
                cudaFree(tupper);
                cudaFree(ulower);
                cudaFree(uupper);
                cudaFree(vlower);
                cudaFree(vupper);
                cudaFree(level);
                cudaFree(qid);
                cudaFree(counters);
                cudaFreeHost(host_counters);
                tlower = tupper = ulower = uupper = vlower = vupper = nullptr;
                level = qid = counters = nullptr;
                host_counters = nullptr;
                cap = 0;
            }
        };

        template <typename T>
        static PersistentDfsStack<T>& persistent_dfs_stack() {
            static PersistentDfsStack<T> stack;
            return stack;
        }


        static inline __device__ int reserve_slots(int* SCCD_RESTRICT counter, int k, int capacity) {
            int old = atomicAdd(counter, 0);
            while (true) {
                if (old + k > capacity) return -1;
                const int prev = atomicCAS(counter, old, old + k);
                if (prev == old) return old;
                old = prev;
            }
        }




    }  // namespace device
}  // namespace sccd

#endif  // SCCD_DEVICE_DFS_STACK_CUH
