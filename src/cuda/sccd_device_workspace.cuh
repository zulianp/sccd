#ifndef SCCD_DEVICE_WORKSPACE_CUH
#define SCCD_DEVICE_WORKSPACE_CUH

#include <cstddef>

#include "sccd_cuda_base.cuh"

namespace sccd {
    namespace device {

        /**
         * \brief A grow-only device scratch buffer that outlives a single call.
         *
         * cudaMalloc and cudaFree are device-wide synchronizing operations, so
         * allocating CUB temp storage (or any other scratch) per call inserts a
         * full pipeline stall on every broadphase/narrowphase invocation. CCD is
         * driven per simulation step with sizes that barely change, so the
         * buffer is kept and only reallocated when it has to grow.
         *
         * Not thread-safe: one workspace instance belongs to one host thread
         * driving one stream, which matches how the pipeline is invoked.
         */
        class DeviceWorkspace {
        public:
            DeviceWorkspace() = default;

            DeviceWorkspace(const DeviceWorkspace&) = delete;
            DeviceWorkspace& operator=(const DeviceWorkspace&) = delete;

            ~DeviceWorkspace() { release(); }

            /**
             * \brief Pointer to at least \p bytes of device scratch.
             *
             * The contents are undefined and are invalidated by any later call
             * that grows the buffer.
             */
            void* get(const size_t bytes) {
                if (bytes == 0) {
                    return nullptr;
                }

                if (bytes > capacity_) {
                    release();
                    // Round up so a slowly growing workload does not reallocate
                    // on every single step.
                    size_t new_capacity = capacity_ ? capacity_ : size_t(1) << 12;
                    while (new_capacity < bytes) {
                        new_capacity *= 2;
                    }
                    SCCD_CHECK_CUDA(cudaMalloc(&data_, new_capacity));
                    capacity_ = new_capacity;
                }

                return data_;
            }

            /** \brief Typed convenience wrapper around get(). */
            template <typename T>
            T* get_as(const size_t count) {
                return static_cast<T*>(get(count * sizeof(T)));
            }

            void release() {
                if (data_) {
                    cudaFree(data_);
                    data_ = nullptr;
                }
                capacity_ = 0;
            }

            size_t capacity() const { return capacity_; }

        private:
            void* data_{nullptr};
            size_t capacity_{0};
        };

        /**
         * \brief Process-wide scratch buffers, one per distinct concurrent use.
         *
         * Separate slots so that two buffers needed at the same time (for
         * example CUB temp storage and the sort's index array) never alias.
         */
        enum class WorkspaceSlot { TempStorage = 0, SortIndex, Scratch, Count };

        inline DeviceWorkspace& workspace(const WorkspaceSlot slot) {
            static DeviceWorkspace slots[static_cast<int>(WorkspaceSlot::Count)];
            return slots[static_cast<int>(slot)];
        }

        /** \brief Pinned host staging for small device-to-host reads. */
        class PinnedHostScalar {
        public:
            PinnedHostScalar() = default;

            PinnedHostScalar(const PinnedHostScalar&) = delete;
            PinnedHostScalar& operator=(const PinnedHostScalar&) = delete;

            ~PinnedHostScalar() {
                if (data_) {
                    cudaFreeHost(data_);
                }
            }

            template <typename T>
            T* get() {
                if (!data_) {
                    SCCD_CHECK_CUDA(cudaHostAlloc(&data_, sizeof(long long), cudaHostAllocDefault));
                }
                return static_cast<T*>(data_);
            }

        private:
            void* data_{nullptr};
        };

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_DEVICE_WORKSPACE_CUH
