// choose_axis must return the axis of largest centre variance, and must return
// it whatever ran on the device before.
//
// Its two kernels reduce per block through a shared array indexed by warp, one
// slot written by lane 0 of each warp. A thread whose element index is past the
// end has no element to contribute, and the obvious way to say so -- return
// before the reduction -- leaves every slot belonging to a wholly out-of-range
// warp unwritten, while the cross-warp pass that follows reads all of them. The
// block's sum then carries whatever that streaming multiprocessor's shared
// memory happened to hold, which is normally the previous launch's partial sums.
// The last block is partial for almost every mesh, so the variance the axis is
// chosen from is contaminated on essentially every call; whether the choice
// changes depends on residue, which is what made it look intermittent. Under
// live device asserts it surfaces as `acc == acc` failing in
// block_reduce_to_gmem when the residue is a NaN bit pattern.
//
// So the test runs a large call with enormous coordinates to fill those slots
// with values far above anything the probe can produce, then asks for the axis
// of a small input whose answer is unambiguous. Sizes are chosen so the probe's
// single block leaves a whole warp out of range.
//
// Deliberately free of smesh: it allocates with cudaMalloc directly, so it runs
// in any CUDA build.

#include "sccd_broadphase.cuh"

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

using scalar_t = double;

namespace {

#define SCCD_CUDA_CHECK(call)                                                                         \
    do {                                                                                              \
        const cudaError_t err_ = (call);                                                              \
        if (err_ != cudaSuccess) {                                                                    \
            std::fprintf(stderr, "CUDA %s at %s:%d\n", cudaGetErrorString(err_), __FILE__, __LINE__); \
            std::exit(2);                                                                             \
        }                                                                                             \
    } while (0)

    // choose_axis reads a device array of six row pointers: the three lower
    // bounds then the three upper bounds, each of length n.
    struct DeviceBoxes {
        scalar_t* rows[6] = {};
        scalar_t** aabbs = nullptr;

        DeviceBoxes(const std::vector<std::vector<scalar_t>>& host) {
            for (int d = 0; d < 6; ++d) {
                SCCD_CUDA_CHECK(cudaMalloc(&rows[d], sizeof(scalar_t) * host[d].size()));
                SCCD_CUDA_CHECK(cudaMemcpy(
                    rows[d], host[d].data(), sizeof(scalar_t) * host[d].size(), cudaMemcpyHostToDevice));
            }
            SCCD_CUDA_CHECK(cudaMalloc(&aabbs, sizeof(scalar_t*) * 6));
            SCCD_CUDA_CHECK(cudaMemcpy(aabbs, rows, sizeof(scalar_t*) * 6, cudaMemcpyHostToDevice));
        }

        ~DeviceBoxes() {
            for (int d = 0; d < 6; ++d) cudaFree(rows[d]);
            cudaFree(aabbs);
        }
    };

    // Degenerate boxes at the given centres, which is all choose_axis looks at.
    std::vector<std::vector<scalar_t>> boxes_at(const std::vector<scalar_t>& cx,
                                                const std::vector<scalar_t>& cy,
                                                const std::vector<scalar_t>& cz) {
        return {cx, cy, cz, cx, cy, cz};
    }

    int expected_axis(const std::vector<std::vector<scalar_t>>& boxes) {
        const size_t n = boxes[0].size();
        int best = 0;
        scalar_t best_var = -1;
        for (int d = 0; d < 3; ++d) {
            scalar_t mean = 0;
            for (size_t i = 0; i < n; ++i) mean += (boxes[d][i] + boxes[3 + d][i]) / 2;
            mean /= (scalar_t)n;
            scalar_t var = 0;
            for (size_t i = 0; i < n; ++i) {
                const scalar_t p = (boxes[d][i] + boxes[3 + d][i]) / 2 - mean;
                var += p * p;
            }
            if (var > best_var) {
                best_var = var;
                best = d;
            }
        }
        return best;
    }

}  // namespace

int main() {
    // The probe. 200 elements against a 256-thread block leaves the last warp
    // wholly out of range, so one shared slot belongs to no active lane. The
    // spread is far wider on z than on x or y, so the answer is 2 and nothing
    // near the rounding of a double can move it.
    const int n_probe = 200;
    std::vector<scalar_t> px(n_probe), py(n_probe), pz(n_probe);
    for (int i = 0; i < n_probe; ++i) {
        px[i] = 0.001 * (scalar_t)(i % 7);
        py[i] = 0.001 * (scalar_t)(i % 5);
        pz[i] = 100.0 * (scalar_t)i;
    }
    const auto probe = boxes_at(px, py, pz);
    const int want = expected_axis(probe);
    if (want != 2) {
        std::fprintf(stderr, "test is malformed: expected axis 2, host says %d\n", want);
        return 2;
    }

    // The dirtying call. Coordinates around 1e7 give per-warp partial sums near
    // 1e16 once squared, which is far beyond anything the probe produces, so a
    // stale slot cannot be mistaken for a plausible variance.
    const int n_dirty = 4096;
    std::vector<scalar_t> dx(n_dirty), dy(n_dirty), dz(n_dirty);
    for (int i = 0; i < n_dirty; ++i) {
        dx[i] = 1.0e7 * (scalar_t)(i + 1);
        dy[i] = 2.0e7 * (scalar_t)(i + 1);
        dz[i] = 3.0e7 * (scalar_t)(i + 1);
    }
    const auto dirty = boxes_at(dx, dy, dz);

    DeviceBoxes d_probe(probe);
    DeviceBoxes d_dirty(dirty);

    const int clean = sccd::device::choose_axis<scalar_t>(n_probe, d_probe.aabbs);
    SCCD_CUDA_CHECK(cudaDeviceSynchronize());
    if (clean != want) {
        std::fprintf(stderr, "choose_axis returned %d on an untouched device, expected %d\n", clean, want);
        return 1;
    }

    // Residue depends on what the scheduler put where, so repeat: one pass that
    // happens to land on a clean multiprocessor proves nothing.
    const int rounds = 64;
    int wrong = 0;
    for (int r = 0; r < rounds; ++r) {
        (void)sccd::device::choose_axis<scalar_t>(n_dirty, d_dirty.aabbs);
        SCCD_CUDA_CHECK(cudaDeviceSynchronize());

        const int got = sccd::device::choose_axis<scalar_t>(n_probe, d_probe.aabbs);
        SCCD_CUDA_CHECK(cudaDeviceSynchronize());
        if (got != want) {
            if (wrong < 5) {
                std::fprintf(stderr, "round %d: choose_axis returned %d, expected %d\n", r, got, want);
            }
            ++wrong;
        }
    }

    if (wrong) {
        std::fprintf(stderr, "choose_axis disagreed with the host in %d of %d rounds\n", wrong, rounds);
        return 1;
    }

    std::printf("choose_axis: axis %d over %d rounds, unaffected by the preceding launch\n", want, rounds);
    return 0;
}
