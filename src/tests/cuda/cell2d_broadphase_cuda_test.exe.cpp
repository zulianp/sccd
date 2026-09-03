// The device cell list must find exactly the pair set the device sweep finds.
//
// The host versions are compared against each other in
// src/tests/cell2d_broadphase_test; this is the same check on the device, because
// the two ports differ in ways that could diverge: the device binning scatters
// with atomics rather than a serial cursor, so cell contents come out in a
// different order, and the grid is sized from extents computed by a reduction
// kernel rather than on the host.
//
// Deliberately free of smesh: it allocates with cudaMalloc directly, so it runs
// in any CUDA build.

#include "sccd_broadphase.cuh"
#include "sccd_cell2d_broadphase.cuh"

#include <cuda_runtime.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <random>
#include <set>
#include <string>
#include <vector>

using scalar_t = double;
using idx_t = std::int32_t;

namespace {

#define CHECK(call)                                                                          \
    do {                                                                                     \
        const cudaError_t err_ = (call);                                                     \
        if (err_ != cudaSuccess) {                                                           \
            std::fprintf(stderr, "CUDA %s at %s:%d\n", cudaGetErrorString(err_), __FILE__, __LINE__); \
            std::exit(2);                                                                    \
        }                                                                                    \
    } while (0)

    template <typename T>
    T* dup(const std::vector<T>& h) {
        T* d = nullptr;
        CHECK(cudaMalloc(&d, sizeof(T) * (h.empty() ? 1 : h.size())));
        if (!h.empty()) CHECK(cudaMemcpy(d, h.data(), sizeof(T) * h.size(), cudaMemcpyHostToDevice));
        return d;
    }

    struct DeviceBoxes {
        std::vector<scalar_t> h[6];
        std::vector<idx_t> h_elem[3];
        scalar_t* rows[6] = {};
        idx_t* elem_rows[3] = {};
        scalar_t** aabbs = nullptr;
        idx_t** elements = nullptr;
        idx_t* idx = nullptr;
        ptrdiff_t n = 0;
        int nxe = 0;

        void upload() {
            for (int d = 0; d < 6; ++d) rows[d] = dup(h[d]);
            aabbs = nullptr;
            CHECK(cudaMalloc(&aabbs, sizeof(scalar_t*) * 6));
            CHECK(cudaMemcpy(aabbs, rows, sizeof(scalar_t*) * 6, cudaMemcpyHostToDevice));

            for (int v = 0; v < nxe; ++v) elem_rows[v] = dup(h_elem[v]);
            CHECK(cudaMalloc(&elements, sizeof(idx_t*) * 3));
            CHECK(cudaMemcpy(elements, elem_rows, sizeof(idx_t*) * (size_t)nxe, cudaMemcpyHostToDevice));

            std::vector<idx_t> ids(n);
            for (ptrdiff_t i = 0; i < n; ++i) ids[i] = (idx_t)i;
            idx = dup(ids);
        }

        void free_all() {
            for (int d = 0; d < 6; ++d) cudaFree(rows[d]);
            for (int v = 0; v < nxe; ++v) cudaFree(elem_rows[v]);
            cudaFree(aabbs);
            cudaFree(elements);
            cudaFree(idx);
        }
    };

    DeviceBoxes make(std::mt19937& rng, const ptrdiff_t n, const int nxe, const double spread, const double size) {
        std::uniform_real_distribution<double> pos(0.0, spread);
        std::uniform_real_distribution<double> ext(0.0, size);
        std::uniform_int_distribution<int> node(0, 999);

        DeviceBoxes b;
        b.n = n;
        b.nxe = nxe;
        for (int d = 0; d < 6; ++d) b.h[d].resize(n);
        for (int v = 0; v < nxe; ++v) b.h_elem[v].resize(n);
        for (ptrdiff_t i = 0; i < n; ++i) {
            for (int d = 0; d < 3; ++d) {
                const double lo = pos(rng);
                b.h[d][i] = lo;
                b.h[3 + d][i] = lo + ext(rng);
            }
            for (int v = 0; v < nxe; ++v) b.h_elem[v][i] = (idx_t)node(rng);
        }
        return b;
    }

    /**
     * \brief Boxes on a 2-manifold, which is where the sweep actually degenerates.
     *
     * Uniformly scattered boxes are the wrong benchmark: with small boxes in a
     * large volume the sort-axis window is short and sweep-and-prune does fine.
     * What breaks it is a surface -- N elements spread over a 2D manifold, so any
     * slab perpendicular to the sort axis holds O(sqrt(N)) of them while only a
     * handful overlap in 3D. That is the ~2,100:1 scan ratio measured on a refined
     * cloth-ball frame, and it is the case the cell list exists for.
     */
    DeviceBoxes make_sheet(std::mt19937& rng, const ptrdiff_t n, const double side) {
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        std::uniform_int_distribution<int> node(0, 999);

        // Spacing of n elements over a side x side sheet, and boxes a few spacings
        // across to stand in for swept AABBs.
        const double h = side / std::sqrt((double)n);
        const double box = 3.0 * h;

        DeviceBoxes b;
        b.n = n;
        b.nxe = 2;
        for (int d = 0; d < 6; ++d) b.h[d].resize(n);
        for (int v = 0; v < 2; ++v) b.h_elem[v].resize(n);

        for (ptrdiff_t i = 0; i < n; ++i) {
            const double u = uni(rng) * side;
            const double v = uni(rng) * side;
            // A gentle fold, so the sheet is not axis-aligned and the third axis
            // carries real but small extent.
            const double w = 0.15 * side * std::sin(3.0 * u / side) * std::cos(3.0 * v / side);
            const double p[3] = {u, v, w};
            for (int d = 0; d < 3; ++d) {
                b.h[d][i] = p[d];
                b.h[3 + d][i] = p[d] + box;
            }
            for (int q = 0; q < 2; ++q) b.h_elem[q][i] = (idx_t)node(rng);
        }
        return b;
    }

    /**
     * \brief Candidate statistics for a strategy heuristic, one pass, host side.
     *
     * lambda[d] = n * (mean extent along d) / (span along d) is the expected
     * number of boxes whose interval covers a given point on axis d -- that is,
     * the length of the window sweep-and-prune walks if it sorts on d. It is
     * computable from quantities the AABB pass already touches.
     */
    struct Stats {
        double lambda[3];
        double lam_min, lam_max;
    };

    Stats stats_of(const DeviceBoxes& b) {
        Stats st{};
        for (int d = 0; d < 3; ++d) {
            double lo = b.h[d][0], hi = b.h[3 + d][0], sum = 0;
            for (ptrdiff_t i = 0; i < b.n; ++i) {
                lo = std::min(lo, b.h[d][i]);
                hi = std::max(hi, b.h[3 + d][i]);
                sum += b.h[3 + d][i] - b.h[d][i];
            }
            const double span = std::max(hi - lo, 1e-300);
            st.lambda[d] = (double)b.n * (sum / (double)b.n) / span;
        }
        st.lam_min = std::min(st.lambda[0], std::min(st.lambda[1], st.lambda[2]));
        st.lam_max = std::max(st.lambda[0], std::max(st.lambda[1], st.lambda[2]));
        return st;
    }

    using PairSet = std::set<std::pair<idx_t, idx_t>>;

    PairSet download(idx_t* a, idx_t* b, const ptrdiff_t n) {
        PairSet out;
        if (n <= 0) return out;
        std::vector<idx_t> ha(n), hb(n);
        CHECK(cudaMemcpy(ha.data(), a, sizeof(idx_t) * n, cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(hb.data(), b, sizeof(idx_t) * n, cudaMemcpyDeviceToHost));
        for (ptrdiff_t i = 0; i < n; ++i) out.insert({ha[i], hb[i]});
        return out;
    }

    PairSet cell2d_self(DeviceBoxes& e) {
        sccd::device::Cell2DGridD<scalar_t> grid;
        ptrdiff_t spans = 0;
        ptrdiff_t* cellptr = nullptr;

        // ncells is not known until setup picks the axes, so size generously and
        // let setup memset what it needs; 4n is the cap setup itself enforces.
        CHECK(cudaMalloc(&cellptr, sizeof(ptrdiff_t) * (size_t)(4 * e.n + 2)));
        sccd::device::cell2d_setup_and_count<scalar_t, idx_t>(3, e.n, e.aabbs, grid, cellptr, &spans);

        idx_t* cellidx = nullptr;
        ptrdiff_t* cursor = nullptr;
        CHECK(cudaMalloc(&cellidx, sizeof(idx_t) * (size_t)(spans > 0 ? spans : 1)));
        CHECK(cudaMalloc(&cursor, sizeof(ptrdiff_t) * (size_t)(grid.ncells() + 1)));
        sccd::device::cell2d_fill<scalar_t, idx_t>(e.n, e.aabbs, grid, cellptr, cellidx, cursor);

        ptrdiff_t* ccdptr = nullptr;
        CHECK(cudaMalloc(&ccdptr, sizeof(ptrdiff_t) * (size_t)(e.n + 1)));
        sccd::device::cell2d_count_self_overlaps<2, scalar_t, idx_t>(
            e.n, e.aabbs, e.idx, 1, e.elements, grid, cellptr, cellidx, ccdptr);

        ptrdiff_t total = 0;
        CHECK(cudaMemcpy(&total, ccdptr + e.n, sizeof(ptrdiff_t), cudaMemcpyDeviceToHost));

        idx_t *o0 = nullptr, *o1 = nullptr;
        CHECK(cudaMalloc(&o0, sizeof(idx_t) * (size_t)(total > 0 ? total : 1)));
        CHECK(cudaMalloc(&o1, sizeof(idx_t) * (size_t)(total > 0 ? total : 1)));
        sccd::device::cell2d_collect_self_overlaps<2, scalar_t, idx_t>(
            e.n, e.aabbs, e.idx, 1, e.elements, grid, cellptr, cellidx, ccdptr, o0, o1);
        CHECK(cudaDeviceSynchronize());

        PairSet out = download(o0, o1, total);
        cudaFree(cellptr);
        cudaFree(cellidx);
        cudaFree(cursor);
        cudaFree(ccdptr);
        cudaFree(o0);
        cudaFree(o1);
        return out;
    }

    PairSet sweep_self(DeviceBoxes& e) {
        scalar_t* scratch = nullptr;
        CHECK(cudaMalloc(&scratch, sizeof(scalar_t) * (size_t)(2 * e.n)));
        const int axis = sccd::device::choose_axis<scalar_t>(3, e.n, e.aabbs);
        sccd::device::sort_along_axis<scalar_t, idx_t>(3, e.n, axis, e.aabbs, e.idx, scratch);

        ptrdiff_t* ccdptr = nullptr;
        CHECK(cudaMalloc(&ccdptr, sizeof(ptrdiff_t) * (size_t)(e.n + 1)));
        sccd::device::count_self_overlaps<2, scalar_t, idx_t>(axis, e.n, e.aabbs, e.idx, 1, e.elements, ccdptr);

        ptrdiff_t total = 0;
        CHECK(cudaMemcpy(&total, ccdptr + e.n, sizeof(ptrdiff_t), cudaMemcpyDeviceToHost));

        idx_t *o0 = nullptr, *o1 = nullptr;
        CHECK(cudaMalloc(&o0, sizeof(idx_t) * (size_t)(total > 0 ? total : 1)));
        CHECK(cudaMalloc(&o1, sizeof(idx_t) * (size_t)(total > 0 ? total : 1)));
        sccd::device::collect_self_overlaps<2, scalar_t, idx_t>(
            axis, e.n, e.aabbs, e.idx, 1, e.elements, ccdptr, o0, o1);
        CHECK(cudaDeviceSynchronize());

        PairSet out = download(o0, o1, total);
        cudaFree(scratch);
        cudaFree(ccdptr);
        cudaFree(o0);
        cudaFree(o1);
        return out;
    }

    int run(const char* name, const ptrdiff_t n, const double spread, const double size) {
        std::mt19937 rng(4242);
        DeviceBoxes a = make(rng, n, 2, spread, size);
        DeviceBoxes b = a;  // the sweep permutes its arrays, so give each its own
        a.upload();
        b.upload();

        const PairSet cell = cell2d_self(a);
        const PairSet sweep = sweep_self(b);

        const bool ok = (cell == sweep);
        std::printf("%-26s n=%-7ld sweep=%-8zu cell=%-8zu  %s\n",
                    name,
                    (long)n,
                    sweep.size(),
                    cell.size(),
                    ok ? "ok" : "MISMATCH");

        a.free_all();
        b.free_all();
        return ok ? 0 : 1;
    }

    // Timing includes the sort for the sweep and the binning for the cell list:
    // removing the sort is a stated reason to prefer the cell list, so charging
    // only the query would flatter it.
    double time_ms(const char* which, DeviceBoxes& e, const int repeats) {
        cudaEvent_t a, b;
        CHECK(cudaEventCreate(&a));
        CHECK(cudaEventCreate(&b));

        double best = 1e30;
        for (int r = 0; r < repeats; ++r) {
            DeviceBoxes copy = e;  // the sweep permutes in place
            copy.upload();
            CHECK(cudaDeviceSynchronize());
            CHECK(cudaEventRecord(a));
            if (std::string(which) == "cell2d") {
                (void)cell2d_self(copy);
            } else {
                (void)sweep_self(copy);
            }
            CHECK(cudaEventRecord(b));
            CHECK(cudaEventSynchronize(b));
            float ms = 0.f;
            CHECK(cudaEventElapsedTime(&ms, a, b));
            best = std::min(best, (double)ms);
            copy.free_all();
        }
        cudaEventDestroy(a);
        cudaEventDestroy(b);
        return best;
    }

}  // namespace

int main(int argc, char** argv) {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count == 0) {
        std::printf("no CUDA device available, skipping\n");
        return 0;
    }

    int bad = 0;
    bad |= run("small boxes", 4000, 100.0, 1.0);
    bad |= run("spanning many cells", 2500, 100.0, 25.0);
    bad |= run("dense overlap", 900, 1.0, 1.0);
    bad |= run("tiny spread", 700, 0.0001, 0.00001);

    std::printf("%s\n", bad ? "FAIL" : "OK: device cell list and device sweep agree");

    if (argc > 1) {
        std::mt19937 rng(7);
        // Sweep the box size at fixed n: this moves the input from "almost no
        // pairs" to "densely overlapping" without changing the shape of the
        // distribution, which is the only way to tell whether the discriminator is
        // the geometry's shape or simply how much work there is to do.
        std::printf("\nvolume, varying density at n = 400,000\n");
        std::printf("%-10s %10s %10s %8s %10s %9s %12s\n",
                    "box", "sweep_ms", "cell2d_ms", "speedup", "lam_min", "aniso", "pairs");
        for (double box : {1.5, 6.0, 25.0, 100.0}) {
            DeviceBoxes e = make(rng, 400000, 2, 1000.0, box);
            const Stats st = stats_of(e);
            DeviceBoxes probe = e;
            probe.upload();
            const size_t pairs = cell2d_self(probe).size();
            probe.free_all();
            const double s_ms = time_ms("sweep", e, 3);
            const double c_ms = time_ms("cell2d", e, 3);
            std::printf("%-10.1f %10.3f %10.3f %7.2fx %10.1f %9.2f %12zu\n",
                        box, s_ms, c_ms, s_ms / c_ms, st.lam_min, st.lam_max / st.lam_min, pairs);
        }

        std::printf("\nscattered boxes in a volume (the sweep's good case)\n");
        std::printf("%-12s %10s %10s %8s %10s %10s %9s\n",
                    "n", "sweep_ms", "cell2d_ms", "speedup", "lam_min", "lam_max", "aniso");
        for (ptrdiff_t n : {ptrdiff_t(100000), ptrdiff_t(400000), ptrdiff_t(1600000)}) {
            DeviceBoxes e = make(rng, n, 2, 1000.0, 1.5);
            const Stats st = stats_of(e);
            const double s_ms = time_ms("sweep", e, 3);
            const double c_ms = time_ms("cell2d", e, 3);
            std::printf("%-12ld %10.3f %10.3f %7.2fx %10.1f %10.1f %9.2f\n",
                        (long)n, s_ms, c_ms, s_ms / c_ms, st.lam_min, st.lam_max, st.lam_max / st.lam_min);
        }

        std::printf("\nboxes on a surface (what a cloth looks like)\n");
        std::printf("%-12s %10s %10s %8s %10s %10s %9s\n",
                    "n", "sweep_ms", "cell2d_ms", "speedup", "lam_min", "lam_max", "aniso");
        for (ptrdiff_t n : {ptrdiff_t(100000), ptrdiff_t(400000), ptrdiff_t(1600000)}) {
            DeviceBoxes e = make_sheet(rng, n, 1000.0);
            const Stats st = stats_of(e);
            const double s_ms = time_ms("sweep", e, 3);
            const double c_ms = time_ms("cell2d", e, 3);
            std::printf("%-12ld %10.3f %10.3f %7.2fx %10.1f %10.1f %9.2f\n",
                        (long)n, s_ms, c_ms, s_ms / c_ms, st.lam_min, st.lam_max, st.lam_max / st.lam_min);
        }
    }

    return bad;
}
