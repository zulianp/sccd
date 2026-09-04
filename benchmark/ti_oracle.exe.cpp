// Phase 0 oracle: compares every narrow-phase mode against TightInclusion on
// the CCD benchmark query sets.
//
// This is the acceptance gate for the TI-accurate vectorized kernel. It answers
// three questions per (dataset, phase, mode):
//   * do we agree with TI on hit/miss, and in which direction do we disagree
//     (a false negative is a missed collision; a false positive is conservative)
//   * how far off is the reported time of impact
//   * how much work did it take
//
// Deliberately free of smesh: it reads the raw rational query CSVs directly, so
// it builds with nothing but SCCD + TightInclusion.
//
//   ti_oracle <dataset-dir> [--phase vf|ee|both] [--max-files N]
//             [--tol T] [--max-depth N] [--csv out.csv]
//
// Built with SCCD_ENABLE_CUDA it also runs the device narrow phase as a fourth
// mode, over the same queries and through the same gate, so the GPU is held to
// the invariant the CPU modes are. Without CUDA that mode is absent and nothing
// else changes.

#include "sccd_narrowphase.hpp"
#include "sccd_math.hpp"
#include "sccd_rootfinder.hpp"

#ifdef SCCD_ENABLE_CUDA
#include "sccd_narrowphase.cuh"

#include <cuda_runtime.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

using scalar_t = double;
using idx_t = std::int32_t;

namespace {

    struct Options {
        fs::path dataset_dir;
        bool do_vf = true;
        bool do_ee = true;
        std::size_t max_files = 0;  // 0 == all
        scalar_t tol = 3e-8;
        int max_depth = 96;
        fs::path csv;
        fs::path violations_csv;
        bool strict = true;
        // Which modes count toward the exit code. Modes 0 and 1 are known not to
        // be conservative, so gating CI on mode 2 alone is the useful default
        // once it is the mode being shipped.
        std::string gate = "all";
        // Run the device narrow phase in single precision. The device kernel is
        // instantiated for both, and the two differ in more than speed: the
        // origin-containment test is padded by the user tolerance rather than by
        // TightInclusion's numerical error bound, and which of the two is larger
        // flips between float and double.
        bool device_float = false;
        // Round the loaded geometry through float once, for *every* mode
        // including the TightInclusion reference. The exact roots in the dataset
        // belong to the original rational geometry, so narrowing the input moves
        // the true root; this makes that effect visible on the host rows too,
        // which is the only way to tell it apart from a kernel defect.
        bool float_geometry = false;
        // > 0 selects throughput mode: merge every query file into one batch and
        // time the device kernels with CUDA events instead of scoring accuracy.
        int bench_repeats = 0;
    };

    // One file's worth of queries, as the flat SoA the narrow phase expects.
    struct QuerySet {
        std::vector<scalar_t> p0[3];
        std::vector<scalar_t> p1[3];
        std::vector<idx_t> q0;
        std::vector<idx_t> q1;
        std::vector<idx_t> prim[3];  // faces (3) or edges (2)
        std::size_t n_queries = 0;

        scalar_t* p0_ptr[3];
        scalar_t* p1_ptr[3];
        idx_t* prim_ptr[3];

        void bind() {
            for (int d = 0; d < 3; ++d) {
                p0_ptr[d] = p0[d].data();
                p1_ptr[d] = p1[d].data();
                prim_ptr[d] = prim[d].empty() ? nullptr : prim[d].data();
            }
        }
    };

#ifdef SCCD_ENABLE_CUDA
#define ORACLE_CUDA_CHECK(call)                                                                   \
    do {                                                                                          \
        const cudaError_t _err = (call);                                                          \
        if (_err != cudaSuccess) {                                                                \
            std::fprintf(stderr,                                                                  \
                         "ti_oracle: CUDA error at %s:%d: %s\n",                                  \
                         __FILE__,                                                                \
                         __LINE__,                                                                \
                         cudaGetErrorString(_err));                                               \
            /* 3, not 2: 2 is the conservativeness gate's failure code, and CI has */              \
            /* to be able to tell "the GPU broke the invariant" from "the GPU was  */              \
            /* not there". */                                                                      \
            std::exit(3);                                                                          \
        }                                                                                         \
    } while (0)

    template <typename T>
    T* device_dup(const T* host, const std::size_t n) {
        if (n == 0) return nullptr;
        T* dev = nullptr;
        ORACLE_CUDA_CHECK(cudaMalloc(&dev, n * sizeof(T)));
        ORACLE_CUDA_CHECK(cudaMemcpy(dev, host, n * sizeof(T), cudaMemcpyHostToDevice));
        return dev;
    }

    /**
     * \brief A QuerySet mirrored into device memory, in the layout the kernels want.
     *
     * The device narrow phase takes the coordinate and connectivity arrays as
     * pointer-to-pointer and dereferences the outer level *on the device*
     * (sp[0][idx] inside load_query_vf), so the array of three pointers has to
     * live in device memory too, not just the arrays it points at. Passing a host
     * array of device pointers compiles and then faults, which is why this is
     * spelled out rather than left to the caller.
     */
    template <typename T>
    struct DeviceQuerySet {
        T* p0[3] = {nullptr, nullptr, nullptr};
        T* p1[3] = {nullptr, nullptr, nullptr};
        idx_t* prim[3] = {nullptr, nullptr, nullptr};
        T** d_p0 = nullptr;
        T** d_p1 = nullptr;
        idx_t** d_prim = nullptr;
        idx_t* d_q0 = nullptr;
        idx_t* d_q1 = nullptr;
        T* d_toi = nullptr;
        std::size_t n_queries = 0;

        explicit DeviceQuerySet(const QuerySet& qs) : n_queries(qs.n_queries) {
            // Narrowed here when T is float. That rounding is part of what a float
            // pipeline does, so a violation found this way is a property of the
            // pipeline, not an artefact of the harness -- but it does mean a float
            // run cannot separate "the kernel rejected a box it should have kept"
            // from "the rounded geometry genuinely moved the root".
            std::vector<T> tmp;
            auto upload = [&](const std::vector<scalar_t>& src) {
                tmp.assign(src.begin(), src.end());
                return device_dup(tmp.data(), tmp.size());
            };
            for (int d = 0; d < 3; ++d) {
                p0[d] = upload(qs.p0[d]);
                p1[d] = upload(qs.p1[d]);
                // Edge-edge leaves the third slot empty; it is never read, but the
                // outer array still has to be three wide.
                prim[d] = qs.prim[d].empty() ? nullptr : device_dup(qs.prim[d].data(), qs.prim[d].size());
            }
            d_p0 = device_dup(p0, 3);
            d_p1 = device_dup(p1, 3);
            d_prim = device_dup(prim, 3);
            d_q0 = device_dup(qs.q0.data(), qs.q0.size());
            d_q1 = device_dup(qs.q1.data(), qs.q1.size());
            ORACLE_CUDA_CHECK(cudaMalloc(&d_toi, sccd::max<std::size_t>(n_queries, 1) * sizeof(T)));
        }

        ~DeviceQuerySet() {
            for (int d = 0; d < 3; ++d) {
                cudaFree(p0[d]);
                cudaFree(p1[d]);
                cudaFree(prim[d]);
            }
            cudaFree(d_p0);
            cudaFree(d_p1);
            cudaFree(d_prim);
            cudaFree(d_q0);
            cudaFree(d_q1);
            cudaFree(d_toi);
        }

        DeviceQuerySet(const DeviceQuerySet&) = delete;
        DeviceQuerySet& operator=(const DeviceQuerySet&) = delete;

        void download(std::vector<scalar_t>& out) const {
            if (n_queries == 0) return;
            std::vector<T> tmp(n_queries);
            ORACLE_CUDA_CHECK(cudaMemcpy(tmp.data(), d_toi, n_queries * sizeof(T), cudaMemcpyDeviceToHost));
            // Widening a float toi is exact, so scoring stays in double throughout.
            out.assign(tmp.begin(), tmp.end());
        }
    };

    /** \brief Run the device narrow phase over one query set, in T precision. */
    template <typename T>
    void run_device_narrow_phase(const QuerySet& qs,
                                 const bool is_vf,
                                 const int max_depth,
                                 const scalar_t tol,
                                 std::vector<scalar_t>& toi_out) {
        if (qs.n_queries == 0) return;
        DeviceQuerySet<T> dev(qs);
        if (is_vf) {
            sccd::device::narrow_phase_vf<3, T, idx_t>(qs.n_queries,
                                                              dev.d_q0,
                                                              dev.d_q1,
                                                              dev.d_p0,
                                                              dev.d_p1,
                                                              /*face_stride=*/1,
                                                              dev.d_prim,
                                                              T(1),
                                                              dev.d_toi,
                                                              max_depth,
                                                              T(tol),
                                                              /*toi_stride=*/1);
        } else {
            sccd::device::narrow_phase_ee<T, idx_t>(qs.n_queries,
                                                           dev.d_q0,
                                                           dev.d_q1,
                                                           dev.d_p0,
                                                           dev.d_p1,
                                                           /*edge_stride=*/1,
                                                           dev.d_prim,
                                                           T(1),
                                                           dev.d_toi,
                                                           max_depth,
                                                           T(tol),
                                                           /*toi_stride=*/1);
        }
        ORACLE_CUDA_CHECK(cudaDeviceSynchronize());
        dev.download(toi_out);
    }
#endif  // SCCD_ENABLE_CUDA

    /**
     * \brief Concatenate one query set onto another, re-basing its indices.
     *
     * Each file owns four consecutive nodes per query, so merging means shifting
     * node ids by the running node count and primitive ids by the running
     * primitive count. Needed because a per-file batch is a few hundred queries,
     * where launch overhead dominates and the measurement says nothing about the
     * kernel -- which is what the "Caveat on the ms column" in
     * benchmark/oracle/README.md is about.
     */
    void append_query_set(QuerySet& dst, const QuerySet& src, const bool is_vf) {
        const idx_t node_offset = static_cast<idx_t>(dst.p0[0].size());
        const idx_t prim_offset = static_cast<idx_t>(dst.prim[0].size());

        for (int d = 0; d < 3; ++d) {
            dst.p0[d].insert(dst.p0[d].end(), src.p0[d].begin(), src.p0[d].end());
            dst.p1[d].insert(dst.p1[d].end(), src.p1[d].begin(), src.p1[d].end());
            for (const idx_t v : src.prim[d]) {
                dst.prim[d].push_back(v + node_offset);
            }
        }
        for (const idx_t v : src.q0) {
            dst.q0.push_back(v + (is_vf ? node_offset : prim_offset));
        }
        for (const idx_t v : src.q1) {
            dst.q1.push_back(v + prim_offset);
        }
        dst.n_queries += src.n_queries;
    }

#ifdef SCCD_ENABLE_CUDA
    /**
     * \brief Time one device mode on a single large batch, uploading once.
     *
     * CUDA events rather than a host clock, and the transfer is outside the timed
     * region: the question is what the kernel costs, not what PCIe costs.
     */
    template <typename T>
    double bench_device(const QuerySet& qs,
                        const bool is_vf,
                        const int max_depth,
                        const scalar_t tol,
                        const int repeats,
                        const int toi_stride) {
        DeviceQuerySet<T> dev(qs);
        cudaEvent_t beg, end;
        ORACLE_CUDA_CHECK(cudaEventCreate(&beg));
        ORACLE_CUDA_CHECK(cudaEventCreate(&end));

        auto once = [&]() {
            if (is_vf) {
                sccd::device::narrow_phase_vf<3, T, idx_t>(qs.n_queries, dev.d_q0, dev.d_q1, dev.d_p0, dev.d_p1, 1,
                                                           dev.d_prim, T(1), dev.d_toi, max_depth, T(tol),
                                                           toi_stride);
            } else {
                sccd::device::narrow_phase_ee<T, idx_t>(qs.n_queries, dev.d_q0, dev.d_q1, dev.d_p0, dev.d_p1, 1,
                                                        dev.d_prim, T(1), dev.d_toi, max_depth, T(tol),
                                                        toi_stride);
            }
        };

        once();  // warm-up: first call also sizes the persistent global stack
        ORACLE_CUDA_CHECK(cudaDeviceSynchronize());

        double best = 1e30;
        for (int r = 0; r < repeats; ++r) {
            ORACLE_CUDA_CHECK(cudaEventRecord(beg));
            once();
            ORACLE_CUDA_CHECK(cudaEventRecord(end));
            ORACLE_CUDA_CHECK(cudaEventSynchronize(end));
            float ms = 0.f;
            ORACLE_CUDA_CHECK(cudaEventElapsedTime(&ms, beg, end));
            best = std::min(best, static_cast<double>(ms));
        }
        cudaEventDestroy(beg);
        cudaEventDestroy(end);
        return best;
    }
#endif

    // Rows are "num,den,num,den,num,den": exact rationals from the benchmark.
    bool parse_query_row(const std::string& line, scalar_t coords[3]) {
        double values[6];
        std::size_t cursor = 0;
        for (int i = 0; i < 6; ++i) {
            const std::size_t end = line.find(',', cursor);
            if (i < 5 && end == std::string::npos) {
                return false;
            }
            try {
                values[i] = std::stod(line.substr(cursor, end - cursor));
            } catch (...) {
                return false;
            }
            cursor = (i < 5) ? end + 1 : end;
        }
        coords[0] = static_cast<scalar_t>(values[0] / values[1]);
        coords[1] = static_cast<scalar_t>(values[2] / values[3]);
        coords[2] = static_cast<scalar_t>(values[4] / values[5]);
        return true;
    }

    bool read_queries(const fs::path& path,
                      const bool is_vf,
                      const bool narrow_geometry_to_float,
                      QuerySet& out) {
        std::ifstream in(path);
        if (!in) {
            return false;
        }

        std::string line;
        std::size_t row = 0;
        while (std::getline(in, line)) {
            if (line.empty()) {
                continue;
            }
            scalar_t coords[3];
            if (!parse_query_row(line, coords)) {
                std::cerr << "error: bad row in " << path << "\n";
                return false;
            }
            auto& points = (row % 8 < 4) ? out.p0 : out.p1;
            for (int d = 0; d < 3; ++d) {
                points[d].push_back(coords[d]);
            }
            ++row;
        }

        if (row == 0 || row % 8 != 0) {
            return false;
        }

        out.n_queries = row / 8;
        out.q0.resize(out.n_queries);
        out.q1.resize(out.n_queries);

        // Same index convention as benchmark/bench.exe.cpp: each query owns four
        // consecutive nodes, laid out (vertex, f0, f1, f2) or (ea0, ea1, eb0, eb1).
        if (is_vf) {
            for (std::size_t i = 0; i < out.n_queries; ++i) {
                const idx_t base = static_cast<idx_t>(4 * i);
                out.q0[i] = base;                        // vertex node
                out.q1[i] = static_cast<idx_t>(i);       // face id
                out.prim[0].push_back(base + 1);
                out.prim[1].push_back(base + 2);
                out.prim[2].push_back(base + 3);
            }
        } else {
            for (std::size_t i = 0; i < out.n_queries; ++i) {
                const idx_t base = static_cast<idx_t>(4 * i);
                out.q0[i] = static_cast<idx_t>(2 * i);
                out.q1[i] = static_cast<idx_t>(2 * i + 1);
                out.prim[0].push_back(base);
                out.prim[1].push_back(base + 1);
                out.prim[0].push_back(base + 2);
                out.prim[1].push_back(base + 3);
            }
        }

        if (narrow_geometry_to_float) {
            for (int d = 0; d < 3; ++d) {
                for (auto& v : out.p0[d]) v = static_cast<scalar_t>(static_cast<float>(v));
                for (auto& v : out.p1[d]) v = static_cast<scalar_t>(static_cast<float>(v));
            }
        }

        out.bind();
        return true;
    }

    /**
     * \brief Exact time of impact per query, from the dataset's own roots.
     *
     * NaN marks a query with no collision. These come from the benchmark's
     * symbolic root finder, so unlike TightInclusion's output they are the true
     * value rather than a certified lower bound on it -- which is what the
     * conservativeness invariant is actually stated against.
     */
    bool read_ground_truth_toi(const fs::path& dataset_dir, const std::string& key,
                               std::vector<double>& out) {
        const fs::path path = dataset_dir / "roots" / key / "toi.float64";
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        if (!in) {
            return false;
        }
        const auto bytes = in.tellg();
        if (bytes <= 0 || (static_cast<std::size_t>(bytes) % sizeof(double)) != 0) {
            return false;
        }
        out.resize(static_cast<std::size_t>(bytes) / sizeof(double));
        in.seekg(0);
        in.read(reinterpret_cast<char*>(out.data()), bytes);
        return static_cast<bool>(in);
    }

    // Optional Mathematica-verified ground truth, one byte per query.
    bool read_ground_truth(const fs::path& dataset_dir, const std::string& key, std::vector<std::uint8_t>& out) {
        const fs::path path = dataset_dir / "mma_bool" / key / "mma_bool.uint8";
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        if (!in) {
            return false;
        }
        const auto bytes = in.tellg();
        if (bytes <= 0) {
            return false;
        }
        out.resize(static_cast<std::size_t>(bytes));
        in.seekg(0);
        in.read(reinterpret_cast<char*>(out.data()), bytes);
        return static_cast<bool>(in);
    }

    // Relative error is meaningless when TI reports toi == 0 (already touching
    // at the start of the step), which happens often enough to dominate a naive
    // average. Relative error is therefore only accumulated for references above
    // REL_ERR_FLOOR; everything else is covered by the absolute error.
    constexpr double REL_ERR_FLOOR = 1e-9;

    // Enough to diagnose without unbounded memory on a bad run.
    constexpr std::size_t MAX_RECORDED_VIOLATIONS = 64;

    // A violation of the conservativeness invariant: either a missed collision
    // or a time of impact reported after the reference's. Both let a simulation
    // step through a contact, so both are failures rather than inaccuracies.
    struct Violation {
        std::string file;
        std::size_t query = 0;
        const char* kind = "";
        double ref_toi = 0;
        double got_toi = 0;
    };

    struct Stats {
        std::size_t n = 0;
        std::size_t hits = 0;
        std::size_t false_negative = 0;  // TI hit, mode missed -- UNSAFE
        std::size_t false_positive = 0;  // mode hit, TI missed -- conservative
        std::size_t late = 0;            // toi later than TI's reported value
        std::size_t gt_late = 0;         // toi later than the TRUE toi -- UNSAFE
        std::size_t gt_missed = 0;       // exact root exists, mode found none -- UNSAFE
        std::size_t gt_checked = 0;      // queries with an exact root available
        double gt_worst_overshoot = 0;   // largest amount by which we ran late
        std::size_t near_zero_ref = 0;   // TI toi below REL_ERR_FLOOR
        std::size_t gt_false_negative = 0;
        std::vector<double> rel_err;
        std::vector<double> abs_err;
        std::vector<Violation> violations;
        double seconds = 0;

        void add(const double reference, const double value) {
            const double err = value - reference;
            abs_err.push_back(std::abs(err));
            // A time of impact reported *after* the true one lets a solver step
            // through the collision; earlier is merely conservative.
            if (err > REL_ERR_FLOOR) {
                late += 1;
            }
            if (std::abs(reference) >= REL_ERR_FLOOR) {
                rel_err.push_back(std::abs(err) / std::abs(reference));
            } else {
                near_zero_ref += 1;
            }
        }

        static double quantile_of(std::vector<double>& v, const double q) {
            if (v.empty()) return 0.0;
            std::sort(v.begin(), v.end());
            const std::size_t i = static_cast<std::size_t>(q * (v.size() - 1) + 0.5);
            return v[std::min(i, v.size() - 1)];
        }
        static double mean_of(const std::vector<double>& v) {
            if (v.empty()) return 0.0;
            double s = 0;
            for (double x : v) s += x;
            return s / v.size();
        }
        static double max_of(const std::vector<double>& v) {
            double m = 0;
            for (double x : v) m = std::max(m, x);
            return m;
        }
    };

    double now_seconds() {
        using namespace std::chrono;
        return duration<double>(steady_clock::now().time_since_epoch()).count();
    }

    // The modes under comparison. TI itself is the reference, not a mode.
    //
    // Device is the CUDA narrow phase. It is deliberately scored by the same
    // gate as the host modes rather than compared against them: the invariant is
    // a property of the answer, not of which processor produced it. Note that it
    // does not honour SCCD_NARROWPHASE_MODE at all -- the device has one kernel
    // -- so unlike the host rows its name describes the hardware, not the
    // algorithm.
    // The names are the shipped mode names: SCCD_NARROWPHASE_MODE 0 is Fast and 2
    // is Tight. There used to be a third host row, "vector", pinned to mode 1;
    // mode 1 no longer exists, so that row was measuring Fast a second time under
    // a different name.
    enum class Mode { Fast, Tight, DeviceFast, DeviceTight };
#ifdef SCCD_ENABLE_CUDA
    constexpr int N_MODES = 4;
#else
    constexpr int N_MODES = 2;
#endif

    const char* mode_name(const Mode m) {
        switch (m) {
            case Mode::Fast: return "fast";
            case Mode::Tight: return "tight";
            case Mode::DeviceFast: return "device-fast";
            case Mode::DeviceTight: return "device-tight";
        }
        return "?";
    }

    Mode mode_of(const int m) { return static_cast<Mode>(m); }

    /**
     * \brief Point the next narrow-phase call at one mode.
     *
     * SCCD_NARROWPHASE_MODE is set for every mode rather than only the host ones.
     * It takes precedence over the legacy SCCD_USE_VNARROW_PHASE, so leaving it
     * unset for some rows would let a previous row's value leak into this one --
     * and now that the device honours the mode too, that leak would silently
     * decide which GPU kernel ran.
     */
    void select_mode(const Mode m) {
        const char* value = "0";
        switch (m) {
            case Mode::Fast: value = "0"; break;
            case Mode::Tight: value = "2"; break;
            case Mode::DeviceFast: value = "0"; break;
            case Mode::DeviceTight: value = "2"; break;
        }
        setenv("SCCD_NARROWPHASE_MODE", value, 1);
        setenv("SCCD_USE_VNARROW_PHASE", value, 1);
        unsetenv("SCCD_VNARROWPHASE_TI_COMPAT");
        unsetenv("SCCD_USE_TI");
    }

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if (a == "--phase") {
            const std::string p = next();
            opt.do_vf = (p == "vf" || p == "both");
            opt.do_ee = (p == "ee" || p == "both");
        } else if (a == "--max-files") {
            opt.max_files = static_cast<std::size_t>(std::stoul(next()));
        } else if (a == "--tol") {
            opt.tol = std::stod(next());
        } else if (a == "--max-depth") {
            opt.max_depth = std::stoi(next());
        } else if (a == "--csv") {
            opt.csv = next();
        } else if (a == "--violations-csv") {
            opt.violations_csv = next();
        } else if (a == "--no-strict") {
            opt.strict = false;
        } else if (a == "--device-float") {
            opt.device_float = true;
        } else if (a == "--float-geometry") {
            opt.float_geometry = true;
        } else if (a == "--bench") {
            const std::string v = next();
            opt.bench_repeats = v.empty() ? 3 : std::max(1, atoi(v.c_str()));
        } else if (a == "--gate") {
            opt.gate = next();
        } else if (!a.empty() && a[0] != '-') {
            opt.dataset_dir = a;
        }
    }

    const bool geometry_is_narrowed = opt.float_geometry || opt.device_float;

    if (opt.dataset_dir.empty()) {
        std::cerr << "usage: ti_oracle <dataset-dir> [--phase vf|ee|both] [--max-files N]\n"
                     "                 [--tol T] [--max-depth N] [--csv out.csv]\n"
                     "                 [--violations-csv out.csv] [--no-strict] [--gate MODE]\n"
#ifdef SCCD_ENABLE_CUDA
                     "                 [--device-float] [--float-geometry]\n"
#endif
                     "\n"
                     "Exits non-zero when a mode misses a collision or reports a time of\n"
                     "impact later than the dataset's exact root; both break\n"
                     "conservativeness. Pass --no-strict to report without failing.\n"
                     "\n"
                     "--gate MODE restricts the exit code to one mode: fast, tight"
#ifdef SCCD_ENABLE_CUDA
                     ", device-fast, device-tight"
#endif
                     ", or all (the default).\n"
#ifdef SCCD_ENABLE_CUDA
                     "\n"
                     "Built with CUDA: the 'device-*' rows are the GPU narrow phase, run\n"
                     "the same queries and held to the same gate. --device-float runs it\n"
                     "in single precision instead of double.\n"
#endif
            ;
        return 1;
    }

    const fs::path queries_dir = opt.dataset_dir / "queries";
    if (!fs::is_directory(queries_dir)) {
        std::cerr << "error: no queries/ under " << opt.dataset_dir << "\n";
        return 1;
    }

    const std::string dataset = opt.dataset_dir.filename().string();

    struct Phase {
        bool is_vf;
        const char* suffix;
        const char* name;
    };
    const Phase phases[2] = {{true, "vf.csv", "VF"}, {false, "ee.csv", "EE"}};

    std::vector<std::string> csv_rows;
    std::vector<std::string> violation_rows;
    const fs::path& violation_csv = opt.violations_csv;
    std::size_t violations_total = 0;

    for (const Phase& phase : phases) {
        if ((phase.is_vf && !opt.do_vf) || (!phase.is_vf && !opt.do_ee)) {
            continue;
        }

        std::vector<fs::path> files;
        for (const auto& entry : fs::directory_iterator(queries_dir)) {
            const std::string name = entry.path().filename().string();
            if (name.size() > 6 && name.rfind(phase.suffix) == name.size() - std::strlen(phase.suffix)) {
                files.push_back(entry.path());
            }
        }
        std::sort(files.begin(), files.end());
        if (opt.max_files && files.size() > opt.max_files) {
            files.resize(opt.max_files);
        }
        if (files.empty()) {
            continue;
        }

        if (opt.bench_repeats > 0) {
            QuerySet batch;
            for (const fs::path& file : files) {
                QuerySet qs;
                if (!read_queries(file, phase.is_vf, opt.float_geometry, qs)) continue;
                append_query_set(batch, qs, phase.is_vf);
            }
            batch.bind();
            if (batch.n_queries == 0) continue;

            const char* omp = getenv("OMP_NUM_THREADS");
            std::printf("\n%s / %s  throughput: %zu queries in one batch, %s, best of %d"
                        "  [host threads: %s]\n",
                        opt.dataset_dir.filename().string().c_str(),
                        phase.name,
                        batch.n_queries,
                        opt.device_float ? "float storage" : "double storage",
                        opt.bench_repeats,
                        omp ? omp : "unset");
            // stride 1 writes a time of impact per query; stride 0 writes one
            // global minimum. They are different questions and different kernels:
            // stride 1 gives a block to each query, stride 0 a thread, and stride
            // 0's single shared toi lets every query prune against every other
            // query's progress. Timed side by side because production callers
            // choose between them.
            std::printf("%-12s %12s %12s %10s\n", "mode", "stride1_ms", "stride0_ms", "ratio");
            for (int m = 0; m < N_MODES; ++m) {
                const Mode mode = mode_of(m);
                const bool is_device = (mode == Mode::DeviceFast || mode == Mode::DeviceTight);
                // Every mode, host and device, at both strides. The host rows are
                // what make a device row readable: `tight` and `device-tight` run
                // the identical search, and `fast` is the host counterpart of
                // `device-fast`. Without a like-for-like pair there is
                // no way to tell "this kernel is inefficient" from "this search is
                // expensive on this geometry".
                select_mode(mode);

                double ms = 0.0;
                (void)ms;
#ifdef SCCD_ENABLE_CUDA
                if (is_device) {
                    double ms_s[2] = {0.0, 0.0};
                    for (int k = 0; k < 2; ++k) {
                        const int stride = (k == 0) ? 1 : 0;
                        ms_s[k] = opt.device_float
                                      ? bench_device<float>(batch, phase.is_vf, opt.max_depth, opt.tol,
                                                            opt.bench_repeats, stride)
                                      : bench_device<double>(batch, phase.is_vf, opt.max_depth, opt.tol,
                                                             opt.bench_repeats, stride);
                    }
                    std::printf("%-12s %12.3f %12.3f %9.2fx\n", mode_name(mode), ms_s[0], ms_s[1],
                                ms_s[1] > 0 ? ms_s[0] / ms_s[1] : 0.0);
                    continue;
                }
#else
                if (is_device) continue;  // no device rows without CUDA
#endif
                {
                    double ms_s[2] = {0.0, 0.0};
                    for (int k = 0; k < 2; ++k) {
                        const int stride = (k == 0) ? 1 : 0;
                        std::vector<scalar_t> toi(stride == 0 ? 1 : batch.n_queries, 1.0);
                        double best = 1e30;
                        for (int r = 0; r < opt.bench_repeats; ++r) {
                            std::fill(toi.begin(), toi.end(), scalar_t(1));
                            const double t0 = now_seconds();
                            if (phase.is_vf) {
                                sccd::narrow_phase_vf<3, scalar_t, idx_t>(
                                    batch.n_queries, batch.q0.data(), batch.q1.data(), batch.p0_ptr, batch.p1_ptr, 1,
                                    batch.prim_ptr, scalar_t(1), toi.data(), opt.max_depth, opt.tol, stride);
                            } else {
                                sccd::narrow_phase_ee<scalar_t, idx_t>(
                                    batch.n_queries, batch.q0.data(), batch.q1.data(), batch.p0_ptr, batch.p1_ptr, 1,
                                    batch.prim_ptr, scalar_t(1), toi.data(), opt.max_depth, opt.tol, stride);
                            }
                            best = std::min(best, (now_seconds() - t0) * 1e3);
                        }
                        ms_s[k] = best;
                    }
                    std::printf("%-12s %12.3f %12.3f %9.2fx\n", mode_name(mode), ms_s[0], ms_s[1],
                                ms_s[1] > 0 ? ms_s[0] / ms_s[1] : 0.0);
                    continue;
                }
            }
            continue;
        }

        Stats stats[N_MODES];
        double ti_seconds = 0;
        std::size_t ti_hits = 0;
        std::size_t total_queries = 0;
        std::size_t gt_available = 0;

        for (const fs::path& file : files) {
            QuerySet qs;
            if (!read_queries(file, phase.is_vf, opt.float_geometry, qs)) {
                std::cerr << "warn: skipping " << file << "\n";
                continue;
            }

            const std::string key = file.stem().string();
            std::vector<std::uint8_t> gt;
            const bool have_gt = read_ground_truth(opt.dataset_dir, key, gt) && gt.size() >= qs.n_queries;

            // The exact roots are the reference the invariant is stated against.
            std::vector<double> gt_toi;
            const bool have_gt_toi =
                read_ground_truth_toi(opt.dataset_dir, key, gt_toi) && gt_toi.size() >= qs.n_queries;

            // --- reference: TightInclusion, one query at a time ---
            std::vector<double> ti_toi(qs.n_queries, 1.0);
            std::vector<std::uint8_t> ti_hit(qs.n_queries, 0);
            const double t_ti0 = now_seconds();
            for (std::size_t i = 0; i < qs.n_queries; ++i) {
                const idx_t b = static_cast<idx_t>(4 * i);
                auto P0 = [&](int k, int d) { return qs.p0[d][b + k]; };
                auto P1 = [&](int k, int d) { return qs.p1[d][b + k]; };
                const double a0[3] = {P0(0, 0), P0(0, 1), P0(0, 2)};
                const double a1[3] = {P0(1, 0), P0(1, 1), P0(1, 2)};
                const double a2[3] = {P0(2, 0), P0(2, 1), P0(2, 2)};
                const double a3[3] = {P0(3, 0), P0(3, 1), P0(3, 2)};
                const double b0[3] = {P1(0, 0), P1(0, 1), P1(0, 2)};
                const double b1[3] = {P1(1, 0), P1(1, 1), P1(1, 2)};
                const double b2[3] = {P1(2, 0), P1(2, 1), P1(2, 2)};
                const double b3[3] = {P1(3, 0), P1(3, 1), P1(3, 2)};

                double t = 1, u = 0, v = 0;
                const bool hit = phase.is_vf
                                     ? sccd::find_root_tight_inclusion_vf<double>(
                                           opt.max_depth, opt.tol, a0, a1, a2, a3, b0, b1, b2, b3, t, u, v)
                                     : sccd::find_root_tight_inclusion_ee<double>(
                                           opt.max_depth, opt.tol, a0, a1, a2, a3, b0, b1, b2, b3, t, u, v);
                ti_hit[i] = hit ? 1 : 0;
                ti_toi[i] = hit ? t : 1.0;
                ti_hits += hit ? 1 : 0;
            }
            ti_seconds += now_seconds() - t_ti0;

            if (have_gt) {
                gt_available += qs.n_queries;
            }
            total_queries += qs.n_queries;

            // --- each mode, batched exactly as production calls it ---
            for (int m = 0; m < N_MODES; ++m) {
                select_mode(mode_of(m));

                std::vector<scalar_t> toi(qs.n_queries, 1.0);
                const double t0 = now_seconds();
#ifdef SCCD_ENABLE_CUDA
                if (mode_of(m) == Mode::DeviceFast || mode_of(m) == Mode::DeviceTight) {
                    if (opt.device_float) {
                        run_device_narrow_phase<float>(qs, phase.is_vf, opt.max_depth, opt.tol, toi);
                    } else {
                        run_device_narrow_phase<double>(qs, phase.is_vf, opt.max_depth, opt.tol, toi);
                    }
                } else
#endif
                if (phase.is_vf) {
                    sccd::narrow_phase_vf<3, scalar_t, idx_t>(qs.n_queries,
                                                              qs.q0.data(),
                                                              qs.q1.data(),
                                                              qs.p0_ptr,
                                                              qs.p1_ptr,
                                                              1,
                                                              qs.prim_ptr,
                                                              scalar_t(1),
                                                              toi.data(),
                                                              opt.max_depth,
                                                              opt.tol,
                                                              1);
                } else {
                    sccd::narrow_phase_ee<scalar_t, idx_t>(qs.n_queries,
                                                           qs.q0.data(),
                                                           qs.q1.data(),
                                                           qs.p0_ptr,
                                                           qs.p1_ptr,
                                                           1,
                                                           qs.prim_ptr,
                                                           scalar_t(1),
                                                           toi.data(),
                                                           opt.max_depth,
                                                           opt.tol,
                                                           1);
                }
                stats[m].seconds += now_seconds() - t0;

                const std::size_t late_before = stats[m].late;
                for (std::size_t i = 0; i < qs.n_queries; ++i) {
                    const bool hit = toi[i] < scalar_t(1);
                    stats[m].n += 1;
                    stats[m].hits += hit ? 1 : 0;
                    if (ti_hit[i] && !hit) {
                        stats[m].false_negative += 1;
                        if (stats[m].violations.size() < MAX_RECORDED_VIOLATIONS) {
                            stats[m].violations.push_back({key, i, "missed", ti_toi[i], double(toi[i])});
                        }
                    }
                    if (!ti_hit[i] && hit) stats[m].false_positive += 1;
                    if (ti_hit[i] && hit) {
                        stats[m].add(ti_toi[i], toi[i]);
                        if (stats[m].late != late_before
                            && stats[m].violations.size() < MAX_RECORDED_VIOLATIONS
                            && double(toi[i]) - ti_toi[i] > REL_ERR_FLOOR) {
                            stats[m].violations.push_back({key, i, "late", ti_toi[i], double(toi[i])});
                        }
                    }
                    if (have_gt && gt[i] && !hit) stats[m].gt_false_negative += 1;

                    // --- the real invariant check, against the exact root ---
                    if (have_gt_toi) {
                        const double truth = gt_toi[i];
                        if (!sccd::is_nan_bits(truth)) {
                            stats[m].gt_checked += 1;
                            if (!hit) {
                                stats[m].gt_missed += 1;
                                if (stats[m].violations.size() < MAX_RECORDED_VIOLATIONS) {
                                    stats[m].violations.push_back(
                                        {key, i, "missed(truth)", truth, double(toi[i])});
                                }
                            } else if (double(toi[i]) > truth) {
                                stats[m].gt_late += 1;
                                const double over = double(toi[i]) - truth;
                                if (over > stats[m].gt_worst_overshoot) {
                                    stats[m].gt_worst_overshoot = over;
                                }
                                if (stats[m].violations.size() < MAX_RECORDED_VIOLATIONS) {
                                    stats[m].violations.push_back(
                                        {key, i, "late(truth)", truth, double(toi[i])});
                                }
                            }
                        }
                    }
                }
            }
        }

        std::printf("\n%s / %s  (%zu files, %zu queries; TI hits %zu)\n",
                    dataset.c_str(), phase.name, files.size(), total_queries, ti_hits);
        std::printf("%-8s %9s %8s %6s %6s %6s %8s %8s   %10s %10s %9s\n",
                    "mode", "queries", "hits", "FN!", "FP", "lateTI",
                    "gtMISS!", "gtLATE!", "relerr_med", "abserr_max", "ms");
        for (int m = 0; m < N_MODES; ++m) {
            Stats& s = stats[m];
            const char* name = mode_name(mode_of(m));
            const double rel_med = Stats::quantile_of(s.rel_err, 0.5);
            const double rel_p95 = Stats::quantile_of(s.rel_err, 0.95);
            const double abs_p95 = Stats::quantile_of(s.abs_err, 0.95);
            const double abs_max = Stats::max_of(s.abs_err);
            std::printf("%-8s %9zu %8zu %6zu %6zu %6zu %8zu %8zu   %10.3e %10.3e %9.1f\n",
                        name, s.n, s.hits, s.false_negative, s.false_positive, s.late,
                        s.gt_missed, s.gt_late, rel_med, abs_max, s.seconds * 1e3);
            (void)rel_p95;
            (void)abs_p95;
            char row[640];
            std::snprintf(row, sizeof(row),
                          "%s,%s,%s,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.6f,"
                          "%zu,%zu,%zu,%.17g",
                          dataset.c_str(), phase.name, name, s.n, s.hits, s.false_negative,
                          s.false_positive, s.late, s.near_zero_ref, s.gt_false_negative,
                          Stats::mean_of(s.rel_err), rel_med, rel_p95, abs_p95, abs_max,
                          s.seconds * 1e3, s.gt_checked, s.gt_missed, s.gt_late,
                          s.gt_worst_overshoot);
            csv_rows.push_back(row);
        }
        std::printf("%-8s %9zu %8zu %6s %6s %6s   %10s %10s   %10s %10s %9.1f\n",
                    "TI(ref)", total_queries, ti_hits, "-", "-", "-", "-", "-", "-", "-",
                    ti_seconds * 1e3);

        // The reference goes in the CSV too: "how fast is this mode" is only
        // meaningful next to the implementation it reproduces.
        {
            char row[640];
            std::snprintf(row, sizeof(row),
                          "%s,%s,%s,%zu,%zu,0,0,0,0,0,0,0,0,0,0,%.6f,0,0,0,0",
                          dataset.c_str(), phase.name, "ti-reference", total_queries, ti_hits,
                          ti_seconds * 1e3);
            csv_rows.push_back(row);
        }
        std::printf("  gtMISS!/gtLATE! are measured against the dataset's exact roots and are the\n"
                    "  real conservativeness test. FN!/lateTI compare with TightInclusion, whose\n"
                    "  own answer is a lower bound on the truth, so lateTI over-reports.\n"
                    "  relerr over the %zu/%zu queries with TI toi >= %g; the rest are covered by abserr.\n",
                    stats[0].rel_err.size(), stats[0].rel_err.size() + stats[0].near_zero_ref, REL_ERR_FLOOR);
#ifdef SCCD_ENABLE_CUDA
        std::printf("  device rows: CUDA narrow phase in %s. 'device-fast' is the\n"
                    "  mode-0 kernel; 'device-tight' is mode 2, with TightInclusion's\n"
                    "  predicate, split rule and numerical error bound.\n",
                    opt.device_float ? "single precision (--device-float)" : "double precision");
#endif
        if (gt_available) {
            std::printf("  (ground truth available for %zu queries; gt_FN: fast=%zu tight=%zu)\n",
                        gt_available, stats[0].gt_false_negative, stats[1].gt_false_negative);
        }

        for (int m = 0; m < N_MODES; ++m) {
            Stats& s = stats[m];
            // Gate on the exact roots when the dataset ships them; fall back to the
            // TightInclusion comparison only where it does not.
            const std::size_t bad = s.gt_checked ? (s.gt_missed + s.gt_late)
                                                 : (s.false_negative + s.late);
            if (bad == 0) {
                continue;
            }
            const bool gated = (opt.gate == "all") || (opt.gate == mode_name(mode_of(m)));
            if (gated) {
                violations_total += bad;
            }
            std::printf("  %s %s violates conservativeness on %zu queries "
                        "(%zu missed, %zu late). First few:\n",
                        gated ? "!!" : "  ", mode_name(mode_of(m)), bad,
                        s.gt_checked ? s.gt_missed : s.false_negative,
                        s.gt_checked ? s.gt_late : s.late);
            std::size_t shown = 0;
            for (const Violation& v : s.violations) {
                if (shown++ >= 5) break;
                std::printf("       %-10s %-8s query %-6zu  ref_toi=%.17g  got=%.17g  (%+.3e)\n",
                            v.file.c_str(), v.kind, v.query, v.ref_toi, v.got_toi, v.got_toi - v.ref_toi);
            }
            if (!violation_csv.empty()) {
                for (const Violation& v : s.violations) {
                    violation_rows.push_back(
                        std::string(dataset) + "," + phase.name + ","
                        + mode_name(mode_of(m)) + "," + v.file + ","
                        + std::to_string(v.query) + "," + v.kind + ","
                        + std::to_string(v.ref_toi) + "," + std::to_string(v.got_toi));
                }
            }
        }
    }

    if (!opt.csv.empty()) {
        std::ofstream out(opt.csv);
        out << "dataset,phase,mode,queries,hits,false_negative,false_positive,late,near_zero_ref,"
               "gt_false_negative,relerr_mean,relerr_median,relerr_p95,abserr_p95,abserr_max,ms,"
               "gt_checked,gt_missed,gt_late,gt_worst_overshoot\n";
        for (const std::string& r : csv_rows) {
            out << r << "\n";
        }
        std::cerr << "wrote " << opt.csv << "\n";
    }

    if (!violation_csv.empty() && !violation_rows.empty()) {
        std::ofstream out(violation_csv);
        out << "dataset,phase,mode,file,query,kind,ref_toi,got_toi\n";
        for (const std::string& r : violation_rows) {
            out << r << "\n";
        }
        std::cerr << "wrote " << violation_csv << "\n";
    }

    if (violations_total != 0) {
        std::printf("\nFAIL: %zu queries break conservativeness "
                    "(a missed collision, or a time of impact after the true one).\n",
                    violations_total);
        if (opt.strict && !geometry_is_narrowed) {
            return 2;
        }
        if (geometry_is_narrowed) {
            std::printf(
                "\nNOT GATED: the geometry was narrowed to float (--device-float or\n"
                "--float-geometry). The dataset's exact roots belong to the original\n"
                "rational geometry, so narrowing the input moves the true root and every\n"
                "correct kernel looks late against them. Control: the host tight kernel,\n"
                "which is bit-identical to TightInclusion and computes in double, reports\n"
                "~4500 'late' queries on armadillo-rollers under --float-geometry and zero\n"
                "without it. Use the double-geometry run to gate.\n");
        }
    } else {
        std::printf("\nOK: no missed collisions and no late times of impact%s.\n",
                    opt.gate == "all" ? "" : (" for mode " + opt.gate).c_str());
    }

    return 0;
}
