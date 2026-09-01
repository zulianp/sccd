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

#include "narrowphase.hpp"
#include "srootfinder.hpp"

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

    bool read_queries(const fs::path& path, const bool is_vf, QuerySet& out) {
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
    enum class Mode { Scalar, Vector, TiVector };
    constexpr int N_MODES = 3;

    const char* mode_name(const Mode m) {
        switch (m) {
            case Mode::Scalar: return "scalar";
            case Mode::Vector: return "vector";
            case Mode::TiVector: return "ti-vec";
        }
        return "?";
    }

    Mode mode_of(const int m) { return static_cast<Mode>(m); }

    void select_mode(const Mode m) {
        // The kernels read this at entry, so the mode is selected per call.
        const char* value = (m == Mode::TiVector) ? "2" : (m == Mode::Vector ? "1" : "0");
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
        } else if (a == "--gate") {
            opt.gate = next();
        } else if (!a.empty() && a[0] != '-') {
            opt.dataset_dir = a;
        }
    }

    if (opt.dataset_dir.empty()) {
        std::cerr << "usage: ti_oracle <dataset-dir> [--phase vf|ee|both] [--max-files N]\n"
                     "                 [--tol T] [--max-depth N] [--csv out.csv]\n"
                     "                 [--violations-csv out.csv] [--no-strict] [--gate MODE]\n"
                     "\n"
                     "Exits non-zero when a mode misses a collision or reports a time of\n"
                     "impact later than TightInclusion's; both break conservativeness.\n"
                     "Pass --no-strict to report without failing.\n";
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

        Stats stats[N_MODES];
        double ti_seconds = 0;
        std::size_t ti_hits = 0;
        std::size_t total_queries = 0;
        std::size_t gt_available = 0;

        for (const fs::path& file : files) {
            QuerySet qs;
            if (!read_queries(file, phase.is_vf, qs)) {
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
                        if (!std::isnan(truth)) {
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
        if (gt_available) {
            std::printf("  (ground truth available for %zu queries; gt_FN: scalar=%zu vector=%zu ti-vec=%zu)\n",
                        gt_available, stats[0].gt_false_negative, stats[1].gt_false_negative,
                        stats[2].gt_false_negative);
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
        if (opt.strict) {
            return 2;
        }
    } else {
        std::printf("\nOK: no missed collisions and no late times of impact%s.\n",
                    opt.gate == "all" ? "" : (" for mode " + opt.gate).c_str());
    }

    return 0;
}
