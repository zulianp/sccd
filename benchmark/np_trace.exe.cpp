// Per-query box counts for the narrow phase, host and device, one query at a time.
//
//     np_trace <dataset-dir> <vf|ee> [options]
//
// Step 0 of wip/CUDA_NARROWPHASE_PLAN.md: establish where the device's 94x box
// count comes from before proposing anything.
//
// The measurement that matters here is **one query in isolation**. Running a
// whole scene mixes two effects that an aggregate cannot separate: how each
// kernel searches a single query, and how the shared time-of-impact bound
// collapses across queries as the run proceeds. This driver removes the second
// entirely -- every query is a separate call with max_toi = 1 and toi_stride = 1,
// so nothing another query found can prune it.
//
// If the device still costs far more than the host on an isolated query, the
// cause is inside the per-query search and a box-by-box trace will find it. If
// the two agree in isolation, the cause is the bound schedule and the search is
// fine.
//
// Requires -DSCCD_NP_COUNT_BOXES for the counts to exist at all; without it the
// driver still runs and reports times of impact, which is a useful cross-check on
// its own.
//
// Options:
//   --file KEY      one query file (e.g. 227vf); default: every file in the set
//   --query N       one query index within the file, printed on its own
//   --top N         print the N costliest queries (default 20)
//   --tol T         distance tolerance (default 3e-8)
//   --max-depth D   depth cap (default 69)
//   --mode M        SCCD_NARROWPHASE_MODE for the host side (default 2)
//   --csv PATH      write every query's counts, for offline analysis
//   --batch         also run the whole file as ONE toi_stride=0 call, which is
//                   how the library is actually used, and report the same counts.
//                   The difference between the two is the value of the collapsing
//                   shared bound, isolated on identical data.

#include "sccd_narrowphase.hpp"
#include "sccd_narrowphase_mode.hpp"

#ifdef SCCD_ENABLE_CUDA
#include "sccd_narrowphase.cuh"
#include <cuda_runtime.h>
#endif

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

using scalar_t = double;
using idx_t = int;

#ifdef SCCD_NP_COUNT_BOXES
extern "C" {
}
#endif

// One query: four points at the start of the step and four at the end. For
// vertex-face they are (vertex, f0, f1, f2); for edge-edge (ea0, ea1, eb0, eb1).
struct Query {
    scalar_t p0[4][3];
    scalar_t p1[4][3];
};

bool parse_row(const std::string& line, scalar_t out[3]) {
    double v[6];
    std::size_t cursor = 0;
    for (int i = 0; i < 6; ++i) {
        const std::size_t end = line.find(',', cursor);
        if (i < 5 && end == std::string::npos) return false;
        try {
            v[i] = std::stod(line.substr(cursor, end - cursor));
        } catch (...) {
            return false;
        }
        cursor = (i < 5) ? end + 1 : end;
    }
    // Coordinates are exact dyadic rationals, numerator and denominator per axis.
    for (int d = 0; d < 3; ++d) out[d] = static_cast<scalar_t>(v[2 * d] / v[2 * d + 1]);
    return true;
}

bool read_queries(const fs::path& path, std::vector<Query>& out) {
    std::ifstream in(path);
    if (!in) return false;
    std::string line;
    std::vector<scalar_t> rows;  // 3 per row, 8 rows per query
    std::size_t n_rows = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        scalar_t c[3];
        if (!parse_row(line, c)) return false;
        rows.insert(rows.end(), c, c + 3);
        ++n_rows;
    }
    if (n_rows == 0 || n_rows % 8 != 0) return false;
    out.resize(n_rows / 8);
    for (std::size_t q = 0; q < out.size(); ++q) {
        for (int k = 0; k < 4; ++k)
            for (int d = 0; d < 3; ++d) {
                out[q].p0[k][d] = rows[(q * 8 + k) * 3 + d];
                out[q].p1[k][d] = rows[(q * 8 + 4 + k) * 3 + d];
            }
    }
    return true;
}

// The host kernels take structure-of-arrays geometry plus a connectivity table.
// One query becomes four nodes and one element.
struct SingleQueryScene {
    std::vector<scalar_t> x0, y0, z0, x1, y1, z1;
    std::vector<idx_t> e[3];
    std::vector<idx_t> c0, c1;
    scalar_t* p0[3];
    scalar_t* p1[3];
    idx_t* elem[3];

    explicit SingleQueryScene(const Query& q, const bool is_vf) {
        for (int k = 0; k < 4; ++k) {
            x0.push_back(q.p0[k][0]); y0.push_back(q.p0[k][1]); z0.push_back(q.p0[k][2]);
            x1.push_back(q.p1[k][0]); y1.push_back(q.p1[k][1]); z1.push_back(q.p1[k][2]);
        }
        if (is_vf) {
            e[0] = {1}; e[1] = {2}; e[2] = {3};
            c0 = {0};   // the vertex node
            c1 = {0};   // the face
        } else {
            e[0] = {0}; e[1] = {1};   // first edge; the second is nodes 2,3
            e[0] = {0, 2};
            e[1] = {1, 3};
            c0 = {0};
            c1 = {1};
        }
        p0[0] = x0.data(); p0[1] = y0.data(); p0[2] = z0.data();
        p1[0] = x1.data(); p1[1] = y1.data(); p1[2] = z1.data();
        for (int d = 0; d < 3; ++d) elem[d] = e[d].empty() ? nullptr : e[d].data();
    }
};

// Every query in a file as one call: 4n nodes, n elements, n candidate pairs.
struct BatchScene {
    std::vector<scalar_t> c0[3], c1[3];
    std::vector<idx_t> e[3];
    std::vector<idx_t> a, b;
    scalar_t* p0[3];
    scalar_t* p1[3];
    idx_t* elem[3];

    BatchScene(const std::vector<Query>& qs, const bool is_vf) {
        for (std::size_t q = 0; q < qs.size(); ++q) {
            for (int k = 0; k < 4; ++k)
                for (int d = 0; d < 3; ++d) {
                    c0[d].push_back(qs[q].p0[k][d]);
                    c1[d].push_back(qs[q].p1[k][d]);
                }
            const idx_t base = static_cast<idx_t>(4 * q);
            if (is_vf) {
                e[0].push_back(base + 1);
                e[1].push_back(base + 2);
                e[2].push_back(base + 3);
                a.push_back(base);                    // vertex node
                b.push_back(static_cast<idx_t>(q));   // face id
            } else {
                e[0].push_back(base);
                e[1].push_back(base + 1);
                e[0].push_back(base + 2);
                e[1].push_back(base + 3);
                a.push_back(static_cast<idx_t>(2 * q));
                b.push_back(static_cast<idx_t>(2 * q + 1));
            }
        }
        for (int d = 0; d < 3; ++d) {
            p0[d] = c0[d].data();
            p1[d] = c1[d].data();
            elem[d] = e[d].empty() ? nullptr : e[d].data();
        }
    }
};

struct Row {
    std::size_t index = 0;
    unsigned long long host_boxes = 0;
    unsigned long long device_boxes = 0;
    scalar_t host_toi = 0;
    scalar_t device_toi = 0;
};

unsigned long long host_box_count() {
#ifdef SCCD_NP_COUNT_BOXES
    return sccd::g_np_host_boxes;
#else
    return 0ull;
#endif
}

void usage(const char* argv0) {
    std::fprintf(stderr,
                 "usage: %s <dataset-dir> <vf|ee> [--file KEY] [--query N] [--top N]\n"
                 "          [--tol T] [--max-depth D] [--mode M] [--csv PATH]\n",
                 argv0);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }
    const fs::path dataset = argv[1];
    const std::string phase = argv[2];
    if (phase != "vf" && phase != "ee") {
        usage(argv[0]);
        return 1;
    }
    const bool is_vf = (phase == "vf");

    std::string only_file;
    long only_query = -1;
    int top = 20;
    scalar_t tol = 3e-8;
    int max_depth = 69;
    int mode = 2;
    std::string csv_path;
    bool batch = false;

    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? argv[++i] : std::string(); };
        if (a == "--file") only_file = next();
        else if (a == "--query") only_query = std::stol(next());
        else if (a == "--top") top = std::stoi(next());
        else if (a == "--tol") tol = std::stod(next());
        else if (a == "--max-depth") max_depth = std::stoi(next());
        else if (a == "--mode") mode = std::stoi(next());
        else if (a == "--csv") csv_path = next();
        else if (a == "--batch") batch = true;
        else { usage(argv[0]); return 1; }
    }

    // The mode is read from the environment on every narrow-phase call, so
    // setting it here covers both kernels.
    setenv("SCCD_NARROWPHASE_MODE", std::to_string(mode).c_str(), 1);

    const fs::path queries_dir = dataset / "queries";
    if (!fs::is_directory(queries_dir)) {
        std::fprintf(stderr, "error: no queries directory at %s\n", queries_dir.c_str());
        return 1;
    }

    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(queries_dir)) {
        const std::string name = entry.path().filename().string();
        if (name.size() < 3 || name.substr(name.size() - 4) != ".csv") continue;
        if (name.find(phase + ".csv") == std::string::npos) continue;
        if (!only_file.empty() && name.find(only_file) != 0) continue;
        files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) {
        std::fprintf(stderr, "error: no %s query files under %s\n", phase.c_str(), queries_dir.c_str());
        return 1;
    }

#ifndef SCCD_NP_COUNT_BOXES
    std::fprintf(stderr,
                 "note: built without -DSCCD_NP_COUNT_BOXES; box counts will read 0.\n"
                 "      Times of impact are still reported.\n");
#endif

    std::printf("# np_trace  phase=%s mode=%d tol=%g max_depth=%d  (one query per call, max_toi=1)\n",
                phase.c_str(), mode, (double)tol, max_depth);

    std::vector<Row> rows;
    std::size_t total = 0;

    for (const auto& file : files) {
        std::vector<Query> queries;
        if (!read_queries(file, queries)) {
            std::fprintf(stderr, "warning: could not read %s\n", file.c_str());
            continue;
        }
        for (std::size_t q = 0; q < queries.size(); ++q) {
            if (only_query >= 0 && (long)q != only_query) continue;

            SingleQueryScene scene(queries[q], is_vf);
            scalar_t toi = 1.0;

            const unsigned long long before = host_box_count();
            if (is_vf) {
                sccd::narrow_phase_vf<3, scalar_t, idx_t>(1, scene.c0.data(), scene.c1.data(),
                                                          scene.p0, scene.p1, 1, scene.elem,
                                                          /*max_toi=*/1.0, &toi, max_depth, tol,
                                                          /*toi_stride=*/1);
            } else {
                sccd::narrow_phase_ee<scalar_t, idx_t>(1, scene.c0.data(), scene.c1.data(),
                                                       scene.p0, scene.p1, 1, scene.elem,
                                                       /*max_toi=*/1.0, &toi, max_depth, tol,
                                                       /*toi_stride=*/1);
            }
            Row r;
            r.index = total;
            r.host_boxes = host_box_count() - before;
            r.host_toi = toi;
            rows.push_back(r);
            ++total;
        }
    }

    std::sort(rows.begin(), rows.end(),
              [](const Row& a, const Row& b) { return a.host_boxes > b.host_boxes; });

    unsigned long long sum = 0;
    for (const Row& r : rows) sum += r.host_boxes;
    std::printf("# %zu queries, %llu host boxes, %.2f per query\n",
                rows.size(), sum, rows.empty() ? 0.0 : (double)sum / (double)rows.size());
    std::printf("%8s %12s %16s\n", "query", "host_boxes", "host_toi");
    for (int i = 0; i < top && i < (int)rows.size(); ++i) {
        std::printf("%8zu %12llu %16.9f\n", rows[i].index, rows[i].host_boxes, (double)rows[i].host_toi);
    }

    if (batch) {
        // The same queries, in one call, the way the library is used. Everything
        // that differs from the per-query total above is the shared bound.
        std::size_t nq = 0;
        unsigned long long batch_boxes = 0;
        scalar_t batch_toi = 1.0;
        for (const auto& file : files) {
            std::vector<Query> queries;
            if (!read_queries(file, queries)) continue;
            BatchScene scene(queries, is_vf);
            nq += queries.size();
            const unsigned long long before = host_box_count();
            if (is_vf) {
                sccd::narrow_phase_vf<3, scalar_t, idx_t>(queries.size(), scene.a.data(), scene.b.data(),
                                                          scene.p0, scene.p1, 1, scene.elem,
                                                          /*max_toi=*/1.0, &batch_toi, max_depth, tol,
                                                          /*toi_stride=*/0);
            } else {
                sccd::narrow_phase_ee<scalar_t, idx_t>(queries.size(), scene.a.data(), scene.b.data(),
                                                       scene.p0, scene.p1, 1, scene.elem,
                                                       /*max_toi=*/1.0, &batch_toi, max_depth, tol,
                                                       /*toi_stride=*/0);
            }
            batch_boxes += host_box_count() - before;
        }
        std::printf("\n# batched: %zu queries in one toi_stride=0 call\n", nq);
        std::printf("#   host boxes %llu, %.2f per query, earliest toi %.9f\n",
                    batch_boxes, nq ? (double)batch_boxes / (double)nq : 0.0, (double)batch_toi);
        std::printf("#   isolated / batched = %.1fx\n",
                    batch_boxes ? (double)sum / (double)batch_boxes : 0.0);
    }

    if (!csv_path.empty()) {
        std::ofstream out(csv_path);
        out << "query,host_boxes,host_toi\n";
        for (const Row& r : rows) out << r.index << ',' << r.host_boxes << ',' << r.host_toi << '\n';
        std::printf("# wrote %s\n", csv_path.c_str());
    }
    return 0;
}
