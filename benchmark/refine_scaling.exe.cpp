// Scaling study: how the narrow and broad phases respond to surface refinement.
//
// The accuracy oracle and bench.exe both run fixed query sets, so neither says
// anything about how cost grows with mesh size. This one takes a closed surface,
// refines it repeatedly -- each level splits every triangle into four -- and runs
// the same collision step at every level.
//
// It reports element counts *and* candidate-pair counts, because those need not
// grow at the same rate: broad-phase pair count is what actually drives the
// narrow phase, and if it grows faster than the element count then a "narrow
// phase is O(n)" claim is about the wrong n. Per-pair time is printed for the
// same reason.
//
//   refine_scaling <levels> [mesh_t0] [mesh_t1]
//
// Without a mesh path it builds a tetrahedral cube and skins it, so the
// benchmark runs with no data dependency.

#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"

#include "sccd_config.hpp"
#include "sccd_narrowphase_mode.hpp"
#include "sccd_smesh_ccd.hpp"

#include <algorithm>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <vector>

using scalar_t = double;

namespace {

    /**
     * \brief Whether the mesh is already the surface we want to collide.
     *
     * The benchmark accepts both a volume mesh, which has to be skinned, and a
     * surface mesh, which must not be -- skinning a surface yields its boundary
     * *edges*, which is not a collidable surface at all. The CCD benchmark frames
     * are surfaces; a generated cube is a volume.
     */
    bool is_surface_mesh(const smesh::ElemType type) {
        return type == smesh::TRI3 || type == smesh::TRISHELL3 || type == smesh::QUAD4;
    }

    double now_ms() { return smesh::time_seconds() * 1e3; }

    /**
     * \brief Displace the surface so it self-intersects at a time set by the
     *        scene rather than by the mesh resolution.
     *
     * Getting this right took three attempts and the toi column caught each
     * failure, which is why it is printed.
     *
     *  - Pushing the two halves through each other reported toi exactly 1, 1/2,
     *    1/4, 1/8, 1/16 across five levels. That is not distant geometry meeting;
     *    it is adjacent elements either side of a discontinuity interpenetrating
     *    at once, so the contact time falls with the element size and every level
     *    measures a different problem.
     *  - A large rigid rotation never collided at all (toi stayed 1) and made the
     *    broad phase degenerate: swept boxes so large that pair counts grew 15x
     *    per level against 4x for the elements.
     *
     * What works is a point reflection through the centre with |scale| < 1: every
     * vertex travels through the middle of the model to the opposite side, so
     * opposite faces cross and the crossing time depends on a vertex's distance
     * from the centre -- a property of the scene -- and not on h. The
     * displacement stays bounded by the model size, so swept boxes stay sane.
     *
     * A convex surface cannot self-collide under any small deformation, so a
     * generated cube needs a motion this drastic. For a representative workload
     * pass a real pair of meshes instead; see the usage text.
     */
    void reflect(const std::vector<std::vector<scalar_t>>& p0,
                 const scalar_t scale,
                 std::vector<std::vector<scalar_t>>& p1) {
        const int dim = static_cast<int>(p0.size());
        const ptrdiff_t n = static_cast<ptrdiff_t>(p0[0].size());

        std::vector<scalar_t> center(dim, 0);
        for (int d = 0; d < dim; ++d) {
            const scalar_t lo = *std::min_element(p0[d].begin(), p0[d].end());
            const scalar_t hi = *std::max_element(p0[d].begin(), p0[d].end());
            center[d] = scalar_t(0.5) * (lo + hi);
        }

        p1 = p0;
        for (int d = 0; d < dim; ++d) {
            for (ptrdiff_t i = 0; i < n; ++i) {
                p1[d][i] = center[d] + scale * (p0[d][i] - center[d]);
            }
        }
    }

    std::vector<std::vector<scalar_t>> read_points(const std::shared_ptr<smesh::Mesh>& mesh) {
        auto pts = smesh::to_host(smesh::astype<scalar_t>(mesh->points()));
        const int dim = mesh->spatial_dimension();
        const ptrdiff_t n = mesh->n_nodes();

        std::vector<std::vector<scalar_t>> out(dim);
        for (int d = 0; d < dim; ++d) {
            out[d].assign(pts->data()[d], pts->data()[d] + n);
        }
        return out;
    }

    smesh::SharedBuffer<scalar_t*> to_2d(std::vector<std::vector<scalar_t>>& values,
                                         const smesh::ExecutionSpace space) {
        std::vector<smesh::SharedBuffer<scalar_t>> buffers;
        buffers.reserve(values.size());
        for (std::vector<scalar_t>& v : values) {
            auto b = smesh::create_host_buffer<scalar_t>(v.size());
            std::copy(v.begin(), v.end(), b->data());
            buffers.push_back(b);
        }
        if (space == smesh::EXECUTION_SPACE_DEVICE) {
            buffers = smesh::to_device(buffers);
        }
        return smesh::create_2d(buffers);
    }

}  // namespace

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);
    auto comm = ctx->communicator();

    if (argc < 2) {
        fprintf(stderr,
                "usage: %s <levels> [mesh_t0] [mesh_t1]\n"
                "\n"
                "Refines a closed surface <levels> times, quadrupling the triangle\n"
                "count at each level, and runs one collision step per level.\n"
                "Without a mesh path a tetrahedral cube is generated and skinned,\n"
                "and the motion is synthesized. Prefer passing a real pair of frames:\n"
                "a synthetic motion violent enough to make a convex mesh self-collide\n"
                "also makes every swept box overlap every other, so the broad phase\n"
                "degenerates to all pairs and the totals scale as n^2 rather than as a\n"
                "solver would see.\n"
                "\n"
                "Env: SCCD_NARROWPHASE_MODE, SCCD_MAX_DEPTH, SCCD_TOL,\n"
                "     SCCD_BENCH_EXECUTION_SPACE=device, SCCD_SCALE,\n"
                "     SCCD_TOPOLOGY=quad (generate a hex cube, so the skin is\n"
                "     QUADSHELL4 rather than TRISHELL3).\n",
                argv[0]);
        return 1;
    }

    const int levels = atoi(argv[1]);

    int SCCD_MAX_DEPTH = 69;
    SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);
    scalar_t SCCD_TOL = scalar_t(3e-8);
    SCCD_READ_ENV(SCCD_TOL, atof);
    // Scale about the centre for the synthesized motion. Negative reflects the
    // model through its centre, which is what makes it self-intersect; see
    // reflect() for why a gentler motion does not work on a convex mesh.
    double SCCD_SCALE = -0.5;
    SCCD_READ_ENV(SCCD_SCALE, atof);

    smesh::ExecutionSpace space = smesh::EXECUTION_SPACE_HOST;
    if (const char* s = getenv("SCCD_BENCH_EXECUTION_SPACE")) {
        const std::string v(s);
        if (v == "device" || v == "cuda" || v == "gpu") {
#if defined(SCCD_ENABLE_CUDA)
            space = smesh::EXECUTION_SPACE_DEVICE;
#else
            fprintf(stderr, "warn: built without CUDA, running on the host\n");
#endif
        }
    }

    // SCCD_TOPOLOGY=quad generates a hexahedral cube instead of a tetrahedral
    // one. Skinning a hex mesh yields QUADSHELL4 where skinning a tet mesh
    // yields TRISHELL3, so this is the whole of what it takes to run the same
    // scaling study against the quad narrow phase -- which needs its own
    // measurements, not triangle numbers assumed to carry over. Ignored when a
    // mesh is supplied, since the file then fixes the topology.
    const char* topology_env = getenv("SCCD_TOPOLOGY");
    const bool want_quad = topology_env && std::string(topology_env) == "quad";

    std::shared_ptr<smesh::Mesh> base;
    std::shared_ptr<smesh::Mesh> base_t1;
    if (argc >= 3) {
        base = smesh::Mesh::create_from_file(comm, smesh::Path(argv[2]));
    } else if (want_quad) {
        base = smesh::Mesh::create_hex8_cube(comm, 2, 2, 2, 0, 0, 0, 1, 1, 1);
    } else {
        base = smesh::Mesh::create_tet4_cube(comm, 2, 2, 2, 0, 0, 0, 1, 1, 1);
    }
    if (argc >= 4) {
        // A real second frame. Much preferable to the synthesized motion: the
        // displacement is then a physical step, so the broad phase sees the pair
        // count a solver would see instead of the all-pairs blowup any drastic
        // synthetic motion produces. Both meshes are refined the same way, which
        // is well defined because they share a connectivity.
        base_t1 = smesh::Mesh::create_from_file(comm, smesh::Path(argv[3]));
        if (base_t1->n_nodes() != base->n_nodes() || base_t1->n_elements() != base->n_elements()) {
            fprintf(stderr, "error: the two meshes must share a topology\n");
            return 1;
        }
    }

    // A volume mesh is refined and re-skinned at each level; a surface mesh is
    // refined directly. Note smesh's refine rejects TRISHELL3, the type skin
    // produces ("Refinement factor not supported for element type 103"), which is
    // why the volume is refined before skinning rather than after.

    const bool base_is_surface = is_surface_mesh(base->block(0)->element_type());

    printf("# mode=%s max_depth=%d tol=%g scale=%g space=%s base_topology=%s\n",
           sccd::narrow_phase_mode_name(sccd::narrow_phase_mode()),
           SCCD_MAX_DEPTH,
           (double)SCCD_TOL,
           SCCD_SCALE,
           space == smesh::EXECUTION_SPACE_DEVICE ? "device" : "host",
           smesh::type_to_string(base->block(0)->element_type()));
    printf("%5s %10s %12s %12s %9s %9s %9s %9s %10s %10s %14s\n",
           "level",
           "faces",
           "vf_pairs",
           "ee_pairs",
           "prep_ms",
           "bp_fv_ms",
           "bp_ee_ms",
           "broad_ms",
           "narrow_ms",
           "ns/pair",
           "toi");

    for (int level = 0; level <= levels; ++level) {
        std::shared_ptr<smesh::Mesh> refined = (level == 0) ? base : smesh::refine(base, level);
        std::shared_ptr<smesh::Mesh> mesh = base_is_surface ? refined : smesh::skin(refined);

        std::vector<std::vector<scalar_t>> h0 = read_points(mesh);
        std::vector<std::vector<scalar_t>> h1;
        if (base_t1) {
            auto refined_t1 = (level == 0) ? base_t1 : smesh::refine(base_t1, level);
            h1 = read_points(base_is_surface ? refined_t1 : smesh::skin(refined_t1));
        } else {
            reflect(h0, static_cast<scalar_t>(SCCD_SCALE), h1);
        }

        auto points0 = to_2d(h0, space);
        auto points1 = to_2d(h1, space);

        auto ccd = sccd::CCD<scalar_t>::create(mesh, space);

        smesh::SharedBuffer<smesh::idx_t> v_overlap, f_overlap, e0_overlap, e1_overlap;

        // Broken out rather than timed as one call: the prep (AABBs + the sort
        // that sweep-and-prune needs) scales with elements, while the two sweeps
        // scale with how much the sorted intervals overlap. Those are different
        // costs with different fixes, and a single broad_ms number hides which
        // one is growing.
        const double t_prep0 = now_ms();
        ccd->broad_phase_prep(points0, points1);
        const double prep_ms = now_ms() - t_prep0;

        const double t_fv0 = now_ms();
        ccd->broad_phase_fv_step(v_overlap, f_overlap);
        const double fv_ms = now_ms() - t_fv0;

        const double t_ee0 = now_ms();
        ccd->broad_phase_ee_step(e0_overlap, e1_overlap);
        const double ee_ms = now_ms() - t_ee0;

        const double broad_ms = prep_ms + fv_ms + ee_ms;

        scalar_t toi = 1;
        smesh::SharedBuffer<scalar_t> vf_tois, ee_tois;
        const double t_narrow0 = now_ms();
        ccd->narrow_phase(toi, vf_tois, ee_tois, SCCD_MAX_DEPTH, SCCD_TOL);
        const double narrow_ms = now_ms() - t_narrow0;

        const ptrdiff_t n_vf = f_overlap ? f_overlap->size() : 0;
        const ptrdiff_t n_ee = e0_overlap ? e0_overlap->size() : 0;
        const ptrdiff_t n_pairs = n_vf + n_ee;

        // The time of impact, not just the timings. Without it the table cannot
        // answer the only question that matters when comparing two execution
        // spaces: do they agree on the answer? Printed at full precision because
        // the interesting differences are in the last digits.
        printf("%5d %10ld %12ld %12ld %9.2f %9.2f %9.2f %9.2f %10.2f %10.1f %14.10f\n",
               level,
               (long)mesh->block(0)->n_elements(),
               (long)n_vf,
               (long)n_ee,
               prep_ms,
               fv_ms,
               ee_ms,
               broad_ms,
               narrow_ms,
               n_pairs > 0 ? (narrow_ms * 1e6 / (double)n_pairs) : 0.0,
               (double)toi);
        fflush(stdout);
    }

    return 0;
}
