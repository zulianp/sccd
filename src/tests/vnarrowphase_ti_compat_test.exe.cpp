#include "sccd_vnarrowphase.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>

int main() {
#ifndef SCCD_ENABLE_TIGHT_INCLUSION
    return 0;
#else
    using I = std::int32_t;
    constexpr int nqueries = 2;
    constexpr int nnodes = 5;
    constexpr double tolerance = 3e-8;

    double x0[nnodes] = {0, 1, 0, 0.25, 2};
    double y0[nnodes] = {0, 0, 1, 0.25, 2};
    double z0[nnodes] = {0, 0, 0, 1, 1};
    double x1[nnodes] = {0, 1, 0, 0.25, 2};
    double y1[nnodes] = {0, 0, 1, 0.25, 2};
    double z1[nnodes] = {0, 0, 0, -1, -1};
    double* start[3] = {x0, y0, z0};
    double* end[3] = {x1, y1, z1};

    I face0[1] = {0};
    I face1[1] = {1};
    I face2[1] = {2};
    I* faces[3] = {face0, face1, face2};
    I vertices[nqueries] = {3, 4};
    I face_ids[nqueries] = {0, 0};
    double vector_toi[nqueries] = {1, 1};

    setenv("SCCD_VNARROWPHASE_TI_COMPAT", "1", 1);
    const int status = sccd::v_narrow_phase_vf<3, double, I>(nqueries,
                                                              vertices,
                                                              face_ids,
                                                              start,
                                                              end,
                                                              1,
                                                              faces,
                                                              1,
                                                              vector_toi,
                                                              96,
                                                              tolerance,
                                                              1);
    unsetenv("SCCD_VNARROWPHASE_TI_COMPAT");
    if (status != 0) {
        return 1;
    }

    double expected_toi[nqueries];
    for (int q = 0; q < nqueries; ++q) {
        const I vertex = vertices[q];
        const double sv[3] = {x0[vertex], y0[vertex], z0[vertex]};
        const double ev[3] = {x1[vertex], y1[vertex], z1[vertex]};
        const double s1[3] = {x0[0], y0[0], z0[0]};
        const double s2[3] = {x0[1], y0[1], z0[1]};
        const double s3[3] = {x0[2], y0[2], z0[2]};
        const double e1[3] = {x1[0], y1[0], z1[0]};
        const double e2[3] = {x1[1], y1[1], z1[1]};
        const double e3[3] = {x1[2], y1[2], z1[2]};
        double ti_t = 1;
        double ti_u = 0;
        double ti_v = 0;
        const bool hit = sccd::find_root_tight_inclusion_vf<double>(
            96, tolerance, sv, s1, s2, s3, ev, e1, e2, e3, ti_t, ti_u, ti_v);
        expected_toi[q] = hit ? ti_t : 1;
        if (vector_toi[q] != expected_toi[q]) {
            return 2 + q;
        }
    }

    double global_toi = 1;
    setenv("SCCD_VNARROWPHASE_TI_COMPAT", "1", 1);
    const int global_status = sccd::v_narrow_phase_vf<3, double, I>(nqueries,
                                                                     vertices,
                                                                     face_ids,
                                                                     start,
                                                                     end,
                                                                     1,
                                                                     faces,
                                                                     1,
                                                                     &global_toi,
                                                                     96,
                                                                     tolerance,
                                                                     0);
    unsetenv("SCCD_VNARROWPHASE_TI_COMPAT");
    if (global_status != 0 || global_toi != std::min(expected_toi[0], expected_toi[1])) {
        return 4;
    }

    return 0;
#endif
}
