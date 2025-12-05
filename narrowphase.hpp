#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "vaabb.h"
#include "assert.h"
#include "roots.hpp"
#include "srootfinder.hpp"

// #define TIGHT_INCLUSION_ORACLE
#ifdef TIGHT_INCLUSION_ORACLE
#include "tight_inclusion/ccd.hpp"
#include "tight_inclusion/interval_root_finder.hpp"
#include <Eigen/Dense>
#endif

namespace sccd {

#ifdef TIGHT_INCLUSION_ORACLE
bool barycentricTriangle3D(
    const ticcd::Vector3& A,
    const ticcd::Vector3& B,
    const ticcd::Vector3& C,
    const ticcd::Vector3& P,
    double& u,
    double& v)
{
    using std::abs;

    ticcd::Vector3 e1 = B - A;
    ticcd::Vector3 e2 = C - A;
    ticcd::Vector3 n = e1.cross(e2).eval();

    ticcd::Vector3 dir = P - A;
    double dist = n.dot(dir);
    if (dist * dist > 1e-5) {
        return false;
    }

    // Compute local coordinates u and v: (P - A) = u * (B - A) + v * (C - A)
    // Solve: dir = u * e1 + v * e2
    // Using dot product method (more numerically stable):
    // dir · e1 = u * (e1 · e1) + v * (e2 · e1)
    // dir · e2 = u * (e1 · e2) + v * (e2 · e2)
    double d00 = e1.dot(e1);
    double d01 = e1.dot(e2);
    double d11 = e2.dot(e2);
    double d20 = dir.dot(e1);
    double d21 = dir.dot(e2);

    double denom = d00 * d11 - d01 * d01;
    if (abs(denom) < 1e-10) {
        // Degenerate triangle
        return false;
    }

    u = (d11 * d20 - d01 * d21) / denom;
    v = (d00 * d21 - d01 * d20) / denom;

    return true;
}

bool isInsideTriangle(
    const ticcd::Vector3& lambda, ticcd::Scalar tol = ticcd::Scalar(1e-6))
{
    return (lambda.array() >= -tol).all()
        && (lambda.array() <= ticcd::Scalar(1) + tol).all()
        && std::abs(lambda.sum() - ticcd::Scalar(1)) <= ticcd::Scalar(1e-6);
}

#endif

#ifdef TIGHT_INCLUSION_ORACLE
template <typename T>
bool find_root_oracle(
    const int max_iter,
    const T atol,
    const T sv[3],
    const T s1[3],
    const T s2[3],
    const T s3[3],
    const T ev[3],
    const T e1[3],
    const T e2[3],
    const T e3[3],
    T& t,
    T& u,
    T& v)
{
    ticcd::Vector3 v_t0(sv[0], sv[1], sv[2]);
    ticcd::Vector3 f0_t0(s1[0], s1[1], s1[2]);
    ticcd::Vector3 f1_t0(s2[0], s2[1], s2[2]);
    ticcd::Vector3 f2_t0(s3[0], s3[1], s3[2]);

    ticcd::Vector3 v_t1(ev[0], ev[1], ev[2]);
    ticcd::Vector3 f0_t1(e1[0], e1[1], e1[2]);
    ticcd::Vector3 f1_t1(e2[0], e2[1], e2[2]);
    ticcd::Vector3 f2_t1(e3[0], e3[1], e3[2]);

    ticcd::Array3 tol(atol, atol, atol);
    ticcd::Array3 err(1e-10, 1e-10, 1e-10);

    ticcd::Scalar ms = 0;
    ticcd::Scalar max_time = 1;
    ticcd::Scalar max_itr = 10000;
    ticcd::Scalar toi = 0;
    ticcd::Scalar output_tolerance = 1e-6;

    bool test_ok = ticcd::vertexFaceCCD(
        v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, err, ms, toi,
        1e-6, 1, 1000, output_tolerance, false,
        ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    double u0 = -1, v0 = -1;
    double discrepancy = -1;
    if (test_ok) {
        auto f0 = f0_t0 * (1 - toi) + toi * (f0_t1);
        auto f1 = f1_t0 * (1 - toi) + toi * (f1_t1);
        auto f2 = f2_t0 * (1 - toi) + toi * (f2_t1);
        auto pt = (1 - toi) * v_t0 + toi * v_t1;

        const bool inplane = barycentricTriangle3D(
            f0.eval(), f1.eval(), f2.eval(), pt.eval(), u0, v0);
        assert(inplane);

        test_ok =
            (u0 >= -1e-8 && v0 >= -1e-8 && u0 + v0 <= 1 + 1e-8 && toi >= -1e-8
             && toi <= 1 + 1e-8);

        auto pt_rec = (1 - u0 - v0) * f0 + u0 * f1 + v0 * f2;
        auto diff = pt_rec - pt;

        discrepancy = diff.dot(diff);
        t = toi;
        u = u0;
        v = v0;
    }
    return test_ok;
}

// #define TEST_ORACLE
#endif

template <int nxe, typename T>
T narrow_phase_vf(
    const size_t noverlaps,
    const idx_t* const SFEM_RESTRICT voveralp,
    const idx_t* const SFEM_RESTRICT foveralp,
    // Geometric data
    T** const SFEM_RESTRICT v0,
    T** const SFEM_RESTRICT v1,
    const size_t face_stride,
    idx_t** const SFEM_RESTRICT faces,
    // Output
    T* const SFEM_RESTRICT toi)
{
    const T infty = 100000;
    T min_t = infty;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, noverlaps),
                      [&](const tbb::blocked_range<size_t> &r) {
                        for (size_t i = r.begin(); i < r.end(); i++) {
    // for (size_t i = 0; i < noverlaps; i++) {
        const idx_t vi = voveralp[i];
        const idx_t fi = foveralp[i];

        idx_t nodes[3] = { faces[0][fi * face_stride],
                           faces[1][fi * face_stride],
                           faces[2][fi * face_stride] };

        const double sv[3] = { v0[0][vi], v0[1][vi], v0[2][vi] };
        const double ev[3] = { v1[0][vi], v1[1][vi], v1[2][vi] };

        const double s1[3] = { v0[0][nodes[0]], v0[1][nodes[0]],
                               v0[2][nodes[0]] };
        const double s2[3] = { v0[0][nodes[1]], v0[1][nodes[1]],
                               v0[2][nodes[1]] };
        const double s3[3] = { v0[0][nodes[2]], v0[1][nodes[2]],
                               v0[2][nodes[2]] };

        const double e1[3] = { v1[0][nodes[0]], v1[1][nodes[0]],
                               v1[2][nodes[0]] };
        const double e2[3] = { v1[0][nodes[1]], v1[1][nodes[1]],
                               v1[2][nodes[1]] };
        const double e3[3] = { v1[0][nodes[2]], v1[1][nodes[2]],
                               v1[2][nodes[2]] };

        // Iteration variables
        double t = 0;
        double u = 0;
        double v = 0;


#ifdef TEST_ORACLE
        double t_oracle = 0;
        double u_oracle = 0;
        double v_oracle = 0;
        bool test_ok = find_root_oracle<double>(
            1000, 1e-6, sv, s1, s2, s3, ev, e1, e2, e3, t_oracle, u_oracle,
            v_oracle);
#endif
        if (find_root_bisection<double>(
                400, 1e-10, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
            toi[i] = t;
            min_t = sccd::min<T>(t, min_t);
#ifdef TEST_ORACLE
            if (!test_ok) {
                printf(
                    "RF(Ours: (%g, %g, %g))\n", double(t), double(u),
                    double(v));
                // assert(false);
            }
#endif

        } else {
            toi[i] = infty;
#ifdef TEST_ORACLE
            if (test_ok) {
                printf(
                    "Missed root: (%g, %g, %g) == (%g, %g, %g)\n", t_oracle,
                    u_oracle, v_oracle, t, u, v);
                assert(false);
            }
#endif
        }
    }

    //
    });

    return min_t;
}

} // namespace sccd

#endif