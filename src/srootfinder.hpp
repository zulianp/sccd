#ifndef S_ROOT_FINDER_HPP
#define S_ROOT_FINDER_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "vaabb.hpp"

#include "roots.hpp"
#include "snumtol.hpp"

// #define SCCD_ENABLE_TIGHT_INCLUSION
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#include <Eigen/Dense>
#include "tight_inclusion/ccd.hpp"
#include "tight_inclusion/interval_root_finder.hpp"
#endif

namespace sccd {

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
    static bool barycentric_triangle_3d(const ticcd::Vector3 &A,
                                 const ticcd::Vector3 &B,
                                 const ticcd::Vector3 &C,
                                 const ticcd::Vector3 &P,
                                 double &u,
                                 double &v) {
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

    static bool isInsideTriangle(const ticcd::Vector3 &lambda, ticcd::Scalar tol = ticcd::Scalar(1e-6)) {
        return (lambda.array() >= -tol).all() && (lambda.array() <= ticcd::Scalar(1) + tol).all() &&
               std::abs(lambda.sum() - ticcd::Scalar(1)) <= ticcd::Scalar(1e-6);
    }

    template <typename T>
    bool find_root_tight_inclusion_vf(const int max_iter,
                                      const T atol,
                                      const T sv[3],
                                      const T s1[3],
                                      const T s2[3],
                                      const T s3[3],
                                      const T ev[3],
                                      const T e1[3],
                                      const T e2[3],
                                      const T e3[3],
                                      T &t,
                                      T &u,
                                      T &v) {
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
        ticcd::Scalar toi = 1;
        ticcd::Scalar output_tolerance = 1e-6;
        bool no_zero_toi = false;

        return ticcd::vertexFaceCCD(v_t0,
                                    f0_t0,
                                    f1_t0,
                                    f2_t0,
                                    v_t1,
                                    f0_t1,
                                    f1_t1,
                                    f2_t1,
                                    err,
                                    ms,
                                    toi,
                                    atol,
                                    1,
                                    max_iter,
                                    output_tolerance,
                                    no_zero_toi,
                                    ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

        // double u0 = -1, v0 = -1;
        // double discrepancy = -1;
        // if (test_ok) {
        //     auto f0 = f0_t0 * (1 - toi) + toi * (f0_t1);
        //     auto f1 = f1_t0 * (1 - toi) + toi * (f1_t1);
        //     auto f2 = f2_t0 * (1 - toi) + toi * (f2_t1);
        //     auto pt = (1 - toi) * v_t0 + toi * v_t1;

        //     const bool inplane = barycentric_triangle_3d(f0.eval(), f1.eval(), f2.eval(), pt.eval(), u0, v0);
        //     assert(inplane);

        //     test_ok = (u0 >= -1e-8 && v0 >= -1e-8 && u0 + v0 <= 1 + 1e-8 && toi >= -1e-8 && toi <= 1 + 1e-8);

        //     auto pt_rec = (1 - u0 - v0) * f0 + u0 * f1 + v0 * f2;
        //     auto diff = pt_rec - pt;

        //     discrepancy = diff.dot(diff);
        //     t = toi;
        //     u = u0;
        //     v = v0;
        // }
    }

    template <typename T>
    bool find_root_tight_inclusion_ee(const int max_iter,
                                      const T atol,
                                      const T s1[3],
                                      const T s2[3],
                                      const T s3[3],
                                      const T s4[3],
                                      const T e1[3],
                                      const T e2[3],
                                      const T e3[3],
                                      const T e4[3],
                                      T &t,
                                      T &u,
                                      T &v) {
        ticcd::Vector3 e1_t0(s1[0], s1[1], s1[2]);
        ticcd::Vector3 e2_t0(s2[0], s2[1], s2[2]);
        ticcd::Vector3 e3_t0(s3[0], s3[1], s3[2]);
        ticcd::Vector3 e4_t0(s4[0], s4[1], s4[2]);

        ticcd::Vector3 e1_t1(e1[0], e1[1], e1[2]);
        ticcd::Vector3 e2_t1(e2[0], e2[1], e2[2]);
        ticcd::Vector3 e3_t1(e3[0], e3[1], e3[2]);
        ticcd::Vector3 e4_t1(e4[0], e4[1], e4[2]);
        ticcd::Array3 tol(atol, atol, atol);
        ticcd::Array3 err(1e-10, 1e-10, 1e-10);

        ticcd::Scalar ms = 0;
        ticcd::Scalar max_time = 1;
        ticcd::Scalar toi = 1;
        ticcd::Scalar output_tolerance = 1e-6;
        bool no_zero_toi = true;
        return ticcd::edgeEdgeCCD(e1_t0,
                                  e2_t0,
                                  e3_t0,
                                  e4_t0,
                                  e1_t1,
                                  e2_t1,
                                  e3_t1,
                                  e4_t1,
                                  err,
                                  ms,
                                  toi,
                                  atol,
                                  1,
                                  max_iter,
                                  output_tolerance,
                                  no_zero_toi,
                                  ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);
    }

#endif

    template <typename T>
    inline void project_uv_simplex(T &u, T &v) {
        u = sccd::max<T>(u, 0);
        v = sccd::max<T>(v, 0);
        const T s = u + v;
        if (s <= static_cast<T>(1)) {
            return;
        }

        T u_proj = static_cast<T>(0.5) * (u - v + 1);
        u_proj = sccd::min<T>(static_cast<T>(1), sccd::max<T>(0, u_proj));
        v = static_cast<T>(1) - u_proj;
        u = u_proj;
    }

    template <typename T>
    inline void diff_vf(const T sv[3],
                        const T s1[3],
                        const T s2[3],
                        const T s3[3],
                        const T ev[3],
                        const T e1[3],
                        const T e2[3],
                        const T e3[3],
                        const T &t,
                        const T &u,
                        const T &v,
                        T *const SCCD_RESTRICT diff) {
        T t0 = (1 - t);
        T t1 = t;
        T o = (1 - u - v);
        for (int d = 0; d < 3; d++) {
            T v_pos = t0 * sv[d] + t1 * ev[d];
            T f0 = t0 * (o * s1[d] + u * s2[d] + v * s3[d]);
            T f1 = t1 * (o * e1[d] + u * e2[d] + v * e3[d]);
            T f = f0 + f1;
            diff[d] = v_pos - f;
        }
    }

    template <typename T>
    inline T norm_diff_vf(const T sv[3],
                          const T s1[3],
                          const T s2[3],
                          const T s3[3],
                          const T ev[3],
                          const T e1[3],
                          const T e2[3],
                          const T e3[3],
                          T &t,
                          T &u,
                          T &v) {
        T diff[3];
        diff_vf(sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, diff);
        return sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
    }

    template <typename T>
    bool find_root_newton(const int max_iter,
                          const T atol,
                          const T sv[3],
                          const T s1[3],
                          const T s2[3],
                          const T s3[3],
                          const T ev[3],
                          const T e1[3],
                          const T e2[3],
                          const T e3[3],
                          T &t,
                          T &u,
                          T &v) {
        project_uv_simplex<T>(u, v);
        t = sccd::min<T>(static_cast<T>(1), sccd::max<T>(0, t));

        T s4[3] = {0, 0, 0};
        T e4[3] = {0, 0, 0};

        T f = 0;
        vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f);

        for (int k = 0; k < max_iter; k++) {
            T p[3] = {0, 0, 0};
            vf_objective_dir<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f, p);

            T best_t = t;
            T best_u = u;
            T best_v = v;
            T best_f = f;

            T alpha = 1;
            bool improved = false;
            for (int j = 0; j < 12; j++) {
                T cand_t = t - alpha * p[0];
                T cand_u = u - alpha * p[1];
                T cand_v = v - alpha * p[2];

                cand_t = sccd::min<T>(static_cast<T>(1), sccd::max<T>(0, cand_t));
                project_uv_simplex<T>(cand_u, cand_v);

                T fnext = 0;
                vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, cand_t, cand_u, cand_v, &fnext);

                if (fnext < best_f) {
                    best_t = cand_t;
                    best_u = cand_u;
                    best_v = cand_v;
                    best_f = fnext;
                    improved = true;
                    break;
                }

                alpha *= static_cast<T>(0.5);
            }

            t = best_t;
            u = best_u;
            v = best_v;
            f = best_f;

            const T norm_diff = norm_diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, t, u, v);
            if (norm_diff < atol) {
                return (u >= -atol && v >= -atol && u + v <= 1 + atol && t >= 0 && t <= 1);
            }

            if (!improved) {
                break;
            }
        }

        return false;
    }

    template <typename T>
    bool sum_less_than_one(const T u, const T v) {
        return u + v <= 1. / (1. - DBL_EPSILON);
    };

    template <typename T>
    struct Interval {
        T lower, upper;
        bool is_terminal() const { return lower >= upper; }
    };

    template <typename T>
    struct Box {
        using Interval = sccd::Interval<T>;
        Interval tuv[3];
        int depth{0};

        friend bool operator<(const Box &l, const Box &r) { return l.tuv[0].lower < r.tuv[0].lower; }

        Box() = default;
        Box(Interval t, Interval u, Interval v, int depth) : tuv{t, u, v}, depth(depth) {}
        bool is_terminal() const { return tuv[0].is_terminal() || tuv[1].is_terminal() || tuv[2].is_terminal(); }
        bool smaller_than_tol(const T tol0, const T tol1, const T tol2) const {
            return tuv[0].upper - tuv[0].lower <= tol0 && tuv[1].upper - tuv[1].lower <= tol1 &&
                   tuv[2].upper - tuv[2].lower <= tol2;
        }

        void print() const {
            std::cout << "Box: t: [" << tuv[0].lower << ", " << tuv[0].upper << "], u: [" << tuv[1].lower << ", "
                      << tuv[1].upper << "], v: [" << tuv[2].lower << ", " << tuv[2].upper << "], depth: " << depth
                      << std::endl;
        }

        bool is_at_depth_limit(const int max_iter) const { return depth >= max_iter; }

        int widest_dimension() const {
            const T dt = tuv[0].upper - tuv[0].lower;
            const T du = tuv[1].upper - tuv[1].lower;
            const T dv = tuv[2].upper - tuv[2].lower;
            if (du > dt && du >= dv) {
                return 1;
            } else if (dv > dt && dv > du) {
                return 2;
            }
            return 0;
        }

        bool bisect_vf(int split_dim, const T toi, std::vector<Box> &stack) const {
            std::pair<Interval, Interval> split_intervals{
                Interval{tuv[split_dim].lower, (tuv[split_dim].lower + tuv[split_dim].upper) * T(0.5)},
                Interval{(tuv[split_dim].lower + tuv[split_dim].upper) * T(0.5), tuv[split_dim].upper}};

            if (split_intervals.first.is_terminal() || split_intervals.second.is_terminal()) {
                return true;
            }

            stack.push_back(*this);
            stack.back().tuv[split_dim] = split_intervals.first;
            stack.back().depth++;

            if (split_dim == 0) {
                if (split_intervals.second.lower < toi) {
                    stack.push_back(*this);
                    stack.back().tuv[split_dim] = split_intervals.second;
                    stack.back().depth++;
                }
            } else {
                if (split_dim == 1) {
                    if (sum_less_than_one(split_intervals.first.lower, tuv[2].lower)) {
                        stack.push_back(*this);
                        stack.back().tuv[split_dim] = split_intervals.second;
                        stack.back().depth++;
                    }
                } else if (split_dim == 2) {
                    if (sum_less_than_one(split_intervals.second.lower, tuv[1].lower)) {
                        stack.push_back(*this);
                        stack.back().tuv[split_dim] = split_intervals.second;
                        stack.back().depth++;
                    }
                }
            }

            return false;
        }

        bool bisect_ee(int split_dim, const T toi, std::vector<Box> &stack) const {
            std::pair<Interval, Interval> split_intervals{
                Interval{tuv[split_dim].lower, (tuv[split_dim].lower + tuv[split_dim].upper) * T(0.5)},
                Interval{(tuv[split_dim].lower + tuv[split_dim].upper) * T(0.5), tuv[split_dim].upper}};

            if (split_intervals.first.is_terminal() || split_intervals.second.is_terminal()) {
                return true;
            }

            stack.push_back(*this);
            stack.back().tuv[split_dim] = split_intervals.first;
            stack.back().depth++;

            if (split_dim == 0) {
                if (split_intervals.second.lower < toi) {
                    stack.push_back(*this);
                    stack.back().tuv[split_dim] = split_intervals.second;
                    stack.back().depth++;
                }
            } else {
                stack.push_back(*this);
                stack.back().tuv[split_dim] = split_intervals.second;
                stack.back().depth++;
            }

            return false;
        }
    };

    template <typename T>
    bool find_root_bisection(const int max_iter,
                             const T tol,
                             const T sv[3],
                             const T s1[3],
                             const T s2[3],
                             const T s3[3],
                             const T ev[3],
                             const T e1[3],
                             const T e2[3],
                             const T e3[3],
                             T &t,
                             T &u,
                             T &v) {
        using Box = sccd::Box<T>;
        using Interval = sccd::Interval<T>;

        auto codomain_box = [=](const Box &domain, Box &codomain) -> void {
            codomain.tuv[0].lower = std::numeric_limits<T>::max();
            codomain.tuv[0].upper = std::numeric_limits<T>::lowest();
            codomain.tuv[1].lower = std::numeric_limits<T>::max();
            codomain.tuv[1].upper = std::numeric_limits<T>::lowest();
            codomain.tuv[2].lower = std::numeric_limits<T>::max();
            codomain.tuv[2].upper = std::numeric_limits<T>::lowest();

            for (int i = 0; i < 8; i++) {
                const T t = (i & 1) ? domain.tuv[0].upper : domain.tuv[0].lower;
                const T u = (i & 2) ? domain.tuv[1].upper : domain.tuv[1].lower;
                const T v = (i & 4) ? domain.tuv[2].upper : domain.tuv[2].lower;

                T F[3];
                diff_vf(sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, F);

                codomain.tuv[0].lower = sccd::min(codomain.tuv[0].lower, F[0]);
                codomain.tuv[0].upper = sccd::max(codomain.tuv[0].upper, F[0]);
                codomain.tuv[1].lower = sccd::min(codomain.tuv[1].lower, F[1]);
                codomain.tuv[1].upper = sccd::max(codomain.tuv[1].upper, F[1]);
                codomain.tuv[2].lower = sccd::min(codomain.tuv[2].lower, F[2]);
                codomain.tuv[2].upper = sccd::max(codomain.tuv[2].upper, F[2]);
            }
        };

        auto contains_origin = [&](const Box &box, T &true_tol, bool &inside_box) -> bool {
            Box codomain;
            codomain_box(box, codomain);

            for (int i = 0; i < 3; i++) {
                if (codomain.tuv[i].lower > tol || codomain.tuv[i].upper < -tol) {
                    return false;
                }
            }

            inside_box = true;
            for (int i = 0; i < 3; i++) {
                if (codomain.tuv[i].lower < tol || codomain.tuv[i].upper > -tol) {
                    inside_box = false;
                }
            }

            true_tol = sccd::max(
                sccd::max(codomain.tuv[0].upper - codomain.tuv[0].lower, codomain.tuv[1].upper - codomain.tuv[1].lower),
                codomain.tuv[2].upper - codomain.tuv[2].lower);
            return true;
        };

        // Compute per-axis tolerances (matching snumtol.hpp signature)
        T tol0 = tol, tol1 = tol, tol2 = tol;
        compute_face_vertex_tolerance_soa<T>(tol,
                                             sv[0],
                                             sv[1],
                                             sv[2],
                                             s1[0],
                                             s1[1],
                                             s1[2],
                                             s2[0],
                                             s2[1],
                                             s2[2],
                                             s3[0],
                                             s3[1],
                                             s3[2],
                                             ev[0],
                                             ev[1],
                                             ev[2],
                                             e1[0],
                                             e1[1],
                                             e1[2],
                                             e2[0],
                                             e2[1],
                                             e2[2],
                                             e3[0],
                                             e3[1],
                                             e3[2],
                                             &tol0,
                                             &tol1,
                                             &tol2);

        // printf("tol %f -> tol0: %f, tol1: %f, tol2: %f\n", tol, tol0, tol1, tol2);

        std::vector<Box> stack;
        stack.reserve(1024);
        stack.push_back(Box(Interval{T(0), T(1)}, Interval{T(0), T(1)}, Interval{T(0), T(1)}, 0));

        T toi = 1;

        bool found_root = false;
        while (!stack.empty()) {
            Box box = stack.back();
            stack.pop_back();

            if (box.tuv[0].lower > toi) {
                continue;
            }

            T min_t = sccd::min(toi, box.tuv[0].lower);

            T true_tol = tol;
            bool inside_box = false;
            if (contains_origin(box, true_tol, inside_box)) {
                // Condition 1: the domain is smaller than the tolerance.
                if (box.smaller_than_tol(tol0, tol1, tol2)) {
                    t = box.tuv[0].lower;
                    u = box.tuv[1].lower;
                    v = box.tuv[2].lower;
                    toi = sccd::min(toi, min_t);
                    found_root = true;
                    continue;
                }

                // Condition 2: the box is inside the epsilon box
                if (inside_box) {
                    t = box.tuv[0].lower;
                    u = box.tuv[1].lower;
                    v = box.tuv[2].lower;
                    toi = sccd::min(toi, min_t);
                    found_root = true;
                    continue;
                }

                // Condition 3: real tolerance is smaller than the int tolerance
                if (true_tol < tol && box.tuv[0].lower > 0) {
                    t = box.tuv[0].lower;
                    u = box.tuv[1].lower;
                    v = box.tuv[2].lower;
                    toi = sccd::min(toi, min_t);
                    found_root = true;
                    continue;
                }

                if (box.is_terminal()) {
                    toi = sccd::min(toi, min_t);
                    found_root = true;
                    continue;
                }

                if (box.depth > max_iter) continue;

                // Split the box along the widest dimension
                int split_dim = box.widest_dimension();
                if (box.bisect_vf(split_dim, toi, stack)) {
                    // Split box too small
                    t = box.tuv[0].lower;
                    u = box.tuv[1].lower;
                    v = box.tuv[2].lower;
                    toi = sccd::min(toi, min_t);
                    found_root = true;
                    continue;
                }
            }
        }

        t = toi;
        return found_root;
    }

    template <int NT, int NU, int NV, typename T>
    inline static void grid_sample_F_vf(const T start_t,
                                        const T start_u,
                                        const T start_v,
                                        const T ht,
                                        const T hu,
                                        const T hv,
                                        const T sv,
                                        const T ev,
                                        const T s1,
                                        const T s2,
                                        const T s3,
                                        const T e1,
                                        const T e2,
                                        const T e3,
                                        T *const SCCD_RESTRICT F) {
        static constexpr int STRIDE_T = (NU + 1) * (NV + 1);
        static constexpr int STRIDE_U = (NV + 1);

        for (int a = 0; a <= NT; a++) {
            for (int b = 0; b <= NU; b++) {
                for (int c = 0; c <= NV; c++) {
                    const int idx = a * STRIDE_T + b * STRIDE_U + c;
                    const T t = start_t + a * ht;
                    const T u = start_u + b * hu;
                    const T v = start_v + c * hv;

                    const T vertex = (ev - sv) * t + sv;
                    const T t0 = (e1 - s1) * t + s1;
                    const T t1 = (e2 - s2) * t + s2;
                    const T t2 = (e3 - s3) * t + s3;

                    const T face = (t1 - t0) * u + (t2 - t0) * v + t0;
                    F[idx] = vertex - face;

                    // T t0 = (1 - t);
                    // T t1 = t;
                    // T o = (1-u-v);
                    // T v_pos = t0 * sv + t1 * ev;
                    // T f0 = t0 * (o * s1 + u * s2 + v * s3);
                    // T f1 = t1 * (o * e1 + u * e2 + v * e3);
                    // T f = f0 + f1;
                    // T diff = v_pos - f;
                    // F[idx] = diff;
                }
            }
        }
    }

    template <int NT, int NU, int NV, typename T>
    inline static void grid_zero_and_accept(const T *const SCCD_RESTRICT F,
                                            const T tol,
                                            const T adaptive_tol,
                                            uint8_t *const SCCD_RESTRICT contains_origin,
                                            uint8_t *const SCCD_RESTRICT accept) {
        static constexpr int STIDE_T = (NU + 1) * (NV + 1);
        static constexpr int STIDE_U = (NV + 1);

        for (int a = 0; a < NT; a++) {
            for (int b = 0; b < NU; b++) {
                const int i0 = a * STIDE_T + b * STIDE_U;
                const int i1 = a * STIDE_T + b * STIDE_U + 1;
                const int i2 = a * STIDE_T + (b + 1) * STIDE_U;
                const int i3 = a * STIDE_T + (b + 1) * STIDE_U + 1;
                const int i4 = (a + 1) * STIDE_T + b * STIDE_U;
                const int i5 = (a + 1) * STIDE_T + b * STIDE_U + 1;
                const int i6 = (a + 1) * STIDE_T + (b + 1) * STIDE_U;
                const int i7 = (a + 1) * STIDE_T + (b + 1) * STIDE_U + 1;

                const T *const SCCD_RESTRICT F000 = &F[i0];
                const T *const SCCD_RESTRICT F001 = &F[i1];
                const T *const SCCD_RESTRICT F010 = &F[i2];
                const T *const SCCD_RESTRICT F011 = &F[i3];
                const T *const SCCD_RESTRICT F100 = &F[i4];
                const T *const SCCD_RESTRICT F101 = &F[i5];
                const T *const SCCD_RESTRICT F110 = &F[i6];
                const T *const SCCD_RESTRICT F111 = &F[i7];

                const int cell_offset = a * NU * NV + b * NV;
                uint8_t *const SCCD_RESTRICT contains_origin_cell = &contains_origin[cell_offset];
                uint8_t *const SCCD_RESTRICT accept_cell = &accept[cell_offset];

                for (int c = 0; c < NV; c++) {
                    const T fmin = sccd::min(sccd::min(sccd::min(F000[c], F001[c]), sccd::min(F010[c], F011[c])),
                                             sccd::min(sccd::min(F100[c], F101[c]), sccd::min(F110[c], F111[c])));

                    const T fmax = sccd::max(sccd::max(sccd::max(F000[c], F001[c]), sccd::max(F010[c], F011[c])),
                                             sccd::max(sccd::max(F100[c], F101[c]), sccd::max(F110[c], F111[c])));

                    contains_origin_cell[c] &= (fmin <= tol) & (fmax >= -tol);  // AND) for all dims
                    bool cond1 = (fmax - fmin <= adaptive_tol);    // AND) The domain is smaller than the tolerance.
                    bool cond2 = !((fmin < tol) | (fmax > -tol));  // AND) The box is inside the epsilon box
                    bool cond3 = (fmax - fmin < tol);  // OR) Real tolerance is smaller than the int tolerance
                    bool cond4 = (fmin >= fmax);       // AND) The interval is terminal

                    uint8_t cond_mask = (cond1 ? (1 & accept_cell[c]) : 0);
                    cond_mask |= (cond2 ? (2 & accept_cell[c]) : 0);
                    cond_mask |= (cond3 ? 4 : 0);
                    cond_mask |= (cond4 ? (8 & accept_cell[c]) : 0);
                    accept_cell[c] = cond_mask & ((fmin <= tol) & ((fmax >= -tol) ? 0xf : 0));
                }
            }
        }
    }


    template <int NT, int NU, int NV, typename T>
    inline bool grid_search_vf(const sccd::Box<T> &domain,
                               const int max_iter,
                               const T tol,
                               const T tols[3],
                               const T sv[3],
                               const T s1[3],
                               const T s2[3],
                               const T s3[3],
                               const T ev[3],
                               const T e1[3],
                               const T e2[3],
                               const T e3[3],
                               T &toi,  // In/Out
                               T &u,
                               T &v,
                               std::vector<sccd::Box<T>> &stack) {
        static constexpr int N_nodes = (NT + 1) * (NU + 1) * (NV + 1);
        static constexpr int N_cells = NT * NU * NV;
        static constexpr int STRIDE_T = (NU + 1) * (NV + 1);
        static constexpr int STRIDE_U = (NV + 1);

        const T t_min = domain.tuv[0].lower;
        const T u_min = domain.tuv[1].lower;
        const T v_min = domain.tuv[2].lower;
        const T t_max = domain.tuv[0].upper;
        const T u_max = domain.tuv[1].upper;
        const T v_max = domain.tuv[2].upper;
        const T t_h = (t_max - t_min) / NT;
        const T u_h = (u_max - u_min) / NU;
        const T v_h = (v_max - v_min) / NV;

        // 1) Generate F_grid
        T F[3][N_nodes];
        grid_sample_F_vf<NT, NU, NV, T>(
            t_min, u_min, v_min, t_h, u_h, v_h, sv[0], ev[0], s1[0], s2[0], s3[0], e1[0], e2[0], e3[0], F[0]);

        grid_sample_F_vf<NT, NU, NV, T>(
            t_min, u_min, v_min, t_h, u_h, v_h, sv[1], ev[1], s1[1], s2[1], s3[1], e1[1], e2[1], e3[1], F[1]);

        grid_sample_F_vf<NT, NU, NV, T>(
            t_min, u_min, v_min, t_h, u_h, v_h, sv[2], ev[2], s1[2], s2[2], s3[2], e1[2], e2[2], e3[2], F[2]);

        // 2) Find cells containing zeros and check for acceptability
        uint8_t contains_zero_and_refine[N_cells];
        uint8_t accept[N_cells];
        for (int i = 0; i < N_cells; i++) {
            contains_zero_and_refine[i] = true;
            accept[i] = 0xf;
        }

        grid_zero_and_accept<NT, NU, NV, T>(F[0], tol, tols[0], contains_zero_and_refine, accept);
        grid_zero_and_accept<NT, NU, NV, T>(F[1], tol, tols[1], contains_zero_and_refine, accept);
        grid_zero_and_accept<NT, NU, NV, T>(F[2], tol, tols[2], contains_zero_and_refine, accept);

        bool found = false;
        // 3) Find earilest toi and schedule for refinement
        for (int a = 0; a < NT; a++) {
            const T t0 = t_min + a * t_h;
            if (t0 > toi) continue;

            for (int b = 0; b < NU; b++) {
                for (int c = 0; c < NV; c++) {
                    const int i = a * NU * NV + b * NV + c;
                    if (accept[i]) {
                        toi = t0;
                        u = u_min + b * u_h;
                        v = v_min + c * v_h;
                        found = true;
                        contains_zero_and_refine[i] = false;
                    }
                }
            }
        }

        // Create new boxes
        for (int a = 0; a < NT; a++) {
            for (int b = 0; b < NU; b++) {
                for (int c = 0; c < NV; c++) {
                    const int i = a * NU * NV + b * NV + c;
                    if (contains_zero_and_refine[i] && !accept[i]) {
                        const T tt_min = t_min + a * t_h;
                        const T uu_min = u_min + b * u_h;
                        const T vv_min = v_min + c * v_h;

                        const T tt_max = t_min + (a + 1) * t_h;
                        const T uu_max = u_min + (b + 1) * u_h;
                        const T vv_max = v_min + (c + 1) * v_h;

                        if (uu_min + vv_min >= 1 + tol || tt_min >= toi) {
                            continue;
                        }

                        Box<T> box({tt_min, tt_max}, {uu_min, uu_max}, {vv_min, vv_max}, domain.depth + 1);
                        if (box.depth > max_iter) {
                            // Conservative approximation
                            const T approx = box.tuv[0].lower;
                            if (approx < toi) {
                                toi = approx;
                                u = box.tuv[1].lower;
                                v = box.tuv[2].lower;
                                found = true;
                            }
                            continue;
                        }

                        int split_dim = box.widest_dimension();
                        box.bisect_vf(split_dim, toi, stack);
                    }
                }
            }
        }

        return found;
    }

    template <typename T>
    bool find_root_grid_vf(const int max_iter,
                           const T tol,
                           const T sv[3],
                           const T s1[3],
                           const T s2[3],
                           const T s3[3],
                           const T ev[3],
                           const T e1[3],
                           const T e2[3],
                           const T e3[3],
                           T &t,
                           T &u,
                           T &v,
                           std::vector<Box<T>> &stack) {
        using Box = sccd::Box<T>;
        using Interval = sccd::Interval<T>;

        // Compute per-axis tolerances (matching snumtol.hpp signature)
        T tols[3];
        compute_face_vertex_tolerance_soa<T>(tol,
                                             sv[0],
                                             sv[1],
                                             sv[2],
                                             s1[0],
                                             s1[1],
                                             s1[2],
                                             s2[0],
                                             s2[1],
                                             s2[2],
                                             s3[0],
                                             s3[1],
                                             s3[2],
                                             ev[0],
                                             ev[1],
                                             ev[2],
                                             e1[0],
                                             e1[1],
                                             e1[2],
                                             e2[0],
                                             e2[1],
                                             e2[2],
                                             e3[0],
                                             e3[1],
                                             e3[2],
                                             &tols[0],
                                             &tols[1],
                                             &tols[2]);

        t = 1.1;
        u = 0;
        v = 0;

        bool found = false;
        stack.clear();
        stack.push_back(Box(Interval{T(0), T(1)}, Interval{T(0), T(1)}, Interval{T(0), T(1)}, 0));
        bool found_root = false;
        while (!stack.empty()) {
            Box box = stack.back();
            stack.pop_back();

            if (box.tuv[0].lower >= t) {
                continue;
            }

            found |= grid_search_vf<4, 4, 4, T>(box, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack);
        }

        return found;
    }

    template <typename T>
    inline bool find_root_grid_rotate_vf(
                               const int max_iter,
                               const T tol,
                               const T sv[3],
                               const T s1[3],
                               const T s2[3],
                               const T s3[3],
                               const T ev[3],
                               const T e1[3],
                               const T e2[3],
                               const T e3[3],
                               T &toi,  // In/Out
                               T &u,
                               T &v,
                               std::vector<sccd::Box<T>> &stack) 
{

        static constexpr T EPS_LEN = static_cast<T>(1e-12);
        static constexpr T EPS_ANG = static_cast<T>(1e-12);

        const T d[3] = {ev[0] - sv[0], ev[1] - sv[1], ev[2] - sv[2]};
        const T dlen2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        const T dlen = std::sqrt(sccd::max<T>(dlen2, T(0)));

        if (dlen < EPS_LEN) {
            return find_root_grid_vf<T>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack);
        }

        T d_hat[3];
        const T inv_len = T(1) / dlen;
        d_hat[0] = d[0] * inv_len;
        d_hat[1] = d[1] * inv_len;
        d_hat[2] = d[2] * inv_len;


        T R[3][3];
        const T sss = d_hat[1] * d_hat[1] + d_hat[2] * d_hat[2];
        const T c = d_hat[0];

        if (sss < EPS_ANG && c < T(0)) {
            const bool use_y = sccd::abs(d_hat[0]) > sccd::abs(d_hat[1]);
            const T axis_seed[3] = {T(0), use_y ? T(1) : T(0), use_y ? T(0) : T(1)};

            T axis[3];
            axis[0] = d_hat[1] * axis_seed[2] - d_hat[2] * axis_seed[1];
            axis[1] = d_hat[2] * axis_seed[0] - d_hat[0] * axis_seed[2];
            axis[2] = d_hat[0] * axis_seed[1] - d_hat[1] * axis_seed[0];

            const T axis_norm2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
            const T inv_axis_norm = T(1) / std::sqrt(sccd::max<T>(axis_norm2, EPS_ANG));

            axis[0] *= inv_axis_norm;
            axis[1] *= inv_axis_norm;
            axis[2] *= inv_axis_norm;

            const T a0a0 = axis[0] * axis[0];
            const T a0a1 = axis[0] * axis[1];
            const T a0a2 = axis[0] * axis[2];
            const T a1a1 = axis[1] * axis[1];
            const T a1a2 = axis[1] * axis[2];
            const T a2a2 = axis[2] * axis[2];

            R[0][0] = T(2) * a0a0 - T(1);
            R[0][1] = T(2) * a0a1;
            R[0][2] = T(2) * a0a2;
            R[1][0] = R[0][1];
            R[1][1] = T(2) * a1a1 - T(1);
            R[1][2] = T(2) * a1a2;
            R[2][0] = R[0][2];
            R[2][1] = R[1][2];
            R[2][2] = T(2) * a2a2 - T(1);
        } else if (sss < EPS_ANG) {
            R[0][0] = T(1);
            R[0][1] = T(0);
            R[0][2] = T(0);
            R[1][0] = T(0);
            R[1][1] = T(1);
            R[1][2] = T(0);
            R[2][0] = T(0);
            R[2][1] = T(0);
            R[2][2] = T(1);
        } else {
            const T factor = (T(1) - c) / sccd::max<T>(sss, EPS_ANG);
            const T s2_factor = sss * factor;
            const T xy = d_hat[1] * d_hat[2] * (-factor);
            const T yy = d_hat[1] * d_hat[1] * factor;
            const T zz = d_hat[2] * d_hat[2] * factor;

            R[0][0] = T(1) - s2_factor;
            R[0][1] = d_hat[1];
            R[0][2] = d_hat[2];

            R[1][0] = -d_hat[1];
            R[1][1] = T(1) - yy;
            R[1][2] = xy;

            R[2][0] = -d_hat[2];
            R[2][1] = xy;
            R[2][2] = T(1) - zz;
        }

        const T scale = T(1) / dlen;

        auto transform_point = [&](const T in[3], T out[3]) {
            const T x = in[0] - sv[0];
            const T y = in[1] - sv[1];
            const T z = in[2] - sv[2];

            out[0] = scale * (R[0][0] * x + R[0][1] * y + R[0][2] * z);
            out[1] = scale * (R[1][0] * x + R[1][1] * y + R[1][2] * z);
            out[2] = scale * (R[2][0] * x + R[2][1] * y + R[2][2] * z);
        };

        T sv_r[3], s1_r[3], s2_r[3], s3_r[3], ev_r[3], e1_r[3], e2_r[3], e3_r[3];
        transform_point(sv, sv_r);
        transform_point(s1, s1_r);
        transform_point(s2, s2_r);
        transform_point(s3, s3_r);
        transform_point(ev, ev_r);
        transform_point(e1, e1_r);
        transform_point(e2, e2_r);
        transform_point(e3, e3_r);

        const T min_y = sccd::min(sccd::min(sccd::min(s1_r[1], s2_r[1]), sccd::min(s3_r[1], e1_r[1])),
                                  sccd::min(e2_r[1], e3_r[1]));
        const T max_y = sccd::max(sccd::max(sccd::max(s1_r[1], s2_r[1]), sccd::max(s3_r[1], e1_r[1])),
                                  sccd::max(e2_r[1], e3_r[1]));
        const T min_z = sccd::min(sccd::min(sccd::min(s1_r[2], s2_r[2]), sccd::min(s3_r[2], e1_r[2])),
                                  sccd::min(e2_r[2], e3_r[2]));
        const T max_z = sccd::max(sccd::max(sccd::max(s1_r[2], s2_r[2]), sccd::max(s3_r[2], e1_r[2])),
                                  sccd::max(e2_r[2], e3_r[2]));

        if ((min_y > T(0) || max_y < T(0)) || (min_z > T(0) || max_z < T(0))) {
            // Early skip good for TTS
            return false;
        }

        return find_root_grid_vf<T>(max_iter, tol, sv_r, s1_r, s2_r, s3_r, ev_r, e1_r, e2_r, e3_r, toi, u, v, stack);
}


    template <int NT, int NU, int NV, typename T>
    inline static void grid_sample_F_ee(const T start_t,
                                        const T start_u,
                                        const T start_v,
                                        const T ht,
                                        const T hu,
                                        const T hv,
                                        const T s1,
                                        const T s2,
                                        const T s3,
                                        const T s4,
                                        const T e1,
                                        const T e2,
                                        const T e3,
                                        const T e4,
                                        T *const SCCD_RESTRICT F) {
        static constexpr int STRIDE_T = (NU + 1) * (NV + 1);
        static constexpr int STRIDE_U = (NV + 1);

        for (int a = 0; a <= NT; a++) {
            for (int b = 0; b <= NU; b++) {
                for (int c = 0; c <= NV; c++) {
                    const int idx = a * STRIDE_T + b * STRIDE_U + c;
                    const T t = start_t + a * ht;
                    const T u = start_u + b * hu;
                    const T v = start_v + c * hv;

                    const T ea0 = (e1 - s1) * t + s1;
                    const T ea1 = (e2 - s2) * t + s2;
                    const T eb0 = (e3 - s3) * t + s3;
                    const T eb1 = (e4 - s4) * t + s4;
                    F[idx] = ((ea1 - ea0) * u + ea0) - ((eb1 - eb0) * v + eb0);
                }
            }
        }
    }

    template <int NT, int NU, int NV, typename T>
    inline bool grid_search_ee(const sccd::Box<T> &domain,
                               const int max_iter,
                               const T tol,
                               const T tols[3],
                               const T s1[3],
                               const T s2[3],
                               const T s3[3],
                               const T s4[3],
                               const T e1[3],
                               const T e2[3],
                               const T e3[3],
                               const T e4[3],
                               T &toi,  // In/Out
                               T &u,
                               T &v,
                               std::vector<sccd::Box<T>> &stack) {
        static constexpr int N_nodes = (NT + 1) * (NU + 1) * (NV + 1);
        static constexpr int N_cells = NT * NU * NV;
        static constexpr int STRIDE_T = (NU + 1) * (NV + 1);
        static constexpr int STRIDE_U = (NV + 1);

        const T t_min = domain.tuv[0].lower;
        const T u_min = domain.tuv[1].lower;
        const T v_min = domain.tuv[2].lower;
        const T t_max = domain.tuv[0].upper;
        const T u_max = domain.tuv[1].upper;
        const T v_max = domain.tuv[2].upper;
        const T t_h = (t_max - t_min) / NT;
        const T u_h = (u_max - u_min) / NU;
        const T v_h = (v_max - v_min) / NV;

        // 1) Generate F_grid
        T F[3][N_nodes];
        grid_sample_F_ee<NT, NU, NV, T>(
            t_min, u_min, v_min, t_h, u_h, v_h, s1[0], s2[0], s3[0], s4[0], e1[0], e2[0], e3[0], e4[0], F[0]);

        grid_sample_F_ee<NT, NU, NV, T>(
            t_min, u_min, v_min, t_h, u_h, v_h, s1[1], s2[1], s3[1], s4[1], e1[1], e2[1], e3[1], e4[1], F[1]);

        grid_sample_F_ee<NT, NU, NV, T>(
            t_min, u_min, v_min, t_h, u_h, v_h, s1[2], s2[2], s3[2], s4[2], e1[2], e2[2], e3[2], e4[2], F[2]);

        // 2) Find cells containing zeros and check for acceptability
        uint8_t contains_zero_and_refine[N_cells];
        uint8_t accept[N_cells];
        for (int i = 0; i < N_cells; i++) {
            contains_zero_and_refine[i] = true;
            accept[i] = 0xf;
        }

        grid_zero_and_accept<NT, NU, NV, T>(F[0], tol, tols[0], contains_zero_and_refine, accept);
        grid_zero_and_accept<NT, NU, NV, T>(F[1], tol, tols[1], contains_zero_and_refine, accept);
        grid_zero_and_accept<NT, NU, NV, T>(F[2], tol, tols[2], contains_zero_and_refine, accept);

        bool found = false;
        // 3) Find earilest toi and schedule for refinement
        for (int a = 0; a < NT; a++) {
            const T t0 = t_min + a * t_h;
            if (t0 > toi) continue;

            for (int b = 0; b < NU; b++) {
                for (int c = 0; c < NV; c++) {
                    const int i = a * NU * NV + b * NV + c;
                    if (accept[i]) {
                        toi = t0;
                        u = u_min + b * u_h;
                        v = v_min + c * v_h;
                        found = true;
                        contains_zero_and_refine[i] = false;
                    }
                }
            }
        }

        // Create new boxes
        for (int a = 0; a < NT; a++) {
            for (int b = 0; b < NU; b++) {
                for (int c = 0; c < NV; c++) {
                    const int i = a * NU * NV + b * NV + c;
                    if (contains_zero_and_refine[i] && !accept[i]) {
                        const T tt_min = t_min + a * t_h;
                        const T uu_min = u_min + b * u_h;
                        const T vv_min = v_min + c * v_h;

                        const T tt_max = t_min + (a + 1) * t_h;
                        const T uu_max = u_min + (b + 1) * u_h;
                        const T vv_max = v_min + (c + 1) * v_h;

                        if (tt_min >= toi) {
                            continue;
                        }

                        Box<T> box({tt_min, tt_max}, {uu_min, uu_max}, {vv_min, vv_max}, domain.depth + 1);
                        if (box.depth > max_iter) {
                            // Conservative approximation
                            const T approx = box.tuv[0].lower;
                            if (approx < toi) {
                                toi = approx;
                                u = box.tuv[1].lower;
                                v = box.tuv[2].lower;
                                found = true;
                            }
                            continue;
                        }

                        int split_dim = box.widest_dimension();
                        box.bisect_ee(split_dim, toi, stack);
                    }
                }
            }
        }

        return found;
    }

    template <typename T>
    bool find_root_grid_ee(const int max_iter,
                           const T tol,
                           const T s1[3],
                           const T s2[3],
                           const T s3[3],
                           const T s4[3],
                           const T e1[3],
                           const T e2[3],
                           const T e3[3],
                           const T e4[3],
                           T &t,
                           T &u,
                           T &v,
                           std::vector<Box<T>> &stack) {
        using Box = sccd::Box<T>;
        using Interval = sccd::Interval<T>;

        // Compute per-axis tolerances (matching snumtol.hpp signature)
        T tols[3];
        compute_edge_edge_tolerance_soa<T>(tol,
                                           s1[0],
                                           s1[1],
                                           s1[2],
                                           s2[0],
                                           s2[1],
                                           s2[2],
                                           s3[0],
                                           s3[1],
                                           s3[2],
                                           s4[0],
                                           s4[1],
                                           s4[2],
                                           e1[0],
                                           e1[1],
                                           e1[2],
                                           e2[0],
                                           e2[1],
                                           e2[2],
                                           e3[0],
                                           e3[1],
                                           e3[2],
                                           e4[0],
                                           e4[1],
                                           e4[2],
                                           &tols[0],
                                           &tols[1],
                                           &tols[2]);

        t = 1.1;
        u = 0;
        v = 0;

        bool found = false;
        stack.clear();
        stack.push_back(Box(Interval{T(0), T(1)}, Interval{T(0), T(1)}, Interval{T(0), T(1)}, 0));
        bool found_root = false;
        while (!stack.empty()) {
            Box box = stack.back();
            stack.pop_back();

            if (box.tuv[0].lower >= t) {
                continue;
            }

            found |= grid_search_ee<4, 4, 4, T>(box, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack);
        }

        return found;
    }
}  // namespace sccd

#endif  // S_ROOT_FINDER_HPP
