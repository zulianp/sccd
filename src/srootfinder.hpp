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

#define ADAPTIVE_NUM_SPLITS 2
#define UNIFORM_NUM_SPLITS 4

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

        ticcd::Array3 err(-1, -1, -1);

        ticcd::Scalar ms = 0;
        ticcd::Scalar max_time = 1;
        // ticcd::Scalar output_tolerance = 1e-6;
        ticcd::Scalar output_tolerance = atol;
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
                                    t,
                                    atol,
                                    1,
                                    max_iter,
                                    output_tolerance,
                                    no_zero_toi,
                                    // ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);
                                    ticcd::CCDRootFindingMethod::DEPTH_FIRST_SEARCH);

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
        ticcd::Array3 err(-1, -1, -1);

        ticcd::Scalar ms = 0;
        ticcd::Scalar max_time = 1;
        // ticcd::Scalar output_tolerance = 1e-6;
        ticcd::Scalar output_tolerance = atol;
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
                                  t,
                                  atol,
                                  1,
                                  max_iter,
                                  output_tolerance,
                                  no_zero_toi,
                                  //   ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);
                                  ticcd::CCDRootFindingMethod::DEPTH_FIRST_SEARCH);
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

            // // NEW
            // if (split_dim == 0) {
            //     split_intervals.first.lower = std::min(split_intervals.first.lower, toi);
            //     split_intervals.second.lower = std::min(split_intervals.second.lower, toi);
            // }

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

            // // NEW
            // if (split_dim == 0) {
            //     split_intervals.first.lower = std::min(split_intervals.first.lower, toi);
            //     split_intervals.second.lower = std::min(split_intervals.second.lower, toi);
            // }

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
    inline Box<T> unit_domain_box() {
        using Interval = sccd::Interval<T>;
        return Box<T>(Interval{T(0), T(1)}, Interval{T(0), T(1)}, Interval{T(0), T(1)}, 0);
    }

    template <int SplitDim, typename T>
    inline Box<T> split_axis_box(const Box<T> &domain, const T sample_min, const T sample_max) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        return Box<T>(Interval<T>{SplitDim == 0 ? sample_min : domain.tuv[0].lower,
                                  SplitDim == 0 ? sample_max : domain.tuv[0].upper},
                      Interval<T>{SplitDim == 1 ? sample_min : domain.tuv[1].lower,
                                  SplitDim == 1 ? sample_max : domain.tuv[1].upper},
                      Interval<T>{SplitDim == 2 ? sample_min : domain.tuv[2].lower,
                                  SplitDim == 2 ? sample_max : domain.tuv[2].upper},
                      domain.depth + 1);
    }

    template <typename T>
    inline void init_codomain_bounds(T fmin[3], T fmax[3]) {
        for (int d = 0; d < 3; ++d) {
            fmin[d] = std::numeric_limits<T>::max();
            fmax[d] = std::numeric_limits<T>::lowest();
        }
    }

    template <typename T>
    inline void update_codomain_bounds(const T F[3], T fmin[3], T fmax[3]) {
        for (int d = 0; d < 3; ++d) {
            fmin[d] = sccd::min<T>(fmin[d], F[d]);
            fmax[d] = sccd::max<T>(fmax[d], F[d]);
        }
    }

    template <typename T>
    inline bool codomain_acceptance(const T fmin[3], const T fmax[3], const T tol, const T tols[3], bool &accept) {
        bool contains_zero = true;
        bool smaller_than_axis_tol = true;
        bool inside_epsilon_box = true;
        bool smaller_than_scalar_tol = true;
        bool degenerate_interval = true;

        const T eps = std::numeric_limits<T>::epsilon();

        for (int d = 0; d < 3; ++d) {
            const T interval_width = fmax[d] - fmin[d];
            contains_zero = contains_zero && (fmin[d] <= eps) && (fmax[d] >= -eps);
            smaller_than_axis_tol = smaller_than_axis_tol && (interval_width <= tols[d]);
            inside_epsilon_box = inside_epsilon_box && (fmin[d] >= -tol) && (fmax[d] <= tol);
            smaller_than_scalar_tol = smaller_than_scalar_tol && (interval_width < tol);
            degenerate_interval = degenerate_interval && (fmin[d] >= fmax[d]);
        }

        accept = contains_zero && (smaller_than_axis_tol || inside_epsilon_box || smaller_than_scalar_tol || degenerate_interval);
        return contains_zero;
    }

    // template <typename T>
    // inline bool codomain_acceptance(const T fmin[3], const T fmax[3], const T tol, const T tols[3], bool &accept) {
    //     // Replicates predicates of TI
    //     accept = true;
    //     bool contains_zero = true;

    //     for (int d = 0; d < 3; ++d) {
    //         contains_zero = contains_zero &&  //
    //                         (fmin[d] <= tols[d]) && (fmax[d] >= -tols[d]);
    //         accept = accept && ((fmin[d] >= -tols[d]) && (fmax[d] <= tols[d]));
    //     }

    //     accept = contains_zero && accept;
    //     return contains_zero;
    // }

    template <typename T>
    inline bool accept_grid_root_vf(const Box<T> &box,
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
                                    T &toi,
                                    T &u,
                                    T &v,
                                    const bool refine) {
        T t_approx = box.tuv[0].lower;
        if (t_approx < toi && box.tuv[1].lower + box.tuv[2].lower < T(1) + tols[1] + tols[2]) {
            T u_approx = box.tuv[1].lower;
            T v_approx = box.tuv[2].lower;

            if (refine) {
                const bool refined =
                    find_root_newton<T>(40, tol, sv, s1, s2, s3, ev, e1, e2, e3, t_approx, u_approx, v_approx);

                if (refined && t_approx < toi) {
                    toi = sccd::min<T>(box.tuv[0].upper, sccd::max<T>(box.tuv[0].lower, T(0.99) * t_approx));
                    u = u_approx;
                    v = v_approx;
                    return true;
                }
            } else {
                toi = t_approx;
                u = u_approx;
                v = v_approx;
                return true;
            }
        }
        return false;
    }

    template <typename T>
    inline bool accept_grid_root_ee(const Box<T> &box, T &toi, T &u, T &v) {
        const T t_approx = box.tuv[0].lower;
        if (t_approx < toi) {
            toi = t_approx;
            u = box.tuv[1].lower;
            v = box.tuv[2].lower;
            return true;
        }
        return false;
    }

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
                             const Box<T> &initial_domain,
                             T &t,
                             T &u,
                             T &v) {
        using Box = sccd::Box<T>;

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
        T axis_tol[3];
        compute_face_vertex_tolerance<T>(tol, sv, s1, s2, s3, ev, e1, e2, e3, axis_tol);
        const T tol0 = axis_tol[0];
        const T tol1 = axis_tol[1];
        const T tol2 = axis_tol[2];

        // printf("tol %f -> tol0: %f, tol1: %f, tol2: %f\n", tol, tol0, tol1, tol2);

        std::vector<Box> stack;
        stack.reserve(1024);
        stack.push_back(initial_domain);

        T toi = t;

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
        return find_root_bisection<T>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, unit_domain_box<T>(), t, u, v);
    }

    template <int SplitDim, int N, typename T>
    inline static void normal_equation_axis_splitters_vf(const Box<T> &domain,
                                                         const T sv[3],
                                                         const T s1[3],
                                                         const T s2[3],
                                                         const T s3[3],
                                                         const T ev[3],
                                                         const T e1[3],
                                                         const T e2[3],
                                                         const T e3[3],
                                                         T *const SCCD_RESTRICT splitters) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);

        const T lo = domain.tuv[SplitDim].lower;
        const T hi = domain.tuv[SplitDim].upper;
        const T h = (hi - lo) / T(N + 1);

        T damping = 0.45;
        if constexpr (N == 1) {
            damping = 0.6;
        }

        const T radius = h * damping;
        const T mid_t = (domain.tuv[0].lower + domain.tuv[0].upper) * T(0.5);
        const T mid_u = (domain.tuv[1].lower + domain.tuv[1].upper) * T(0.5);
        const T mid_v = (domain.tuv[2].lower + domain.tuv[2].upper) * T(0.5);
        const T eps = std::numeric_limits<T>::epsilon();

        T F_base[3];
        T J_axis[3];
        T H_axis = T(0);

        const T base_t = SplitDim == 0 ? T(0) : mid_t;
        const T base_u = SplitDim == 1 ? T(0) : mid_u;
        const T base_v = SplitDim == 2 ? T(0) : mid_v;
        const T base_omt = T(1) - base_t;
        const T base_o = T(1) - base_u - base_v;

        for (int d = 0; d < 3; ++d) {
            const T vertex = base_omt * sv[d] + base_t * ev[d];
            const T face = base_omt * (base_o * s1[d] + base_u * s2[d] + base_v * s3[d]) +
                           base_t * (base_o * e1[d] + base_u * e2[d] + base_v * e3[d]);
            F_base[d] = vertex - face;

            if constexpr (SplitDim == 0) {
                const T o = T(1) - mid_u - mid_v;
                J_axis[d] = (ev[d] - sv[d]) - (o * (e1[d] - s1[d]) + mid_u * (e2[d] - s2[d]) + mid_v * (e3[d] - s3[d]));
            } else if constexpr (SplitDim == 1) {
                J_axis[d] = -((T(1) - mid_t) * (s2[d] - s1[d]) + mid_t * (e2[d] - e1[d]));
            } else {
                J_axis[d] = -((T(1) - mid_t) * (s3[d] - s1[d]) + mid_t * (e3[d] - e1[d]));
            }
            H_axis += J_axis[d] * J_axis[d];
        }
        const T step_scale = H_axis > eps ? T(1) / H_axis : T(0.00001);

#pragma omp simd aligned(splitters)
        for (int i = 0; i < N; ++i) {
            const T x0 = lo + h * T(i + 1);
            auto xmin = sccd::max<T>(lo, x0 - radius);
            auto xmax = sccd::min<T>(hi, x0 + radius);

            T g = T(0);
            for (int d = 0; d < 3; ++d) {
                const T J = J_axis[d];
                g += (F_base[d] + x0 * J) * J;
            }

            const T step = g * step_scale;
            splitters[i] = sccd::min<T>(xmax, sccd::max<T>(xmin, x0 - step));
        }
    }

    template <int SplitDim, int N, typename T>
    inline bool grid_search_adaptive_split_vf_axis(const sccd::Box<T> &domain,
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
                                                   T &toi,
                                                   T &u,
                                                   T &v,
                                                   std::vector<sccd::Box<T>> &stack,
                                                   const bool refine) {
        const T lo = domain.tuv[SplitDim].lower;
        const T hi = domain.tuv[SplitDim].upper;

        alignas(64) T splitters[N];
        normal_equation_axis_splitters_vf<SplitDim, N, T>(domain, sv, s1, s2, s3, ev, e1, e2, e3, splitters);

        alignas(64) T samples[N + 2];
        samples[0] = lo;
        samples[N + 1] = hi;
#pragma omp simd aligned(splitters, samples)
        for (int i = 0; i < N; ++i) {
            samples[i + 1] = splitters[i];
        }

        auto stack_size = stack.size();

        bool found = false;
        for (int i = 0; i < N + 1; ++i) {
            // for (int i = N; i >= 0; --i) {
            const T sample_min = samples[i];
            const T sample_max = samples[i + 1];
            const T tt_min = SplitDim == 0 ? sample_min : domain.tuv[0].lower;
            const T tt_max = SplitDim == 0 ? sample_max : domain.tuv[0].upper;
            const T uu_min = SplitDim == 1 ? sample_min : domain.tuv[1].lower;
            const T uu_max = SplitDim == 1 ? sample_max : domain.tuv[1].upper;
            const T vv_min = SplitDim == 2 ? sample_min : domain.tuv[2].lower;
            const T vv_max = SplitDim == 2 ? sample_max : domain.tuv[2].upper;

            if (tt_min >= toi || uu_min + vv_min >= T(1) + tol) {
                continue;
            }

            T fmin[3];
            T fmax[3];
            init_codomain_bounds<T>(fmin, fmax);

            for (int mask = 0; mask < 8; ++mask) {
                const T ct = (mask & 1) ? tt_max : tt_min;
                const T cu = (mask & 2) ? uu_max : uu_min;
                const T cv = (mask & 4) ? vv_max : vv_min;
                T F[3];
                diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, ct, cu, cv, F);
                update_codomain_bounds<T>(F, fmin, fmax);
            }

            bool accepted = false;
            if (!codomain_acceptance<T>(fmin, fmax, tol, tols, accepted)) {
                // Does not contain origin
                continue;
            }

            accepted = accepted && (tt_min > 0);

            Box<T> box = split_axis_box<SplitDim, T>(domain, sample_min, sample_max);
            if (accepted || box.depth > max_iter) {
                found |= accept_grid_root_vf<T>(box, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, refine);
                continue;
            }

            stack.push_back(box);
        }

        // Make sure the tmin is on top of the stack
        if constexpr (SplitDim == 0) {
            std::reverse(stack.begin() + stack_size, stack.end());
        }

        return found;
    }

    template <int N, typename T>
    inline bool grid_search_adaptive_split_vf(const sccd::Box<T> &domain,
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
                                              T &toi,
                                              T &u,
                                              T &v,
                                              std::vector<sccd::Box<T>> &stack,
                                              const bool refine) {
        const int split_dim = domain.widest_dimension();
        if (split_dim == 0) {
            return grid_search_adaptive_split_vf_axis<0, N, T>(
                domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_adaptive_split_vf_axis<1, N, T>(
                domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        return grid_search_adaptive_split_vf_axis<2, N, T>(
            domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vf(const int max_iter,
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
                                          const Box<T> &initial_domain,
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        if (initial_domain.tuv[0].lower >= t) {
            return false;
        }

        return grid_search_adaptive_split_vf<ADAPTIVE_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vf(const int max_iter,
                                          const T tol,
                                          const T sv[3],
                                          const T s1[3],
                                          const T s2[3],
                                          const T s3[3],
                                          const T ev[3],
                                          const T e1[3],
                                          const T e2[3],
                                          const T e3[3],
                                          const Box<T> &initial_domain,
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        T tols[3];
        compute_face_vertex_tolerance<T>(tol, sv, s1, s2, s3, ev, e1, e2, e3, tols);
        return find_root_grid_adaptive_split_vf<T>(
            max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vf(const int max_iter,
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
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        return find_root_grid_adaptive_split_vf<T>(
            max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, unit_domain_box<T>(), t, u, v, stack, refine);
    }

    template <int SplitDim, int N, typename T>
    inline bool grid_search_uniform_split_vf_axis(const sccd::Box<T> &domain,
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
                                                  T &toi,
                                                  T &u,
                                                  T &v,
                                                  std::vector<sccd::Box<T>> &stack,
                                                  const bool refine) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);

        const T lo = domain.tuv[SplitDim].lower;
        const T h = (domain.tuv[SplitDim].upper - lo) / T(N);

        alignas(64) T sample_min[N];
        alignas(64) T sample_max[N];
        alignas(64) T fmin[3][N];
        alignas(64) T fmax[3][N];
        alignas(64) bool contains_zero[N];
        alignas(64) bool accept[N];

#pragma omp simd aligned(sample_min, sample_max)
        for (int i = 0; i < N; ++i) {
            sample_min[i] = lo + h * T(i);
            sample_max[i] = lo + h * T(i + 1);
        }

        for (int d = 0; d < 3; ++d) {
#pragma omp simd aligned(fmin, fmax)
            for (int i = 0; i < N; ++i) {
                fmin[d][i] = std::numeric_limits<T>::max();
                fmax[d][i] = std::numeric_limits<T>::lowest();
            }
        }

        for (int mask = 0; mask < 8; ++mask) {
            const bool mt = mask & 1;
            const bool mu = mask & 2;
            const bool mv = mask & 4;
#pragma omp simd aligned(sample_min, sample_max, fmin, fmax)
            for (int i = 0; i < N; ++i) {
                const T tt = SplitDim == 0 ? (mt ? sample_max[i] : sample_min[i])
                                           : (mt ? domain.tuv[0].upper : domain.tuv[0].lower);
                const T uu = SplitDim == 1 ? (mu ? sample_max[i] : sample_min[i])
                                           : (mu ? domain.tuv[1].upper : domain.tuv[1].lower);
                const T vv = SplitDim == 2 ? (mv ? sample_max[i] : sample_min[i])
                                           : (mv ? domain.tuv[2].upper : domain.tuv[2].lower);
                T F[3];
                diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, tt, uu, vv, F);
                for (int d = 0; d < 3; ++d) {
                    fmin[d][i] = sccd::min<T>(fmin[d][i], F[d]);
                    fmax[d][i] = sccd::max<T>(fmax[d][i], F[d]);
                }
            }
        }

#pragma omp simd aligned(fmin, fmax, contains_zero, accept)
        for (int i = 0; i < N; ++i) {
            const T fmin_i[3] = {fmin[0][i], fmin[1][i], fmin[2][i]};
            const T fmax_i[3] = {fmax[0][i], fmax[1][i], fmax[2][i]};
            contains_zero[i] = codomain_acceptance<T>(fmin_i, fmax_i, tol, tols, accept[i]);
        }

        auto stack_size = stack.size();

        bool found = false;
        for (int i = 0; i < N; ++i) {
            if (!contains_zero[i]) {
                continue;
            }

            const T sample_lo = sample_min[i];
            const T sample_hi = sample_max[i];
            const T tt_min = SplitDim == 0 ? sample_lo : domain.tuv[0].lower;
            const T tt_max = SplitDim == 0 ? sample_hi : domain.tuv[0].upper;
            const T uu_min = SplitDim == 1 ? sample_lo : domain.tuv[1].lower;
            const T uu_max = SplitDim == 1 ? sample_hi : domain.tuv[1].upper;
            const T vv_min = SplitDim == 2 ? sample_lo : domain.tuv[2].lower;
            const T vv_max = SplitDim == 2 ? sample_hi : domain.tuv[2].upper;

            if (tt_min >= toi || uu_min + vv_min >= T(1) + tol) {
                continue;
            }

            accept[i] = accept[i] && (tt_min > 0);

            Box<T> box = split_axis_box<SplitDim, T>(domain, sample_lo, sample_hi);
            if (accept[i] || box.depth > max_iter) {
                found |= accept_grid_root_vf<T>(box, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, refine);
                continue;
            }

            // box.tuv[0].lower = std::min(box.tuv[0].lower, toi);
            stack.push_back(box);
        }

        // Make sure the tmin is on top of the stack
        if constexpr (SplitDim == 0) {
            std::reverse(stack.begin() + stack_size, stack.end());
        }

        return found;
    }

    template <int N, typename T>
    inline bool grid_search_uniform_split_vf(const sccd::Box<T> &domain,
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
                                             T &toi,
                                             T &u,
                                             T &v,
                                             std::vector<sccd::Box<T>> &stack,
                                             const bool refine) {
        const int split_dim = domain.widest_dimension();
        if (split_dim == 0) {
            return grid_search_uniform_split_vf_axis<0, N, T>(
                domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_uniform_split_vf_axis<1, N, T>(
                domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        return grid_search_uniform_split_vf_axis<2, N, T>(
            domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_vf(const int max_iter,
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
                                         const Box<T> &initial_domain,
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        if (initial_domain.tuv[0].lower >= t) {
            return false;
        }

        return grid_search_uniform_split_vf<UNIFORM_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_vf(const int max_iter,
                                         const T tol,
                                         const T sv[3],
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T ev[3],
                                         const T e1[3],
                                         const T e2[3],
                                         const T e3[3],
                                         const Box<T> &initial_domain,
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        T tols[3];
        compute_face_vertex_tolerance<T>(tol, sv, s1, s2, s3, ev, e1, e2, e3, tols);
        return find_root_grid_uniform_split_vf<T>(
            max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_vf(const int max_iter,
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
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        return find_root_grid_uniform_split_vf<T>(
            max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, unit_domain_box<T>(), t, u, v, stack, refine);
    }

    template <typename T>
    inline void diff_ee(const T s1[3],
                        const T s2[3],
                        const T s3[3],
                        const T s4[3],
                        const T e1[3],
                        const T e2[3],
                        const T e3[3],
                        const T e4[3],
                        const T t,
                        const T u,
                        const T v,
                        T *const SCCD_RESTRICT diff) {
        for (int d = 0; d < 3; ++d) {
            const T ea0 = (e1[d] - s1[d]) * t + s1[d];
            const T ea1 = (e2[d] - s2[d]) * t + s2[d];
            const T eb0 = (e3[d] - s3[d]) * t + s3[d];
            const T eb1 = (e4[d] - s4[d]) * t + s4[d];
            diff[d] = ((ea1 - ea0) * u + ea0) - ((eb1 - eb0) * v + eb0);
        }
    }

    template <int SplitDim, int N, typename T>
    inline static void normal_equation_axis_splitters_ee(const Box<T> &domain,
                                                         const T s1[3],
                                                         const T s2[3],
                                                         const T s3[3],
                                                         const T s4[3],
                                                         const T e1[3],
                                                         const T e2[3],
                                                         const T e3[3],
                                                         const T e4[3],
                                                         T *const SCCD_RESTRICT splitters) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);

        const T lo = domain.tuv[SplitDim].lower;
        const T hi = domain.tuv[SplitDim].upper;
        const T h = (hi - lo) / T(N + 1);
        const T radius = h * T(0.45);
        const T mid_t = (domain.tuv[0].lower + domain.tuv[0].upper) * T(0.5);
        const T mid_u = (domain.tuv[1].lower + domain.tuv[1].upper) * T(0.5);
        const T mid_v = (domain.tuv[2].lower + domain.tuv[2].upper) * T(0.5);
        const T eps = std::numeric_limits<T>::epsilon();

        T F_base[3];
        T J_axis[3];
        T H_axis = T(0);

        const T base_t = SplitDim == 0 ? T(0) : mid_t;
        const T base_u = SplitDim == 1 ? T(0) : mid_u;
        const T base_v = SplitDim == 2 ? T(0) : mid_v;
        diff_ee<T>(s1, s2, s3, s4, e1, e2, e3, e4, base_t, base_u, base_v, F_base);

        for (int d = 0; d < 3; ++d) {
            if constexpr (SplitDim == 0) {
                J_axis[d] = (T(1) - mid_u) * (e1[d] - s1[d]) + mid_u * (e2[d] - s2[d]) -
                            (T(1) - mid_v) * (e3[d] - s3[d]) - mid_v * (e4[d] - s4[d]);
            } else if constexpr (SplitDim == 1) {
                J_axis[d] = (T(1) - mid_t) * (s2[d] - s1[d]) + mid_t * (e2[d] - e1[d]);
            } else {
                J_axis[d] = -((T(1) - mid_t) * (s4[d] - s3[d]) + mid_t * (e4[d] - e3[d]));
            }
            H_axis += J_axis[d] * J_axis[d];
        }
        const T step_scale = H_axis > eps ? T(1) / H_axis : T(0.00001);

#pragma omp simd
        for (int i = 0; i < N; ++i) {
            const T x0 = lo + h * T(i + 1);
            auto xmin = sccd::max<T>(lo, x0 - radius);
            auto xmax = sccd::min<T>(hi, x0 + radius);

            T g = T(0);
            for (int d = 0; d < 3; ++d) {
                const T J = J_axis[d];
                g += (F_base[d] + x0 * J) * J;
            }

            const T step = g * step_scale;
            splitters[i] = sccd::min<T>(xmax, sccd::max<T>(xmin, x0 - step));
        }
    }

    template <int SplitDim, int N, typename T>
    inline bool grid_search_adaptive_split_ee_axis(const sccd::Box<T> &domain,
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
                                                   T &toi,
                                                   T &u,
                                                   T &v,
                                                   std::vector<sccd::Box<T>> &stack,
                                                   const bool refine) {
        (void)refine;

        const T lo = domain.tuv[SplitDim].lower;
        const T hi = domain.tuv[SplitDim].upper;

        alignas(64) T splitters[N];
        normal_equation_axis_splitters_ee<SplitDim, N, T>(domain, s1, s2, s3, s4, e1, e2, e3, e4, splitters);

        alignas(64) T samples[N + 2];
        samples[0] = lo;
        samples[N + 1] = hi;
#pragma omp simd
        for (int i = 0; i < N; ++i) {
            samples[i + 1] = splitters[i];
        }

        bool found = false;
        for (int i = 0; i < N + 1; ++i) {
            const T sample_min = samples[i];
            const T sample_max = samples[i + 1];
            const T tt_min = SplitDim == 0 ? sample_min : domain.tuv[0].lower;
            const T tt_max = SplitDim == 0 ? sample_max : domain.tuv[0].upper;
            const T uu_min = SplitDim == 1 ? sample_min : domain.tuv[1].lower;
            const T uu_max = SplitDim == 1 ? sample_max : domain.tuv[1].upper;
            const T vv_min = SplitDim == 2 ? sample_min : domain.tuv[2].lower;
            const T vv_max = SplitDim == 2 ? sample_max : domain.tuv[2].upper;

            if (tt_min >= toi) {
                continue;
            }

            T fmin[3];
            T fmax[3];
            init_codomain_bounds<T>(fmin, fmax);

            for (int mask = 0; mask < 8; ++mask) {
                const T ct = (mask & 1) ? tt_max : tt_min;
                const T cu = (mask & 2) ? uu_max : uu_min;
                const T cv = (mask & 4) ? vv_max : vv_min;
                T F[3];
                diff_ee<T>(s1, s2, s3, s4, e1, e2, e3, e4, ct, cu, cv, F);
                update_codomain_bounds<T>(F, fmin, fmax);
            }

            bool accepted = false;
            if (!codomain_acceptance<T>(fmin, fmax, tol, tols, accepted)) {
                continue;
            }

            accepted = accepted && (tt_min > 0);

            Box<T> box = split_axis_box<SplitDim, T>(domain, sample_min, sample_max);
            if (accepted || box.depth > max_iter) {
                found |= accept_grid_root_ee<T>(box, toi, u, v);
                continue;
            }

            // box.tuv[0].lower = std::min(box.tuv[0].lower, toi);
            stack.push_back(box);
        }

        return found;
    }

    template <int N, typename T>
    inline bool grid_search_adaptive_split_ee(const sccd::Box<T> &domain,
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
                                              T &toi,
                                              T &u,
                                              T &v,
                                              std::vector<sccd::Box<T>> &stack,
                                              const bool refine) {
        const int split_dim = domain.widest_dimension();
        if (split_dim == 0) {
            return grid_search_adaptive_split_ee_axis<0, N, T>(
                domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_adaptive_split_ee_axis<1, N, T>(
                domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        return grid_search_adaptive_split_ee_axis<2, N, T>(
            domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_ee(const int max_iter,
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
                                          const Box<T> &initial_domain,
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        if (initial_domain.tuv[0].lower >= t) {
            return false;
        }

        return grid_search_adaptive_split_ee<ADAPTIVE_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_ee(const int max_iter,
                                          const T tol,
                                          const T s1[3],
                                          const T s2[3],
                                          const T s3[3],
                                          const T s4[3],
                                          const T e1[3],
                                          const T e2[3],
                                          const T e3[3],
                                          const T e4[3],
                                          const Box<T> &initial_domain,
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        T tols[3];
        compute_edge_edge_tolerance<T>(tol, s1, s2, s3, s4, e1, e2, e3, e4, tols);
        return find_root_grid_adaptive_split_ee<T>(
            max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_ee(const int max_iter,
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
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        return find_root_grid_adaptive_split_ee<T>(
            max_iter, tol, s1, s2, s3, s4, e1, e2, e3, e4, unit_domain_box<T>(), t, u, v, stack, refine);
    }

    template <int SplitDim, int N, typename T>
    inline bool grid_search_uniform_split_ee_axis(const sccd::Box<T> &domain,
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
                                                  T &toi,
                                                  T &u,
                                                  T &v,
                                                  std::vector<sccd::Box<T>> &stack,
                                                  const bool refine) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);
        (void)refine;

        const T lo = domain.tuv[SplitDim].lower;
        const T h = (domain.tuv[SplitDim].upper - lo) / T(N);

        alignas(64) T sample_min[N];
        alignas(64) T sample_max[N];
        alignas(64) T fmin[3][N];
        alignas(64) T fmax[3][N];
        alignas(64) bool contains_zero[N];
        alignas(64) bool accept[N];

#pragma omp simd
        for (int i = 0; i < N; ++i) {
            sample_min[i] = lo + h * T(i);
            sample_max[i] = lo + h * T(i + 1);
        }

        for (int d = 0; d < 3; ++d) {
#pragma omp simd
            for (int i = 0; i < N; ++i) {
                fmin[d][i] = std::numeric_limits<T>::max();
                fmax[d][i] = std::numeric_limits<T>::lowest();
            }
        }

        for (int mask = 0; mask < 8; ++mask) {
            const bool mt = mask & 1;
            const bool mu = mask & 2;
            const bool mv = mask & 4;
#pragma omp simd
            for (int i = 0; i < N; ++i) {
                const T tt = SplitDim == 0 ? (mt ? sample_max[i] : sample_min[i])
                                           : (mt ? domain.tuv[0].upper : domain.tuv[0].lower);
                const T uu = SplitDim == 1 ? (mu ? sample_max[i] : sample_min[i])
                                           : (mu ? domain.tuv[1].upper : domain.tuv[1].lower);
                const T vv = SplitDim == 2 ? (mv ? sample_max[i] : sample_min[i])
                                           : (mv ? domain.tuv[2].upper : domain.tuv[2].lower);
                T F[3];
                diff_ee<T>(s1, s2, s3, s4, e1, e2, e3, e4, tt, uu, vv, F);
                for (int d = 0; d < 3; ++d) {
                    fmin[d][i] = sccd::min<T>(fmin[d][i], F[d]);
                    fmax[d][i] = sccd::max<T>(fmax[d][i], F[d]);
                }
            }
        }

#pragma omp simd
        for (int i = 0; i < N; ++i) {
            const T fmin_i[3] = {fmin[0][i], fmin[1][i], fmin[2][i]};
            const T fmax_i[3] = {fmax[0][i], fmax[1][i], fmax[2][i]};
            contains_zero[i] = codomain_acceptance<T>(fmin_i, fmax_i, tol, tols, accept[i]);
        }

        bool found = false;
        for (int i = 0; i < N; ++i) {
            if (!contains_zero[i]) {
                continue;
            }

            const T sample_lo = sample_min[i];
            const T sample_hi = sample_max[i];
            const T tt_min = SplitDim == 0 ? sample_lo : domain.tuv[0].lower;
            const T tt_max = SplitDim == 0 ? sample_hi : domain.tuv[0].upper;
            const T uu_min = SplitDim == 1 ? sample_lo : domain.tuv[1].lower;
            const T uu_max = SplitDim == 1 ? sample_hi : domain.tuv[1].upper;
            const T vv_min = SplitDim == 2 ? sample_lo : domain.tuv[2].lower;
            const T vv_max = SplitDim == 2 ? sample_hi : domain.tuv[2].upper;

            if (tt_min >= toi) {
                continue;
            }

            accept[i] = accept[i] && (tt_min > 0);

            Box<T> box = split_axis_box<SplitDim, T>(domain, sample_lo, sample_hi);
            if (accept[i] || box.depth > max_iter) {
                found |= accept_grid_root_ee<T>(box, toi, u, v);
                continue;
            }

            // box.tuv[0].lower = std::min(box.tuv[0].lower, toi);
            stack.push_back(box);
        }

        return found;
    }

    template <int N, typename T>
    inline bool grid_search_uniform_split_ee(const sccd::Box<T> &domain,
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
                                             T &toi,
                                             T &u,
                                             T &v,
                                             std::vector<sccd::Box<T>> &stack,
                                             const bool refine) {
        const int split_dim = domain.widest_dimension();
        if (split_dim == 0) {
            return grid_search_uniform_split_ee_axis<0, N, T>(
                domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_uniform_split_ee_axis<1, N, T>(
                domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        return grid_search_uniform_split_ee_axis<2, N, T>(
            domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_ee(const int max_iter,
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
                                         const Box<T> &initial_domain,
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        if (initial_domain.tuv[0].lower >= t) {
            return false;
        }

        return grid_search_uniform_split_ee<UNIFORM_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_ee(const int max_iter,
                                         const T tol,
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T s4[3],
                                         const T e1[3],
                                         const T e2[3],
                                         const T e3[3],
                                         const T e4[3],
                                         const Box<T> &initial_domain,
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        T tols[3];
        compute_edge_edge_tolerance<T>(tol, s1, s2, s3, s4, e1, e2, e3, e4, tols);
        return find_root_grid_uniform_split_ee<T>(
            max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_ee(const int max_iter,
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
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        return find_root_grid_uniform_split_ee<T>(
            max_iter, tol, s1, s2, s3, s4, e1, e2, e3, e4, unit_domain_box<T>(), t, u, v, stack, refine);
    }

}  // namespace sccd

#endif  // S_ROOT_FINDER_HPP
