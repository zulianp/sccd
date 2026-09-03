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
#include "snumerical_error.hpp"
#include "snumtol.hpp"

// #define SCCD_ENABLE_TIGHT_INCLUSION
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#include <Eigen/Dense>
#include "tight_inclusion/ccd.hpp"
#include "tight_inclusion/interval_root_finder.hpp"
#endif

#ifndef ADAPTIVE_NUM_SPLITS
#define ADAPTIVE_NUM_SPLITS 2
#endif
#ifndef UNIFORM_NUM_SPLITS
#define UNIFORM_NUM_SPLITS 4
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

        ticcd::Array3 err(-1, -1, -1);

        ticcd::Scalar ms = 0;
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

        int widest_dimension(const T scale[3]) const {
            const T dt = (tuv[0].upper - tuv[0].lower) * scale[0];
            const T du = (tuv[1].upper - tuv[1].lower) * scale[1];
            const T dv = (tuv[2].upper - tuv[2].lower) * scale[2];
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
    inline bool codomain_acceptance(const T fmin[3],
                                    const T fmax[3],
                                    const T tol,
                                    const T tols[3],
                                    const T numerical_error[3],
                                    bool &accept) {
        bool contains_zero = true;
        bool smaller_than_axis_tol = true;
        bool inside_error_box = true;
        bool smaller_than_scalar_tol = true;
        bool degenerate_interval = true;

        for (int d = 0; d < 3; ++d) {
            const T interval_width = fmax[d] - fmin[d];
            contains_zero = contains_zero && (fmin[d] <= numerical_error[d]) && (fmax[d] >= -numerical_error[d]);
            smaller_than_axis_tol = smaller_than_axis_tol && (interval_width <= tols[d]);
            inside_error_box =
                inside_error_box && (fmin[d] >= -numerical_error[d]) && (fmax[d] <= numerical_error[d]);
            smaller_than_scalar_tol = smaller_than_scalar_tol && (interval_width < tol);
            degenerate_interval = degenerate_interval && (fmin[d] >= fmax[d]);
        }

        accept = contains_zero && (smaller_than_axis_tol || inside_error_box || smaller_than_scalar_tol || degenerate_interval);
        return contains_zero;
    }

    /**
     * \brief Acceptance for the vertex-quad search.
     *
     * `numerical_error` is the certified bound, and it is what makes the
     * rejection sound. This used to pad the origin-containment test with
     * std::numeric_limits<T>::epsilon() instead -- about 30x smaller than the
     * bound for unit-scale geometry -- which is an *unsound rejection*: it can
     * discard a box that contains a root, and a discarded box is a missed
     * collision. Every other acceptance condition here is free to be as loose as
     * it likes, because accepting early only reports an earlier time of impact;
     * only this one test can lose a root.
     */
    template <typename T>
    inline bool codomain_acceptance_vq(const T fmin[3],
                                       const T fmax[3],
                                       const T tol,
                                       const T tols[3],
                                       const T numerical_error[3],
                                       bool &accept) {
        bool contains_zero = true;
        bool smaller_than_axis_tol = true;
        bool inside_epsilon_box = true;
        bool smaller_than_scalar_tol = true;
        bool degenerate_interval = true;

        for (int d = 0; d < 3; ++d) {
            const T interval_width = fmax[d] - fmin[d];
            contains_zero =
                contains_zero && (fmin[d] <= numerical_error[d]) && (fmax[d] >= -numerical_error[d]);
            smaller_than_axis_tol = smaller_than_axis_tol && (interval_width <= tols[d]);
            inside_epsilon_box =
                inside_epsilon_box && (fmin[d] >= -numerical_error[d]) && (fmax[d] <= numerical_error[d]);
            smaller_than_scalar_tol = smaller_than_scalar_tol && (interval_width < tol);
            degenerate_interval = degenerate_interval && (fmin[d] >= fmax[d]);
        }

        accept = contains_zero &&
                 (smaller_than_axis_tol || inside_epsilon_box || smaller_than_scalar_tol || degenerate_interval);
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

        T numerical_error[3];
        numerical_error_bound<true, T>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);

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
                if (codomain.tuv[i].lower > numerical_error[i] ||
                    codomain.tuv[i].upper < -numerical_error[i]) {
                    return false;
                }
            }

            inside_box = true;
            for (int i = 0; i < 3; i++) {
                if (codomain.tuv[i].lower < -numerical_error[i] ||
                    codomain.tuv[i].upper > numerical_error[i]) {
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

    /** \brief Componentwise codomain bounds of one face of a sub-box. */
    template <typename T>
    struct FaceBounds {
        T fmin[3];
        T fmax[3];
    };

    /**
     * \brief Codomain bounds of the 4 corners on the face where the split axis
     *        equals \p level.
     *
     * Adjacent sub-intervals of a split share exactly this face, so evaluating
     * it once and carrying it forward halves the corner evaluations: a split
     * into N+1 sub-intervals needs N+2 faces (4 evaluations each) rather than
     * N+1 boxes at 8 evaluations each.
     *
     * Returned by value rather than through out-pointers: the caller's bounds
     * would otherwise be assumed to alias the geometry the evaluator reads,
     * which blocks the cross-corner CSE that makes these evaluations cheap.
     *
     * Reassociating the min/max this way is exact -- min and max introduce no
     * rounding, so the resulting bounds are bit-identical to reducing over all
     * 8 corners in sequence.
     */
    template <int SplitDim, typename T, typename Eval>
    SCCD_ALWAYS_INLINE static FaceBounds<T> axis_face_bounds(const T level,
                                                             const T free_lo[3],
                                                             const T free_hi[3],
                                                             Eval eval) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        constexpr int A = (SplitDim + 1) % 3;
        constexpr int B = (SplitDim + 2) % 3;

        const T a_lo = free_lo[A];
        const T a_hi = free_hi[A];
        const T b_lo = free_lo[B];
        const T b_hi = free_hi[B];

        T corner[4][3];
        {
            T c[3];
            c[SplitDim] = level;
            c[A] = a_lo;
            c[B] = b_lo;
            eval(c[0], c[1], c[2], corner[0]);
            c[A] = a_hi;
            c[B] = b_lo;
            eval(c[0], c[1], c[2], corner[1]);
            c[A] = a_lo;
            c[B] = b_hi;
            eval(c[0], c[1], c[2], corner[2]);
            c[A] = a_hi;
            c[B] = b_hi;
            eval(c[0], c[1], c[2], corner[3]);
        }

        FaceBounds<T> out;
        for (int d = 0; d < 3; ++d) {
            const T m01 = sccd::min<T>(corner[0][d], corner[1][d]);
            const T m23 = sccd::min<T>(corner[2][d], corner[3][d]);
            const T x01 = sccd::max<T>(corner[0][d], corner[1][d]);
            const T x23 = sccd::max<T>(corner[2][d], corner[3][d]);
            out.fmin[d] = sccd::min<T>(m01, m23);
            out.fmax[d] = sccd::max<T>(x01, x23);
        }
        return out;
    }

    // ---------------------------------------------------------------------
    // On the adaptive splitter below, from a measured study of it:
    //
    //  * It earns its keep. Against a uniform grid it is ~1.6x faster in the
    //    scalar search, using three sub-intervals where uniform uses five.
    //
    //  * It is not, however, adapting much. The x0 term cancels exactly:
    //
    //        x_new = x0 - Jt(F + x0 J)/(JtJ) = -JtF/(JtJ)
    //
    //    (confirmed numerically to 4e-15 over 200k boxes, on the H_axis > eps
    //    branch; the degenerate fallback does not cancel). So the loop below
    //    does not nudge each grid point toward a nearby root -- every iteration
    //    solves for the *same* least-squares point, and the N splitters differ
    //    only in which window each is then clamped into. That is why the clamp
    //    saturates: 95.9% of splitters are pinned at it here (N = 2), and 97.6%
    //    in the vectorized kernel's copy of this function (N = 4), whose
    //    narrower windows pin it harder still.
    //
    //    A consequence worth taking: the loop recomputes one quantity N times
    //    and could hoist -JtF/(JtJ) out, bit-for-bit unchanged.
    //
    //  * The clamp is load-bearing, not cosmetic. It is what guarantees each
    //    child is a bounded fraction of its parent. A splitter free to place a
    //    cut anywhere has no guaranteed progress; tried directly, the search
    //    failed to terminate.
    //
    //  * Placement cannot affect correctness. The sub-intervals tile the parent
    //    (verified: no monotonicity or range violations in 600k boxes, which the
    //    damping < 0.5 guarantees) and every child is still tested exactly. The
    //    splitter is purely a performance knob.
    //
    //  * Placing the cuts analytically was tried properly and lost. A chord
    //    across the split axis's two faces locates 68.3% of the 68.4% that is
    //    genuinely root-free, soundly (min-of-affine is concave, max-of-affine
    //    convex, so a chord can only understate it -- 0 over-claims in 4002
    //    boxes). A splitter built on that, cutting at the dead prefix and suffix
    //    and reserving one cut to subdivide what survives, ran 2.7x SLOWER than
    //    this one on cloth-funnel (137 ms against 51 ms) and isolated LESS dead
    //    volume: 56.7% against 61.1% at N = 2, and 58.4% against 64.5% at N = 4.
    //
    //    The reason is worth keeping, because it is easy to get wrong: that
    //    68.4% is dead *pointwise*, as a union over the three dimensions -- one
    //    dimension rules out one stretch of the axis and another dimension a
    //    different stretch. Rejecting a child box needs a single dimension to
    //    rule out the whole child. So a cut derived from the union cannot
    //    harvest it, and only about 58% of the axis is rejectable by one
    //    dimension over one contiguous stretch. A multi-way split does better
    //    precisely because each child may then die by a different dimension --
    //    which is what this splitter, crude as it is, already gets.
    //
    //    Two cheaper variants also lost: a bare crossing cut has no progress
    //    guarantee and failed to terminate, and a fixed asymmetric cut is a
    //    workload-dependent speed-versus-precision dial (0.05 ran 3x faster than
    //    bisection on cloth-funnel and did not finish on armadillo-rollers).
    //
    // The conservative kernel therefore bisects, which is the only variant that
    // is uniformly well behaved. See benchmark/oracle/README.md.
    // ---------------------------------------------------------------------
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
                                                   const T numerical_error[3],
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

        // The two axes that are not split keep the parent box's extents for
        // every sub-interval, so they are hoisted out of the loop.
        const T free_lo[3] = {domain.tuv[0].lower, domain.tuv[1].lower, domain.tuv[2].lower};
        const T free_hi[3] = {domain.tuv[0].upper, domain.tuv[1].upper, domain.tuv[2].upper};
        const auto eval_vf = [&](const T ct, const T cu, const T cv, T F[3]) {
            diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, ct, cu, cv, F);
        };

        // Face shared with the previous sub-interval; `carried_level` is the
        // sample index it was evaluated at, or -1 when nothing is carried.
        FaceBounds<T> lower_face;
        int carried_level = -1;

        bool found = false;
        for (int i = 0; i < N + 1; ++i) {
            // for (int i = N; i >= 0; --i) {
            const T sample_min = samples[i];
            const T sample_max = samples[i + 1];
            const T tt_min = SplitDim == 0 ? sample_min : domain.tuv[0].lower;
            const T uu_min = SplitDim == 1 ? sample_min : domain.tuv[1].lower;
            const T vv_min = SplitDim == 2 ? sample_min : domain.tuv[2].lower;

            if (tt_min >= toi || uu_min + vv_min >= T(1) + tol) {
                continue;
            }

            if (carried_level != i) {
                lower_face = axis_face_bounds<SplitDim, T>(sample_min, free_lo, free_hi, eval_vf);
            }

            const FaceBounds<T> upper_face = axis_face_bounds<SplitDim, T>(sample_max, free_lo, free_hi, eval_vf);

            T fmin[3];
            T fmax[3];
            for (int d = 0; d < 3; ++d) {
                fmin[d] = sccd::min<T>(lower_face.fmin[d], upper_face.fmin[d]);
                fmax[d] = sccd::max<T>(lower_face.fmax[d], upper_face.fmax[d]);
            }
            lower_face = upper_face;
            carried_level = i + 1;

            bool accepted = false;
            if (!codomain_acceptance<T>(fmin, fmax, tol, tols, numerical_error, accepted)) {
                // Does not contain origin
                continue;
            }

            // Reject a box whose t lower bound is zero. This is TightInclusion's
            // no_zero_toi policy: a contact reported at exactly t == 0 stalls a
            // solver, which cannot advance the step at all.
            //
            // It does not cost a collision. Measured against the datasets' exact
            // roots, every mode reports gt_missed == 0 and gt_late == 0 with this
            // in place -- see benchmark/oracle/README.md. TightInclusion does
            // report a hit at t == 0 on some of these queries, but its answer is
            // a conservative lower bound and over-reports; the exact roots agree
            // with rejecting them.
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
                                              const T numerical_error[3],
                                              const T codomain_widths[3],
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
        const int split_dim = domain.widest_dimension(codomain_widths);
        if (split_dim == 0) {
            return grid_search_adaptive_split_vf_axis<0, N, T>(
                domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_adaptive_split_vf_axis<1, N, T>(
                domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        return grid_search_adaptive_split_vf_axis<2, N, T>(
            domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
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
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        numerical_error_bound<true, T>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);
        return grid_search_adaptive_split_vf<N, T>(
            domain, max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vf(const int max_iter,
                                          const T tol,
                                          const T tols[3],
                                          const T numerical_error[3],
                                          const T codomain_widths[3],
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
            initial_domain, max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack, refine);
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
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        numerical_error_bound<true, T>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);
        return find_root_grid_adaptive_split_vf<T>(
            max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, initial_domain, t, u, v, stack, refine);
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

    template <typename T>
    inline void compute_edge_edge_codomain_widths(const T s1[3],
                                                  const T s2[3],
                                                  const T s3[3],
                                                  const T s4[3],
                                                  const T e1[3],
                                                  const T e2[3],
                                                  const T e3[3],
                                                  const T e4[3],
                                                  T widths[3]) {
        T wt = T(0);
        T wu = T(0);
        T wv = T(0);
        for (int d = 0; d < 3; ++d) {
            const T a0 = e1[d] - s1[d];
            const T a1 = e2[d] - s2[d];
            const T b0 = e3[d] - s3[d];
            const T b1 = e4[d] - s4[d];
            wt = sccd::max<T>(wt,
                              sccd::max<T>(sccd::max<T>(sccd::abs<T>(a0 - b0), sccd::abs<T>(a0 - b1)),
                                           sccd::max<T>(sccd::abs<T>(a1 - b0), sccd::abs<T>(a1 - b1))));
            wu = sccd::max<T>(wu, sccd::max<T>(sccd::abs<T>(s2[d] - s1[d]), sccd::abs<T>(e2[d] - e1[d])));
            wv = sccd::max<T>(wv, sccd::max<T>(sccd::abs<T>(s4[d] - s3[d]), sccd::abs<T>(e4[d] - e3[d])));
        }
        widths[0] = wt;
        widths[1] = wu;
        widths[2] = wv;
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
                                                   const T numerical_error[3],
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

        const T free_lo[3] = {domain.tuv[0].lower, domain.tuv[1].lower, domain.tuv[2].lower};
        const T free_hi[3] = {domain.tuv[0].upper, domain.tuv[1].upper, domain.tuv[2].upper};
        const auto eval_ee = [&](const T ct, const T cu, const T cv, T F[3]) {
            diff_ee<T>(s1, s2, s3, s4, e1, e2, e3, e4, ct, cu, cv, F);
        };

        FaceBounds<T> lower_face;
        int carried_level = -1;

        bool found = false;
        for (int i = 0; i < N + 1; ++i) {
            const T sample_min = samples[i];
            const T sample_max = samples[i + 1];
            const T tt_min = SplitDim == 0 ? sample_min : domain.tuv[0].lower;

            if (tt_min >= toi) {
                continue;
            }

            if (carried_level != i) {
                lower_face = axis_face_bounds<SplitDim, T>(sample_min, free_lo, free_hi, eval_ee);
            }

            const FaceBounds<T> upper_face = axis_face_bounds<SplitDim, T>(sample_max, free_lo, free_hi, eval_ee);

            T fmin[3];
            T fmax[3];
            for (int d = 0; d < 3; ++d) {
                fmin[d] = sccd::min<T>(lower_face.fmin[d], upper_face.fmin[d]);
                fmax[d] = sccd::max<T>(lower_face.fmax[d], upper_face.fmax[d]);
            }
            lower_face = upper_face;
            carried_level = i + 1;

            bool accepted = false;
            if (!codomain_acceptance<T>(fmin, fmax, tol, tols, numerical_error, accepted)) {
                continue;
            }

            // Reject a box whose t lower bound is zero. This is TightInclusion's
            // no_zero_toi policy: a contact reported at exactly t == 0 stalls a
            // solver, which cannot advance the step at all.
            //
            // It does not cost a collision. Measured against the datasets' exact
            // roots, every mode reports gt_missed == 0 and gt_late == 0 with this
            // in place -- see benchmark/oracle/README.md. TightInclusion does
            // report a hit at t == 0 on some of these queries, but its answer is
            // a conservative lower bound and over-reports; the exact roots agree
            // with rejecting them.
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
                                              const T numerical_error[3],
                                              const T codomain_widths[3],
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
        const int split_dim = domain.widest_dimension(codomain_widths);
        if (split_dim == 0) {
            return grid_search_adaptive_split_ee_axis<0, N, T>(
                domain, max_iter, tol, tols, numerical_error, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_adaptive_split_ee_axis<1, N, T>(
                domain, max_iter, tol, tols, numerical_error, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        return grid_search_adaptive_split_ee_axis<2, N, T>(
            domain, max_iter, tol, tols, numerical_error, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
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
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        numerical_error_bound<false, T>(s1, s2, s3, s4, e1, e2, e3, e4, numerical_error);
        return grid_search_adaptive_split_ee<N, T>(
            domain, max_iter, tol, tols, numerical_error, codomain_widths, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_ee(const int max_iter,
                                          const T tol,
                                          const T tols[3],
                                          const T numerical_error[3],
                                          const T codomain_widths[3],
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
            initial_domain, max_iter, tol, tols, numerical_error, codomain_widths, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack, refine);
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
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        numerical_error_bound<false, T>(s1, s2, s3, s4, e1, e2, e3, e4, numerical_error);
        return find_root_grid_adaptive_split_ee<T>(
            max_iter, tol, tols, numerical_error, codomain_widths, s1, s2, s3, s4, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
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


}  // namespace sccd

#endif  // S_ROOT_FINDER_HPP
