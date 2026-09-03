#ifndef SCCD_BROADPHASE_STRATEGY_HPP
#define SCCD_BROADPHASE_STRATEGY_HPP

#include "sccd_base.hpp"
#include "smath.hpp"

#include <cstdlib>
#include <cstring>

/**
 * \file
 * \brief Choosing between sweep-and-prune and the 2D cell list.
 *
 * Both are first class: neither wins everywhere and the difference is large in
 * both directions, so the choice has to be explicit.
 *
 * ## What the measurements say
 *
 * GH200, self-overlap, 400k boxes, isotropic distribution, varying only the box
 * size so the *shape* is fixed and only the amount of work changes:
 *
 * | box | sweep | cell2d | speedup | pairs |
 * |---|---|---|---|---|
 * | 1.5 | 1.3 ms | 3.1 ms | 0.41x | 264 |
 * | 6 | 6.2 ms | 5.1 ms | 1.23x | 17k |
 * | 25 | 484 ms | 126 ms | 3.83x | 1.2M |
 * | 100 | 62,842 ms | 9,058 ms | 6.94x | 72M |
 *
 * The crossover is between the first two rows, and **the cell list wins whenever
 * there is a non-trivial amount of work to do**. Where it loses, it loses on an
 * input whose broad phase costs a millisecond either way.
 *
 * ## Why there is no shape heuristic here
 *
 * Several cheap statistics were tried as a selector and none survived contact
 * with the data. They are recorded so they are not retried:
 *
 *  - **Anisotropy** of the box distribution (how far from filling a volume it is)
 *    looked decisive on two synthetic families -- 1.00 where the sweep won, 3.3
 *    where the cell list won. It is refuted by the table above, where anisotropy
 *    is 1.00 on every row and the answer flips anyway. It also misclassifies the
 *    real target: a cloth-ball frame measures 1.96, which a threshold fitted to
 *    those families would call "sweep", while the cell list is 4.2x faster on it.
 *  - **The expected sweep window** `lambda_min = n * mean_extent / span` does not
 *    separate either: 1198 with the sweep winning, 940 with it losing by 3.3x.
 *    The sweep walks its window over sorted contiguous memory with a 32-wide
 *    bitmask, so it absorbs a long window far better than a candidate count
 *    suggests.
 *  - **Estimated pair density** from the same statistics is accurate for uniform
 *    distributions (predicting 270, 16,940, 1.158M against measured 264, 16,992,
 *    1.217M) and then underestimates a real cloth by four orders of magnitude,
 *    because the uniform assumption is exactly what a surface violates.
 *
 * ## So: default to the cell list
 *
 * Not for want of a heuristic, but because the payoff is asymmetric and the
 * measurements bound both sides. The worst observed loss is 1.8 ms on a broad
 * phase that costs 1.3 ms; the best observed gain is 54 seconds. A selector that
 * could read the geometry would be worth having, but guessing wrong in the cheap
 * direction costs milliseconds and guessing wrong in the expensive direction
 * costs minutes.
 *
 * `SCCD_BROADPHASE=sweep` forces the old path, and `broadphase_stats` is exposed
 * so a caller with a scene the defaults suit badly can measure and decide.
 */

namespace sccd {

    enum class BroadPhaseStrategy : int { Auto = 0, Sweep = 1, Cell2D = 2 };

    /**
     * \brief Cheap shape and density statistics for the AABB set.
     *
     * `lambda[d] = sum of extents on axis d / span on axis d` is the expected
     * number of boxes whose interval on `d` covers a given point -- the window
     * sweep-and-prune would walk if it sorted on `d`. `anisotropy` is the spread
     * of that across the three axes, near 1 when the boxes fill a volume and
     * larger when the geometry is flattened.
     *
     * **These describe the input; they do not select a strategy.** Both were
     * tried as selectors and both are refuted -- see the file comment. They are
     * kept because they are one cheap pass, they are the right things to log when
     * a scene behaves unexpectedly, and a caller with its own evidence can build
     * on them.
     */
    template <typename T>
    struct BroadPhaseStats {
        T lambda[3];
        T lam_min;
        T lam_max;

        double anisotropy() const {
            return (double)lam_max / (double)(lam_min > T(0) ? lam_min : T(1));
        }
    };

    template <typename T>
    static BroadPhaseStats<T> broadphase_stats(const ptrdiff_t n, T** const SCCD_RESTRICT aabb) {
        BroadPhaseStats<T> st{};
        if (n <= 0) {
            for (int d = 0; d < 3; ++d) st.lambda[d] = T(0);
            st.lam_min = st.lam_max = T(0);
            return st;
        }

        for (int d = 0; d < 3; ++d) {
            T lo = aabb[d][0];
            T hi = aabb[3 + d][0];
            T sum = T(0);
            for (ptrdiff_t i = 0; i < n; ++i) {
                lo = sccd::min<T>(lo, aabb[d][i]);
                hi = sccd::max<T>(hi, aabb[3 + d][i]);
                sum += aabb[3 + d][i] - aabb[d][i];
            }
            const T span = hi - lo;
            st.lambda[d] = span > T(0) ? (sum / span) : (T)n;
        }

        st.lam_min = sccd::min<T>(st.lambda[0], sccd::min<T>(st.lambda[1], st.lambda[2]));
        st.lam_max = sccd::max<T>(st.lambda[0], sccd::max<T>(st.lambda[1], st.lambda[2]));
        return st;
    }

    /** \brief The strategy the user asked for, or Auto. */
    static inline BroadPhaseStrategy broadphase_strategy_setting() {
        const char* v = getenv("SCCD_BROADPHASE");
        if (!v) return BroadPhaseStrategy::Auto;
        if (std::strcmp(v, "sweep") == 0) return BroadPhaseStrategy::Sweep;
        if (std::strcmp(v, "cell2d") == 0) return BroadPhaseStrategy::Cell2D;
        return BroadPhaseStrategy::Auto;
    }

    /**
     * \brief Resolve Auto; a forced setting passes through.
     *
     * Auto is the cell list, for the reasons in the file comment. \p out_stats is
     * filled when asked so a caller can log or override on its own evidence;
     * computing it costs one pass over the AABBs, so it is skipped when nobody
     * wants it and the strategy is already decided.
     */
    template <typename T>
    static BroadPhaseStrategy choose_broadphase_strategy(const ptrdiff_t n,
                                                         T** const SCCD_RESTRICT aabb,
                                                         BroadPhaseStats<T>* out_stats = nullptr) {
        const BroadPhaseStrategy forced = broadphase_strategy_setting();
        if (out_stats) *out_stats = broadphase_stats<T>(n, aabb);
        if (forced != BroadPhaseStrategy::Auto) return forced;
        return BroadPhaseStrategy::Cell2D;
    }

    static inline const char* broadphase_strategy_name(const BroadPhaseStrategy s) {
        switch (s) {
            case BroadPhaseStrategy::Sweep: return "sweep";
            case BroadPhaseStrategy::Cell2D: return "cell2d";
            case BroadPhaseStrategy::Auto: return "auto";
        }
        return "?";
    }

}  // namespace sccd

#endif  // SCCD_BROADPHASE_STRATEGY_HPP
