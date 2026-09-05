#ifndef SCCD_BROADPHASE_STRATEGY_HPP
#define SCCD_BROADPHASE_STRATEGY_HPP

#include "sccd_base.hpp"
#include "sccd_math.hpp"

#include <chrono>
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
 *  - **Mesh size.** The cell list's large wins were measured on refined meshes,
 *    which suggested a crossover in element count. The end-to-end assessment
 *    refutes that too, and in the wrong direction: across three real scenes the
 *    cell list wins the *smallest* (cloth-funnel) and loses the *largest*
 *    (cloth-ball, where the sweep is 1.36x faster) and the middle one
 *    (armadillo-rollers, 1.59x). See `wip/ASSESSMENT.md`.
 *
 * ## So: Auto measures, it does not guess
 *
 * Four heuristics are refuted above, and the end-to-end assessment refuted the
 * constant that stood in for them. `Auto` used to resolve to the cell list
 * unconditionally, on the strength of synthetic box-list benchmarks where it
 * wins by 4-7x and where the worst case for choosing it is 1.8 ms on a broad
 * phase costing 1.3 ms. On real scenes that argument does not hold: the sweep
 * wins two of three, by 1.36x on cloth-ball and 1.59x on armadillo-rollers,
 * which is 13 ms and 78 ms rather than 1.8 ms (`wip/ASSESSMENT.md`).
 *
 * Neither constant is right and no cheap statistic has separated the cases. But
 * the two produce **identical pair sets**, so they can simply be raced: run one
 * on a step, the other on the next, then keep the winner. A broad phase runs
 * every step of a simulation, so the cost is two probe steps out of thousands,
 * and unlike a heuristic it cannot be wrong about a scene nobody tested.
 *
 * `BroadPhaseAutoTuner` does that. It re-probes periodically, because a
 * simulation's geometry changes -- cloth that starts flat and ends crumpled is
 * not one workload -- and a verdict reached on frame one should not bind frame
 * ten thousand.
 *
 * `SCCD_BROADPHASE=sweep` or `=cell2d` forces one and skips the race entirely,
 * and `broadphase_stats` is still exposed for a caller who wants to look at the
 * geometry itself.
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


    /** \brief Wall clock in milliseconds, for the race below. */
    static inline double broadphase_now_ms() {
        return std::chrono::duration<double, std::milli>(
                   std::chrono::steady_clock::now().time_since_epoch())
            .count();
    }

    /**
     * \brief Picks a broad phase by racing the two, not by predicting a winner.
     *
     * The two implementations return identical pair sets, so whichever is faster
     * on this geometry is simply the right one, and measuring that is cheaper and
     * more honest than any statistic tried so far -- see the file comment for the
     * four that failed and the constant that failed after them.
     *
     * One probe per strategy, on consecutive steps, then the winner until the
     * next re-probe. Consecutive frames of a simulation are similar enough for
     * that to be a fair race; frames far apart are not, which is why it re-probes
     * rather than deciding once and for all.
     *
     * A forced `SCCD_BROADPHASE` short-circuits it.
     */
    class BroadPhaseAutoTuner {
    public:
        /** \brief The strategy to use for the next broad phase. */
        BroadPhaseStrategy next() {
            const BroadPhaseStrategy forced = broadphase_strategy_setting();
            if (forced != BroadPhaseStrategy::Auto) return forced;

            if (reprobe_now_) {
                // Start a fresh race rather than defending the old verdict, so a
                // scene that has changed character can change the answer. The
                // warm-up is not repeated: by now everything is allocated.
                sweep_ms_ = -1.0;
                cell2d_ms_ = -1.0;
                reprobe_now_ = false;
            }

            // Each strategy runs once unmeasured before the race, because
            // otherwise it is not a race. A strategy's first execution pays for
            // allocating its own scratch -- the cell list's grid, the sweep's sort
            // buffers -- and for first touch of everything the two share. Measured
            // cold on a small mesh that showed up as 1.195 ms against 0.152 ms,
            // an eight-fold difference that was entirely allocation and would flip
            // any close race.
            //
            // Four steps before a verdict rather than two. Out of a simulation
            // that is nothing, and a caller with fewer steps than that has nothing
            // to gain from racing anyway.
            if (warmups_done_ < 2) {
                warmup_pending_ = true;
                return warmups_done_ == 0 ? BroadPhaseStrategy::Cell2D : BroadPhaseStrategy::Sweep;
            }
            // Cell2D is probed first, and which one goes first matters more than it
            // looks. A caller that runs a single broad phase on a CCD object -- one
            // query, or a driver that rebuilds the object per frame -- never
            // completes the race and only ever gets the first probe. So the first
            // probe has to be the choice that bounds the worst case, and that is
            // the cell list: the sweep's bad case on dense synthetic input was
            // measured at 54 seconds against 9, while the cell list's bad case on a
            // real scene is 78 ms against 22. Racing improves the steady state; it
            // must not make the one-shot case worse than the constant it replaced.
            if (cell2d_ms_ < 0.0) return BroadPhaseStrategy::Cell2D;
            if (sweep_ms_ < 0.0) return BroadPhaseStrategy::Sweep;
            return sweep_ms_ <= cell2d_ms_ ? BroadPhaseStrategy::Sweep : BroadPhaseStrategy::Cell2D;
        }

        /** \brief Report what the strategy `next()` handed out actually cost. */
        void record(const BroadPhaseStrategy used, const double elapsed_ms) {
            if (warmup_pending_) {
                warmup_pending_ = false;
                ++warmups_done_;
                return;
            }
            if (used == BroadPhaseStrategy::Sweep && sweep_ms_ < 0.0) {
                sweep_ms_ = elapsed_ms;
                return;
            }
            if (used == BroadPhaseStrategy::Cell2D && cell2d_ms_ < 0.0) {
                cell2d_ms_ = elapsed_ms;
                return;
            }
            // Both timed: this was a step running the winner. Count it down to the
            // next race.
            if (++steps_since_race_ >= kReprobeInterval) {
                steps_since_race_ = 0;
                reprobe_now_ = true;
            }
        }

        /** \brief The verdict, or Auto while a race is still in progress. */
        BroadPhaseStrategy decided() const {
            if (sweep_ms_ < 0.0 || cell2d_ms_ < 0.0) return BroadPhaseStrategy::Auto;
            return sweep_ms_ <= cell2d_ms_ ? BroadPhaseStrategy::Sweep : BroadPhaseStrategy::Cell2D;
        }

        double sweep_ms() const { return sweep_ms_; }
        double cell2d_ms() const { return cell2d_ms_; }

    private:
        // Long enough that two probe steps disappear into the noise, short enough
        // to follow a scene whose character changes over a run.
        static constexpr int kReprobeInterval = 64;

        double sweep_ms_ = -1.0;
        double cell2d_ms_ = -1.0;
        int steps_since_race_ = 0;
        bool reprobe_now_ = false;
        int warmups_done_ = 0;
        bool warmup_pending_ = false;
    };

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
