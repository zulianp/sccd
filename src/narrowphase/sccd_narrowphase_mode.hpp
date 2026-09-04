#ifndef SCCD_NARROWPHASE_MODE_HPP
#define SCCD_NARROWPHASE_MODE_HPP

#include "sccd_base.hpp"

#include <cstdio>
#include <cstdlib>

namespace sccd {

    /**
     * \brief Which narrow-phase kernel to run.
     *
     * **Every mode here is conservative, and that is a property of the
     * algorithm rather than of any one kernel.** All of them are branch and
     * bound over boxes in (t, u, v) with an inclusion function, and for the
     * multilinear F these queries use, the hull of the eight corner values is
     * the exact range over a box. That makes acceptance safe however loose the
     * test -- accept reports the box's t *lower* bound, so a looser test accepts
     * earlier and reports an earlier time of impact, costing accuracy and false
     * positives but never safety. Refinement is safe for the same reason, and
     * exhaustion accepts rather than drops.
     *
     * The only way to lose a root is an unsound *rejection*: discarding a box
     * that contains one, which happens when the origin-containment test is
     * padded by less than the certified numerical error bound in
     * snumerical_error.hpp. Every mode below pads by that bound, and the oracle
     * measures it: all report gt_missed == 0 and gt_late == 0 on every dataset
     * that ships exact roots.
     *
     * So the modes differ in accuracy and speed, not in safety. If one of them
     * ever measures as non-conservative, that is a defect in that kernel -- look
     * at its rejection test, at what it reports on accept, and at what it does
     * on depth, resolution or stack exhaustion -- and not a property to document
     * and live with. See benchmark/oracle/README.md.
     */
    enum class NarrowPhaseMode : int {
        /// Scalar search with the looser acceptance test: it compares codomain
        /// widths against domain tolerances, so it accepts sooner and reports a
        /// time of impact further before the true one. Median earliness 1.18e-01
        /// on cloth-funnel. Ships.
        Fast = 0,
        /// The Fast predicate, lane-packed. Vertex-face only; edge-edge falls
        /// back to the scalar search. Validation only -- it lost every scene the
        /// assessment measured.
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
        FastVectorized = 1,
#endif
        /// Lane-packed search that compares *domain* widths against domain
        /// tolerances, which is what makes it tight rather than merely
        /// conservative: median earliness 1.72e-03 on cloth-funnel, 69x tighter
        /// than Fast, and 22x tighter on armadillo-rollers. Ships.
        Tight = 2,
        /// Runs the FastVectorized kernel and then corrects every result with
        /// the external TightInclusion library. An oracle, not a code path to
        /// ship -- and the only mode here that touches that library, which is
        /// why it is the only one still named after it.
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
        TightInclusionCorrected = 3,
#endif

        // Values 1 and 3 exist only in a validation build. They are not
        // omissions: the kernel behind them is installed only with
        // SCCD_ENABLE_TIGHT_INCLUSION, so exposing the names unconditionally
        // meant a caller could write NarrowPhaseMode::TightInclusionCorrected,
        // compile, and silently run Fast. An enum value that cannot do what it
        // says is worse than one that does not compile. The numbers stay 1 and 3
        // so SCCD_NARROWPHASE_MODE keeps meaning the same thing either way.

        // Reserved, deliberately unused: RationalMode. Nothing in SCCD computes
        // with rational arithmetic -- the dependency's rational mode is forced
        // OFF in cmake/SCCDDependencies.cmake, and the datasets' num/den
        // coordinates are parsed straight to double. If an exact-arithmetic
        // kernel is ever added it takes that name, so that no approximate kernel
        // is ever called "exact" again.
    };

    static inline const char* narrow_phase_mode_name(const NarrowPhaseMode mode) {
        switch (mode) {
            case NarrowPhaseMode::Fast: return "fast";
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
            case NarrowPhaseMode::FastVectorized: return "fast-vectorized";
#endif
            case NarrowPhaseMode::Tight: return "tight";
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
            case NarrowPhaseMode::TightInclusionCorrected: return "tight-inclusion-corrected";
#endif
        }
        return "unknown";
    }

    /**
     * \brief True when the mode runs the tight predicate and split.
     *
     * This selects an *algorithm*, not a safety property -- every mode is
     * conservative, see the note on the enum. It is what the CUDA dispatch keys
     * on to choose between its two kernels.
     */
    static inline bool narrow_phase_mode_is_tight(const NarrowPhaseMode mode) {
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
        return mode == NarrowPhaseMode::Tight || mode == NarrowPhaseMode::TightInclusionCorrected;
#else
        return mode == NarrowPhaseMode::Tight;
#endif
    }

    /**
     * \brief The mode selected for this call.
     *
     * Read per call rather than cached, because harnesses switch modes between
     * calls in the same process.
     *
     * SCCD_NARROWPHASE_MODE is the one knob. The older SCCD_USE_VNARROW_PHASE
     * and SCCD_VNARROWPHASE_TI_COMPAT still work and are consulted when it is
     * unset, so existing scripts keep behaving as they did.
     */
    /**
     * \brief Modes that exist to be compared against TightInclusion.
     *
     * `FastVector` and `TightInclusionCompat` are the two entries into the
     * vectorised kernel, and neither is a mode to ship:
     *
     *  - `TightInclusionCompat` says so itself -- it runs the vectorised kernel
     *    and then corrects its answers with TightInclusion, which is an oracle
     *    and not a code path. Without TightInclusion it returned -1 at the point
     *    of use, so it was already unavailable, just later and less clearly.
     *  - `FastVector` is a duplicate that loses. The assessment measured it
     *    behind the scalar reference on all three real scenes and behind it by
     *    5.2x on armadillo-rollers, and it is vertex-face only. What it is still
     *    good for is reproducing that comparison, which needs TightInclusion in
     *    the build anyway.
     *
     * So the vectorised kernel is validation-only: available when the build has
     * TightInclusion, and otherwise not selectable at all. Asking for one of
     * these in a build without TightInclusion gets the scalar reference rather
     * than an error, because the caller asked for a time of impact and the
     * scalar path is the supported way to get one.
     */
    /// True for the mode numbers that only a validation build provides, 1 and 3.
    /// Takes an int rather than the enum so it can be asked about a value the
    /// current build does not name.
    static inline bool narrow_phase_mode_is_validation_only(const int mode) {
        return mode == 1 || mode == 3;
    }

    /**
     * \brief Resolve the requested narrow-phase kernel, saying so when it cannot
     *        be honoured.
     *
     * The mode is read from the environment on **every call** rather than cached,
     * because callers switch it between calls -- `sccd_bench` and
     * `sccd_narrowphase_cuda_test` both do -- so the warnings below are one-shot flags
     * rather than a resolved-once result.
     *
     * Two requests used to be swallowed in silence, which is the thing worth
     * fixing here. A value that is not 0-3, or not a number at all, fell straight
     * through to the legacy variables, so `SCCD_NARROWPHASE_MODE=banana` selected
     * a kernel via `atoi` returning 0 and a typo'd `SCCD_NARROWPHASE_MODE=20` ran
     * whatever the legacy path decided. And asking for a validation-only mode in a
     * build without TightInclusion silently downgraded to the scalar reference.
     * Both are still handled the same way -- the caller asked for a time of impact
     * and gets one -- but they now say so once on stderr, because a caller who
     * measures the wrong kernel and concludes something about it is a worse
     * outcome than a line of output.
     */
    static inline NarrowPhaseMode narrow_phase_mode() {
        if (const char* explicit_mode = getenv("SCCD_NARROWPHASE_MODE")) {
            char* end = nullptr;
            const long value = strtol(explicit_mode, &end, 10);
            const bool parsed = (end != explicit_mode) && (end != nullptr) && (*end == '\0');

            if (parsed && value >= 0 && value <= 3) {
#ifndef SCCD_ENABLE_TIGHT_INCLUSION
                if (narrow_phase_mode_is_validation_only((int)value)) {
                    static bool warned_unavailable = false;
                    if (!warned_unavailable) {
                        warned_unavailable = true;
                        fprintf(stderr,
                                "SCCD: SCCD_NARROWPHASE_MODE=%ld is validation-only and needs a "
                                "build with SCCD_ENABLE_TIGHT_INCLUSION=ON; running mode 0 "
                                "(Fast) instead.\n",
                                value);
                    }
                    return NarrowPhaseMode::Fast;
                }
#endif
                return static_cast<NarrowPhaseMode>(value);
            }

            static bool warned_invalid = false;
            if (!warned_invalid) {
                warned_invalid = true;
                fprintf(stderr,
                        "SCCD: ignoring SCCD_NARROWPHASE_MODE=\"%s\" -- expected 0 (Fast), "
                        "1, 2 (Tight) or 3.\n",
                        explicit_mode);
            }
        }

        int SCCD_VNARROWPHASE_TI_COMPAT = SCCD_VNARROWPHASE_TI_COMPAT_DEFAULT;
        SCCD_READ_ENV(SCCD_VNARROWPHASE_TI_COMPAT, atoi);
        if (SCCD_VNARROWPHASE_TI_COMPAT) {
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
            return NarrowPhaseMode::TightInclusionCorrected;
#else
            return NarrowPhaseMode::Fast;
#endif
        }

        int SCCD_USE_VNARROW_PHASE = SCCD_USE_VNARROW_PHASE_DEFAULT;
        SCCD_READ_ENV(SCCD_USE_VNARROW_PHASE, atoi);
        if (SCCD_USE_VNARROW_PHASE == 2) {
            return NarrowPhaseMode::Tight;
        }
#ifdef SCCD_ENABLE_TIGHT_INCLUSION
        return SCCD_USE_VNARROW_PHASE ? NarrowPhaseMode::FastVectorized : NarrowPhaseMode::Fast;
#else
        return NarrowPhaseMode::Fast;
#endif
    }

    /**
     * \brief Say once that a requested mode does not reach the quad path.
     *
     * The quad root finder has one variant and never consults the mode enum, so
     * `SCCD_NARROWPHASE_MODE` is ignored there. That is a defensible state -- there
     * is nothing for the enum to select between -- but it was silent, and a caller
     * comparing modes on a quad mesh would have measured the same kernel twice and
     * concluded the modes were equivalent.
     *
     * Only a request for a non-default mode is worth a line: asking for mode 0 and
     * getting the one kernel there is is not a surprise. A malformed value is left
     * alone here, since `narrow_phase_mode` has already reported it.
     */
    static inline void narrow_phase_mode_note_quads_ignore() {
        const char* explicit_mode = getenv("SCCD_NARROWPHASE_MODE");
        if (explicit_mode == nullptr) return;

        char* end = nullptr;
        const long value = strtol(explicit_mode, &end, 10);
        if (end == explicit_mode || end == nullptr || *end != '\0') return;
        if (value <= 0 || value > 3) return;

        static bool warned = false;
        if (warned) return;
        warned = true;
        fprintf(stderr,
                "SCCD: SCCD_NARROWPHASE_MODE=%ld does not apply to quads -- the vertex-quad "
                "root finder has one variant and does not consult the mode.\n",
                value);
    }

}  // namespace sccd

#endif  // SCCD_NARROWPHASE_MODE_HPP
