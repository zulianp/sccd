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
        /// on cloth-funnel.
        Relaxed = 0,
        /// Lane-packed search comparing *domain* widths against domain
        /// tolerances, which is what makes it tight rather than merely
        /// conservative: median earliness 1.72e-03 on cloth-funnel, 69x tighter
        /// than Relaxed, and 22x tighter on armadillo-rollers.
        Tight = 2

        // 1 and 3 are not defined, and the gap is deliberate.
        //
        // Mode 1 was the vectorised form of the Relaxed predicate. It was the
        // slowest kernel in the library on every scene measured -- 15.9 ms
        // against 6.4 and 8.1 on cloth-funnel, 90.9 against 19.4 and 17.5 on
        // armadillo-rollers -- so its name promised the opposite of what it did.
        // It is in spikes/src/sccd_narrowphase_fast_vector.hpp.
        //
        // Mode 3 ran mode 1 and then corrected each answer with TightInclusion.
        // If you want TightInclusion's answer, call TightInclusion: set
        // SCCD_USE_TI=1, which dispatches straight to the library, or use the
        // sccd_find_root_tight_inclusion_* entry points. A hybrid that is
        // neither this library's kernel nor the reference is not worth shipping
        // or maintaining.
        //
        // Reserved, still unused: RationalMode. Nothing here computes with
        // rational arithmetic; the name is held so no approximate kernel is ever
        // called "exact".
    };


    static inline const char* narrow_phase_mode_name(const NarrowPhaseMode mode) {
        switch (mode) {
            case NarrowPhaseMode::Relaxed: return "relaxed";
            case NarrowPhaseMode::Tight: return "tight";
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
        return mode == NarrowPhaseMode::Tight;
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

            if (parsed && (value == 0 || value == 2)) {
                return static_cast<NarrowPhaseMode>(value);
            }

            if (parsed && (value == 1 || value == 3)) {
                static bool warned_removed = false;
                if (!warned_removed) {
                    warned_removed = true;
                    fprintf(stderr,
                            "SCCD: SCCD_NARROWPHASE_MODE=%ld no longer exists; running 0 (Relaxed). "
                            "Mode 1 was the slowest kernel here and is now in spikes/. Mode 3 "
                            "corrected mode 1 with TightInclusion -- set SCCD_USE_TI=1 to call "
                            "that library directly instead.\n",
                            value);
                }
                return NarrowPhaseMode::Relaxed;
            }

            static bool warned_invalid = false;
            if (!warned_invalid) {
                warned_invalid = true;
                fprintf(stderr,
                        "SCCD: ignoring SCCD_NARROWPHASE_MODE=\"%s\" -- expected 0 (Relaxed) or "
                        "2 (Tight).\n",
                        explicit_mode);
            }
        }

        // The two pre-enum switches. Both selected kernels that no longer exist
        // here, so both are dead letters now; they warn rather than silently
        // doing something else, and SCCD_USE_VNARROW_PHASE=2 still means Tight
        // because that mode is still real.
        int SCCD_USE_VNARROW_PHASE = SCCD_USE_VNARROW_PHASE_DEFAULT;
        SCCD_READ_ENV(SCCD_USE_VNARROW_PHASE, atoi);
        if (SCCD_USE_VNARROW_PHASE == 2) {
            return NarrowPhaseMode::Tight;
        }

        int SCCD_VNARROWPHASE_TI_COMPAT = SCCD_VNARROWPHASE_TI_COMPAT_DEFAULT;
        SCCD_READ_ENV(SCCD_VNARROWPHASE_TI_COMPAT, atoi);
        if (SCCD_VNARROWPHASE_TI_COMPAT || SCCD_USE_VNARROW_PHASE) {
            static bool warned_legacy = false;
            if (!warned_legacy) {
                warned_legacy = true;
                fprintf(stderr,
                        "SCCD: SCCD_VNARROWPHASE_TI_COMPAT and SCCD_USE_VNARROW_PHASE selected "
                        "kernels that have moved to spikes/; running 0 (Relaxed). Use "
                        "SCCD_NARROWPHASE_MODE=2 for the tight kernel, or SCCD_USE_TI=1 to call "
                        "TightInclusion directly.\n");
            }
        }
        return NarrowPhaseMode::Relaxed;
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
