#ifndef SCCD_NARROWPHASE_MODE_HPP
#define SCCD_NARROWPHASE_MODE_HPP

#include "sccd_base.hpp"

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
        /// Scalar reference search. Accepts on several conditions beyond
        /// TightInclusion's, so it stops earlier and reports a less accurate
        /// (earlier) time of impact.
        ScalarReference = 0,
        /// Lane-packed vectorized kernel, vertex-face only (edge-edge falls back
        /// to the scalar search).
        FastVector = 1,
        /// Lane-packed vectorized kernel reproducing TightInclusion exactly, for
        /// both vertex-face and edge-edge. The most accurate of the three, and
        /// the reference the others are compared against.
        TightInclusionExact = 2,
        /// Runs the fast kernel and then corrects every result with
        /// TightInclusion's own search. An oracle, not a code path to ship.
        TightInclusionCompat = 3
    };

    /// Kept so existing callers and env values keep compiling and meaning the
    /// same thing. The old name implied the other modes were not conservative.
    static constexpr NarrowPhaseMode NarrowPhaseMode_Conservative = NarrowPhaseMode::TightInclusionExact;

    static inline const char* narrow_phase_mode_name(const NarrowPhaseMode mode) {
        switch (mode) {
            case NarrowPhaseMode::ScalarReference: return "scalar";
            case NarrowPhaseMode::FastVector: return "fast-vector";
            case NarrowPhaseMode::TightInclusionExact: return "ti-exact";
            case NarrowPhaseMode::TightInclusionCompat: return "ti-compat";
        }
        return "unknown";
    }

    /**
     * \brief True when the mode runs TightInclusion's exact predicate and split.
     *
     * This selects an *algorithm*, not a safety property -- every mode is
     * conservative, see the note on the enum. It is what the CUDA dispatch keys
     * on to choose between its two kernels.
     */
    static inline bool narrow_phase_mode_is_ti_exact(const NarrowPhaseMode mode) {
        return mode == NarrowPhaseMode::TightInclusionExact || mode == NarrowPhaseMode::TightInclusionCompat;
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
    static inline NarrowPhaseMode narrow_phase_mode() {
        if (const char* explicit_mode = getenv("SCCD_NARROWPHASE_MODE")) {
            const int value = atoi(explicit_mode);
            if (value >= 0 && value <= 3) {
                return static_cast<NarrowPhaseMode>(value);
            }
        }

        int SCCD_VNARROWPHASE_TI_COMPAT = SCCD_VNARROWPHASE_TI_COMPAT_DEFAULT;
        SCCD_READ_ENV(SCCD_VNARROWPHASE_TI_COMPAT, atoi);
        if (SCCD_VNARROWPHASE_TI_COMPAT) {
            return NarrowPhaseMode::TightInclusionCompat;
        }

        int SCCD_USE_VNARROW_PHASE = SCCD_USE_VNARROW_PHASE_DEFAULT;
        SCCD_READ_ENV(SCCD_USE_VNARROW_PHASE, atoi);
        if (SCCD_USE_VNARROW_PHASE == 2) {
            return NarrowPhaseMode::TightInclusionExact;
        }
        return SCCD_USE_VNARROW_PHASE ? NarrowPhaseMode::FastVector : NarrowPhaseMode::ScalarReference;
    }

}  // namespace sccd

#endif  // SCCD_NARROWPHASE_MODE_HPP
