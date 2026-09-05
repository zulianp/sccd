#ifndef SCCD_BASE_HPP
#define SCCD_BASE_HPP

#include "sccd_config.hpp"

#ifndef SCCD_READ_ENV
#define SCCD_READ_ENV(name, conversion) \
    do {                                \
        char* var = getenv(#name);      \
        if (var) {                      \
            name = conversion(var);     \
        }                               \
    } while (0)
#endif

#ifndef SCCD_RESTRICT
#ifndef _WIN32
#define SCCD_RESTRICT __restrict__
#else
#define SCCD_RESTRICT __restrict
#endif
#endif

#ifndef SCCD_ALWAYS_INLINE
#if defined(_MSC_VER)
#define SCCD_ALWAYS_INLINE __forceinline
#else
#define SCCD_ALWAYS_INLINE inline __attribute__((always_inline))
#endif
#endif

// Pin IEEE semantics across a region of function definitions.
//
// TightInclusion's numerical error bound is derived for one specific sequence of
// floating-point operations, and the acceptance predicate is only sound because
// it pads by that bound. -ffast-math licenses the compiler to reassociate that
// sequence, which voids the derivation, and -ffinite-math-only lets it assume
// the infinities that snumtol.hpp's tolerance clamp exists to catch cannot
// occur. Neither assumption holds here.
//
// The macros below remove the licence; the #warning catches the case where a
// consumer's flags reintroduce it. Note that putting -Ofast in CMAKE_CXX_FLAGS
// does NOT reintroduce it, because CMAKE_CXX_FLAGS_RELEASE appends -O3 after it
// and the last -O wins -- which is why the warning stays silent on a build that
// looks like it should trip it. Setting -ffast-math explicitly, or putting it in
// CMAKE_CXX_FLAGS_RELEASE, does trip it.
//
// clang applies float_control per instruction, so a strict function keeps its
// semantics even when inlined into a fast-math caller -- which is what makes
// this usable on SCCD_ALWAYS_INLINE helpers. GCC models the equivalent as a
// function optimization attribute, and mixing those with always_inline makes it
// refuse the inline outright, so there the macros stay empty and the warning
// below is the whole defence.
#if defined(__clang__)
#define SCCD_FP_STRICT_BEGIN         \
    _Pragma("float_control(push)")   \
    _Pragma("float_control(precise, on)")
#define SCCD_FP_STRICT_END _Pragma("float_control(pop)")
#else
#define SCCD_FP_STRICT_BEGIN
#define SCCD_FP_STRICT_END
#endif

#if defined(__FAST_MATH__) && !defined(SCCD_ALLOW_FAST_MATH)
#if defined(__clang__) || defined(__GNUC__)
#warning \
    "SCCD: built with -ffast-math (-Ofast). The narrow phase's conservativeness rests on an error bound derived for an exact operation order and on infinities reaching the tolerance clamp; fast-math voids both. Prefer -O3. Define SCCD_ALLOW_FAST_MATH to silence."
#endif
#endif

// Force the next loop through the vectorizer. Needed where a short, known-length
// lane loop gets fully unrolled and the cost model then declines to vectorize it,
// which is how this library's hottest kernel ended up emitting no vector
// instructions at all. #pragma omp simd cannot serve here: it is inert unless the
// consumer compiles with OpenMP, and correctness must not depend on that.
#ifndef SCCD_VECTORIZE_LOOP
#if defined(__clang__)
#define SCCD_VECTORIZE_LOOP _Pragma("clang loop vectorize(enable)")
#elif defined(__GNUC__)
#define SCCD_VECTORIZE_LOOP _Pragma("GCC ivdep")
#else
#define SCCD_VECTORIZE_LOOP
#endif
#endif

/*
 * Return convention for the C++ entry points.
 *
 * They return an int STATUS, zero on success -- the convention
 * narrow_phase_tight_vf follows, and the one every narrow-phase entry point in
 * the library uses. The smesh CCD<T> methods report the same thing through an
 * `int& err` compared against SCCD_SUCCESS.
 *
 * The C ABI in sccd.h is deliberately different and says so at its own
 * declarations: there the int is a HIT, non-zero meaning a collision was found.
 * Do not read one as the other.
 */
#define SCCD_SUCCESS 0
#define SCCD_FAILURE 1

#endif
