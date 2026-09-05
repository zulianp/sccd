#ifndef SCCD_H
#define SCCD_H

/*
 * SCCD's C ABI: single-query continuous collision detection root finders.
 *
 * This header is the one declaration of these symbols. Before it existed they
 * were declared three times -- in the definitions, in SCCD's own ABI test, and
 * again through ctypes in python/sccd.py -- with nothing checking any copy
 * against another, so a changed argument broke callers at run time rather than
 * at compile time.
 *
 * Scope: one query per call. There is no broad phase here and no batching; those
 * live in the C++ headers (sccd_narrowphase.hpp, sccd_broadphase_sweep.hpp).
 *
 * Link against `sccd`. sccd_config.hpp is generated at configure time and says
 * which optional entry points the library was built with.
 */

#include "sccd_config.hpp"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * ## The shape every entry point shares
 *
 *     int sccd_find_root_<kind>(int max_iter, T tol,
 *                               const T a0[3], a1[3], a2[3], a3[3],   // step start
 *                               const T b0[3], b1[3], b2[3], b3[3],   // step end
 *                               T *t, T *u, T *v);
 *
 * The first four points are the configuration at the start of the step and the
 * next four the configuration at the end.
 *
 * For vertex-face the points are the vertex and the face's three corners, and
 * (u, v) are barycentric on the face. For edge-edge they are the two endpoints
 * of each edge, and u and v parameterise the first and second edge.
 *
 * ## The return value is a hit, not a status
 *
 * Non-zero means a collision was found; zero means none exists in the step.
 * There is no error code: these functions do not validate their arguments, and
 * a null pointer or a non-positive max_iter is undefined behaviour. Note this
 * differs from the C++ entry points, whose int is a status.
 *
 * ## *t is in AND out, which is easy to miss
 *
 * On entry *t is the upper bound of the search: the root finder ignores any
 * contact at or after it, and a box whose t lower bound is at or past *t is
 * discarded without being opened. **Pass a value greater than the latest time
 * you care about** -- 1.1 for a full step, which is what python/sccd.py sends.
 * Passing 0, or uninitialised memory, searches nothing and returns 0.
 *
 * On return, if the call reported a hit, *t holds the time of impact and
 * *u, *v the parameter coordinates of the contact. *u and *v are likewise read
 * on entry and should be initialised, conventionally to 0.
 *
 * ## What the time of impact means
 *
 * *t lies in [0, 1] and is **conservative**: at or before the true time of
 * impact, never after. Reporting early costs a shorter step; reporting late
 * lets a simulation pass through the contact, which is what the whole search is
 * built to prevent. A returned "no collision" means no collision exists --
 * misses are not traded against speed.
 *
 * max_iter caps the subdivision depth and tol sets the codomain tolerance.
 * Tightening either sharpens *t and costs time; neither can make the answer
 * unsafe, because the rejection test pads by a certified numerical error bound
 * that neither touches.
 */

/* Vertex-face, the shipped adaptive-split search. */
int sccd_find_root_vf_f(int max_iter,
                        float tol,
                        const float sv[3],
                        const float s1[3],
                        const float s2[3],
                        const float s3[3],
                        const float ev[3],
                        const float e1[3],
                        const float e2[3],
                        const float e3[3],
                        float* t,
                        float* u,
                        float* v);

int sccd_find_root_vf_d(int max_iter,
                        double tol,
                        const double sv[3],
                        const double s1[3],
                        const double s2[3],
                        const double s3[3],
                        const double ev[3],
                        const double e1[3],
                        const double e2[3],
                        const double e3[3],
                        double* t,
                        double* u,
                        double* v);

/* Edge-edge. */
int sccd_find_root_ee_d(int max_iter,
                        double tol,
                        const double s0[3],
                        const double s1[3],
                        const double s2[3],
                        const double s3[3],
                        const double e0[3],
                        const double e1[3],
                        const double e2[3],
                        const double e3[3],
                        double* t,
                        double* u,
                        double* v);

/*
 * Vertex-face with a bisecting splitter instead of the adaptive one. A second
 * splitter, reachable only from here: no C++ entry point and no environment
 * variable selects it.
 */
int sccd_find_root_bisection_vf_f(int max_iter,
                                  float tol,
                                  const float sv[3],
                                  const float s1[3],
                                  const float s2[3],
                                  const float s3[3],
                                  const float ev[3],
                                  const float e1[3],
                                  const float e2[3],
                                  const float e3[3],
                                  float* t,
                                  float* u,
                                  float* v);

int sccd_find_root_bisection_vf_d(int max_iter,
                                  double tol,
                                  const double sv[3],
                                  const double s1[3],
                                  const double s2[3],
                                  const double s3[3],
                                  const double ev[3],
                                  const double e1[3],
                                  const double e2[3],
                                  const double e3[3],
                                  double* t,
                                  double* u,
                                  double* v);

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
/*
 * The external TightInclusion library's answer, for validation rather than use.
 * Present only in a build configured with SCCD_ENABLE_TIGHT_INCLUSION=ON, which
 * is why these are behind the macro: a caller gets a compile error naming the
 * missing option rather than an undefined symbol at link time.
 */
int sccd_find_root_tight_inclusion_vf_d(int max_iter,
                                        double tol,
                                        const double sv[3],
                                        const double s1[3],
                                        const double s2[3],
                                        const double s3[3],
                                        const double ev[3],
                                        const double e1[3],
                                        const double e2[3],
                                        const double e3[3],
                                        double* t,
                                        double* u,
                                        double* v);

int sccd_find_root_tight_inclusion_ee_d(int max_iter,
                                        double tol,
                                        const double s0[3],
                                        const double s1[3],
                                        const double s2[3],
                                        const double s3[3],
                                        const double e0[3],
                                        const double e1[3],
                                        const double e2[3],
                                        const double e3[3],
                                        double* t,
                                        double* u,
                                        double* v);
#endif /* SCCD_ENABLE_TIGHT_INCLUSION */

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* SCCD_H */
