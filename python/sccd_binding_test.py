#!/usr/bin/env python3
"""Check that python/sccd.py agrees with the C ABI it wraps.

The binding had no test at all, and it had drifted: its docstring advertised
three functions where six existed, and it covered six of the C ABI's seven
exports. Neither would have been caught by src/tests/sccd_c_abi_test.exe.cpp,
which calls the same exports from C++ and so shares none of this file's code.

The cases here are the ones that test uses, deliberately: same geometry, same
expected times of impact. If the two disagree, the binding is the thing that
changed.

Run directly, or through ctest as `sccd_binding_test`:

    SCCD_LIB_PATH=build python3 python/sccd_binding_test.py

Not registered in a sanitizer build: the runtime has to be loaded before main and
ctypes loads it after, so an uninstrumented interpreter cannot dlopen the library
at all. That is a property of the sanitizer, not of the binding.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sccd  # noqa: E402

failures = []


def check(condition, message):
    if condition:
        print(f"  ok    {message}")
    else:
        print(f"  FAIL  {message}")
        failures.append(message)


# A vertex dropping through the triangle (0,0,0), (1,0,0), (0,1,0), which lies in
# z = 0 and does not move. The crossing time is where the vertex's z reaches
# zero, so it can be written down.
TRIANGLE = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))

VF_CASES = (
    # x,    y,    z0,   z1,   exact, hit,   name
    (0.25, 0.25,  1.0, -1.0,  0.500, True,  "crosses at t=1/2, inside"),
    (0.20, 0.30,  1.0, -3.0,  0.250, True,  "crosses at t=1/4, inside"),
    (0.90, 0.90,  1.0, -1.0,  0.000, False, "crosses outside the triangle"),
    (0.25, 0.25,  1.0,  0.5,  0.000, False, "never reaches the plane"),
)

# Two segments: one static along x at z = 0, one along y that starts above it and
# drops through, optionally offset in y so it passes to the side.
EE_STATIC = ((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0))

EE_CASES = (
    # z0,   z1,  y offset, exact, hit,   name
    ( 1.0, -1.0, 0.0,      0.500, True,  "crosses at t=1/2"),
    ( 1.0, -3.0, 0.0,      0.250, True,  "crosses at t=1/4"),
    ( 1.0, -1.0, 5.0,      0.000, False, "passes far to the side"),
    ( 1.0,  0.5, 0.0,      0.000, False, "never reaches the plane"),
)


def conservative(root, exact, hit, label, slack):
    """A reported time of impact must be at or before the true one.

    Late is a correctness failure. Early is safe, but a wildly early answer means
    the search is not converging, so that is flagged as well.
    """
    if bool(root) != hit:
        check(False, f"{label}: expected {'a hit' if hit else 'a miss'}, got the other")
        return
    if not hit:
        check(True, f"{label}: miss, as expected")
        return
    err = root.t - exact
    if err > slack:
        check(False, f"{label}: LATE, {root.t:.9f} against a true {exact:.6f}")
    elif err < -0.05:
        check(False, f"{label}: far too early, {root.t:.9f} against a true {exact:.6f}")
    else:
        check(True, f"{label}: toi={root.t:.9f}  err={err:+.2e}")


def main():
    print("loading the library")
    lib = sccd.load_library()
    check(lib is not None, f"loaded {sccd.library_path()}")

    # Every export the C ABI documents is either present, or absent only because
    # this build has no TightInclusion. A missing export outside that pair means
    # the binding and the library have parted company.
    have = sccd.available()
    print("\nexports")
    for name, present in have.items():
        optional = "tight_inclusion" in name
        if present:
            check(True, f"{name}")
        else:
            check(optional, f"{name} absent"
                            + (" (no SCCD_ENABLE_TIGHT_INCLUSION, as allowed)" if optional
                               else " -- the library does not export it"))

    print("\nvertex-face, double")
    for x, y, z0, z1, exact, hit, name in VF_CASES:
        root = sccd.find_root_vf(69, 1e-8, (x, y, z0), *TRIANGLE, (x, y, z1), *TRIANGLE)
        conservative(root, exact, hit, name, slack=1e-5)

    print("\nvertex-face, float")
    for x, y, z0, z1, exact, hit, name in VF_CASES:
        root = sccd.find_root_vf(69, 1e-6, (x, y, z0), *TRIANGLE, (x, y, z1), *TRIANGLE,
                                 dtype="float")
        conservative(root, exact, hit, name, slack=1e-5)

    print("\nvertex-face, bisecting splitter")
    for x, y, z0, z1, exact, hit, name in VF_CASES:
        for dtype, tol in (("double", 1e-8), ("float", 1e-6)):
            root = sccd.find_root_bisection_vf(69, tol, (x, y, z0), *TRIANGLE,
                                               (x, y, z1), *TRIANGLE, dtype=dtype)
            conservative(root, exact, hit, f"{name} [{dtype}]", slack=1e-5)

    print("\nedge-edge, double")
    a0, a1 = EE_STATIC
    for z0, z1, yoff, exact, hit, name in EE_CASES:
        b0s, b1s = (0.0, -1.0 + yoff, z0), (0.0, 1.0 + yoff, z0)
        b0e, b1e = (0.0, -1.0 + yoff, z1), (0.0, 1.0 + yoff, z1)
        root = sccd.find_root_ee(69, 1e-8, a0, a1, b0s, b1s, a0, a1, b0e, b1e)
        conservative(root, exact, hit, name, slack=1e-5)

    print("\nthe binding itself")
    miss = sccd.find_root_vf(69, 1e-8, (0.25, 0.25, 1.0), *TRIANGLE,
                             (0.25, 0.25, 0.5), *TRIANGLE)
    check(not miss, "a miss is falsey, so `if find_root_vf(...)` reads correctly")
    hit_root = sccd.find_root_vf(69, 1e-8, (0.25, 0.25, 1.0), *TRIANGLE,
                                 (0.25, 0.25, -1.0), *TRIANGLE)
    check(bool(hit_root), "a hit is truthy")
    check(tuple(hit_root)[1:] == (hit_root.t, hit_root.u, hit_root.v),
          "Root unpacks as (hit, t, u, v)")
    try:
        sccd.find_root_vf(69, 1e-8, (0.0, 0.0), *TRIANGLE, (0.0, 0.0, 1.0), *TRIANGLE)
        check(False, "a two-coordinate point is rejected")
    except ValueError:
        check(True, "a two-coordinate point is rejected")
    try:
        sccd.find_root_vf(69, 1e-8, (0.0, 0.0, 0.0), *TRIANGLE, (0.0, 0.0, 1.0), *TRIANGLE,
                          dtype="half")
        check(False, "an unknown dtype is rejected")
    except ValueError:
        check(True, "an unknown dtype is rejected")

    print()
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        return 1
    print("OK: the binding agrees with the C ABI, and is conservative on every case")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
