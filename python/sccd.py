#!/usr/bin/env python3
"""ctypes binding for SCCD's C ABI.

SCCD's only compiled translation unit is ``src/api/sccd_c_api.cpp``, and this is
a thin wrapper over it. Each entry point takes a single query -- eight points,
the configuration at the start of the step and at the end -- runs the
branch-and-bound search on it, and returns the time of impact together with the
parameter coordinates of the contact. There is no broad phase here and no
batching; for those, call the C++ headers.

    import sccd
    hit = sccd.find_root_vf(max_iter=100, tol=1e-10,
                            sv=(0.25, -0.2, 0.0),
                            s1=(0.0, 0.0, 0.0), s2=(1.0, 0.0, 0.0), s3=(0.0, 1.0, 0.0),
                            ev=(0.25, 0.8, 0.0),
                            e1=(0.0, 0.0, 0.0), e2=(1.0, 0.0, 0.0), e3=(0.0, 1.0, 0.0))
    if hit:
        print(hit.t, hit.u, hit.v)

The reported time of impact is conservative: at or before the true one, never
after. ``hit`` is falsey when no collision exists in the step.

Finding the shared library
--------------------------
``SCCD_LIB_PATH`` wins if set, and may name either the library file or a
directory holding it. Otherwise this module looks beside itself, then falls back
to the system loader. The library is loaded on first use rather than on import,
so importing this module never fails for want of a build.
"""

from __future__ import annotations

import ctypes as ct
import os
import sys
from typing import Iterable, NamedTuple, Optional, Sequence

__all__ = [
    "Root",
    "find_root_vf",
    "find_root_ee",
    "find_root_bisection_vf",
    "find_root_tight_inclusion_vf",
    "find_root_tight_inclusion_ee",
    "available",
    "library_path",
    "load_library",
]


class Root(NamedTuple):
    """The result of one query.

    ``hit`` is False when no collision exists in the step, in which case ``t``,
    ``u`` and ``v`` hold whatever was passed in as the initial guess and mean
    nothing. The tuple is falsey in that case, so ``if find_root_vf(...)`` reads
    the way it should.
    """

    hit: bool
    t: float
    u: float
    v: float

    def __bool__(self) -> bool:
        return self.hit


# --- the library -----------------------------------------------------------

_lib: Optional[ct.CDLL] = None
_lib_path: Optional[str] = None
_prototyped: set = set()


def _candidate_names() -> Sequence[str]:
    if sys.platform == "darwin":
        return ("libsccd.dylib", "sccd.dylib", "libsccd.so")
    if os.name == "nt":
        return ("sccd.dll",)
    return ("libsccd.so", "sccd.so")


def _try_open(path: str) -> Optional[ct.CDLL]:
    try:
        return ct.CDLL(path)
    except OSError:
        return None


def load_library(path: Optional[str] = None) -> ct.CDLL:
    """Load the shared library, caching the handle.

    Called automatically by every entry point below; call it directly only to
    force a particular build, or to fail early with a clear message.
    """
    global _lib, _lib_path
    if _lib is not None and path is None:
        return _lib

    tried = []
    search = []
    for root in (path, os.environ.get("SCCD_LIB_PATH")):
        if not root:
            continue
        root = os.path.expanduser(os.path.expandvars(root))
        if os.path.isdir(root):
            search.extend(os.path.join(root, name) for name in _candidate_names())
        else:
            search.append(root)
        break  # an explicit path is not a hint; do not fall through past it

    here = os.path.abspath(os.path.dirname(__file__))
    search.extend(os.path.join(here, name) for name in _candidate_names())
    search.extend(_candidate_names())  # the system loader's own search path

    for candidate in search:
        tried.append(candidate)
        lib = _try_open(candidate)
        if lib is not None:
            _lib, _lib_path = lib, candidate
            _prototyped.clear()
            return lib

    raise OSError(
        "could not load the SCCD shared library. Build it with CMake, then set "
        "SCCD_LIB_PATH to the build directory or to the library file. Tried:\n  "
        + "\n  ".join(tried)
    )


def library_path() -> Optional[str]:
    """Path the library was loaded from, or None if it has not been loaded."""
    return _lib_path


# --- calling one entry point ------------------------------------------------
#
# Every export has the same shape, so there is one implementation rather than one
# per export:
#
#   int sccd_find_root_<kind>(int max_iter, T tol,
#                             const T a0[3], ... const T b3[3],
#                             T* t, T* u, T* v)


def _prototype(name: str, scalar):
    if name in _prototyped:
        return getattr(_lib, name)
    try:
        func = getattr(_lib, name)
    except AttributeError as exc:
        raise AttributeError(
            f"{name} is not exported by {_lib_path}. The TightInclusion entry "
            f"points need a build with SCCD_ENABLE_TIGHT_INCLUSION=ON."
        ) from exc
    func.argtypes = [ct.c_int, scalar] + [ct.POINTER(scalar)] * 11
    func.restype = ct.c_int
    _prototyped.add(name)
    return func


def _point(scalar, xs: Iterable[float]):
    values = list(xs)
    if len(values) != 3:
        raise ValueError(f"expected 3 coordinates, got {len(values)}")
    return (scalar * 3)(*values)


def _call(name, scalar, max_iter, tol, points, t0, u0, v0) -> Root:
    load_library()
    func = _prototype(name, scalar)
    args = [_point(scalar, p) for p in points]
    t, u, v = scalar(t0), scalar(u0), scalar(v0)
    hit = func(int(max_iter), scalar(tol), *args, ct.byref(t), ct.byref(u), ct.byref(v))
    return Root(bool(hit), float(t.value), float(u.value), float(v.value))


def _scalar_for(dtype: str):
    if dtype in ("d", "double", "float64"):
        return ct.c_double, "d"
    if dtype in ("f", "float", "float32"):
        return ct.c_float, "f"
    raise ValueError(f"dtype must be 'double' or 'float', got {dtype!r}")


# The initial (t, u, v) is what comes back when the search finds no root. 1.1 is
# outside the unit step on purpose, so a caller that ignores `hit` gets an
# obviously invalid time rather than a plausible one.
_NO_ROOT = (1.1, 0.0, 0.0)


def find_root_vf(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3,
                 dtype="double", t0=_NO_ROOT[0], u0=_NO_ROOT[1], v0=_NO_ROOT[2]) -> Root:
    """Vertex against triangle.

    ``sv`` and ``ev`` are the vertex at the start and end of the step; ``s1..s3``
    and ``e1..e3`` are the triangle's corners. ``u`` and ``v`` are barycentric on
    the triangle. ``dtype`` selects the ``_d`` or ``_f`` export.
    """
    scalar, suffix = _scalar_for(dtype)
    return _call(f"sccd_find_root_vf_{suffix}", scalar, max_iter, tol,
                 (sv, s1, s2, s3, ev, e1, e2, e3), t0, u0, v0)


def find_root_ee(max_iter, tol, s0, s1, s2, s3, e0, e1, e2, e3,
                 dtype="double", t0=_NO_ROOT[0], u0=_NO_ROOT[1], v0=_NO_ROOT[2]) -> Root:
    """Edge against edge.

    The first four points are the two edges' endpoints at the start of the step
    and the last four at the end; ``u`` and ``v`` parameterise the first and the
    second edge. Only the ``double`` export exists.
    """
    scalar, suffix = _scalar_for(dtype)
    return _call(f"sccd_find_root_ee_{suffix}", scalar, max_iter, tol,
                 (s0, s1, s2, s3, e0, e1, e2, e3), t0, u0, v0)


def find_root_bisection_vf(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3,
                           dtype="double", t0=_NO_ROOT[0], u0=_NO_ROOT[1], v0=_NO_ROOT[2]) -> Root:
    """Vertex against triangle, with the bisecting splitter instead of the adaptive one."""
    scalar, suffix = _scalar_for(dtype)
    return _call(f"sccd_find_root_bisection_vf_{suffix}", scalar, max_iter, tol,
                 (sv, s1, s2, s3, ev, e1, e2, e3), t0, u0, v0)


def find_root_tight_inclusion_vf(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3,
                                 t0=_NO_ROOT[0], u0=_NO_ROOT[1], v0=_NO_ROOT[2]) -> Root:
    """Vertex against triangle, answered by the external TightInclusion library.

    Present only in a build with ``SCCD_ENABLE_TIGHT_INCLUSION=ON``. This is the
    accuracy oracle, not a code path to use.
    """
    return _call("sccd_find_root_tight_inclusion_vf_d", ct.c_double, max_iter, tol,
                 (sv, s1, s2, s3, ev, e1, e2, e3), t0, u0, v0)


def find_root_tight_inclusion_ee(max_iter, tol, s0, s1, s2, s3, e0, e1, e2, e3,
                                 t0=_NO_ROOT[0], u0=_NO_ROOT[1], v0=_NO_ROOT[2]) -> Root:
    """Edge against edge, answered by the external TightInclusion library.

    Present only in a build with ``SCCD_ENABLE_TIGHT_INCLUSION=ON``.
    """
    return _call("sccd_find_root_tight_inclusion_ee_d", ct.c_double, max_iter, tol,
                 (s0, s1, s2, s3, e0, e1, e2, e3), t0, u0, v0)


# The seven exports of the C ABI, in the order docs/API.md lists them.
_EXPORTS = (
    "sccd_find_root_vf_f",
    "sccd_find_root_vf_d",
    "sccd_find_root_ee_d",
    "sccd_find_root_bisection_vf_f",
    "sccd_find_root_bisection_vf_d",
    "sccd_find_root_tight_inclusion_vf_d",
    "sccd_find_root_tight_inclusion_ee_d",
)


def available() -> dict:
    """Which of the C ABI's exports the loaded library actually has.

    The two TightInclusion entry points are absent unless the build had
    ``SCCD_ENABLE_TIGHT_INCLUSION=ON``, so check here rather than catching an
    AttributeError from the call.
    """
    lib = load_library()
    return {name: hasattr(lib, name) for name in _EXPORTS}


if __name__ == "__main__":
    lib = load_library()
    print(f"loaded {library_path()}")
    for name, present in available().items():
        print(f"  {'yes' if present else ' no'}  {name}")

    # A vertex crossing the plane of a stationary triangle. It starts at
    # y = -0.2 and ends at y = 0.8, so it reaches y = 0 one fifth of the way
    # through the step, at (0.25, 0) -- inside the triangle.
    tri = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    hit = find_root_vf(100, 1e-10, (0.25, -0.2, 0.0), *tri, (0.25, 0.8, 0.0), *tri)
    print(f"\nvertex-face: hit={hit.hit} t={hit.t:.9f} u={hit.u:.6f} v={hit.v:.6f}")
    print(f"exact t is 0.2; reported {'at or before' if hit.t <= 0.2 else 'AFTER'} it")
    raise SystemExit(0 if hit and hit.t <= 0.2 else 1)
