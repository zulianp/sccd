#!/usr/bin/env python3
"""
ctypes wrapper for the sccd shared library.
Exports:
  - find_root_vf_f(...)
  - find_root_vf_d(...)
"""
from __future__ import annotations

import os
import sys
import ctypes as ct
from typing import Iterable, Tuple


def _load_library() -> ct.CDLL:
    here = os.path.abspath(os.path.dirname(__file__))
    # 1) Environment override
    env_path = os.environ.get("SCCD_LIB_PATH")
    if env_path:
        env_path = os.path.expanduser(os.path.expandvars(env_path))
        # If a directory is provided, try platform-specific names inside it
        if os.path.isdir(env_path):
            if sys.platform == "darwin":
                env_candidates = ["libsccd.dylib", "sccd.dylib", "libsccd.so"]
            elif os.name == "nt":
                env_candidates = ["sccd.dll"]
            else:
                env_candidates = ["libsccd.so", "sccd.so"]
            for name in env_candidates:
                p = os.path.join(env_path, name)
                if os.path.exists(p):
                    return ct.CDLL(p)
        else:
            # Full path to a library file provided
            try:
                return ct.CDLL(env_path)
            except OSError:
                pass

    # 2) Same-directory search
    candidates = []
    if sys.platform == "darwin":
        candidates = ["libsccd.dylib", "sccd.dylib", "libsccd.so"]
    elif os.name == "nt":
        candidates = ["sccd.dll"]
    else:
        candidates = ["libsccd.so", "sccd.so"]
    for name in candidates:
        path = os.path.join(here, name)
        if os.path.exists(path):
            return ct.CDLL(path)
    # As a fallback, rely on system loader PATH
    for name in candidates:
        try:
            return ct.CDLL(name)
        except OSError:
            pass
    raise OSError("Could not locate the sccd shared library. Build it with CMake first.")


_lib = _load_library()

# Prototypes
_lib.sccd_find_root_vf_f.argtypes = [
    ct.c_int, ct.c_float,
    ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
    ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
    ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.POINTER(ct.c_float)
]
_lib.sccd_find_root_vf_f.restype = ct.c_int

_lib.sccd_find_root_vf_d.argtypes = [
    ct.c_int, ct.c_double,
    ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
    ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)
]
_lib.sccd_find_root_vf_d.restype = ct.c_int


def _as_arr3_f(xs: Iterable[float]) -> Tuple[ct.Array,]:
    a = (ct.c_float * 3)(*list(xs))
    return a,


def _as_arr3_d(xs: Iterable[float]) -> Tuple[ct.Array,]:
    a = (ct.c_double * 3)(*list(xs))
    return a,


def find_root_vf_f(
    max_iter: int,
    tol: float,
    sv, s1, s2, s3,
    ev, e1, e2, e3,
    t0: float = 0.0, u0: float = 0.0, v0: float = 0.0
) -> Tuple[bool, float, float, float]:
    svp = _as_arr3_f(sv)[0]
    s1p = _as_arr3_f(s1)[0]
    s2p = _as_arr3_f(s2)[0]
    s3p = _as_arr3_f(s3)[0]
    evp = _as_arr3_f(ev)[0]
    e1p = _as_arr3_f(e1)[0]
    e2p = _as_arr3_f(e2)[0]
    e3p = _as_arr3_f(e3)[0]
    t = ct.c_float(t0)
    u = ct.c_float(u0)
    v = ct.c_float(v0)
    ok = _lib.sccd_find_root_vf_f(
        int(max_iter), float(tol),
        svp, s1p, s2p, s3p, evp, e1p, e2p, e3p,
        ct.byref(t), ct.byref(u), ct.byref(v)
    )
    return (bool(ok), float(t.value), float(u.value), float(v.value))


def find_root_vf_d(
    max_iter: int,
    tol: float,
    sv, s1, s2, s3,
    ev, e1, e2, e3,
    t0: float = 0.0, u0: float = 0.0, v0: float = 0.0
) -> Tuple[bool, float, float, float]:
    svp = _as_arr3_d(sv)[0]
    s1p = _as_arr3_d(s1)[0]
    s2p = _as_arr3_d(s2)[0]
    s3p = _as_arr3_d(s3)[0]
    evp = _as_arr3_d(ev)[0]
    e1p = _as_arr3_d(e1)[0]
    e2p = _as_arr3_d(e2)[0]
    e3p = _as_arr3_d(e3)[0]
    t = ct.c_double(t0)
    u = ct.c_double(u0)
    v = ct.c_double(v0)
    ok = _lib.sccd_find_root_vf_d(
        int(max_iter), float(tol),
        svp, s1p, s2p, s3p, evp, e1p, e2p, e3p,
        ct.byref(t), ct.byref(u), ct.byref(v)
    )
    return (bool(ok), float(t.value), float(u.value), float(v.value))


if __name__ == "__main__":
    # Quick smoke test (will not guarantee correctness)
    sv = (0.25, -0.2, 0.0)
    s1 = (0.0, 0.0, 0.0)
    s2 = (1.0, 0.0, 0.0)
    s3 = (0.0, 1.0, 0.0)
    ev = (0.25, 0.8, 0.0)
    e1 = (0.0, 0.0, 0.0)
    e2 = (1.0, 0.0, 0.0)
    e3 = (0.0, 1.0, 0.0)
    ok, t, u, v = find_root_vf_d(100, 1e-10, sv, s1, s2, s3, ev, e1, e2, e3)
    print("ok=", ok, "tuv=", (t, u, v))


