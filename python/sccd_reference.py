#!/usr/bin/env python3
"""A NumPy reference search, used to cross-check the C ABI's answers.

`find_root_dfs_3D` is a depth-first branch and bound over boxes in (t, u, v),
written to be read rather than to be fast: it is the independent implementation
that python/sccd_crosscheck.py compares the library against, so it deliberately
shares no code with it. `vf_F_3d` and `ee_F_3d` are the residuals themselves,
also used to print the residual at a disputed root.

Numba is optional. Without it `myjit` is the identity and everything still runs,
about two orders of magnitude slower.

The 2D helpers and the 3D visualisations that used to live here had no caller in
the shipped tree -- 958 of 1344 lines, reachable only from two demoted scripts --
and now live with those scripts in spikes/python/.
"""

import numpy as np
from math import floor
from typing import Tuple, Optional

try:
    from numba import njit

    def myjit(f):
        return njit(f, fastmath=True, boundscheck=False, nogil=True)
except ImportError:
    def myjit(f):
        return f

# Keep the chunk size similar to the C++ heuristic.
ROOT_FINDING_CHUNK_SIZE = 4 * 4096


@myjit
def vf_F_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, tt, uu, vv):
    """
    3D vertex-face residual:
      F(t,u,v) = V(t) - ((1-u-v) * V1(t) + u * V2(t) + v * V3(t))
    Inputs are 3D sequences (len==3). Returns tuple (Fx, Fy, Fz).
    """
    t0 = 1.0 - tt
    t1 = tt
    o = 1.0 - uu - vv
    vx = t0 * sv_3d[0] + t1 * ev_3d[0]
    vy = t0 * sv_3d[1] + t1 * ev_3d[1]
    vz = t0 * sv_3d[2] + t1 * ev_3d[2]

    f1x = t0 * s1_3d[0] + t1 * e1_3d[0]
    f1y = t0 * s1_3d[1] + t1 * e1_3d[1]
    f1z = t0 * s1_3d[2] + t1 * e1_3d[2]
    f2x = t0 * s2_3d[0] + t1 * e2_3d[0]
    f2y = t0 * s2_3d[1] + t1 * e2_3d[1]
    f2z = t0 * s2_3d[2] + t1 * e2_3d[2]
    f3x = t0 * s3_3d[0] + t1 * e3_3d[0]
    f3y = t0 * s3_3d[1] + t1 * e3_3d[1]
    f3z = t0 * s3_3d[2] + t1 * e3_3d[2]

    fx = o * f1x + uu * f2x + vv * f3x
    fy = o * f1y + uu * f2y + vv * f3y
    fz = o * f1z + uu * f2z + vv * f3z
    return (vx - fx, vy - fy, vz - fz)

def ee_F_3d(s1, s2, s3, s4, e1, e2, e3, e4, tt, uu, vv):
    s1 = np.array(s1)
    s2 = np.array(s2)
    s3 = np.array(s3)
    s4 = np.array(s4)
    e1 = np.array(e1)
    e2 = np.array(e2)
    e3 = np.array(e3)
    e4 = np.array(e4)

    ea0 = (e1 - s1) * tt + s1;
    ea1 = (e2 - s2) * tt + s2;
    eb0 = (e3 - s3) * tt + s3;
    eb1 = (e4 - s4) * tt + s4;
    return((ea1 - ea0) * uu + ea0) - ((eb1 - eb0) * vv + eb0);

@myjit
def _sample_Fvf_component_3d(
    n_t, n_u, n_v,
    t_stride, u_stride,
    t_start, u_start, v_start,
    t_step, u_step, v_step,
    comp,  # 0:x, 1:y, 2:z
    sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
):
    """
    Sample a single component (x/y/z) of F_vf on a regular 3D grid over (t,u,v).
    Returns a flat list/array of length (n_t+1)*(n_u+1)*(n_v+1).
    """
    total = (n_t + 1) * (n_u + 1) * (n_v + 1)
    F = [0.0] * total
    for ti in range(n_t + 1):
        t_val = t_start + ti * t_step
        base_t = ti * t_stride
        for ui in range(n_u + 1):
            u_val = u_start + ui * u_step
            base_u = base_t + ui * u_stride
            for vi in range(n_v + 1):
                v_val = v_start + vi * v_step
                idx = base_u + vi
                Fx, Fy, Fz = vf_F_3d(
                    sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, t_val, u_val, v_val
                )
                F[idx] = Fx if comp == 0 else (Fy if comp == 1 else Fz)
    return F

@myjit
def _detect_zero_cells_3d(n_t, n_u, n_v, t_stride, u_stride, Fx, Fy, Fz, tol, contains_zero):
    """
    For each 3D cell (8 corners) check whether 0 is contained in the interval
    of Fx, Fy and Fz components: [min(corners), max(corners)] overlapping 0.
    contains_zero is a list sized n_t*n_u*n_v (one per cell) updated in-place.
    """
    for ti in range(n_t):
        for ui in range(n_u):
            for vi in range(n_v):
                cell_idx = ti * n_u * n_v + ui * n_v + vi
                i0 = ti * t_stride + ui * u_stride + vi
                i1 = i0 + 1
                i2 = i0 + u_stride
                i3 = i2 + 1
                i4 = i0 + t_stride
                i5 = i1 + t_stride
                i6 = i2 + t_stride
                i7 = i3 + t_stride

                xs = (Fx[i0], Fx[i1], Fx[i2], Fx[i3], Fx[i4], Fx[i5], Fx[i6], Fx[i7])
                ys = (Fy[i0], Fy[i1], Fy[i2], Fy[i3], Fy[i4], Fy[i5], Fy[i6], Fy[i7])
                zs = (Fz[i0], Fz[i1], Fz[i2], Fz[i3], Fz[i4], Fz[i5], Fz[i6], Fz[i7])

                has_zero_x = (min(xs) <= tol) and (max(xs) >= -tol)
                has_zero_y = (min(ys) <= tol) and (max(ys) >= -tol)
                has_zero_z = (min(zs) <= tol) and (max(zs) >= -tol)
                contains_zero[cell_idx] &= int(has_zero_x and has_zero_y and has_zero_z)

@myjit
def _project_uv_simplex(u: float, v: float):
    """
    Project (u,v) onto the 2-simplex S = {u>=0, v>=0, u+v<=1} in L2 sense.
    """
    # Clamp negatives
    u = max(u, 0.0)
    v = max(v, 0.0)
    s = u + v
    if s <= 1.0:
        return u, v
    # Project onto the line u+v=1 with u,v>=0
    u_proj = 0.5 * (u - v + 1.0)
    u_proj = min(max(u_proj, 0.0), 1.0)
    v_proj = 1.0 - u_proj
    return u_proj, v_proj

def _refine_vf_root_3d(
    sv_3d, s1_3d, s2_3d, s3_3d,
    ev_3d, e1_3d, e2_3d, e3_3d,
    t_init: float, u_init: float, v_init: float,
    max_iter: int = 25, tol_f: float = 1e-8
):
    """
    Gauss-Newton refinement with simple backtracking and projection to [0,1]xS.
    Returns (ok, t,u,v, fnorm).
    """
    import numpy as np

    def F_eval(t, u, v):
        Fx, Fy, Fz = vf_F_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, t, u, v)
        return np.array([Fx, Fy, Fz], dtype=float)

    x = np.array([t_init, u_init, v_init], dtype=float)
    # Ensure within domain
    x[0] = min(max(x[0], 0.0), 1.0)
    x[1], x[2] = _project_uv_simplex(x[1], x[2])

    for _ in range(max_iter):
        F0 = F_eval(x[0], x[1], x[2])
        fn = float(np.linalg.norm(F0))
        if fn <= tol_f:
            return True, float(x[0]), float(x[1]), float(x[2]), fn

        # Numerical Jacobian (forward diff; do not project for Jacobian eval)
        J = np.zeros((3, 3), dtype=float)
        eps_t = 1e-7
        eps_u = 1e-7
        eps_v = 1e-7
        F_t = F_eval(min(1.0, x[0] + eps_t), x[1], x[2])
        F_u = F_eval(x[0], x[1] + eps_u, x[2])
        F_v = F_eval(x[0], x[1], x[2] + eps_v)
        J[:, 0] = (F_t - F0) / eps_t
        J[:, 1] = (F_u - F0) / eps_u
        J[:, 2] = (F_v - F0) / eps_v

        # Solve J p = -F (least-squares for robustness)
        try:
            p, *_ = np.linalg.lstsq(J, -F0, rcond=None)
        except Exception:
            # Fallback to gradient descent step
            p = -J.T @ F0

        # Backtracking with projection
        alpha = 1.0
        improved = False
        for _ls in range(12):
            xn = x + alpha * p
            xn[0] = min(max(xn[0], 0.0), 1.0)
            xn[1], xn[2] = _project_uv_simplex(xn[1], xn[2])
            F1 = F_eval(xn[0], xn[1], xn[2])
            if float(np.linalg.norm(F1)) < fn:
                x = xn
                improved = True
                break
            alpha *= 0.5
        if not improved:
            # Small step - stop
            break

    F_final = np.linalg.norm(F_eval(x[0], x[1], x[2]))
    ok = F_final <= tol_f
    return ok, float(x[0]), float(x[1]), float(x[2]), float(F_final)


def find_root_dfs_3D(
    _max_iter: int,
    tol: float,
    sv_3d, s1_3d, s2_3d, s3_3d,
    ev_3d, e1_3d, e2_3d, e3_3d
) -> Tuple[bool, Optional[float], Optional[float], Optional[float]]:
    """
    Narrow-phase root search for 3D vertex-face CCD using interval tests on Fx,Fy,Fz.
    Returns (found, t, u, v). If not found, remaining values are None.
    """
    import math

    # Choose a cubic grid sized by ROOT_FINDING_CHUNK_SIZE
    side = max(4, int(round(ROOT_FINDING_CHUNK_SIZE ** (1.0 / 3.0))))
    while (side ** 3) > ROOT_FINDING_CHUNK_SIZE and side > 2:
        side -= 1

    t_n = side
    u_n = side
    v_n = side
    t_min = 0.0
    u_min = 0.0
    v_min = 0.0
    t_max = 1.0
    u_max = 1.0
    v_max = 1.0
    t_h = (t_max - t_min) / t_n
    u_h = (u_max - u_min) / u_n
    v_h = (v_max - v_min) / v_n

    t_stride = (u_n + 1) * (v_n + 1)
    u_stride = (v_n + 1)

    Fx = _sample_Fvf_component_3d(
        t_n, u_n, v_n, t_stride, u_stride, t_min, u_min, v_min, t_h, u_h, v_h,
        0, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
    )
    Fy = _sample_Fvf_component_3d(
        t_n, u_n, v_n, t_stride, u_stride, t_min, u_min, v_min, t_h, u_h, v_h,
        1, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
    )
    Fz = _sample_Fvf_component_3d(
        t_n, u_n, v_n, t_stride, u_stride, t_min, u_min, v_min, t_h, u_h, v_h,
        2, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
    )

    contains_zero = [1] * (t_n * u_n * v_n)
    _detect_zero_cells_3d(t_n, u_n, v_n, t_stride, u_stride, Fx, Fy, Fz, tol, contains_zero)

    def _cell_has_zero(t0, t1, u0, u1, v0, v1):
        if (max(0.0, u0) + max(0.0, v0)) > 1.0 + 1e-8:
            return False
        fx_min = fy_min = fz_min = 1e30
        fx_max = fy_max = fz_max = -1e30
        for tt in (t0, t1):
            for uu in (u0, u1):
                for vv in (v0, v1):
                    fx, fy, fz = vf_F_3d(
                        sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, tt, uu, vv
                    )
                    fx_min = fx if fx < fx_min else fx_min
                    fy_min = fy if fy < fy_min else fy_min
                    fz_min = fz if fz < fz_min else fz_min
                    fx_max = fx if fx > fx_max else fx_max
                    fy_max = fy if fy > fy_max else fy_max
                    fz_max = fz if fz > fz_max else fz_max
        return (fx_min <= tol and fx_max >= -tol) and (fy_min <= tol and fy_max >= -tol) and (fz_min <= tol and fz_max >= -tol)

    seeds = []
    for ti in range(t_n):
        for ui in range(u_n):
            u0 = u_min + ui * u_h
            for vi in range(v_n):
                cell_idx = ti * u_n * v_n + ui * v_n + vi
                if not contains_zero[cell_idx]:
                    continue
                v0 = v_min + vi * v_h
                if (u0 + v0) > 1.0 + 1e-8:
                    continue

                seeds.append((t_min + ti * t_h, t_min + (ti + 1) * t_h, u0, u0 + u_h, v0, v0 + v_h, 0))

    if not seeds:
        return (False, None, None, None)

    seeds.sort(key=lambda s: s[0])
    stack = seeds[::-1]  # pop prefers small t first
    min_dim = max(1e-6, tol * 10.0)

    max_depth = 100

    while stack:
        t0, t1, u0, u1, v0, v1, depth = stack.pop()
        if (max(0.0, u0) + max(0.0, v0)) > 1.0 + 1e-8:
            continue

        tc = 0.5 * (t0 + t1)
        uc = 0.5 * (u0 + u1)
        vc = 0.5 * (v0 + v1)

        fn_center = None
        if uc >= -1e-8 and vc >= -1e-8 and (uc + vc) <= 1.0 + 1e-8:
            fx_c, fy_c, fz_c = vf_F_3d(
                sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, tc, uc, vc
            )
            fn_center = math.sqrt(fx_c * fx_c + fy_c * fy_c + fz_c * fz_c)
            if fn_center <= tol:
                ok_ref, tr, ur, vr, fnr = _refine_vf_root_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, tc, uc, vc, max_iter=100, tol_f=max(1e-10, tol * 1e-2))
                if ok_ref:
                    return (True, tr, ur, vr)
                else:
                    continue
                # return (True, tc, uc, vc)

        size = max(t1 - t0, u1 - u0, v1 - v0)
        if size <= min_dim or depth >= _max_iter:
            continue

        tm = 0.5 * (t0 + t1)
        um = 0.5 * (u0 + u1)
        vm = 0.5 * (v0 + v1)

        subcells = (
            (t0, tm, u0, um, v0, vm),
            (t0, tm, u0, um, vm, v1),
            (t0, tm, um, u1, v0, vm),
            (t0, tm, um, u1, vm, v1),
            (tm, t1, u0, um, v0, vm),
            (tm, t1, u0, um, vm, v1),
            (tm, t1, um, u1, v0, vm),
            (tm, t1, um, u1, vm, v1),
        )
        next_depth = depth + 1
        if next_depth >= max_depth:
            break

        for t0s, t1s, u0s, u1s, v0s, v1s in subcells:
            if _cell_has_zero(t0s, t1s, u0s, u1s, v0s, v1s):
                stack.append((t0s, t1s, u0s, u1s, v0s, v1s, next_depth))

    return (False, None, None, None)


