#!/usr/bin/env python3
"""Demoted to spikes/ -- no caller anywhere in the repository.

Two-dimensional continuous collision detection: point against segment, with the
residual plots that go with it. SCCD is three-dimensional, and nothing in the
shipped tree calls any of this.

The implementations used to live in python/sccd_reference.py, where they were
582 of its 1344 lines and reachable only from here. They moved with the script
that uses them, so the shipped module holds only what the shipped tree calls.

Needs python/ on sys.path for `myjit`, and Matplotlib for the plots:

    PYTHONPATH=python python3 spikes/python/ccd2d.py
"""

from math import floor
from typing import Tuple, Optional

import numpy as np

from sccd_reference import ROOT_FINDING_CHUNK_SIZE, myjit


@myjit
def _vf_F_2d(sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d, tt, uu, vv):
    """
    2D version of the vertex-face residual:
      F(t,u,v) = V(t) - BaryFace(t,u,v)
    All inputs are 2D sequences (len==2). Returns tuple (Fx, Fy).
    """
    t0 = 1.0 - tt
    t1 = tt
    one_minus_uv = 1.0 - uu - vv

    vx = t0 * sv_2d[0] + t1 * ev_2d[0]
    vy = t0 * sv_2d[1] + t1 * ev_2d[1]

    fx = t0 * (one_minus_uv * s1_2d[0] + uu * s2_2d[0] + vv * s3_2d[0]) \
       + t1 * (one_minus_uv * e1_2d[0] + uu * e2_2d[0] + vv * e3_2d[0])
    fy = t0 * (one_minus_uv * s1_2d[1] + uu * s2_2d[1] + vv * s3_2d[1]) \
       + t1 * (one_minus_uv * e1_2d[1] + uu * e2_2d[1] + vv * e3_2d[1])

    return (vx - fx, vy - fy)


@myjit
def _sample_Fvf_component_2d(
    n_t, n_u, n_v,
    t_stride, u_stride,
    t_start, u_start, v_start,
    t_step, u_step, v_step,
    comp,  # 0 for x, 1 for y
    sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d
):
    """
    Sample a single component (x or y) of F_vf on a regular 3D grid over (t,u,v).
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
                Fx, Fy = _vf_F_2d(sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d, t_val, u_val, v_val)
                F[idx] = Fx if comp == 0 else Fy
    return F


@myjit
def _detect_zero_cells_2d(n_t, n_u, n_v, t_stride, u_stride, Fx, Fy, tol, contains_zero):
    """
    For each 3D cell (8 corners) check whether 0 is contained in the interval
    of both Fx and Fy components: [min(corners), max(corners)] overlapping 0.
    contains_zero is a list sized n_t*n_u*n_v (one per cell) that gets updated in-place.
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

                # Fx mins/maxs over 8 corners
                fxs = (Fx[i0], Fx[i1], Fx[i2], Fx[i3], Fx[i4], Fx[i5], Fx[i6], Fx[i7])
                fys = (Fy[i0], Fy[i1], Fy[i2], Fy[i3], Fy[i4], Fy[i5], Fy[i6], Fy[i7])

                fx_min, fx_max = min(fxs), max(fxs)
                fy_min, fy_max = min(fys), max(fys)

                # Check interval containment of zero with tolerance
                has_zero_fx = (fx_min <= tol) and (fx_max >= -tol)
                has_zero_fy = (fy_min <= tol) and (fy_max >= -tol)
                contains_zero[cell_idx] &= int(has_zero_fx and has_zero_fy)


@myjit
def find_root_2d(
    _max_iter: int,
    tol: float,
    sv_2d, s1_2d, s2_2d, s3_2d,
    ev_2d, e1_2d, e2_2d, e3_2d
) -> Tuple[bool, Optional[float], Optional[float], Optional[float]]:
    """
    Broad-phase root search for 2D vertex-face CCD using interval tests on Fx,Fy.
    Returns (found, t, u, v). If found=False, the remaining values are None.
    """
    # Coarse 2x2x2 sampling to estimate ranges
    Fx_coarse = _sample_Fvf_component_2d(
        1, 1, 1,
        t_stride=4, u_stride=2,
        t_start=0.0, u_start=0.0, v_start=0.0,
        t_step=1.0, u_step=1.0, v_step=1.0,
        comp=0, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
    )
    Fy_coarse = _sample_Fvf_component_2d(
        1, 1, 1,
        t_stride=4, u_stride=2,
        t_start=0.0, u_start=0.0, v_start=0.0,
        t_step=1.0, u_step=1.0, v_step=1.0,
        comp=1, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
    )

    fx_min, fx_max = min(Fx_coarse), max(Fx_coarse)
    fy_min, fy_max = min(Fy_coarse), max(Fy_coarse)

    # Quick AABB reject
    fx_intersect = (fx_min <= tol) and (fx_max >= -tol)
    fy_intersect = (fy_min <= tol) and (fy_max >= -tol)
    if not (fx_intersect and fy_intersect):
        return (False, None, None, None)

    # Allocate more resolution to dimensions that have smaller value ranges
    def safe_inv_range(vmin, vmax):
        rng = max(1e-5, (vmax - vmin))
        return max(2.0, 1.0 / rng)

    inv_x = safe_inv_range(fx_min, fx_max)
    inv_y = safe_inv_range(fy_min, fy_max)
    total_inv = inv_x * inv_y

    base = ROOT_FINDING_CHUNK_SIZE
    bias_t, bias_u, bias_v = 1.2, 1.0, 1.0
    weight_sum = bias_t + bias_u + bias_v
    Nt = max(2, int(floor((base ** (1/3)) * (bias_t / weight_sum) * total_inv ** (1/3))))
    Nu = max(2, int(floor((base ** (1/3)) * (bias_u / weight_sum) * total_inv ** (1/3))))
    Nv = max(2, int(floor((base ** (1/3)) * (bias_v / weight_sum) * total_inv ** (1/3))))
    while (Nt * Nu * Nv) > ROOT_FINDING_CHUNK_SIZE:
        if Nt >= Nu and Nt >= Nv and Nt > 2:
            Nt -= 1
        elif Nu >= Nt and Nu >= Nv and Nu > 2:
            Nu -= 1
        elif Nv > 2:
            Nv -= 1
        else:
            break

    # Convert to cell counts
    t_n = Nt - 1
    u_n = Nu - 1
    v_n = Nv - 1
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

    Fx = _sample_Fvf_component_2d(
        t_n, u_n, v_n,
        t_stride, u_stride,
        t_min, u_min, v_min,
        t_h, u_h, v_h,
        comp=0, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
    )
    Fy = _sample_Fvf_component_2d(
        t_n, u_n, v_n,
        t_stride, u_stride,
        t_min, u_min, v_min,
        t_h, u_h, v_h,
        comp=1, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
    )

    contains_zero = [1] * (t_n * u_n * v_n)
    _detect_zero_cells_2d(t_n, u_n, v_n, t_stride, u_stride, Fx, Fy, tol, contains_zero)

    found = False
    best_t = None
    best_u = None
    best_v = None
    for ti in range(t_n):
        if found:
            break
        for ui in range(u_n):
            for vi in range(v_n):
                cell_idx = ti * u_n * v_n + ui * v_n + vi
                if contains_zero[cell_idx]:
                    t_current = t_min + ti * t_h
                    u_current = u_min + ui * u_h
                    v_current = v_min + vi * v_h
                    if (u_h + v_h) > 1.0 + tol:
                        continue
                    best_t, best_u, best_v = t_current, u_current, v_current
                    found = True
                    break
    if not found:
        return (False, None, None, None)
    return (True, best_t, best_u, best_v)


def visualize_Fvf_xy_2d(
    sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d,
    t_value=0.5,
    resolution=256,
    out_path="vf2d_Fxy.png"
):
    """
    Save two heatmaps (Fx and Fy) of F_vf(t,u,v) evaluated at a fixed t over (u,v) in [0,1], u+v<=1.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sv_arr = np.asarray(sv_2d, dtype=float)
    s1_arr = np.asarray(s1_2d, dtype=float)
    s2_arr = np.asarray(s2_2d, dtype=float)
    s3_arr = np.asarray(s3_2d, dtype=float)
    ev_arr = np.asarray(ev_2d, dtype=float)
    e1_arr = np.asarray(e1_2d, dtype=float)
    e2_arr = np.asarray(e2_2d, dtype=float)
    e3_arr = np.asarray(e3_2d, dtype=float)

    U, V = np.meshgrid(
        np.linspace(0.0, 1.0, num=resolution),
        np.linspace(0.0, 1.0, num=resolution),
        indexing="xy"
    )
    mask = (U + V) <= 1.0
    O = 1.0 - U - V
    t0 = 1.0 - t_value
    t1 = t_value

    vx = t0 * sv_arr[0] + t1 * ev_arr[0]
    vy = t0 * sv_arr[1] + t1 * ev_arr[1]

    fx = t0 * (O * s1_arr[0] + U * s2_arr[0] + V * s3_arr[0]) + t1 * (O * e1_arr[0] + U * e2_arr[0] + V * e3_arr[0])
    fy = t0 * (O * s1_arr[1] + U * s2_arr[1] + V * s3_arr[1]) + t1 * (O * e1_arr[1] + U * e2_arr[1] + V * e3_arr[1])

    Fx = vx - fx
    Fy = vy - fy

    Fx_plot = np.full_like(Fx, np.nan)
    Fy_plot = np.full_like(Fy, np.nan)
    Fx_plot[mask] = Fx[mask]
    Fy_plot[mask] = Fy[mask]

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    im0 = axs[0].imshow(
        Fx_plot.T, origin="lower", extent=(0, 1, 0, 1), cmap="coolwarm"
    )
    axs[0].plot([0, 1], [1, 0], 'k--', lw=0.8)
    axs[0].set_title(f"F_vf.x at t={t_value:.3f}")
    axs[0].set_xlabel("u"); axs[0].set_ylabel("v")
    fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

    im1 = axs[1].imshow(
        Fy_plot.T, origin="lower", extent=(0, 1, 0, 1), cmap="coolwarm"
    )
    axs[1].plot([0, 1], [1, 0], 'k--', lw=0.8)
    axs[1].set_title(f"F_vf.y at t={t_value:.3f}")
    axs[1].set_xlabel("u"); axs[1].set_ylabel("v")
    fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

    fig.suptitle("2D Vertex-Face Residual Components")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


@myjit
def _vl_F_2d(sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d, tt, uu):
    """
    2D vertex-line residual:
      F(t,u) = V(t) - ((1-u) * E1(t) + u * E2(t))
    Inputs are 2D sequences (len==2). Returns tuple (Fx, Fy).
    """
    t0 = 1.0 - tt
    t1 = tt
    vx = t0 * sv_2d[0] + t1 * ev_2d[0]
    vy = t0 * sv_2d[1] + t1 * ev_2d[1]

    e1x = t0 * s1_2d[0] + t1 * e1_2d[0]
    e1y = t0 * s1_2d[1] + t1 * e1_2d[1]
    e2x = t0 * s2_2d[0] + t1 * e2_2d[0]
    e2y = t0 * s2_2d[1] + t1 * e2_2d[1]

    ex = (1.0 - uu) * e1x + uu * e2x
    ey = (1.0 - uu) * e1y + uu * e2y
    return (vx - ex, vy - ey)


@myjit
def _sample_Fvl_component_2d(
    n_t,
    n_u,
    t_stride,
    t_start,
    u_start,
    t_step,
    u_step,
    comp,  # 0 for x, 1 for y
    sv_2d,
    s1_2d,
    s2_2d,
    ev_2d,
    e1_2d,
    e2_2d,
):
    """
    Sample a single component (x or y) of F_vl on a regular 2D grid over (t,u).
    Returns a flat list/array of length (n_t+1)*(n_u+1).
    """
    total = (n_t + 1) * (n_u + 1)
    F = [0.0] * total
    for ti in range(n_t + 1):
        t_val = t_start + ti * t_step
        base_t = ti * t_stride
        for ui in range(n_u + 1):
            u_val = u_start + ui * u_step
            idx = base_t + ui
            Fx, Fy = _vl_F_2d(sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d, t_val, u_val)
            F[idx] = Fx if comp == 0 else Fy
    return F


def _detect_zero_cells_2d_line(n_t, n_u, t_stride, Fx, Fy, tol, contains_zero):
    """
    For each 2D cell (4 corners) check whether 0 is contained in the interval
    of both Fx and Fy components: [min(corners), max(corners)] overlapping 0.
    contains_zero is a list sized n_t*n_u (one per cell) that gets updated in-place.
    """
    for ti in range(n_t):
        for ui in range(n_u):
            cell_idx = ti * n_u + ui
            i00 = ti * t_stride + ui
            i01 = i00 + 1
            i10 = (ti + 1) * t_stride + ui
            i11 = i10 + 1

            fxs = (Fx[i00], Fx[i01], Fx[i10], Fx[i11])
            fys = (Fy[i00], Fy[i01], Fy[i10], Fy[i11])
            fx_min, fx_max = min(fxs), max(fxs)
            fy_min, fy_max = min(fys), max(fys)

            has_zero_fx = (fx_min <= tol) and (fx_max >= -tol)
            has_zero_fy = (fy_min <= tol) and (fy_max >= -tol)
            contains_zero[cell_idx] &= int(has_zero_fx and has_zero_fy)


@myjit
def find_root_line_2d(
    _max_iter: int,
    tol: float,
    sv_2d,
    s1_2d,
    s2_2d,
    ev_2d,
    e1_2d,
    e2_2d,
):
    """
    Broad-phase root search for 2D vertex-line CCD using interval tests on Fx,Fy.
    Returns (found, t, u). If found=False, the remaining values are None.
    """
    # Coarse 2x2 sampling to estimate ranges
    Fx_coarse = _sample_Fvl_component_2d(
        1, 1, 2, 0.0, 0.0, 1.0, 1.0, 0, sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d
    )
    Fy_coarse = _sample_Fvl_component_2d(
        1, 1, 2, 0.0, 0.0, 1.0, 1.0, 1, sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d
    )

    fx_min, fx_max = min(Fx_coarse), max(Fx_coarse)
    fy_min, fy_max = min(Fy_coarse), max(Fy_coarse)

    fx_intersect = (fx_min <= tol) and (fx_max >= -tol)
    fy_intersect = (fy_min <= tol) and (fy_max >= -tol)
    if not (fx_intersect and fy_intersect):
        return (False, None, None)

    def safe_inv_range(vmin, vmax):
        rng = max(1e-5, (vmax - vmin))
        return max(2.0, 1.0 / rng)

    inv_x = safe_inv_range(fx_min, fx_max)
    inv_y = safe_inv_range(fy_min, fy_max)
    total_inv = inv_x * inv_y

    base = ROOT_FINDING_CHUNK_SIZE
    bias_t, bias_u = 1.2, 1.0
    weight_sum = bias_t + bias_u
    Nt = max(2, int((base ** 0.5) * (bias_t / weight_sum) * (total_inv ** 0.5)))
    Nu = max(2, int((base ** 0.5) * (bias_u / weight_sum) * (total_inv ** 0.5)))
    while (Nt * Nu) > ROOT_FINDING_CHUNK_SIZE:
        if Nt >= Nu and Nt > 2:
            Nt -= 1
        elif Nu > 2:
            Nu -= 1
        else:
            break

    t_n = Nt - 1
    u_n = Nu - 1
    t_min = 0.0
    u_min = 0.0
    t_max = 1.0
    u_max = 1.0
    t_h = (t_max - t_min) / t_n
    u_h = (u_max - u_min) / u_n

    t_stride = (u_n + 1)

    Fx = _sample_Fvl_component_2d(
        t_n, u_n, t_stride, t_min, u_min, t_h, u_h, 0, sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d
    )
    Fy = _sample_Fvl_component_2d(
        t_n, u_n, t_stride, t_min, u_min, t_h, u_h, 1, sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d
    )

    contains_zero = [1] * (t_n * u_n)
    _detect_zero_cells_2d_line(t_n, u_n, t_stride, Fx, Fy, tol, contains_zero)

    found = False
    best_t = None
    best_u = None
    for ti in range(t_n):
        if found:
            break
        for ui in range(u_n):
            cell_idx = ti * u_n + ui
            if contains_zero[cell_idx]:
                t_current = t_min + ti * t_h
                u_current = u_min + ui * u_h
                best_t, best_u = t_current, u_current
                found = True
                break
    if not found:
        return (False, None, None)
    return (True, best_t, best_u)


def visualize_Fvl_tu_2d(
    sv_2d,
    s1_2d,
    s2_2d,
    ev_2d,
    e1_2d,
    e2_2d,
    resolution=256,
    out_path="vl2d_Ftu.png",
):
    """
    Save two heatmaps (Fx and Fy) of F_vl(t,u) over t in [0,1], u in [0,1].
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sv = np.asarray(sv_2d, dtype=float)
    s1 = np.asarray(s1_2d, dtype=float)
    s2 = np.asarray(s2_2d, dtype=float)
    ev = np.asarray(ev_2d, dtype=float)
    e1 = np.asarray(e1_2d, dtype=float)
    e2 = np.asarray(e2_2d, dtype=float)

    T, U = np.meshgrid(
        np.linspace(0.0, 1.0, num=resolution),
        np.linspace(0.0, 1.0, num=resolution),
        indexing="xy",
    )

    t0 = 1.0 - T
    t1 = T

    vx = t0 * sv[0] + t1 * ev[0]
    vy = t0 * sv[1] + t1 * ev[1]

    e1x = t0 * s1[0] + t1 * e1[0]
    e1y = t0 * s1[1] + t1 * e1[1]
    e2x = t0 * s2[0] + t1 * e2[0]
    e2y = t0 * s2[1] + t1 * e2[1]

    ex = (1.0 - U) * e1x + U * e2x
    ey = (1.0 - U) * e1y + U * e2y

    Fx = vx - ex
    Fy = vy - ey

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    im0 = axs[0].imshow(Fx.T, origin="lower", extent=(0, 1, 0, 1), cmap="coolwarm")
    axs[0].set_title("F_vl.x over (t,u)")
    axs[0].set_xlabel("t"); axs[0].set_ylabel("u")
    fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

    im1 = axs[1].imshow(Fy.T, origin="lower", extent=(0, 1, 0, 1), cmap="coolwarm")
    axs[1].set_title("F_vl.y over (t,u)")
    axs[1].set_xlabel("t"); axs[1].set_ylabel("u")
    fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

    fig.suptitle("2D Vertex-Line Residual Components")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def visualize_Fvl_norm_tu_2d(
    sv_2d,
    s1_2d,
    s2_2d,
    ev_2d,
    e1_2d,
    e2_2d,
    resolution=256,
    out_path="vl2d_Fnorm.png",
    log_scale=False,
    zero_contour=1e-6,
):
    """
    Save a heatmap of ||F_vl(t,u)|| over t in [0,1], u in [0,1].
    Optionally plot on log scale and overlay a near-zero contour.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sv = np.asarray(sv_2d, dtype=float)
    s1 = np.asarray(s1_2d, dtype=float)
    s2 = np.asarray(s2_2d, dtype=float)
    ev = np.asarray(ev_2d, dtype=float)
    e1 = np.asarray(e1_2d, dtype=float)
    e2 = np.asarray(e2_2d, dtype=float)

    T, U = np.meshgrid(
        np.linspace(0.0, 1.0, num=resolution),
        np.linspace(0.0, 1.0, num=resolution),
        indexing="xy",
    )

    t0 = 1.0 - T
    t1 = T

    vx = t0 * sv[0] + t1 * ev[0]
    vy = t0 * sv[1] + t1 * ev[1]

    e1x = t0 * s1[0] + t1 * e1[0]
    e1y = t0 * s1[1] + t1 * e1[1]
    e2x = t0 * s2[0] + t1 * e2[0]
    e2y = t0 * s2[1] + t1 * e2[1]

    ex = (1.0 - U) * e1x + U * e2x
    ey = (1.0 - U) * e1y + U * e2y

    Fx = vx - ex
    Fy = vy - ey
    Fn = np.sqrt(Fx * Fx + Fy * Fy)

    data = np.log10(Fn + 1e-16) if log_scale else Fn

    fig, ax = plt.subplots(1, 1, figsize=(6, 5), constrained_layout=True)
    im = ax.imshow(data.T, origin="lower", extent=(0, 1, 0, 1), cmap="magma")
    ax.set_title("||F_vl|| over (t,u)" + (" (log10)" if log_scale else ""))
    ax.set_xlabel("t"); ax.set_ylabel("u")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    if zero_contour is not None and zero_contour > 0:
        cs = ax.contour(
            T, U, Fn, levels=[zero_contour], colors="cyan", linewidths=1.0
        )
        # Add a proxy artist for legend instead of accessing cs.collections (API changes across versions)
        try:
            has_segments = len(cs.allsegs[0]) > 0
        except Exception:
            has_segments = False
        if has_segments:
            from matplotlib.lines import Line2D
            proxy = Line2D([], [], color="cyan", lw=1.0, label=f"||F||={zero_contour:g}")
            ax.legend(handles=[proxy], loc="upper right")

    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def find_root_line(
    max_iter: int,
    tol: float,
    sv_2d, s1_2d, s2_2d,  # 2D point-like (len==2)
    ev_2d, e1_2d, e2_2d,  # 2D point-like (len==2)
) -> Tuple[bool, Optional[float], Optional[float]]:
    """
    2D CCD root-finding for point vs line segment (barycentric u on segment; t in [0,1]).
    Returns (found, t, u).
    """
    return find_root_line_2d(max_iter, tol, sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d)


def visualize_F_tu_line(
    sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d,
    resolution: int = 256,
    out_path: str = "vl2d_Ftu.png",
) -> None:
    """
    Visualize the components (Fx,Fy) of the 2D vertex-line residual over (t,u).
    Saves a side-by-side heatmap image to out_path.
    """
    visualize_Fvl_tu_2d(
        sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d,
        resolution=resolution,
        out_path=out_path,
    )


def visualize_scene_line(
    sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d,
    out_path: str = "vl2d_scene.png",
    t_impact: Optional[float] = None,
    u_impact: Optional[float] = None,
) -> None:
    """
    Draw start/end point and line segment positions with motion arrows in one timestep.
    Optionally overlay the impact point (at t_impact, u_impact).
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sv = np.asarray(sv_2d, dtype=float)
    s1 = np.asarray(s1_2d, dtype=float)
    s2 = np.asarray(s2_2d, dtype=float)
    ev = np.asarray(ev_2d, dtype=float)
    e1 = np.asarray(e1_2d, dtype=float)
    e2 = np.asarray(e2_2d, dtype=float)

    seg0 = np.vstack([s1, s2])
    seg1 = np.vstack([e1, e2])

    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)

    # Start segment
    ax.plot(seg0[:, 0], seg0[:, 1], color="#555555", lw=2, label="segment t=0")
    # End segment
    ax.plot(seg1[:, 0], seg1[:, 1], color="#1f77b4", lw=2, label="segment t=1")

    # Start and end point
    ax.scatter([sv[0]], [sv[1]], color="#d62728", s=40, zorder=3, label="point t=0")
    ax.scatter([ev[0]], [ev[1]], color="#d62728", marker="x", s=60, zorder=3, label="point t=1")
    ax.arrow(sv[0], sv[1], ev[0]-sv[0], ev[1]-sv[1],
             length_includes_head=True, head_width=0.02, head_length=0.03,
             color="#d62728", alpha=0.9, lw=1.5)

    # Vertex motion arrows
    for a, b, c in [(s1, e1, "#2ca02c"), (s2, e2, "#2ca02c")]:
        ax.arrow(a[0], a[1], b[0]-a[0], b[1]-a[1],
                 length_includes_head=True, head_width=0.02, head_length=0.03,
                 color=c, alpha=0.8, lw=1.2)

    # Bounds
    pts = np.vstack([seg0, seg1, sv[None, :], ev[None, :]])
    xmin, ymin = pts.min(axis=0) - 0.1
    xmax, ymax = pts.max(axis=0) + 0.1
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, lw=0.3, alpha=0.4)
    ax.legend(loc="upper right")
    ax.set_title("2D CCD (point-line): one timestep (start/end configs and motions)")

    # Optional impact overlay
    if t_impact is not None and u_impact is not None:
        t0 = 1.0 - float(t_impact)
        t1 = float(t_impact)
        v_imp = t0 * sv + t1 * ev
        e1_imp = t0 * s1 + t1 * e1
        e2_imp = t0 * s2 + t1 * e2
        u = float(u_impact)
        e_imp = (1.0 - u) * e1_imp + u * e2_imp
        ax.scatter([v_imp[0]], [v_imp[1]], color="#ff9900", s=60, zorder=4, label="impact V(t)")
        # Draw the segment at the impact time instead of a single point on the segment
        seg_imp = np.vstack([e1_imp, e2_imp])
        ax.plot(seg_imp[:, 0], seg_imp[:, 1], color="#ff9900", lw=2.0, alpha=0.9, label="segment at impact t")
        ax.plot([v_imp[0], e_imp[0]], [v_imp[1], e_imp[1]], linestyle=":", color="#ff9900", alpha=0.8)
        ax.legend(loc="upper right")

    fig.savefig(out_path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    # Minimal demo: moving point against a moving line segment
    sv_demo = (0.25, -0.2)
    s1_demo = (0.0, 0.0)
    s2_demo = (1.0, 0.0)

    ev_demo = (0.6, 0.4)
    # Move segment slightly to a nearby location at t=1
    tx, ty = 0.06, -0.04
    e1_demo = (s1_demo[0] + tx, s1_demo[1] + ty)
    e2_demo = (s2_demo[0] + tx, s2_demo[1] + ty)

    ok, t, u = find_root_line(20, 1e-10, sv_demo, s1_demo, s2_demo, ev_demo, e1_demo, e2_demo)
    print(f"found={ok}, t={t}, u={u}")

    visualize_F_tu_line(sv_demo, s1_demo, s2_demo, ev_demo, e1_demo, e2_demo, resolution=256, out_path="vl2d_Ftu.png")
    visualize_scene_line(
        sv_demo, s1_demo, s2_demo, ev_demo, e1_demo, e2_demo,
        out_path="vl2d_scene.png",
        t_impact=t if ok else None,
        u_impact=u if ok else None,
    )
    visualize_Fvl_norm_tu_2d(
        sv_demo, s1_demo, s2_demo, ev_demo, e1_demo, e2_demo,
        resolution=256, out_path="vl2d_Fnorm.png", log_scale=False, zero_contour=1e-2
    )