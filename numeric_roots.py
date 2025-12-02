#!/usr/bin/env python3

import numpy as np
from math import floor
from typing import Tuple, Optional

try:
    from numba import njit
    def myjit(f): return njit(f, fastmath=True, boundscheck=False, nogil=True)
except:
    print("Could not find numba")
    def myjit(f): return f

# Keep chunk size similar to C++ heuristic
ROOT_FINDING_CHUNK_SIZE = 4*4096

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

############################################################
# 2D CCD (point vs line segment) utilities and visualization
############################################################

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


############################################################
# 3D CCD (point vs triangle) utilities and visualization
############################################################
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

@myjit
def find_root_vf_3d(
    _max_iter: int,
    tol: float,
    sv_3d, s1_3d, s2_3d, s3_3d,
    ev_3d, e1_3d, e2_3d, e3_3d
) -> Tuple[bool, Optional[float], Optional[float], Optional[float]]:
    """
    Broad-phase root search for 3D vertex-face CCD using interval tests on Fx,Fy,Fz.
    Returns (found, t, u, v). If not found, remaining values are None.
    """
    # Coarse 2x2x2 sampling
    Fx_coarse = _sample_Fvf_component_3d(
        1, 1, 1, 4, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        0, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
    )
    Fy_coarse = _sample_Fvf_component_3d(
        1, 1, 1, 4, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        1, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
    )
    Fz_coarse = _sample_Fvf_component_3d(
        1, 1, 1, 4, 2, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        2, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
    )

    def aabb(arr):
        return (min(arr), max(arr))

    ax0, ax1 = aabb(Fx_coarse)
    ay0, ay1 = aabb(Fy_coarse)
    az0, az1 = aabb(Fz_coarse)

    if not ((ax0 <= tol and ax1 >= -tol) and (ay0 <= tol and ay1 >= -tol) and (az0 <= tol and az1 >= -tol)):
        return (False, None, None, None)

    def safe_inv_range(vmin, vmax):
        rng = max(1e-5, (vmax - vmin))
        return max(2.0, 1.0 / rng)

    invx = safe_inv_range(ax0, ax1)
    invy = safe_inv_range(ay0, ay1)
    invz = safe_inv_range(az0, az1)

    # Allocate grid resolution bounded by ROOT_FINDING_CHUNK_SIZE
    base = ROOT_FINDING_CHUNK_SIZE
    invs = [invx, invy, invz]
    order = sorted(range(3), key=lambda i: invs[i])  # increasing
    N = [0, 0, 0]
    total = 1
    tot_inv = invx * invy * invz
    for idx in order:
        N[idx] = max(8 if idx == 0 else 2, int(floor((base / float(total)) * (invs[idx] / max(1e-12, tot_inv)) ** (1.0 / (3 - order.index(idx))))))
        total *= N[idx]
        tot_inv = max(1e-12, tot_inv / invs[idx])
    while total > ROOT_FINDING_CHUNK_SIZE:
        # reduce the largest dimension
        i = N.index(max(N))
        if N[i] > 2:
            total //= N[i]
            N[i] -= 1
            total *= N[i]
        else:
            break

    Nt, Nu, Nv = N
    t_n, u_n, v_n = Nt - 1, Nu - 1, Nv - 1
    t_min = 0.0; u_min = 0.0; v_min = 0.0
    t_max = 1.0; u_max = 1.0; v_max = 1.0
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

    # print(f'Nt = {Nt}, Nu = {Nu}, Nv = {Nv}')
    # print(f'potential roots {np.sum(contains_zero)}')

    # Build seed list from all cells:
    # - Priority 0: cells that pass interval zero test
    # - Priority 1: remaining cells, sorted by center ||F|| value
    seeds = []
    import math as _math_np_seed
    # Precompute center values and priorities
    for ti in range(t_n):
        for ui in range(u_n):
            for vi in range(v_n):
                cell_idx = ti * u_n * v_n + ui * v_n + vi
                if not contains_zero[cell_idx]:
                    continue
                t_c = t_min + (ti + 0.5) * t_h
                u_c = u_min + (ui + 0.5) * u_h
                v_c = v_min + (vi + 0.5) * v_h

                # print(f'{t_c}, {u_c}, {v_c}')
                
                if (u_c + v_c) > 1.0 + 1e-8:
                    u_c, v_c = _project_uv_simplex(u_c, v_c)
                
                Fx_c, Fy_c, Fz_c = vf_F_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, t_c, u_c, v_c)
                fn_c = _math_np_seed.sqrt(Fx_c * Fx_c + Fy_c * Fy_c + Fz_c * Fz_c)

                seeds.append((t_c, fn_c, ti, ui, vi, t_c, u_c, v_c))

    # Sort seeds: priority first, then by smallest ||F|| at center
    seeds.sort(key=lambda x: (x[0], x[1]))

    # Refine a limited number of seeds to keep compute bounded
    max_seeds = min(1024, len(seeds))
    found_any = False
    best_t = None; best_u = None; best_v = None; best_fn = None
    for k in range(max_seeds):
        _prio, _fnc, ti, ui, vi, t_current, u_current, v_current = seeds[k]
        ok_ref, tr, ur, vr, fnr = _refine_vf_root_3d(
            sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
            t_current, u_current, v_current, max_iter=100, tol_f=max(1e-10, tol * 1e-2)
        )
        if ok_ref and 0.0 <= tr <= 1.0 and ur >= -1e-8 and vr >= -1e-8 and (ur + vr) <= 1.0 + 1e-8:
            if (not found_any) or (tr < best_t):
                found_any = True
                best_t, best_u, best_v, best_fn = tr, ur, vr, fnr
    if not found_any:
        # Fallback: 1D sweep over t with per-t projection of (u,v)
        ok_sw, ts, us, vs = find_root_vf_3d_sweep(
            tol, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
        )
        if ok_sw:
            return (True, ts, us, vs)
        return (False, None, None, None)
    return (True, best_t, best_u, best_v)


def _closest_point_on_triangle_barycentric(p, a, b, c):
    """
    Return barycentric coords (w0,w1,w2) of closest point on triangle ABC to P
    and the closest point itself. Based on Ericson's Real-Time Collision Detection.
    """
    import numpy as np
    p = np.asarray(p, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    ab = b - a
    ac = c - a
    ap = p - a

    d1 = np.dot(ab, ap)
    d2 = np.dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
        return 1.0, 0.0, 0.0, a

    bp = p - b
    d3 = np.dot(ab, bp)
    d4 = np.dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
        return 0.0, 1.0, 0.0, b

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        v = d1 / (d1 - d3)
        closest = a + v * ab
        return 1.0 - v, v, 0.0, closest

    cp = p - c
    d5 = np.dot(ab, cp)
    d6 = np.dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
        return 0.0, 0.0, 1.0, c

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        w = d2 / (d2 - d6)
        closest = a + w * ac
        return 1.0 - w, 0.0, w, closest

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        closest = b + w * (c - b)
        return 0.0, 1.0 - w, w, closest

    denom = 1.0 / (va + vb + vc)
    v = vb * denom
    w = vc * denom
    u = 1.0 - v - w
    closest = u * a + v * b + w * c
    return u, v, w, closest


def find_root_vf_3d_sweep(
    tol: float,
    sv_3d, s1_3d, s2_3d, s3_3d,
    ev_3d, e1_3d, e2_3d, e3_3d,
    samples: int = 201
) -> Tuple[bool, Optional[float], Optional[float], Optional[float]]:
    """
    Fallback 1D search over t. For each t, compute triangle at t and closest
    (u,v) on triangle to point V(t). If min ||F|| across t is <= tol_g, accept.
    """
    import numpy as np
    tol_g = max(1e-6, tol * 1e4)
    best = (None, None, None, None)  # fn, t, u, v
    for i in range(samples):
        t = i / float(samples - 1)
        t0 = 1.0 - t; t1 = t
        v = np.array([
            t0 * sv_3d[0] + t1 * ev_3d[0],
            t0 * sv_3d[1] + t1 * ev_3d[1],
            t0 * sv_3d[2] + t1 * ev_3d[2],
        ], dtype=float)
        f1 = np.array([
            t0 * s1_3d[0] + t1 * e1_3d[0],
            t0 * s1_3d[1] + t1 * e1_3d[1],
            t0 * s1_3d[2] + t1 * e1_3d[2],
        ], dtype=float)
        f2 = np.array([
            t0 * s2_3d[0] + t1 * e2_3d[0],
            t0 * s2_3d[1] + t1 * e2_3d[1],
            t0 * s2_3d[2] + t1 * e2_3d[2],
        ], dtype=float)
        f3 = np.array([
            t0 * s3_3d[0] + t1 * e3_3d[0],
            t0 * s3_3d[1] + t1 * e3_3d[1],
            t0 * s3_3d[2] + t1 * e3_3d[2],
        ], dtype=float)

        w0, w1, w2, closest = _closest_point_on_triangle_barycentric(v, f1, f2, f3)
        uv = (w1, w2)
        res = v - closest
        fn = float(np.linalg.norm(res))
        if best[0] is None or fn < best[0]:
            best = (fn, t, uv[0], uv[1])

    fn, t_b, u_b, v_b = best
    if fn is None or fn > tol_g:
        return (False, None, None, None)

    # Refine with Gauss-Newton from this seed
    ok_ref, tr, ur, vr, fnr = _refine_vf_root_3d(
        sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
        t_b, u_b, v_b, max_iter=60, tol_f=max(1e-10, tol * 1e-2)
    )
    if ok_ref and 0.0 <= tr <= 1.0 and ur >= -1e-8 and vr >= -1e-8 and (ur + vr) <= 1.0 + 1e-8:
        return (True, tr, ur, vr)
    return (False, None, None, None)


def visualize_Fvf_uv_norm_3d(
    sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
    t_value=0.5,
    resolution=256,
    out_path="vf3d_uv_norm.png",
    log_scale=False
):
    """
    Save a heatmap of ||F_vf(t_value,u,v)|| over (u,v) in [0,1], u+v<=1.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sv = np.asarray(sv_3d, dtype=float)
    s1 = np.asarray(s1_3d, dtype=float)
    s2 = np.asarray(s2_3d, dtype=float)
    s3 = np.asarray(s3_3d, dtype=float)
    ev = np.asarray(ev_3d, dtype=float)
    e1 = np.asarray(e1_3d, dtype=float)
    e2 = np.asarray(e2_3d, dtype=float)
    e3 = np.asarray(e3_3d, dtype=float)

    U, V = np.meshgrid(
        np.linspace(0.0, 1.0, num=resolution),
        np.linspace(0.0, 1.0, num=resolution),
        indexing="xy"
    )
    mask = (U + V) <= 1.0
    O = 1.0 - U - V
    t0 = 1.0 - t_value
    t1 = t_value

    vx = t0 * sv[0] + t1 * ev[0]
    vy = t0 * sv[1] + t1 * ev[1]
    vz = t0 * sv[2] + t1 * ev[2]

    f1x = t0 * s1[0] + t1 * e1[0]; f1y = t0 * s1[1] + t1 * e1[1]; f1z = t0 * s1[2] + t1 * e1[2]
    f2x = t0 * s2[0] + t1 * e2[0]; f2y = t0 * s2[1] + t1 * e2[1]; f2z = t0 * s2[2] + t1 * e2[2]
    f3x = t0 * s3[0] + t1 * e3[0]; f3y = t0 * s3[1] + t1 * e3[1]; f3z = t0 * s3[2] + t1 * e3[2]

    fx = O * f1x + U * f2x + V * f3x
    fy = O * f1y + U * f2y + V * f3y
    fz = O * f1z + U * f2z + V * f3z

    Fx = vx - fx
    Fy = vy - fy
    Fz = vz - fz
    Fn = np.sqrt(Fx * Fx + Fy * Fy + Fz * Fz)

    data = np.log10(Fn + 1e-16) if log_scale else Fn

    Fn_plot = np.full_like(data, np.nan)
    Fn_plot[mask] = data[mask]

    fig, ax = plt.subplots(1, 1, figsize=(6, 5), constrained_layout=True)
    im = ax.imshow(Fn_plot.T, origin="lower", extent=(0, 1, 0, 1), cmap="magma")
    ax.plot([0, 1], [1, 0], 'k--', lw=0.8)
    ax.set_title("||F_vf|| over (u,v) at fixed t" + (" (log10)" if log_scale else ""))
    ax.set_xlabel("u"); ax.set_ylabel("v")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def visualize_scene_vf_3d(
    sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
    out_path="vf3d_scene.png",
    t_impact: Optional[float] = None,
    u_impact: Optional[float] = None,
    v_impact: Optional[float] = None
):
    """
    3D scene: draw start/end triangle meshes and point, with motion arrows, and optional impact overlay.
    """
    import numpy as np
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sv = np.asarray(sv_3d, dtype=float)
    s1 = np.asarray(s1_3d, dtype=float)
    s2 = np.asarray(s2_3d, dtype=float)
    s3 = np.asarray(s3_3d, dtype=float)
    ev = np.asarray(ev_3d, dtype=float)
    e1 = np.asarray(e1_3d, dtype=float)
    e2 = np.asarray(e2_3d, dtype=float)
    e3 = np.asarray(e3_3d, dtype=float)

    tri0 = np.vstack([s1, s2, s3])
    tri1 = np.vstack([e1, e2, e3])

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Start and end triangles
    poly0 = Poly3DCollection([tri0], alpha=0.25, facecolor="#bbbbbb", edgecolor="#555555")
    poly1 = Poly3DCollection([tri1], alpha=0.15, facecolor="#1f77b4", edgecolor="#1f77b4")
    ax.add_collection3d(poly0)
    ax.add_collection3d(poly1)

    # Start and end point
    ax.scatter([sv[0]], [sv[1]], [sv[2]], color="#d62728", s=30, label="point t=0")
    ax.scatter([ev[0]], [ev[1]], [ev[2]], color="#d62728", marker="x", s=40, label="point t=1")
    ax.plot([sv[0], ev[0]], [sv[1], ev[1]], [sv[2], ev[2]], color="#d62728", lw=1.5, alpha=0.9)

    # Vertex motion arrows (as lines)
    for a, b in [(s1, e1), (s2, e2), (s3, e3)]:
        ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], color="#2ca02c", lw=1.2, alpha=0.8)

    # Impact overlay
    if t_impact is not None and u_impact is not None and v_impact is not None:
        t0 = 1.0 - float(t_impact)
        t1 = float(t_impact)
        v_imp = t0 * sv + t1 * ev
        f1_imp = t0 * s1 + t1 * e1
        f2_imp = t0 * s2 + t1 * e2
        f3_imp = t0 * s3 + t1 * e3

        # Triangle at impact
        tri_imp = np.vstack([f1_imp, f2_imp, f3_imp])
        poly_imp = Poly3DCollection([tri_imp], alpha=0.35, facecolor="#ff9900", edgecolor="#ff9900")
        ax.add_collection3d(poly_imp)

        u = float(u_impact); v = float(v_impact)
        e_imp = (1.0 - u - v) * f1_imp + u * f2_imp + v * f3_imp
        ax.scatter([v_imp[0]], [v_imp[1]], [v_imp[2]], color="#ff9900", s=40, label="impact V(t)")
        ax.plot([v_imp[0], e_imp[0]], [v_imp[1], e_imp[1]], [v_imp[2], e_imp[2]], linestyle=":", color="#ff9900", alpha=0.8)

    # Bounds/equal-ish aspect
    pts = np.vstack([tri0, tri1, sv[None, :], ev[None, :]])
    xyz_min = pts.min(axis=0) - 0.1
    xyz_max = pts.max(axis=0) + 0.1
    ax.set_xlim(xyz_min[0], xyz_max[0])
    ax.set_ylim(xyz_min[1], xyz_max[1])
    ax.set_zlim(xyz_min[2], xyz_max[2])

    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title("3D CCD (point-triangle): start/end configs and impact")
    ax.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
