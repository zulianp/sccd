#!/usr/bin/env python3
"""Demoted to spikes/ -- no caller anywhere in the repository.

A grid-and-refine vertex-face search, and the two visualisations that go with it.
Its only importer, python/sccd_crosscheck.py, bound `find_root_vf` from here and
then never called it: that path goes through the C ABI binding instead.

The implementations used to live in python/sccd_reference.py, where they were 376
of its 1344 lines and reachable only from here. They moved with the script that
uses them; the residual and the shared helpers stayed, because the cross-check
still calls those.

Needs python/ on sys.path, and Matplotlib for the plots:

    PYTHONPATH=python python3 spikes/python/ccd3d.py
"""

from math import floor
from typing import Tuple, Optional

import numpy as np

from sccd_reference import (
    ROOT_FINDING_CHUNK_SIZE,
    _detect_zero_cells_3d,
    _project_uv_simplex,
    _refine_vf_root_3d,
    _sample_Fvf_component_3d,
    myjit,
    vf_F_3d,
)


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


def find_root_vf(
    max_iter: int,
    tol: float,
    sv_3d, s1_3d, s2_3d, s3_3d,
    ev_3d, e1_3d, e2_3d, e3_3d,
) -> Tuple[bool, Optional[float], Optional[float], Optional[float]]:
    """
    3D CCD root-finding for point vs triangle (barycentric u,v on triangle; t in [0,1]).
    Returns (found, t, u, v).
    """
    return find_root_vf_3d(
        max_iter, tol, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d
    )


def visualize_uv_norm(
    sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
    t_value: float = 0.5,
    resolution: int = 256,
    out_path: str = "vf3d_uv_norm.png",
    log_scale: bool = False,
) -> None:
    """
    Visualize ||F|| over (u,v) at a fixed t slice.
    """
    visualize_Fvf_uv_norm_3d(
        sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
        t_value=t_value, resolution=resolution, out_path=out_path, log_scale=log_scale
    )


def visualize_scene(
    sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
    out_path: str = "vf3d_scene.png",
    t_impact: Optional[float] = None,
    u_impact: Optional[float] = None,
    v_impact: Optional[float] = None,
) -> None:
    """
    Visualize 3D scene (start/end triangles and point) with optional impact overlay.
    """
    visualize_scene_vf_3d(
        sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d,
        out_path=out_path, t_impact=t_impact, u_impact=u_impact, v_impact=v_impact
    )


if __name__ == "__main__":
    # Sample demo: point and triangle moving over one timestep
    sv_demo = (0.2, 0.2, 0.2)
    ev_demo = (0.5, 0.6, 0.2)

    s1_demo = (0.0, 0.0, 0.0)
    s2_demo = (1.0, 0.0, 0.1)
    s3_demo = (0.0, 1.0, -0.1)

    
    tx, ty, tz = 0.4, 0.4, 0.4
    e1_demo = (s1_demo[0] + tx, s1_demo[1] + ty, s1_demo[2] + tz)
    e2_demo = (s2_demo[0] + tx, s2_demo[1] + ty, s2_demo[2] + tz)
    e3_demo = (s3_demo[0] + tx, s3_demo[1] + ty, s3_demo[2] + tz)

    ok, t, u, v = find_root_vf(1000, 1e-6, sv_demo, s1_demo, s2_demo, s3_demo, ev_demo, e1_demo, e2_demo, e3_demo)
    print(f"found={ok}, t={t}, u={u}, v={v}")

    visualize_uv_norm(
        sv_demo, s1_demo, s2_demo, s3_demo, ev_demo, e1_demo, e2_demo, e3_demo,
        t_value=0.5, resolution=256, out_path="vf3d_uv_norm.png", log_scale=False
    )
    visualize_scene(
        sv_demo, s1_demo, s2_demo, s3_demo, ev_demo, e1_demo, e2_demo, e3_demo,
        out_path="vf3d_scene.png",
        t_impact=t if ok else None,
        u_impact=u if ok else None,
        v_impact=v if ok else None,
    )


