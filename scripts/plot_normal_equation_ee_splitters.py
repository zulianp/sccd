#!/usr/bin/env python3
"""
Reimplementation of sccd::srootfinder::normal_equation_axis_splitters_ee
from src/srootfinder.hpp (lines ~1117–1178) plus diff_ee.

Plots adaptive bisection splitters along t, u, and v for a simple edge–edge
example on [0,1]³ with contact at **t = 0.25** (u = v = 0.5), plus a 3D figure:
swept quad for pair A, fixed segment for pair B, trajectories, and gap.
"""

from __future__ import annotations

import argparse
import sys

try:
    import numpy as np
    import matplotlib.pyplot as plt
except ImportError as exc:
    print(
        "This script needs NumPy and Matplotlib.\n"
        "  python3 -m venv .venv && source .venv/bin/activate\n"
        "  pip install numpy matplotlib\n",
        file=sys.stderr,
    )
    raise SystemExit(1) from exc


def diff_ee(
    s1: np.ndarray,
    s2: np.ndarray,
    s3: np.ndarray,
    s4: np.ndarray,
    e1: np.ndarray,
    e2: np.ndarray,
    e3: np.ndarray,
    e4: np.ndarray,
    t: float,
    u: float,
    v: float,
) -> np.ndarray:
    """Matches srootfinder.hpp diff_ee (vector length 3)."""
    diff = np.empty(3, dtype=np.float64)
    for d in range(3):
        ea0 = (e1[d] - s1[d]) * t + s1[d]
        ea1 = (e2[d] - s2[d]) * t + s2[d]
        eb0 = (e3[d] - s3[d]) * t + s3[d]
        eb1 = (e4[d] - s4[d]) * t + s4[d]
        diff[d] = ((ea1 - ea0) * u + ea0) - ((eb1 - eb0) * v + eb0)
    return diff


def normal_equation_axis_splitters_ee(
    split_dim: int,
    n: int,
    t_lo: float,
    t_hi: float,
    u_lo: float,
    u_hi: float,
    v_lo: float,
    v_hi: float,
    s1: np.ndarray,
    s2: np.ndarray,
    s3: np.ndarray,
    s4: np.ndarray,
    e1: np.ndarray,
    e2: np.ndarray,
    e3: np.ndarray,
    e4: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    split_dim: 0 -> t axis, 1 -> u, 2 -> v (same as C++ SplitDim).

    Returns (x0_uniform, splitters) each of shape (n,).
    """
    assert split_dim in (0, 1, 2)
    assert n > 0

    if split_dim == 0:
        lo, hi = t_lo, t_hi
    elif split_dim == 1:
        lo, hi = u_lo, u_hi
    else:
        lo, hi = v_lo, v_hi

    mid_t = 0.5 * (t_lo + t_hi)
    mid_u = 0.5 * (u_lo + u_hi)
    mid_v = 0.5 * (v_lo + v_hi)

    h = (hi - lo) / (n + 1)
    radius = h * 0.45
    eps = np.finfo(np.float64).eps

    base_t = 0.0 if split_dim == 0 else mid_t
    base_u = 0.0 if split_dim == 1 else mid_u
    base_v = 0.0 if split_dim == 2 else mid_v
    f_base = diff_ee(s1, s2, s3, s4, e1, e2, e3, e4, base_t, base_u, base_v)

    j_axis = np.empty(3, dtype=np.float64)
    for d in range(3):
        if split_dim == 0:
            j_axis[d] = (
                (1.0 - mid_u) * (e1[d] - s1[d])
                + mid_u * (e2[d] - s2[d])
                - (1.0 - mid_v) * (e3[d] - s3[d])
                - mid_v * (e4[d] - s4[d])
            )
        elif split_dim == 1:
            j_axis[d] = (1.0 - mid_t) * (s2[d] - s1[d]) + mid_t * (e2[d] - e1[d])
        else:
            j_axis[d] = -((1.0 - mid_t) * (s4[d] - s3[d]) + mid_t * (e4[d] - e3[d]))

    h_axis = float(np.dot(j_axis, j_axis))
    step_scale = 1.0 / h_axis if h_axis > eps else 0.00001

    x0 = np.empty(n, dtype=np.float64)
    splitters = np.empty(n, dtype=np.float64)
    for i in range(n):
        x0[i] = lo + h * (i + 1)
        xmin = max(lo, x0[i] - radius)
        xmax = min(hi, x0[i] + radius)
        g = 0.0
        for d in range(3):
            j = j_axis[d]
            g += (f_base[d] + x0[i] * j) * j
        step = g * step_scale
        splitters[i] = min(xmax, max(xmin, x0[i] - step))

    return x0, splitters


def endpoints_a(
    s1: np.ndarray,
    e1: np.ndarray,
    s2: np.ndarray,
    e2: np.ndarray,
    t: float | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Moving endpoints of the first edge pair at time(s) t."""
    t = np.asarray(t, dtype=np.float64)
    a0 = (e1 - s1) * t[..., None] + s1
    a1 = (e2 - s2) * t[..., None] + s2
    return a0, a1


def endpoints_b(
    s3: np.ndarray,
    e3: np.ndarray,
    s4: np.ndarray,
    e4: np.ndarray,
    t: float | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Moving endpoints of the second edge pair at time(s) t."""
    t = np.asarray(t, dtype=np.float64)
    b0 = (e3 - s3) * t[..., None] + s3
    b1 = (e4 - s4) * t[..., None] + s4
    return b0, b1


def point_on_a(
    s1: np.ndarray,
    e1: np.ndarray,
    s2: np.ndarray,
    e2: np.ndarray,
    t: float,
    u: float,
) -> np.ndarray:
    a0, a1 = endpoints_a(s1, e1, s2, e2, t)
    return (1.0 - u) * a0 + u * a1


def point_on_b(
    s3: np.ndarray,
    e3: np.ndarray,
    s4: np.ndarray,
    e4: np.ndarray,
    t: float,
    v: float,
) -> np.ndarray:
    b0, b1 = endpoints_b(s3, e3, s4, e4, t)
    return (1.0 - v) * b0 + v * b1


def plot_ee_geometry_3d(
    fig: plt.Figure,
    s1: np.ndarray,
    s2: np.ndarray,
    s3: np.ndarray,
    s4: np.ndarray,
    e1: np.ndarray,
    e2: np.ndarray,
    e3: np.ndarray,
    e4: np.ndarray,
    _t_lo: float,
    _t_hi: float,
    u_lo: float,
    u_hi: float,
    v_lo: float,
    v_hi: float,
    t_result: float,
) -> None:
    """
    3D view of the EE configuration:

    * Two translucent swept quads (each drawn as two triangles): motion of
      each moving edge segment from t=0 to t=1 (bilinear patch boundary).
    * Line trajectories of P_a(t) and P_b(t) with fixed u=mid_u, v=mid_v
      (same mids as in normal_equation_axis_splitters_ee).
    * Markers at P_a(t_result, mid_u), P_b(t_result, mid_v) and the segment
      between them (the 3D gap whose vector is diff_ee at that (t,u,v)).
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    mid_u = 0.5 * (u_lo + u_hi)
    mid_v = 0.5 * (v_lo + v_hi)

    ax = fig.add_subplot(111, projection="3d")

    # Swept strip boundary = quad (s1,s2,e2,e1): segment endpoints at t=0
    # and t=1, connected by motion of each endpoint.  Two tris = fill only.
    tri_a = np.array([[s1, s2, e2], [s1, e2, e1]], dtype=np.float64)
    pa = Poly3DCollection(
        tri_a, alpha=0.35, facecolor="C0", edgecolor="C0", linewidths=0.5
    )
    ax.add_collection3d(pa)

    # Pair B: here both edges are time-invariant (s3=e3, s4=e4) → no area
    # sweep; draw a thick translucent rod along the fixed EE segment.
    if np.allclose(s3, e3) and np.allclose(s4, e4):
        ax.plot(
            [s3[0], s4[0]],
            [s3[1], s4[1]],
            [s3[2], s4[2]],
            color="C2",
            linewidth=14.0,
            alpha=0.35,
            solid_capstyle="round",
            label="pair B fixed segment (no sweep)",
        )
    else:
        tri_b = np.array([[s3, s4, e4], [s3, e4, e3]], dtype=np.float64)
        pb = Poly3DCollection(
            tri_b, alpha=0.35, facecolor="C2", edgecolor="C2", linewidths=0.5
        )
        ax.add_collection3d(pb)

    # Highlight quad boundaries (4 edges of each strip)
    # loop_a = np.vstack([s1, s2, e2, e1, s1])
    # loop_b = np.vstack([s3, s4, e4, e3, s3])
    # ax.plot(
    #     loop_a[:, 0],
    #     loop_a[:, 1],
    #     loop_a[:, 2],
    #     color="C0",
    #     linewidth=1.8,
    #     label="quad boundary (pair A)",
    # )
    # ax.plot(
    #     loop_b[:, 0],
    #     loop_b[:, 1],
    #     loop_b[:, 2],
    #     color="C2",
    #     linewidth=1.8,
    #     label="quad boundary (pair B)",
    # )

    ts = np.linspace(0.0, 1.0, 80)
    a0, a1 = endpoints_a(s1, e1, s2, e2, ts)
    b0, b1 = endpoints_b(s3, e3, s4, e4, ts)
    traj_a = (1.0 - mid_u) * a0 + mid_u * a1
    traj_b = (1.0 - mid_v) * b0 + mid_v * b1
    ax.plot(
        traj_a[:, 0],
        traj_a[:, 1],
        traj_a[:, 2],
        color="C1",
        linewidth=2.0,
        label="Pa(t) at fixed mid_u",
    )
    if np.allclose(traj_b[0], traj_b[-1], atol=1e-12):
        ax.scatter(
            traj_b[0, 0],
            traj_b[0, 1],
            traj_b[0, 2],
            color="C3",
            s=120,
            depthshade=True,
            marker="D",
            label="Pb(t) at fixed mid_v (constant — pair B static)",
        )
    else:
        ax.plot(
            traj_b[:, 0],
            traj_b[:, 1],
            traj_b[:, 2],
            color="C3",
            linewidth=2.0,
            label="Pb(t) at fixed mid_v",
        )

    pa_mid = point_on_a(s1, e1, s2, e2, t_result, mid_u)
    pb_mid = point_on_b(s3, e3, s4, e4, t_result, mid_v)
    if np.linalg.norm(pa_mid - pb_mid) <= 1e-9:
        ax.scatter(
            *pa_mid,
            color="gold",
            s=260,
            depthshade=True,
            marker="*",
            edgecolors="k",
            linewidths=1.0,
            label=f"contact Pa = Pb (t={t_result:g}, mid_u, mid_v)",
        )
    else:
        ax.scatter(
            *pa_mid,
            color="C1",
            s=100,
            depthshade=True,
            marker="*",
            label=f"Pa(t={t_result:g}, mid_u)",
        )
        ax.scatter(
            *pb_mid,
            color="C3",
            s=100,
            depthshade=True,
            marker="*",
            label=f"Pb(t={t_result:g}, mid_v)",
        )
        gap = np.stack([pa_mid, pb_mid], axis=0)
        ax.plot(
            gap[:, 0],
            gap[:, 1],
            gap[:, 2],
            "k-",
            linewidth=2.0,
            label="gap Pa - Pb (diff_ee)",
        )

    # Initial (t=0) and final (t=1) EE segments: thick lines — these are the
    # two segments Pa lives on at t=0 / t=1 (pair A) and same for Pb (pair B).
    lw_pose = 6.0
    ax.plot(
        [s1[0], s2[0]],
        [s1[1], s2[1]],
        [s1[2], s2[2]],
        color="C0",
        linewidth=lw_pose,
        solid_capstyle="round",
        label="pair A segment (t=0)",
    )
    ax.plot(
        [e1[0], e2[0]],
        [e1[1], e2[1]],
        [e1[2], e2[2]],
        color="C0",
        linewidth=lw_pose,
        linestyle="--",
        solid_capstyle="round",
        alpha=0.95,
        label="pair A segment (t=1)",
    )
    ax.plot(
        [s3[0], s4[0]],
        [s3[1], s4[1]],
        [s3[2], s4[2]],
        color="C2",
        linewidth=lw_pose,
        solid_capstyle="round",
        label="pair B segment (t=0)",
    )
    ax.plot(
        [e3[0], e4[0]],
        [e3[1], e4[1]],
        [e3[2], e4[2]],
        color="C2",
        linewidth=lw_pose,
        linestyle="--",
        solid_capstyle="round",
        alpha=0.95,
        label="pair B segment (t=1)",
    )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title(
        "EE: pair A sweeps a quad; pair B is fixed. Contact when Pa meets Pb.\n"
        f"Highlighted pose uses t = {t_result:g} (diff_ee = 0 at u = v = 0.5); "
        "thick solid = t=0, thick dashed = t=1.",
    )
    ax.legend(loc="upper left", fontsize=7)
    _set_axes_equal(ax)


def _set_axes_equal(ax) -> None:
    """Roughly equal scale on x,y,z for mplot3d."""
    limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    center = limits.mean(axis=1)
    span = (limits[:, 1] - limits[:, 0]).max() * 0.55
    ax.set_xlim3d(center[0] - span, center[0] + span)
    ax.set_ylim3d(center[1] - span, center[1] + span)
    ax.set_zlim3d(center[2] - span, center[2] + span)
    try:
        ax.set_box_aspect((1, 1, 1))
    except AttributeError:
        pass


def simple_ee_example() -> tuple[
    tuple[float, float, float, float, float, float],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
]:
    """
    Edge–edge configuration with a **single** collision time t = 0.25 (and
    u = v = 0.5 in the search domain [0,1]³).

    Pair A: both edges rise in z with t; at fixed u = 0.5, Pa = (0.5, 0.5, t).

    Pair B: both constituent edges are **static** (s3 = e3, s4 = e4), so the
    segment between B0 and B1 is fixed from (0,1,0) to (1,0,0.5).  At v = 0.5,
    Pb = (0.5, 0.5, 0.25) for all t.  Then diff_ee = (0, 0, t − 0.25), vanishing
    only at t = 0.25.

    Returns t_collision = 0.25 for marking the contact in the 3D plot.
    """
    s1 = np.array([0.0, 0.0, 0.0])
    e1 = np.array([0.0, 0.0, 1.0])
    s2 = np.array([1.0, 1.0, 0.0])
    e2 = np.array([1.0, 1.0, 1.0])
    s3 = np.array([0.0, 1.0, 0.0])
    e3 = np.array([0.0, 1.0, 0.0])  # static first edge of pair B
    s4 = np.array([1.0, 0.0, 0.5])
    e4 = np.array([1.0, 0.0, 0.5])  # static second edge of pair B
    bounds = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    t_collision = 0.25
    return bounds, s1, s2, s3, s4, e1, e2, e3, e4, t_collision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-n", type=int, default=4, help="number of uniform candidates per axis"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="normal_equation_ee_splitters.png",
        help="output PNG path for 2D splitter plots (relative to cwd unless absolute)",
    )
    parser.add_argument(
        "--output-3d",
        default=None,
        metavar="PATH",
        help="output PNG for 3D geometry plot (default: same stem as -o with _geometry.png)",
    )
    parser.add_argument(
        "--no-3d", action="store_true", help="skip the 3D geometry figure"
    )
    parser.add_argument(
        "--show", action="store_true", help="open an interactive window after saving"
    )
    args = parser.parse_args()

    n = args.n
    (t_lo, t_hi, u_lo, u_hi, v_lo, v_hi), s1, s2, s3, s4, e1, e2, e3, e4, t_collision = (
        simple_ee_example()
    )

    axis_names = ("t", "u", "v")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=False)
    titles = tuple(f"SplitDim {d} ({axis_names[d]} axis)" for d in (0, 1, 2))

    for ax, split_dim, title in zip(axes, (0, 1, 2), titles):
        x0, sp = normal_equation_axis_splitters_ee(
            split_dim,
            n,
            t_lo,
            t_hi,
            u_lo,
            u_hi,
            v_lo,
            v_hi,
            s1,
            s2,
            s3,
            s4,
            e1,
            e2,
            e3,
            e4,
        )
        ax.axhline(0, color="k", linewidth=0.5, alpha=0.3)
        ax.scatter(np.arange(n), x0, c="C0", s=22, label="uniform x₀", zorder=3)
        ax.scatter(
            np.arange(n), sp, c="C1", s=22, marker="s", label="splitter", zorder=3
        )
        for i in range(n):
            ax.plot([i, i], [x0[i], sp[i]], "k-", alpha=0.25, linewidth=1)
        ax.set_xlabel("candidate index i")
        ax.set_ylabel(f"{axis_names[split_dim]} (split position)")
        ax.set_title(title)
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "normal_equation_axis_splitters_ee: uniform candidates vs Newton-like step\n"
        f"(EE contact at t = {t_collision:g}, u = v = 0.5; domain [0,1]³)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")
    if args.show:
        plt.show()
    plt.close(fig)

    if not args.no_3d:
        import os

        if args.output_3d:
            out3d = args.output_3d
        else:
            root, ext = os.path.splitext(args.output)
            ext = ext if ext else ".png"
            out3d = f"{root}_geometry{ext}"

        mid_u = 0.5 * (u_lo + u_hi)
        mid_v = 0.5 * (v_lo + v_hi)
        _, sp_t = normal_equation_axis_splitters_ee(
            0, n, t_lo, t_hi, u_lo, u_hi, v_lo, v_hi, s1, s2, s3, s4, e1, e2, e3, e4
        )
        t_star = float(sp_t[n // 2])

        fig3d = plt.figure(figsize=(8, 7))
        plot_ee_geometry_3d(
            fig3d,
            s1,
            s2,
            s3,
            s4,
            e1,
            e2,
            e3,
            e4,
            t_lo,
            t_hi,
            u_lo,
            u_hi,
            v_lo,
            v_hi,
            t_collision,
        )
        fig3d.tight_layout()
        fig3d.savefig(out3d, dpi=150)
        print(
            f"Wrote {out3d} (contact t = {t_collision:g}; mid adaptive t-splitter ≈ {t_star:.4f}; mids u={mid_u:.3f}, v={mid_v:.3f})"
        )
        if args.show:
            plt.show()
        plt.close(fig3d)


if __name__ == "__main__":
    main()
