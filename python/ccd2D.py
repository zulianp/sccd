#!/usr/bin/env python3

from typing import Tuple, Optional

# Reuse the 2D helpers implemented in numeric_roots.py
from numeric_roots import (
    find_root_line_2d as _find_root_line_2d_impl,
    visualize_Fvl_tu_2d as _visualize_Fvl_tu_2d_impl,
    visualize_Fvl_norm_tu_2d as _visualize_Fvl_norm_tu_2d_impl,
)


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
    return _find_root_line_2d_impl(max_iter, tol, sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d)


def visualize_F_tu_line(
    sv_2d, s1_2d, s2_2d, ev_2d, e1_2d, e2_2d,
    resolution: int = 256,
    out_path: str = "vl2d_Ftu.png",
) -> None:
    """
    Visualize the components (Fx,Fy) of the 2D vertex-line residual over (t,u).
    Saves a side-by-side heatmap image to out_path.
    """
    _visualize_Fvl_tu_2d_impl(
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
    _visualize_Fvl_norm_tu_2d_impl(
        sv_demo, s1_demo, s2_demo, ev_demo, e1_demo, e2_demo,
        resolution=256, out_path="vl2d_Fnorm.png", log_scale=False, zero_contour=1e-2
    )