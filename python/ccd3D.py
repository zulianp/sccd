#!/usr/bin/env python3

from typing import Tuple, Optional

from numeric_roots import (
    find_root_vf_3d as _find_root_vf_3d_impl,
    visualize_Fvf_uv_norm_3d as _visualize_Fvf_uv_norm_3d_impl,
    visualize_scene_vf_3d as _visualize_scene_vf_3d_impl,
)


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
    return _find_root_vf_3d_impl(
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
    _visualize_Fvf_uv_norm_3d_impl(
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
    _visualize_scene_vf_3d_impl(
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


