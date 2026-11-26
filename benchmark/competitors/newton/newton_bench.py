"""
Does Newton's self-contact detection find the collisions in the CCD datasets?

    python newton_bench.py <data-dir> <scene> [--csv out.csv] [--repeats N] [--max-steps N]

Newton (github.com/newton-physics/newton) detects cloth self-contact by
proximity, not by swept volumes: SolverVBD refits two Warp BVHs at the current
positions and asks, per vertex and per edge, for the triangles and edges within
a query radius (`particle_self_contact_margin`, default 0.2), then keeps the
simulation penetration-free by truncating each vertex's displacement until the
next detection. The detector is used here exactly as SolverVBD uses it -- the
solver is constructed with self-contact on and its defaults, and its own
detection step is called -- and nothing of Newton is changed.

The datasets' trajectories were not truncated by Newton, so the question this
answers is whether detection at the start of a step returns every pair the
dataset says collides during that step. For each case `<step><vf|ee>` of a
scene, the ground truth is the dataset's candidate list (boxes/<key>/c0,c1)
restricted to the queries with an exact root (roots/<key>/toi.float64 not NaN)
or, without roots, a true mma_bool. A colliding pair Newton does not return is
a miss, and is attributed to:

  filtered   Newton excludes the pair as topologically close
             (particle_topological_contact_filter_threshold, default 2 rings);
  overflow   the primitive's contact buffer was full (fixed sizes, 32 per
             vertex and 64 per edge by default, not resized at runtime);
  distance   otherwise: the primitives were farther apart than the radius at
             the start of the step.

A pair checked in several sub-steps (see --substep) and returned in none takes
the strongest cause any of them gave, in the order filtered, overflow, distance:
a pair that is far apart in the first sub-steps and lost to a full buffer in a
later one is an overflow.

The radius is in scene units, and the scenes are not to one scale, so each row
also carries the scene's bounding-box diagonal at the start of the step.

The radius is a SolverVBD parameter. By default it is left at Newton's 0.2.
--margin sets it to a value, --margin edge:<f> to f times the scene's mean edge
length (a radius on the mesh's own scale, for a sweep from very small upwards,
whatever it misses), and --margin auto to the smallest one under which
VBD would not have had to truncate the dataset's motion: between detections VBD
caps each vertex's displacement at margin * relaxation * 0.5 (relaxation 0.85 by
default, apply_truncation_ts), so the margin is the scene's largest per-step
vertex displacement divided by 0.425. It comes from the motion alone, not from
the ground truth. The contact radius is then set to half the margin, inside the
1.5-2x ratio Newton's documentation advises.

The contact buffers are SolverVBD parameters too (32 entries per vertex, 64 per
edge), and Newton counts every contact even when a buffer overflows. With
--buffers auto, a detection that overflows rebuilds the solver with buffers 1.25x
the largest count seen and is run again, as a user would size them; the rebuild
is setup, not detection. Growth is clamped to --buffer-budget-gb of buffer memory and to
Warp's limit of 2^31 entries per array, and whatever still overflows is
reported as overflow.

The edge-edge distance comes from Warp's closest_point_edge_edge, which treats
an edge as a point at its first vertex when its squared length is at most
SolverVBD's particle_edge_parallel_epsilon (1e-5, an absolute value in scene
units). On a mesh with edges shorter than 3.2e-3 that overestimates the
distance and drops pairs within the radius. --edge-epsilon auto keeps Newton's
1e-5 when the shortest edge h_min of the first curated step is longer than
3.2e-3, and otherwise sets it to (1e-2 h_min)^2, so that no edge counts as a point.

With --substep the detector is given the guarantee VBD relies on. Its solver
keeps a simulation free of penetration by truncation, not by finding every
collision of a motion: each stored pair gets a separating plane, and no vertex
moves more than 0.425 r between two detections (apply_truncation_ts), so two
primitives farther apart than r cannot meet before the next one. A dataset step
moves farther than that, so it is covered the way the cap requires: in
K = ceil(d / (0.425 r)) sub-steps, d the step's largest vertex displacement, with
one detection at the start of each, at positions on the dataset's own linear
trajectory. A colliding pair with an exact root t is then checked in the sub-step
that contains t, where it is within 0.85 r of contact at the sub-step's start; a
pair the dataset marks colliding without a root counts as returned if any
sub-step returns it; the curated cases ship their roots as archives of
Mathematica files (roots/<key>_roots.tar.gz), not as toi.float64, so every pair
of them is of this kind. A miss for
distance is impossible by construction, and the run reports it so that it is
checked. --margin auto/<s> sets the radius to r_m / s,
which needs up to s sub-steps per step.

Timings are the median over --repeats of SolverVBD's detection step (BVH refit
plus both queries), after one warm-up, with the device synchronised. The step is
made of the calls SolverVBD's _collision_detection_penetration_free makes, in its
order, with the device synchronised after each, so that it is also timed in three
parts: the refit (with the reset of the two arrays that precedes it), the
vertex-triangle query and the edge-edge query; detect_ms is their sum. Building
the model and the solver (BVH construction, filter lists) is timed once per
topology and reported separately. Newton stores positions in float32; that is
part of what is measured. With --substep a step's time is the sum of its K
detections, each timed once, after one warm-up per topology.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import warp as wp
import newton
from newton.solvers import SolverVBD

FIELDS = [
    "dataset", "case", "type", "queries", "colliding", "returned_pairs",
    "missed", "missed_filtered", "missed_overflow", "missed_distance",
    "detect_ms", "refit_ms", "vt_ms", "ee_ms", "setup_ms", "radius", "scene_diag", "vt_buffer", "ee_buffer", "substeps", "edge_epsilon",
]
# A pair checked in several sub-steps and returned in none is attributed to the
# strongest cause any of them gave.
CAUSE_RANK = {"distance": 0, "overflow": 1, "filtered": 2}


def frame_stem(scene: str, step: int) -> str:
    """Same naming as frame_path() in benchmark/bench.exe.cpp."""
    if scene == "cloth-ball":
        return f"cloth_ball{step}"
    if scene == "n-body-simulation":
        return f"balls16_{step}"
    if scene == "cloth-funnel":
        return f"{step:03d}"
    return str(step)


def read_raw_mesh(folder: Path):
    xyz = []
    for c in "xyz":
        f64, f32 = folder / f"{c}.float64", folder / f"{c}.float32"
        if f64.exists():
            xyz.append(np.fromfile(f64, dtype=np.float64))
        elif f32.exists():
            xyz.append(np.fromfile(f32, dtype=np.float32).astype(np.float64))
        else:
            raise FileNotFoundError(f"no {c} coordinates in {folder}")
    faces = np.stack([np.fromfile(folder / f"i{k}.int32", dtype=np.int32) for k in range(3)], axis=1)
    return np.stack(xyz, axis=1), faces


def benchmark_edges(faces: np.ndarray) -> np.ndarray:
    """The dataset's edge numbering: unique edges in face order, sorted by
    (larger vertex, smaller vertex) -- benchmark_ordered_edges() in bench.exe.cpp."""
    seen, keys = set(), []
    for tri in faces:
        for e in range(3):
            a, b = sorted((int(tri[e]), int(tri[(e + 1) % 3])))
            if (a, b) not in seen:
                seen.add((a, b))
                keys.append((a, b))
    keys.sort(key=lambda ab: (ab[1], ab[0]))
    return np.array(keys, dtype=np.int64)


def ground_truth(case_dir: Path, roots_dir: Path, mma_dir: Path):
    c0 = np.fromfile(case_dir / "c0.int32", dtype=np.int32).astype(np.int64)
    c1 = np.fromfile(case_dir / "c1.int32", dtype=np.int32).astype(np.int64)
    toi = roots_dir / "toi.float64"
    if toi.exists():
        roots = np.fromfile(toi, dtype=np.float64)
        collides = ~np.isnan(roots)
    else:
        collides = np.fromfile(mma_dir / "mma_bool.uint8", dtype=np.uint8).astype(bool)
        roots = np.full(collides.size, np.nan)
    if collides.size != c0.size:
        raise ValueError(f"{case_dir.name}: {c0.size} queries but {collides.size} answers")
    return c0, c1, collides, roots


class NewtonScene:
    """Newton's model and solver for one mesh topology, built as a user would."""

    def __init__(self, vertices: np.ndarray, faces: np.ndarray, margin: float | None = None,
                 vt_buffer: int = 32, ee_buffer: int = 64, edge_epsilon: float = 1e-5):
        t = time.perf_counter()
        builder = newton.ModelBuilder()
        builder.add_cloth_mesh(
            pos=wp.vec3(0.0, 0.0, 0.0), rot=wp.quat_identity(), scale=1.0, vel=wp.vec3(0.0, 0.0, 0.0),
            vertices=[wp.vec3(*p) for p in vertices], indices=faces.reshape(-1).tolist(), density=1.0,
        )
        builder.color()
        self.model = builder.finalize(device="cuda")
        tri = self.model.tri_indices.numpy()
        if tri.shape != faces.shape or not np.array_equal(tri, faces):
            raise RuntimeError("Newton reordered or dropped triangles; ids would not match the dataset")
        kw = dict(particle_vertex_contact_buffer_size=vt_buffer, particle_edge_contact_buffer_size=ee_buffer,
                  particle_edge_parallel_epsilon=edge_epsilon)
        if margin is not None:
            kw.update(particle_self_contact_margin=margin, particle_self_contact_radius=0.5 * margin)
        self.solver = SolverVBD(self.model, particle_enable_self_contact=True, **kw)
        self.vt_buffer, self.ee_buffer, self.edge_epsilon = vt_buffer, ee_buffer, edge_epsilon
        self.state = self.model.state()
        wp.synchronize()
        self.setup_ms = (time.perf_counter() - t) * 1e3

        self.radius = float(self.solver.particle_self_contact_margin)
        edges = self.model.edge_indices.numpy()[:, 2:4].astype(np.int64)
        self.edges = edges
        self.edge_id = {(min(a, b), max(a, b)): i for i, (a, b) in enumerate(edges)}
        s = self.solver
        self.vt_filter = self._filter_sets(s.particle_vertex_triangle_contact_filtering_list,
                                           s.particle_vertex_triangle_contact_filtering_list_offsets)
        self.ee_filter = self._filter_sets(s.particle_edge_edge_contact_filtering_list,
                                           s.particle_edge_edge_contact_filtering_list_offsets)

    @staticmethod
    def _filter_sets(values, offsets):
        if values is None or offsets is None:
            return None
        v, o = values.numpy(), offsets.numpy()
        return [set(v[o[i]:o[i + 1]].tolist()) for i in range(len(o) - 1)]

    def _detection_step(self):
        """SolverVBD._collision_detection_penetration_free (Newton 1.3.0), call for
        call, timed in three parts: refit, vertex-triangle query, edge-edge query."""
        s, q = self.solver, self.state.particle_q
        det = s.trimesh_collision_detector
        t0 = time.perf_counter()
        s.pos_prev_collision_detection.assign(q)
        s.particle_displacements.zero_()
        det.refit(q)
        wp.synchronize()
        t1 = time.perf_counter()
        det.vertex_triangle_collision_detection(
            s.particle_self_contact_margin,
            min_query_radius=s.particle_rest_shape_contact_exclusion_radius,
            min_distance_filtering_ref_pos=s.particle_q_rest)
        wp.synchronize()
        t2 = time.perf_counter()
        det.edge_edge_collision_detection(
            s.particle_self_contact_margin,
            min_query_radius=s.particle_rest_shape_contact_exclusion_radius,
            min_distance_filtering_ref_pos=s.particle_q_rest)
        wp.synchronize()
        t3 = time.perf_counter()
        return 1e3 * (t1 - t0), 1e3 * (t2 - t1), 1e3 * (t3 - t2)

    def detect(self, positions: np.ndarray, repeats: int, warmup: bool = True):
        """The median over repeats of each part, (refit, vt, ee) in ms, and the lists."""
        self.state.particle_q.assign(positions.astype(np.float32))
        wp.synchronize()
        if warmup:
            self._detection_step()
        parts = [self._detection_step() for _ in range(repeats)]
        info = self.solver.trimesh_collision_detector.collision_info
        return tuple(statistics.median(p[i] for p in parts) for i in range(3)), {
            "vt": [a.numpy() for a in (info.vertex_colliding_triangles, info.vertex_colliding_triangles_offsets,
                                       info.vertex_colliding_triangles_buffer_sizes,
                                       info.vertex_colliding_triangles_count)],
            "ee": [a.numpy() for a in (info.edge_colliding_edges, info.edge_colliding_edges_offsets,
                                       info.edge_colliding_edges_buffer_sizes, info.edge_colliding_edges_count)],
        }


def returned(lists, primitive: int):
    """What Newton returned for one primitive, and whether its buffer overflowed.
    Records are pairs (primitive, other) laid out from 2*offset."""
    flat, offsets, sizes, counts = lists
    n = min(int(counts[primitive]), int(sizes[primitive]))
    start = 2 * int(offsets[primitive])
    others = set(flat[start + 1 : start + 2 * n : 2].tolist())
    return others, int(counts[primitive]) > int(sizes[primitive])


def point_triangle_distance(p, a, b, c):
    """Exact distance from p to triangle abc (Ericson, Real-Time Collision Detection 5.1.5)."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = ab @ ap, ac @ ap
    if d1 <= 0 and d2 <= 0:
        return np.linalg.norm(p - a)
    bp = p - b
    d3, d4 = ab @ bp, ac @ bp
    if d3 >= 0 and d4 <= d3:
        return np.linalg.norm(p - b)
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        return np.linalg.norm(p - (a + d1 / (d1 - d3) * ab))
    cp = p - c
    d5, d6 = ab @ cp, ac @ cp
    if d6 >= 0 and d5 <= d6:
        return np.linalg.norm(p - c)
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        return np.linalg.norm(p - (a + d2 / (d2 - d6) * ac))
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        return np.linalg.norm(p - (b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)))
    denom = 1.0 / (va + vb + vc)
    return np.linalg.norm(p - (a + ab * vb * denom + ac * vc * denom))


def validate(vertices, faces, lists, radius, n_sample=2000, seed=0):
    """Positive control: the vertex-triangle pairs Newton returned must lie
    within the radius. Returns (pairs checked, largest distance)."""
    flat, offsets, sizes, counts = lists
    rng = np.random.default_rng(seed)
    vs = np.flatnonzero(np.minimum(counts, sizes) > 0)
    worst, checked = 0.0, 0
    for v in rng.choice(vs, size=min(n_sample, vs.size), replace=False):
        others, _ = returned(lists, int(v))
        t = int(rng.choice(sorted(others)))
        a, b, c = vertices[faces[t]]
        worst = max(worst, point_triangle_distance(vertices[v], a, b, c))
        checked += 1
    return checked, worst


def segment_distance(p1, q1, p2, q2):
    """Exact distance between segments p1q1 and p2q2 (Ericson 5.1.9)."""
    d1, d2, r = q1 - p1, q2 - p2, p1 - p2
    a, e, f = d1 @ d1, d2 @ d2, d2 @ r
    if a <= 1e-30 and e <= 1e-30:
        return float(np.linalg.norm(r))
    if a <= 1e-30:
        s, t = 0.0, min(max(f / e, 0.0), 1.0)
    else:
        c = d1 @ r
        if e <= 1e-30:
            t, s = 0.0, min(max(-c / a, 0.0), 1.0)
        else:
            b = d1 @ d2
            denom = a * e - b * b
            s = min(max((b * f - c * e) / denom, 0.0), 1.0) if denom > 1e-30 else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t, s = 0.0, min(max(-c / a, 0.0), 1.0)
            elif t > 1.0:
                t, s = 1.0, min(max((b - c) / a, 0.0), 1.0)
    return float(np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t)))


def validate_ee(vertices, edges, lists, n_sample=2000, seed=0):
    """Positive control for edges: the edge-edge pairs Newton returned must lie
    within the radius. `edges` is Newton's edge list (vertex pairs)."""
    flat, offsets, sizes, counts = lists
    rng = np.random.default_rng(seed)
    es = np.flatnonzero(np.minimum(counts, sizes) > 0)
    worst, checked = 0.0, 0
    for e in rng.choice(es, size=min(n_sample, es.size), replace=False):
        others, _ = returned(lists, int(e))
        o = int(rng.choice(sorted(others)))
        (a, b), (c, d) = edges[int(e)], edges[o]
        worst = max(worst, segment_distance(vertices[a], vertices[b], vertices[c], vertices[d]))
        checked += 1
    return checked, worst


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("scene")
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--max-steps", type=int, default=0, help="stop after this many steps (0: all)")
    ap.add_argument("--margin", default=None,
                    help="SolverVBD's particle_self_contact_margin: a number, 'edge:<f>' (f times the mean edge "
                         "length), or 'auto' (see above); "
                         "default: Newton's own")
    ap.add_argument("--buffers", choices=("default", "auto"), default="default",
                    help="contact buffer sizes: Newton's defaults, or grown on overflow (see above)")
    ap.add_argument("--buffer-budget-gb", type=float, default=40.0,
                    help="largest contact-buffer memory --buffers auto may allocate")
    ap.add_argument("--edge-epsilon", default="auto",
                    help="SolverVBD's particle_edge_parallel_epsilon: a number, 'default' (1e-5), or 'auto' "
                         "(see above)")
    ap.add_argument("--substep", action="store_true",
                    help="cover each step in as many detections as VBD's displacement cap requires (see above)")
    ap.add_argument("--validate", action="store_true",
                    help="check that sampled returned vertex-triangle and edge-edge pairs lie within the radius")
    args = ap.parse_args(argv[1:])

    wp.init()
    if newton.__version__ != "1.3.0":
        # _detection_step repeats SolverVBD's detection call for call, as 1.3.0 has it.
        print(f"newton {newton.__version__}: the detection step is transcribed from 1.3.0", file=sys.stderr)
        return 2
    scene_dir = args.data_dir / args.scene
    # The cases with curated queries, as scan_cases() in the C++ harnesses picks
    # them: the dataset also has a box directory, with empty pair lists, for
    # every step it curated nothing for.
    cases = sorted((d for d in (scene_dir / "boxes").iterdir()
                    if d.is_dir() and d.name[-2:] in ("vf", "ee") and (d / "c0.int32").exists()
                    and (scene_dir / "queries" / f"{d.name}.csv").exists()),
                   key=lambda d: (int(d.name[:-2]), d.name[-2:]))
    by_step: dict[int, list[Path]] = {}
    for d in cases:
        by_step.setdefault(int(d.name[:-2]), []).append(d)

    relaxation = 0.85  # SolverVBD's default particle_conservative_bound_relaxation
    cap = relaxation * 0.5  # apply_truncation_ts: displacement <= margin * cap

    def step_displacement(step):
        x0, _ = read_raw_mesh(scene_dir / "frames_raw" / frame_stem(args.scene, step))
        x1, _ = read_raw_mesh(scene_dir / "frames_raw" / frame_stem(args.scene, step + 1))
        return float(np.linalg.norm(x1 - x0, axis=1).max())

    margin = None
    if args.margin is not None and (args.margin == "auto" or args.margin.startswith("auto/")):
        # VBD's isotropic cap on the displacement between detections, solved for
        # the margin that leaves every vertex of every curated step untruncated;
        # auto/<s> divides it by s.
        max_disp = max(step_displacement(step) for step in sorted(by_step))
        r_m = max_disp / cap
        margin = r_m / (float(args.margin[5:]) if args.margin.startswith("auto/") else 1.0)
        print(f"{args.scene}: largest per-step vertex displacement {max_disp:.6g}; r_m {r_m:.6g}; "
              f"margin {margin:.6g}", file=sys.stderr)
    elif args.margin is not None and args.margin.startswith("edge:"):
        # A radius relative to the mesh, so that "small" means the same on scenes
        # of different scale: a fraction of the mean edge length at the first step.
        x0, f0 = read_raw_mesh(scene_dir / "frames_raw" / frame_stem(args.scene, min(by_step)))
        e = benchmark_edges(f0)
        mean_edge = float(np.linalg.norm(x0[e[:, 0]] - x0[e[:, 1]], axis=1).mean())
        margin = float(args.margin[5:]) * mean_edge
        print(f"{args.scene}: mean edge length {mean_edge:.6g}; margin {margin:.6g}", file=sys.stderr)
    elif args.margin is not None:
        margin = float(args.margin)

    if args.edge_epsilon == "auto":
        x0, f0 = read_raw_mesh(scene_dir / "frames_raw" / frame_stem(args.scene, min(by_step)))
        e0 = benchmark_edges(f0)
        h_min = float(np.linalg.norm(x0[e0[:, 0]] - x0[e0[:, 1]], axis=1).min())
        edge_epsilon = 1e-5 if h_min ** 2 > 1e-5 else (1e-2 * h_min) ** 2
        print(f"{args.scene}: shortest edge {h_min:.6g}; edge epsilon {edge_epsilon:.6g}", file=sys.stderr)
    else:
        edge_epsilon = 1e-5 if args.edge_epsilon == "default" else float(args.edge_epsilon)

    out = open(args.csv, "w", newline="") if args.csv else sys.stdout
    writer = csv.DictWriter(out, fieldnames=FIELDS)
    writer.writeheader()
    totals = dict(colliding=0, missed=0, filtered=0, overflow=0, distance=0)
    scene, scene_faces = None, None

    def detect_sized(positions, faces, step, repeats, warmup):
        """One detection, grown and re-run on overflow under --buffers auto.
        Returns the detection time and the lists; rebuilds count as setup."""
        nonlocal scene
        setup = 0.0
        ms, lists = scene.detect(positions, repeats, warmup)
        while args.buffers == "auto":
            need_vt, need_ee = int(lists["vt"][3].max(initial=0)), int(lists["ee"][3].max(initial=0))
            if need_vt <= scene.vt_buffer and need_ee <= scene.ee_buffer:
                break
            n_v, n_e = len(positions), len(lists["ee"][1]) - 1
            vt = max(scene.vt_buffer, math.ceil(1.25 * need_vt))
            ee = max(scene.ee_buffer, math.ceil(1.25 * need_ee))
            # Two limits, and the buffers are clamped to them rather than left
            # where they were, so Newton gets the largest buffers it can have:
            # Warp indexes an array with int32 and stores each contact as a pair
            # of entries, so a buffer holds fewer than 2^30 contacts in all; and
            # the buffers together stay within --buffer-budget-gb.
            vt = min(vt, (2**31 - 1) // (2 * n_v))
            ee = min(ee, (2**31 - 1) // (2 * max(n_e, 1)))
            gb = 8.0 * (n_v * vt + n_e * ee) / 1e9  # a pair of int32 per entry
            if gb > args.buffer_budget_gb:
                scale = args.buffer_budget_gb / gb
                vt, ee = int(vt * scale), int(ee * scale)
                gb = 8.0 * (n_v * vt + n_e * ee) / 1e9
            clamped = vt < math.ceil(1.25 * need_vt) or ee < math.ceil(1.25 * need_ee)
            if vt <= scene.vt_buffer and ee <= scene.ee_buffer:
                print(f"step {step}: {need_vt} vertex and {need_ee} edge contacts; buffers already at the "
                      f"limit, {scene.vt_buffer}/{scene.ee_buffer}", file=sys.stderr)
                break
            vt, ee = max(vt, scene.vt_buffer), max(ee, scene.ee_buffer)
            print(f"step {step}: growing buffers to {vt} per vertex, {ee} per edge ({gb:.2f} GB)"
                  f"{', clamped to the limits' if clamped else ''}", file=sys.stderr)
            scene = NewtonScene(positions, faces, margin, vt, ee, edge_epsilon)
            setup += scene.setup_ms
            ms, lists = scene.detect(positions, repeats, True)
        return ms, lists, setup

    for n_step, (step, step_cases) in enumerate(sorted(by_step.items())):
        if args.max_steps and n_step >= args.max_steps:
            break
        vertices, faces = read_raw_mesh(scene_dir / "frames_raw" / frame_stem(args.scene, step))
        setup_ms = 0.0
        if scene is None or not np.array_equal(faces, scene_faces):
            sizes = (32, 64) if scene is None or scene_faces is None else (scene.vt_buffer, scene.ee_buffer)
            scene, scene_faces = NewtonScene(vertices, faces, margin, *sizes, edge_epsilon), faces
            setup_ms = scene.setup_ms
            if args.substep:
                scene.detect(vertices, 1, True)  # warm-up per topology, untimed
        n_nodes = len(vertices)
        dataset_edges = benchmark_edges(faces)
        n_edges = len(dataset_edges)
        diag = float(np.linalg.norm(vertices.max(axis=0) - vertices.min(axis=0)))

        # The step's sub-steps: one detection at the start of each.
        if args.substep:
            end, _ = read_raw_mesh(scene_dir / "frames_raw" / frame_stem(args.scene, step + 1))
            d = float(np.linalg.norm(end - vertices, axis=1).max())
            # A relative slack keeps a step that needs exactly s sub-steps at s,
            # where rounding in d / (cap * r) would otherwise make it s + 1.
            K = max(1, math.ceil(d / (cap * scene.radius) * (1 - 1e-9)))
        else:
            end, K = vertices, 1

        # Every colliding pair of the step, keyed by the sub-step it must be
        # returned in: the one containing its exact root, or -1 without a root.
        pending = []  # (case index, is_vf, key, x, y, sub-step)
        case_info = []
        for ci, case in enumerate(step_cases):
            is_vf = case.name.endswith("vf")
            c0, c1, collides, roots = ground_truth(case, scene_dir / "roots" / case.name,
                                                   scene_dir / "mma_bool" / case.name)
            case_info.append((case, is_vf, c0.size, int(collides.sum())))
            for x, y, t in zip(c0[collides], c1[collides], roots[collides]):
                k = -1 if np.isnan(t) else min(K - 1, int(math.floor(float(t) * K)))
                pending.append([ci, is_vf, int(x), int(y), k])
        found = [False] * len(pending)
        cause = [None] * len(pending)
        by_sub = {}
        for i, pr in enumerate(pending):
            by_sub.setdefault(pr[4], []).append(i)

        def check(lists, idx):
            ci, is_vf, x, y, _k = pending[idx]
            if is_vf:
                v, f = (x, y - n_nodes - n_edges) if x < n_nodes else (y, x - n_nodes - n_edges)
                others, full = returned(lists["vt"], v)
                if f in others:
                    return True, None
                filt = scene.vt_filter is not None and f in scene.vt_filter[v]
            else:
                ea, eb = dataset_edges[x - n_nodes], dataset_edges[y - n_nodes]
                a, b = scene.edge_id[tuple(ea)], scene.edge_id[tuple(eb)]
                oa, fa = returned(lists["ee"], a)
                ob, fb = returned(lists["ee"], b)
                if b in oa or a in ob:
                    return True, None
                full = fa or fb
                filt = scene.ee_filter is not None and (b in scene.ee_filter[a] or a in scene.ee_filter[b])
            return False, "filtered" if filt else "overflow" if full else "distance"

        parts_ms, returned_pairs = [0.0, 0.0, 0.0], 0
        for k in range(K):
            positions = vertices + (k / K) * (end - vertices)
            ms, lists, grow = detect_sized(positions, faces, step,
                                           1 if args.substep else args.repeats, not args.substep)
            parts_ms = [a + b for a, b in zip(parts_ms, ms)]
            setup_ms += grow
            returned_pairs += int(np.minimum(lists["vt"][3], lists["vt"][2]).sum()
                                  + np.minimum(lists["ee"][3], lists["ee"][2]).sum())
            if args.validate and k == 0:
                cv, wv = validate(positions, faces, lists["vt"], scene.radius)
                ce, we = validate_ee(positions, scene.edges, lists["ee"])
                print(f"validate step {step}: radius {scene.radius:.6g}; {cv} vertex-triangle pairs, largest "
                      f"distance {wv:.6g}; {ce} edge-edge pairs, largest distance {we:.6g}; K={K}",
                      file=sys.stderr)
            for i in by_sub.get(k, []) + by_sub.get(-1, []):
                if found[i]:
                    continue
                ok, why = check(lists, i)
                if ok:
                    found[i], cause[i] = True, None
                elif cause[i] is None or CAUSE_RANK[why] > CAUSE_RANK[cause[i]]:
                    # A pair without a root is checked in every sub-step and is
                    # far apart in most; a full buffer or the filter in any one
                    # of them is what kept it out.
                    cause[i] = why

        for ci, (case, is_vf, n_queries, n_colliding) in enumerate(case_info):
            miss = dict(filtered=0, overflow=0, distance=0)
            for i, pr in enumerate(pending):
                if pr[0] == ci and not found[i]:
                    miss[cause[i] or "distance"] += 1
            missed = sum(miss.values())
            writer.writerow(dict(
                dataset=args.scene, case=case.name, type="vf" if is_vf else "ee", queries=n_queries,
                colliding=n_colliding, returned_pairs=returned_pairs, missed=missed,
                missed_filtered=miss["filtered"], missed_overflow=miss["overflow"],
                missed_distance=miss["distance"], detect_ms=f"{sum(parts_ms):.4f}",
                refit_ms=f"{parts_ms[0]:.4f}", vt_ms=f"{parts_ms[1]:.4f}", ee_ms=f"{parts_ms[2]:.4f}",
                setup_ms=f"{setup_ms:.2f}",
                radius=scene.radius, scene_diag=f"{diag:.6g}",
                vt_buffer=scene.vt_buffer, ee_buffer=scene.ee_buffer, substeps=K,
                edge_epsilon=f"{scene.edge_epsilon:.6g}"))
            totals["colliding"] += n_colliding
            totals["missed"] += missed
            for key in miss:
                totals[key] += miss[key]
            setup_ms = 0.0
        out.flush()

    print(f"{args.scene}: {totals['missed']} of {totals['colliding']} colliding pairs missed "
          f"(filtered {totals['filtered']}, overflow {totals['overflow']}, distance {totals['distance']})",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
