"""
The refinement scaling study: how cost responds to element count.

`sccd_refine_scaling` answers what neither other driver can -- it refines one
surface repeatedly, quadrupling the triangle count at each level, and runs a
collision step on each -- and docs/BENCHMARKS.md did not mention it once. Its
output is a fixed-width table rather than a CSV, so it is parsed here instead of
being rewritten: the driver's format is the one every recorded run on disk is
already in.

    # mode=relaxed max_depth=69 tol=3e-08 scale=1 space=host base_topology=TRISHELL3
    level      faces     vf_pairs     ee_pairs   prep_ms  bp_fv_ms ... narrow_ms  ns/pair  toi
        0       1536         2841         4210      0.42      0.31 ...      1.94    123.4  0.5
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from .style import (COLUMN_WIDTH_IN, SERIES, apply_rcparams, figsize)
from .tables import TAGS

_HEADER_RE = re.compile(r"^\s*level\s+faces\s+")
_META_RE = re.compile(r"(\w+)=(\S+)")


@dataclass
class ScalingRun:
    meta: dict[str, str] = field(default_factory=dict)
    levels: list[int] = field(default_factory=list)
    faces: list[int] = field(default_factory=list)
    vf_pairs: list[int] = field(default_factory=list)
    ee_pairs: list[int] = field(default_factory=list)
    prep_ms: list[float] = field(default_factory=list)
    step_ms: list[float] = field(default_factory=list)
    broad_ms: list[float] = field(default_factory=list)
    narrow_ms: list[float] = field(default_factory=list)
    ns_per_pair: list[float] = field(default_factory=list)

    @property
    def label(self) -> str:
        mode = self.meta.get("mode", "?")
        space = self.meta.get("space", "host")
        topology = self.meta.get("base_topology", "")
        kind = "quad" if "QUAD" in topology.upper() else "tri"
        # The broad phase is part of the identity of the run, not a footnote:
        # the two strategies differ by 1.8x at the largest size measured.
        bp = self.meta.get("broadphase")
        base = f"{mode} / {space} / {kind}"
        return f"{base} / {bp}" if bp else base


def parse(path: Path) -> ScalingRun:
    run = ScalingRun()
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            run.meta.update(dict(_META_RE.findall(line)))
            continue
        if _HEADER_RE.match(line):
            continue
        fields = line.split()
        # level faces vf_pairs ee_pairs prep bp_fv bp_ee broad narrow ns/pair toi
        if len(fields) < 11:
            continue
        try:
            run.levels.append(int(fields[0]))
            run.faces.append(int(fields[1]))
            run.vf_pairs.append(int(fields[2]))
            run.ee_pairs.append(int(fields[3]))
            run.prep_ms.append(float(fields[4]))
            run.step_ms.append(float(fields[5]) + float(fields[6]))
            run.broad_ms.append(float(fields[7]))
            run.narrow_ms.append(float(fields[8]))
            run.ns_per_pair.append(float(fields[9]))
        except ValueError:
            continue
    return run


def fitted_exponent(x: list[float], y: list[float]) -> float:
    """
    Least-squares slope of log y against log x: the empirical exponent p in
    y ~ x^p. Stated rather than asserted -- a claim that cost is linear in the
    element count should come with the number that was measured.
    """
    pairs = [(math.log(a), math.log(b))
             for a, b in zip(x, y) if a > 0 and b > 0]
    if len(pairs) < 2:
        return math.nan
    n = len(pairs)
    sx = sum(a for a, _ in pairs)
    sy = sum(b for _, b in pairs)
    sxx = sum(a * a for a, _ in pairs)
    sxy = sum(a * b for a, b in pairs)
    denom = n * sxx - sx * sx
    if denom == 0:
        return math.nan
    return (n * sxy - sx * sy) / denom


def _processor(run: ScalingRun) -> str:
    return "GPU" if run.meta.get("space", "host") == "device" else "CPU"


def _strategy(run: ScalingRun) -> str:
    """
    The broad phase a run actually ran, or nothing when it did not choose one.

    A run made without SCCD_BROADPHASE records `broadphase=auto`, and `auto` is
    not a strategy: it is the tuner, which races the two and keeps the winner.
    Printing it beside `sweep` and `cell2d` puts a mechanism in a list of
    algorithms, and on this driver it is worse than uninformative -- the driver
    builds one CCD object per refinement level, so each level starts a fresh
    tuner that has not raced anything yet. Such a series is left unnamed here
    rather than given a label it did not earn.
    """
    bp = run.meta.get("broadphase", "")
    return "" if bp == "auto" else bp


def series_labels(runs: list[ScalingRun]) -> list[str]:
    """
    Name each run by what distinguishes it from the others, and nothing else.

    Every series in one figure or table shares most of its configuration --
    the same mode, the same topology -- and repeating those on every entry is
    noise around the one axis the reader is being shown. So a field is printed
    only where the runs disagree about it.
    """
    def fields(r: ScalingRun) -> tuple[str, str, str, str]:
        return (r.meta.get("mode", "?"), _processor(r),
                "quad" if "QUAD" in r.meta.get("base_topology", "").upper() else "tri",
                _strategy(r))

    varying = [i for i in range(4) if len({fields(r)[i] for r in runs}) > 1]
    out = []
    for r in runs:
        f = fields(r)
        parts = [f[i] for i in (varying or [0, 1, 2, 3]) if f[i]]
        out.append(" / ".join(parts) if parts else r.label)
    return out


def table(runs: list[ScalingRun], source: str):
    """Cost per refinement level, with the fitted exponent per series."""
    from .tables import Column, Table

    labels = series_labels(runs)

    def _fields0(r):
        return (r.meta.get("mode", "?"), _processor(r),
                "quad" if "QUAD" in r.meta.get("base_topology", "").upper() else "tri",
                _strategy(r))
    _varying0 = [i for i in range(4) if len({_fields0(r)[i] for r in runs}) > 1]
    _first_header = {0: "mode", 1: "processor", 2: "topology", 3: "broad phase"}.get(
        _varying0[0], "series") if len(_varying0) == 1 else "series"

    t = Table(
        label="tab:scaling",
        caption=("Cost against element count on a repeatedly refined surface, "
                 "each level quadrupling the element count. Every time is one "
                 "collision step at that level, in milliseconds. " + TAGS +
                 "BP prep here is the acceleration structure and BP queries the "
                 "two overlap queries over it. "
                 "$p$ is the least-squares exponent in $t \\sim n^{p}$ fitted "
                 "over all levels of that series."),
        columns=[Column(_first_header, "l"), Column("level"), Column("elements"),
                 Column("candidate pairs"),
                 Column("structure ms", tex_header="BP prep (ms)"),
                 Column("traversal ms", tex_header="BP queries (ms)"),
                 Column("broad ms", tex_header="BP full (ms)"),
                 Column("narrow ms", tex_header="NP EToI (ms)"),
                 Column("p")],
        source=source,
        notes=("The two frames used here do not come into contact, so the narrow "
               "phase has almost no work to do and its column is dominated by "
               "noise rather than by element count; what this measures is the "
               "broad phase and the preparation that feeds it. Narrow-phase cost "
               "against problem size is in the per-case figure, over cases that "
               "do collide. Where the exponent is below 1 it is because the fixed "
               "cost visible at the smallest size is amortised as the mesh "
               "grows. `prep` builds the acceleration structure -- the cell "
               "list's grid or the sweep's sorted intervals -- and `step` is "
               "the traversal that reports pairs; the two strategies divide "
               "the work between those columns quite differently."),
    )
    for run, label in zip(runs, labels):
        faces = [float(f) for f in run.faces]
        totals = [b + n for b, n in zip(run.broad_ms, run.narrow_ms)]
        p_fit = fitted_exponent(faces, totals)
        for i, level in enumerate(run.levels):
            t.add(label if i == 0 else "", f"{level}", f"{run.faces[i]:,}",
                  f"{run.vf_pairs[i] + run.ee_pairs[i]:,}",
                  f"{run.prep_ms[i]:.1f}", f"{run.step_ms[i]:.1f}",
                  f"{run.broad_ms[i]:.1f}", f"{run.narrow_ms[i]:.1f}",
                  f"{p_fit:.2f}" if i == 0 else "")
    return t


def figure(runs: list[ScalingRun], out_dir: Path):
    """Cost against element count, with the fitted exponent on each series."""
    import matplotlib.pyplot as plt

    from .figures import Figure

    apply_rcparams()
    fig, ax = plt.subplots(figsize=figsize(COLUMN_WIDTH_IN, 0.70))

    labels = series_labels(runs)
    for i, (run, label) in enumerate(zip(runs, labels)):
        if not run.faces:
            continue
        total = [b + n for b, n in zip(run.broad_ms, run.narrow_ms)]
        p = fitted_exponent([float(f) for f in run.faces], total)
        colour = SERIES[i % len(SERIES)]
        ax.plot(run.faces, total, marker="o", color=colour,
                label=f"{label} ($p={p:.2f}$)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("surface elements")
    ax.set_ylabel("BP full + NP EToI (ms)")
    ax.grid(True, which="both", axis="both")
    ax.legend(loc="upper left")
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "refine-scaling.pdf")
    fig.savefig(out_dir / "refine-scaling.png")
    plt.close(fig)

    return Figure(
        "refine-scaling", "fig:scaling",
        "Cost against element count on a repeatedly refined surface, each level "
        "quadrupling the element count. Both axes are logarithmic, so the slope "
        "is the exponent $p$ in $t \\sim n^{p}$; the fitted value is given in "
        "the legend of each series.",
        COLUMN_WIDTH_IN)
