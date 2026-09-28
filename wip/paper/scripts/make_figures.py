#!/usr/bin/env python3
"""
Figures the library's own report does not draw, from the same committed data.

    <python-with-matplotlib> scripts/make_figures.py

Three plots, each answering a question the reference study of Belgrod et al.
answers for their pipeline and ours did not answer for this one:

  strong-scaling   speedup of each host phase against thread count, with the
                   perfect line, from benchmark/results/profile/strong-*.csv
  per-frame        cost through a simulation rather than aggregated over it
  broad-per-frame  the same for the broad phase alone, with both strategies
                   on both processors
  broad-vs-scalable  the device broad phase against Scalable CCD's, per frame

Written into figures/ as PDF. They are committed, so building the article needs
no Python; this is only for regenerating them after new data.
"""

from __future__ import annotations

import collections
import csv
import statistics as st
import re
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
REPO = PAPER.parent.parent
sys.path.insert(0, str(REPO / "benchmark"))

from report import style  # noqa: E402  (needs the path above)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SWEEP = REPO / "benchmark" / "assessment" / "broadphase-cell2dmin.csv"
# The strategy the library runs by default, and the one every figure that is
# not explicitly comparing strategies should draw. Named once: it was hard-coded
# in four places, and changing three of them left one figure silently empty.
SHIPPED_BP = "cell2dmin"
COMPARE = REPO / "benchmark" / "competitors" / "results"
PROF = REPO / "benchmark" / "results" / "profile"
OUT = PAPER / "figures"

SCENES = ["armadillo-rollers", "cloth-ball", "cloth-funnel",
          "n-body-simulation", "puffer-ball", "rod-twist"]
LABEL = {s: style.SCENE_LABEL.get(s, s) for s in SCENES}
HOST, DEV = "tight", "device-tight"


def sweep_rows():
    with SWEEP.open() as f:
        for r in csv.DictReader(f):
            yield r


def frame_of(case: str):
    m = re.match(r"(\d+)", case or "")
    return int(m.group(1)) if m else None


def lighten(hex_color, amount):
    """The same hue at a lower intensity, blended that far towards white.

    Intensity carries the broad-phase strategy where a dash pattern would break
    up curves that already change direction every frame.
    """
    h = hex_color.lstrip("#")
    rgb = [int(h[i:i + 2], 16) for i in (0, 2, 4)]
    return "#" + "".join(f"{round(c + (255 - c) * amount):02x}" for c in rgb)


def log_y(ax):
    """A log y axis a reader can read on a panel spanning under a decade.

    A decade-only locator leaves such a panel with one tick or none, which
    happened to four of these six. Put majors at 1, 2 and 5 times each power of
    ten and label them as plain numbers.
    """
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(
        base=10.0, subs=(1.0, 2.0, 5.0), numticks=12))
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def check_ticks(fig, axes, scenes, name):
    """Fail rather than publish an axis with no labelled tick.

    A log panel spanning less than a decade can end up with no labelled tick at
    all, which is how the per-frame figure first shipped.
    """
    fig.canvas.draw()
    for ax, scene in zip(axes, scenes):
        labelled = [t for t in ax.get_yticklabels()
                    if t.get_text() and ax.get_ylim()[0] <= t.get_position()[1]
                    <= ax.get_ylim()[1]]
        if len(labelled) < 2:
            raise SystemExit(
                f"{name}: {scene} has {len(labelled)} y tick labels; "
                "the locator is not covering its range")



def scene_grid(rows=2, cols=3, height=0.62):
    """A panel per scene on `rows` x `cols`, which is what six scenes want.

    One row of six leaves each panel about an inch wide, and these panels carry
    a few thousand frames of a noisy curve. Two rows of three give each of them
    twice the width and enough height to separate curves that sit within a
    factor of two of each other.
    """
    fig, axes = plt.subplots(rows, cols, figsize=style.figsize(
        style.FULL_WIDTH_IN, height), squeeze=False)
    return fig, [a for row in axes for a in row]


# ---------------------------------------------------------------- scaling ---
def strong_scaling():
    """Speedup per phase against thread count, against the perfect line."""
    files = sorted(PROF.glob("strong-*.csv"))
    if not files:
        print("no strong-*.csv; skipping strong scaling", file=sys.stderr)
        return None
    fig, axes = plt.subplots(1, len(files), figsize=style.figsize(
        style.FULL_WIDTH_IN, 0.42), sharey=True)
    axes = [axes] if len(files) == 1 else list(axes)

    for ax, path in zip(axes, files):
        scene = path.stem.replace("strong-", "")
        # threads -> phase -> best total over repeats
        acc = collections.defaultdict(lambda: collections.defaultdict(
            lambda: collections.defaultdict(float)))
        for line in path.read_text().splitlines():
            f = line.split(",")
            if len(f) < 11:
                continue
            try:
                t, rep = int(f[0]), int(f[1])
                # threads,rep then the benchmark's own row. The broadphase
                # column sits between mode and case, so prep is field 8.
                prep, broad, narrow = float(f[8]), float(f[9]), float(f[10])
            except ValueError:
                continue
            acc[t][rep]["prep"] += prep
            acc[t][rep]["broad"] += broad
            acc[t][rep]["narrow"] += narrow
            acc[t][rep]["total"] += prep + broad + narrow

        threads = sorted(acc)
        best = {t: {p: min(acc[t][r][p] for r in acc[t]) for p in
                    ("prep", "broad", "narrow", "total")} for t in threads}
        one = best[threads[0]]

        ax.plot(threads, threads, ls="--", lw=0.9, color=style.REFERENCE_INK,
                label="perfect", zorder=1)
        for i, phase in enumerate(("total", "prep", "broad", "narrow")):
            # The tags the tables and captions use, so one vocabulary covers both.
            tag = {"total": "BP full + NP EToI", "prep": "BP prep",
                   "broad": "BP queries", "narrow": "NP EToI"}[phase]
            ax.plot(threads, [one[phase] / best[t][phase] for t in threads],
                    marker="os^D"[i], ms=3.4, lw=1.2,
                    color=style.SERIES[i], label=tag, zorder=2)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        shown = [t for t in threads if t in (1, 2, 4, 8, 16, 32, 72)]
        ax.set_xticks(shown)
        ax.set_xticklabels([str(t) for t in shown], fontsize=6)
        ax.set_xticks(threads, minor=True)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_title(LABEL.get(scene, scene))
        ax.set_xlabel("threads")
        ax.grid(True, which="both", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("speedup")
    axes[-1].legend(frameon=False, fontsize=6, loc="upper left")
    fig.tight_layout()
    p = OUT / "strong-scaling.pdf"
    fig.savefig(p); plt.close(fig)
    return p


# -------------------------------------------------------------- per frame ---
def per_frame():
    """Cost through a simulation, which every aggregate in the paper hides."""
    # scene -> mode -> frame -> summed ms over the vf and ee case of that frame
    acc = collections.defaultdict(lambda: collections.defaultdict(
        lambda: collections.defaultdict(float)))
    for r in sweep_rows():
        if r["mode"] not in (HOST, DEV) or (r.get("broadphase") or SHIPPED_BP) != SHIPPED_BP:
            continue
        fr = frame_of(r["case"])
        if fr is None:
            continue
        try:
            acc[r["dataset"]][r["mode"]][fr] += (
                float(r["prep_ms"]) + float(r["broad_ms"]) + float(r["narrow_ms"]))
        except ValueError:
            continue

    fig, axes = scene_grid()
    for i, (ax, scene) in enumerate(zip(axes, SCENES)):
        drawn = 0
        for mode, lab in ((HOST, "CPU"), (DEV, "GPU")):
            d = acc[scene][mode]
            if not d:
                continue
            xs = sorted(d)
            ax.plot(xs, [d[x] for x in xs], lw=0.7,
                    color=style.MODE_COLOR[mode], label=lab)
            drawn += 1
        # An empty panel is a figure that compiles, publishes and says nothing.
        # This one shipped blank when the broad-phase filter above stopped
        # matching any row, so it is checked here rather than by eye.
        if not drawn:
            raise SystemExit(
                f"per-frame: {scene} drew no series; no row of {SWEEP.name} "
                f"has broadphase={SHIPPED_BP!r} on {HOST} or {DEV}")
        log_y(ax)
        ax.set_title(LABEL.get(scene, scene), fontsize=7)
        if i >= 3:
            ax.set_xlabel("frame", fontsize=7)
        if i % 3 == 0:
            ax.set_ylabel("BP full + NP EToI, ms per step", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)

    check_ticks(fig, axes, SCENES, "per-frame")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    p = OUT / "per-frame.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


# -------------------------------------------------- broad phase per frame ---
# The three series: colour carries the processor, dash carries the strategy.
# The device broad phase does not implement the choice -- it always sorts, and
# `sccd::device::cell2d_*` has no caller outside its own unit test -- so its two
# settings are one code path measured twice and are pooled as repeats rather
# than drawn as a comparison that does not exist.
# Colour carries the processor and intensity the strategy: the cell list at full
# strength, the sweep at the same hue lightened.
SWEEP_TINT = 0.45
SERIES_BP = (("cell list, CPU", HOST, (SHIPPED_BP,), 0.0),
             ("sweep, CPU", HOST, ("sweep",), SWEEP_TINT),
             ("cell list, GPU", DEV, (SHIPPED_BP,), 0.0),
             ("sweep, GPU", DEV, ("sweep",), SWEEP_TINT))


def broad_per_frame():
    """Broad-phase cost through a simulation, both strategies on both processors.

    The broad phase here is the whole of it, the acceleration structure and the
    traversal over it, summed over the vertex-face and edge-edge case of the
    frame: what one step of a simulation pays before any root finding.
    """
    # scene -> (mode, strategy) -> frame -> summed ms over the vf and ee case
    acc = collections.defaultdict(lambda: collections.defaultdict(
        lambda: collections.defaultdict(float)))
    for r in sweep_rows():
        if r["mode"] not in (HOST, DEV):
            continue
        fr = frame_of(r["case"])
        if fr is None:
            continue
        try:
            acc[r["dataset"]][(r["mode"], r.get("broadphase") or SHIPPED_BP)][fr] += (
                float(r["prep_ms"]) + float(r["broad_ms"]))
        except ValueError:
            continue

    fig, axes = scene_grid()
    for i, (ax, scene) in enumerate(zip(axes, SCENES)):
        for lab, mode, strategies, tint in SERIES_BP:
            got = [acc[scene][(mode, s)] for s in strategies
                   if acc[scene][(mode, s)]]
            if not got:
                continue
            frames = sorted(set().union(*(set(d) for d in got)))
            ys = [st.median([d[f] for d in got if f in d]) for f in frames]
            ax.plot(frames, ys, lw=0.7,
                    color=lighten(style.MODE_COLOR[mode], tint), label=lab)
        log_y(ax)
        ax.set_title(LABEL.get(scene, scene), fontsize=7)
        # The frame axis is only labelled on the bottom row, and the value axis
        # only on the left column; every panel keeps its own ticks, because the
        # scenes do not share a range.
        if i >= 3:
            ax.set_xlabel("frame", fontsize=7)
        if i % 3 == 0:
            ax.set_ylabel("BP full, ms per step", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)

    check_ticks(fig, axes, SCENES, "broad-per-frame")

    # Every panel must carry all three series, or the figure claims a
    # comparison it does not draw.
    for ax, scene in zip(axes, SCENES):
        drawn = len(ax.get_legend_handles_labels()[0])
        if drawn != len(SERIES_BP):
            raise SystemExit(
                f"broad-per-frame: {scene} draws {drawn} of {len(SERIES_BP)} series")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    p = OUT / "broad-per-frame.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


# ------------------------------------------------ device against the competitor ---
def compare_rows():
    """Every per-case row of the newest competitor comparison, gzipped or not."""
    import gzip
    cands = sorted(COMPARE.glob("compare-gh200-full-*.csv.gz")) + \
        sorted(COMPARE.glob("compare-gh200-full-*.csv"))
    if not cands:
        return None
    path = max(cands, key=lambda p: p.name.replace(".csv.gz", "").replace(".csv", ""))
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return [r for r in csv.DictReader(fh) if r.get("type") in ("vf", "ee")]


# Our two strategies and theirs. Colour separates the library, intensity the
# strategy, matching broad-per-frame.
SERIES_VS = (("cell list (ours)", "device-tight", SHIPPED_BP, DEV, 0.0),
             ("sweep (ours)", "device-tight", "sweep", DEV, SWEEP_TINT),
             ("Scalable CCD", "scalable-ccd-device", None, "relaxed", 0.0))


def broad_vs_scalable():
    """The device broad phase against Scalable CCD's, step by step.

    All three come from one comparison run, because between allocations this
    harness varies by about 40% -- far more than the differences drawn here.
    """
    data = compare_rows()
    if not data:
        print("no comparison CSV; skipping broad-vs-scalable", file=sys.stderr)
        return None

    # scene -> series -> frame -> summed ms over the vf and ee case of the frame
    acc = collections.defaultdict(lambda: collections.defaultdict(
        lambda: collections.defaultdict(list)))
    for r in data:
        fr = frame_of(r["case"])
        if fr is None:
            continue
        for lab, mode, strategy, _, _ in SERIES_VS:
            if r["mode"] != mode:
                continue
            if strategy is not None and (r.get("broadphase") or "") != strategy:
                continue
            try:
                acc[r["dataset"]][lab][fr].append(
                    float(r["prep_ms"]) + float(r["broad_ms"]))
            except ValueError:
                pass

    present = [s for s in SCENES if acc[s]]
    if not present:
        print("comparison CSV has no usable rows; skipping", file=sys.stderr)
        return None

    fig, axes = scene_grid() if len(present) == 6 else (
        lambda f, a: (f, list(a) if len(present) > 1 else [a]))(
            *plt.subplots(1, len(present),
                          figsize=style.figsize(style.FULL_WIDTH_IN, 0.30)))
    for i, (ax, scene) in enumerate(zip(axes, present)):
        for lab, _, _, colour_mode, tint in SERIES_VS:
            d = acc[scene][lab]
            if not d:
                continue
            frames = sorted(d)
            # The two cases of a frame are summed per repeat, then the repeats
            # reduced, so one slow pass cannot move the curve.
            ys = [st.median(d[f]) for f in frames]
            ax.plot(frames, ys, lw=0.7,
                    color=lighten(style.MODE_COLOR[colour_mode], tint), label=lab)
        log_y(ax)
        ax.set_title(LABEL.get(scene, scene), fontsize=7)
        if len(present) != 6 or i >= 3:
            ax.set_xlabel("frame", fontsize=7)
        if len(present) != 6 or i % 3 == 0:
            ax.set_ylabel("BP full, ms per step", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)

    check_ticks(fig, axes, present, "broad-vs-scalable")
    for ax, scene in zip(axes, present):
        drawn = len(ax.get_legend_handles_labels()[0])
        if drawn != len(SERIES_VS):
            raise SystemExit(
                f"broad-vs-scalable: {scene} draws {drawn} of {len(SERIES_VS)} series")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=3,
               loc="lower center", bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    p = OUT / "broad-vs-scalable.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


COUNTER_TEX = PAPER / "generated" / "tables"


def _tex_rows(name):
    """The body rows of a generated table, as lists of cells."""
    f = COUNTER_TEX / f"{name}.tex"
    if not f.is_file():
        return []
    body = f.read_text()
    try:
        body = body[body.index("\\midrule"):body.index("\\bottomrule")]
    except ValueError:
        return []
    out = []
    for line in body.splitlines():
        if "&" not in line:
            continue
        out.append([c.strip() for c in line.replace("\\\\", "").split("&")])
    return out


# ------------------------------------------------ competitors, per scene ---
def competitors():
    """Both competitor comparisons as whole-scene bars, one panel each.

    Read from the two generated tables and not recomputed, so the figure and the
    appendix tables it summarises cannot disagree. Those tables are 19 rows
    apiece; this is the form the body of the paper reads.
    """
    PANELS = (
        ("tab-competitor-earliest", 6, "earliest time of impact, per case (ms)",
         ("SCCD (CPU)", "SCCD (GPU)", "Scalable CCD (GPU)"),
         ("SCCD (CPU)", "SCCD (GPU)", "Scalable CCD"), "log"),
        # Linear: this panel spans well under two decades, where a log axis
        # costs readability and buys nothing.
        ("tab-competitor-pair", 7, "per pair, narrow phase (ns/pair)",
         ("SCCD (CPU)", "SCCD (GPU)", "Additive CCD (CPU)"),
         ("SCCD (CPU)", "SCCD (GPU)", "Additive CCD"), "linear"),
    )
    panels = []
    for name, col, ylab, keys, labels, scale in PANELS:
        rows, scene, d = _tex_rows(name), None, collections.defaultdict(dict)
        for r in rows:
            if len(r) <= col:
                continue
            scene = r[0] or scene
            try:
                d[scene][r[1]] = float(r[col].replace(",", ""))
            except ValueError:
                pass
        if not d:
            print(f"no {name}; skipping competitors", file=sys.stderr)
            return None
        panels.append((ylab, list(d), keys, labels, d, scale))

    fig, axes = plt.subplots(1, 2, figsize=style.figsize(style.FULL_WIDTH_IN, 0.34))
    for ax, (ylab, names, keys, labels, d, scale) in zip(axes, panels):
        x, w = range(len(names)), 0.8 / len(keys)
        for i, (k, lab) in enumerate(zip(keys, labels)):
            ax.bar([j + i * w - 0.4 + w / 2 for j in x],
                   [d[n].get(k, float("nan")) for n in names], w,
                   color=style.SERIES[i], label=lab, zorder=2)
        ax.set_yscale(scale)
        ax.set_ylabel(ylab, fontsize=7)
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=6)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", axis="y", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=6)
        drawn = sum(1 for n in names for k in keys if k in d[n])
        if drawn != len(names) * len(keys):
            raise SystemExit(
                f"competitors: {ylab} drew {drawn} of {len(names) * len(keys)} bars")
    fig.tight_layout()
    p = OUT / "competitors.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


# ------------------------------------------------- device counters ---
def device_counters():
    """Pipe utilisation and warp stall reasons, the two counter tables as bars."""
    pipes, stalls = _tex_rows("tab-pipes"), _tex_rows("tab-counters")
    if not pipes or not stalls:
        print("no counter tables; skipping device-counters", file=sys.stderr)
        return None

    def col(rows, idx):
        out = {}
        for r in rows:
            try:
                out[r[0]] = float(r[idx])
            except (ValueError, IndexError):
                pass
        return out

    PANELS = (
        (r"pipe, \% of peak while active",
         [("FP64", col(pipes, 1)), ("ALU", col(pipes, 2)), ("LSU", col(pipes, 3))]),
        ("warps stalled, cycles per issued instruction",
         [("dependency", col(stalls, 7)), ("barrier", col(stalls, 8)),
          ("memory", col(stalls, 9))]),
    )
    scenes = [s for s in pipes]
    names = [r[0] for r in pipes]
    fig, axes = plt.subplots(1, 2, figsize=style.figsize(style.FULL_WIDTH_IN, 0.34))
    x = range(len(names))
    for ax, (ylab, series) in zip(axes, PANELS):
        w = 0.8 / len(series)
        for i, (lab, d) in enumerate(series):
            ax.bar([j + i * w - 0.4 + w / 2 for j in x],
                   [d.get(n, float("nan")) for n in names], w,
                   color=style.SERIES[i], label=lab, zorder=2)
        ax.set_ylabel(ylab, fontsize=7)
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=6)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", axis="y", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    p = OUT / "device-counters.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


# ------------------------------------------------ earliness distributions ---
def earliness():
    """Per-case earliness, SCCD above and additive CCD below, same bins.

    Drawn from the competitor run, so both rows describe the same cases measured
    in one allocation; that is what makes the two rows comparable at all. The
    quantity is the largest earliness in a case, one-sided by construction.
    """
    import numpy as np

    data = compare_rows()
    if not data:
        print("no comparison CSV; skipping earliness", file=sys.stderr)
        return None

    lo, hi = 1e-9, 1.0
    edges = np.logspace(np.log10(lo), np.log10(hi), 26)
    ROWS = (("SCCD", (("CPU", HOST, "-"), ("GPU", DEV, "--"))),
            ("Additive CCD", (("Additive CCD (CPU)", "accd", "-"),)))

    # scene -> mode -> one value per case
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in data:
        if r["mode"].startswith("device-") and r.get("broadphase") != SHIPPED_BP:
            continue
        try:
            v = float(r["toi_max_early"])
        except (ValueError, KeyError):
            continue
        if v > 0:
            acc[r["dataset"]][r["mode"]].append(min(max(v, lo), hi))

    present = [s_ for s_ in SCENES if acc[s_]]
    if not present:
        print("comparison CSV has no earliness; skipping", file=sys.stderr)
        return None

    fig, axes = plt.subplots(len(ROWS), len(present), squeeze=False, sharex=True,
                             figsize=style.figsize(style.FULL_WIDTH_IN, 0.40))
    for row, (rowlab, series) in enumerate(ROWS):
        for col, scene in enumerate(present):
            ax = axes[row][col]
            for i, (lab, mode, ls) in enumerate(series):
                vals = acc[scene].get(mode)
                if not vals:
                    continue
                ax.hist(vals, bins=edges, histtype="step", linewidth=1.0,
                        linestyle=ls, color=style.SERIES[i if row == 0 else 2],
                        label=lab)
            ax.set_xscale("log")
            # Two decade ticks over nine decades leave the peaks unplaceable,
            # which is the whole point of putting the two rows on one axis.
            ax.set_xticks([1e-9, 1e-6, 1e-3, 1e0])
            ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            ax.set_xlim(lo, hi)
            ax.tick_params(labelsize=6, length=2, pad=1)
            ax.grid(True, linewidth=0.3, color=style.GRID_INK)
            ax.set_axisbelow(True)
            if row == 0:
                ax.set_title(LABEL.get(scene, scene), fontsize=7, pad=3)
            if row == len(ROWS) - 1:
                ax.set_xlabel("earliness", fontsize=6.5)
            if col == 0:
                ax.set_ylabel(f"{rowlab}\ncases", fontsize=6.5)

    handles, labels = [], []
    for row in range(len(ROWS)):
        h, l = axes[row][0].get_legend_handles_labels()
        handles += h; labels += l
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), fontsize=6.5,
               frameon=False, handletextpad=0.4, columnspacing=1.4,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    p = OUT / "earliness.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


def main() -> int:
    style.apply_rcparams()
    OUT.mkdir(parents=True, exist_ok=True)
    for fn in (strong_scaling, per_frame, broad_per_frame, broad_vs_scalable,
               competitors, device_counters, earliness):
        p = fn()
        if p:
            print(f"wrote {p.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
