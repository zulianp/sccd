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

SWEEP = REPO / "benchmark" / "results" / "sweep-gh200-bp.csv"
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
            tag = {"total": "BP full + NP", "prep": "BP prep",
                   "broad": "BP queries", "narrow": "NP"}[phase]
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
        if r["mode"] not in (HOST, DEV) or (r.get("broadphase") or "cell2d") != "cell2d":
            continue
        fr = frame_of(r["case"])
        if fr is None:
            continue
        try:
            acc[r["dataset"]][r["mode"]][fr] += (
                float(r["prep_ms"]) + float(r["broad_ms"]) + float(r["narrow_ms"]))
        except ValueError:
            continue

    fig, axes = plt.subplots(1, len(SCENES), figsize=style.figsize(
        style.FULL_WIDTH_IN, 0.30))
    for ax, scene in zip(axes, SCENES):
        for mode, lab in ((HOST, "CPU"), (DEV, "GPU")):
            d = acc[scene][mode]
            if not d:
                continue
            xs = sorted(d)
            ax.plot(xs, [d[x] for x in xs], lw=0.7,
                    color=style.MODE_COLOR[mode], label=lab)
        log_y(ax)
        ax.set_title(LABEL.get(scene, scene), fontsize=7)
        ax.set_xlabel("frame", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("BP full + NP, ms per step", fontsize=7)

    check_ticks(fig, axes, SCENES, "per-frame")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout()
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
SERIES_BP = (("cell list, CPU", HOST, ("cell2d",), 0.0),
             ("sweep, CPU", HOST, ("sweep",), SWEEP_TINT),
             ("cell list, GPU", DEV, ("cell2d",), 0.0),
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
            acc[r["dataset"]][(r["mode"], r.get("broadphase") or "cell2d")][fr] += (
                float(r["prep_ms"]) + float(r["broad_ms"]))
        except ValueError:
            continue

    fig, axes = plt.subplots(1, len(SCENES), figsize=style.figsize(
        style.FULL_WIDTH_IN, 0.30))
    for ax, scene in zip(axes, SCENES):
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
        ax.set_xlabel("frame", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("BP full, ms per step", fontsize=7)

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
               loc="lower center", bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
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
SERIES_VS = (("cell list (ours)", "device-tight", "cell2d", DEV, 0.0),
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

    fig, axes = plt.subplots(1, len(present), figsize=style.figsize(
        style.FULL_WIDTH_IN, 0.30))
    axes = [axes] if len(present) == 1 else list(axes)
    for ax, scene in zip(axes, present):
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
        ax.set_xlabel("frame", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("BP full, ms per step", fontsize=7)

    check_ticks(fig, axes, present, "broad-vs-scalable")
    for ax, scene in zip(axes, present):
        drawn = len(ax.get_legend_handles_labels()[0])
        if drawn != len(SERIES_VS):
            raise SystemExit(
                f"broad-vs-scalable: {scene} draws {drawn} of {len(SERIES_VS)} series")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=3,
               loc="lower center", bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    p = OUT / "broad-vs-scalable.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


def main() -> int:
    style.apply_rcparams()
    OUT.mkdir(parents=True, exist_ok=True)
    for fn in (strong_scaling, per_frame, broad_per_frame, broad_vs_scalable):
        p = fn()
        if p:
            print(f"wrote {p.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
