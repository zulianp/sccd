#!/usr/bin/env python3
"""
Figures the library's own report does not draw, from the same committed data.

    <python-with-matplotlib> scripts/make_figures.py

Three plots, each answering a question the reference study of Belgrod et al.
answers for their pipeline and ours did not answer for this one:

  strong-scaling   speedup of each host phase against thread count, with the
                   perfect line, from benchmark/results/profile/strong-*.csv
  per-frame        cost through a simulation rather than aggregated over it

Written into figures/ as PDF. They are committed, so building the article needs
no Python; this is only for regenerating them after new data.
"""

from __future__ import annotations

import collections
import csv
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
            ax.plot(threads, [one[phase] / best[t][phase] for t in threads],
                    marker="os^D"[i], ms=3.4, lw=1.2,
                    color=style.SERIES[i], label=phase, zorder=2)
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
        ax.set_yscale("log")
        # A decade-only locator leaves a panel spanning less than a decade with
        # one tick or none, which happened to four of these six. Put majors at
        # 1, 2 and 5 times each power of ten and label them as plain numbers.
        ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(
            base=10.0, subs=(1.0, 2.0, 5.0), numticks=12))
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_title(LABEL.get(scene, scene), fontsize=7)
        ax.set_xlabel("frame", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, which="major", lw=0.4, color=style.GRID_INK)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("ms per step", fontsize=7)

    # A log panel spanning less than a decade can end up with no labelled tick
    # at all, which is how this figure first shipped. Fail rather than publish
    # an axis a reader cannot read.
    fig.canvas.draw()
    for ax, scene in zip(axes, SCENES):
        labelled = [t for t in ax.get_yticklabels()
                    if t.get_text() and ax.get_ylim()[0] <= t.get_position()[1]
                    <= ax.get_ylim()[1]]
        if len(labelled) < 2:
            raise SystemExit(
                f"per-frame: {scene} has {len(labelled)} y tick labels; "
                "the locator is not covering its range")

    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout()
    p = OUT / "per-frame.pdf"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    return p


def main() -> int:
    style.apply_rcparams()
    OUT.mkdir(parents=True, exist_ok=True)
    for fn in (strong_scaling, per_frame):
        p = fn()
        if p:
            print(f"wrote {p.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
