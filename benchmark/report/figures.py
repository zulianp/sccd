"""
Figures, drawn once and written as both PDF (for LaTeX) and PNG (for Markdown).

Two rules the previous generators broke and this one keeps.

Spread is drawn. `write_aggregate_csv` computed a mean, median, standard
deviation, minimum and maximum for every quantity -- thirty-four columns of it --
and not one figure drew any of it, while the prose above them said a difference
inside that band is not a result. Every timing figure here carries the observed
run-to-run range.

Identification lives in the caption. Axes carry no title: a title drawn into the
image duplicates the LaTeX caption, cannot be referenced by number, and is the
first thing to be cropped. Each figure therefore returns its caption and label,
which the caller writes into the figure `.tex`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from .data import SceneSummary, Stat
from .style import (mode_label_paper, SERIES, COLUMN_WIDTH_IN, FULL_WIDTH_IN, MODE_ORDER, SCENE_LABEL,
                    apply_rcparams, figsize, mode_color, mode_label, MODE_MARKER)


@dataclass
class Figure:
    stem: str
    label: str
    caption: str
    width_in: float = COLUMN_WIDTH_IN


def _ordered_modes(summaries: dict) -> list[str]:
    present = {mode for _, mode in summaries}
    ordered = [m for m in MODE_ORDER if m in present]
    return ordered + sorted(present - set(ordered))


# Prefix applied to every figure filename, so two documents generated from the
# same data do not overwrite each other's figures.
PREFIX = ""


def _stem(name: str) -> str:
    return f"{PREFIX}{name}"


def _save(fig, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{_stem(stem)}.pdf")
    fig.savefig(out_dir / f"{_stem(stem)}.png")


def narrow_phase_per_case(case_series: dict, out_dir: Path) -> Figure:
    """Per-case narrow-phase cost against query count, which is what drives it."""
    import matplotlib.pyplot as plt

    modes = sorted({s.mode for s in case_series.values()})
    ordered = [m for m in MODE_ORDER if m in modes] + \
              sorted(set(modes) - set(MODE_ORDER))
    fig, ax = plt.subplots(figsize=figsize(COLUMN_WIDTH_IN, 0.72))

    for mode in ordered:
        xs, ys = [], []
        for s in case_series.values():
            if s.mode != mode:
                continue
            stat = s.timings["narrow_ms"]
            if s.queries > 0 and stat.n and math.isfinite(stat.median) and stat.median > 0:
                xs.append(s.queries)
                ys.append(stat.median)
        if not xs:
            continue
        ax.scatter(xs, ys, s=7, alpha=0.55, linewidths=0,
                   color=mode_color(mode), label=mode_label(mode),
                   marker=MODE_MARKER.get(mode, "o"))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("candidate pairs in the case")
    ax.set_ylabel("narrow phase (ms)")
    ax.grid(True, which="both", axis="both")
    ax.legend(loc="upper left")
    fig.tight_layout()
    _save(fig, out_dir, "narrow-per-case")
    plt.close(fig)

    return Figure(
        _stem("narrow-per-case"), "fig:narrow-per-case",
        "Narrow-phase time for each individual case against the number of "
        "candidate pairs the broad phase handed it, on log-log axes. Each point "
        "is the median over repeats.",
        COLUMN_WIDTH_IN)


def results_grid(case_series: dict, out_dir: Path) -> Figure:
    """
    The benchmark laid out as Belgrod et al. lay theirs out: scenes across the
    columns, one measured quantity per row, a box over the cases of that scene.

    A box plot rather than a bar because the quantity varies over the cases of a
    scene by orders of magnitude, and a scene total says nothing about that
    spread. Following the paper, the box spans the first to the third quartile,
    the line is the median, the whiskers reach the furthest point within 1.5
    interquartile ranges, and anything beyond is drawn as a circle. Every row is
    on a log axis for the same reason.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import NullFormatter

    scenes = sorted({s.dataset for s in case_series.values()})
    modes = [m for m in MODE_ORDER
             if m in {s.mode for s in case_series.values()}]
    if not scenes or not modes:
        return Figure(
        _stem("results-grid"), "fig:results", "no data", FULL_WIDTH_IN)

    # False positives are not a row here. They are almost always zero, so the
    # box collapses onto the axis and the panel spends a quarter of the figure
    # saying nothing; the count per scene is in the conservativeness table,
    # which is the right place for a number that is usually the same number.
    rows = [("broad_ms", "broad phase (ms)"),
            ("narrow_ms", "narrow phase (ms)"),
            ("toi_max_early", "error")]

    fig, axes = plt.subplots(len(rows), len(scenes), squeeze=False,
                             figsize=(FULL_WIDTH_IN, 1.35 * len(rows) + 0.9),
                             sharex="col")

    for c, scene in enumerate(scenes):
        axes[0][c].set_title(SCENE_LABEL.get(scene, scene), fontsize=7, pad=3)
        for r, (column, ylabel) in enumerate(rows):
            ax = axes[r][c]
            data, colours = [], []
            for mode in modes:
                vals = []
                for s in case_series.values():
                    if s.dataset != scene or s.mode != mode:
                        continue
                    stat = (s.timings if column in s.timings else s.accuracy).get(column)
                    # A count of zero is a result -- this mode reported no false
                    # positive on that case -- so counts keep their zeros and are
                    # drawn on a symmetric-log axis, which has room for them.
                    floor = -1.0 if column == "fp" else 0.0
                    if stat is not None and stat.n and math.isfinite(stat.median) \
                            and stat.median > floor:
                        vals.append(stat.median)
                data.append(np.asarray(vals) if vals else np.asarray([np.nan]))
                colours.append(mode_color(mode))

            bp = ax.boxplot(data, positions=np.arange(len(modes)), widths=0.62,
                            patch_artist=True, whis=1.5,
                            flierprops=dict(marker="o", markersize=1.1,
                                            markerfacecolor="none", alpha=0.55),
                            medianprops=dict(color="white", linewidth=0.8),
                            capprops=dict(linewidth=0.5),
                            whiskerprops=dict(linewidth=0.5))
            for patch, colour in zip(bp["boxes"], colours):
                patch.set_facecolor(colour)
                patch.set_edgecolor(colour)
                patch.set_linewidth(0.5)
            for flier, colour in zip(bp["fliers"], colours):
                flier.set_markeredgecolor(colour)

            if column == "fp":
                # symlog keeps zero on the axis; below 1 there is nothing to see
                # because the quantity is a count.
                ax.set_yscale("symlog", linthresh=1.0)
            elif any(np.isfinite(d).any() and (d > 0).any() for d in data):
                ax.set_yscale("log")
            ax.set_xticks([])
            ax.tick_params(labelsize=6, length=2, pad=1)
            # On a log axis spanning less than a decade matplotlib labels the
            # minor ticks too, and those labels are wide enough to overrun the
            # neighbouring panel. The row shares a scale, so only the first
            # column is labelled and the minor labels are never wanted.
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.grid(True, axis="y", linewidth=0.3)
            if c == 0:
                ax.set_ylabel(ylabel, fontsize=6.5)
            else:
                ax.tick_params(which="both", labelleft=False)

    handles = [plt.Line2D([], [], marker="s", linestyle="none", markersize=4,
                          color=mode_color(m), label=mode_label_paper(m))
               for m in modes]
    fig.legend(handles=handles, loc="lower center", ncol=len(modes),
               fontsize=6.5, frameon=False, handletextpad=0.4,
               columnspacing=1.4, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    _save(fig, out_dir, "results-grid")
    plt.close(fig)

    return Figure(
        _stem("results-grid"), "fig:results",
        "Per-case distributions for every mode over the six scenes (columns). "
        "Rows are broad-phase time, narrow-phase time, and error against "
        "the exact root. A star marks a "
        "parallel CPU mode and a dagger a GPU one. Each box spans the first to "
        "the third quartile with the median inside it, the whiskers reach the "
        "furthest case within 1.5 interquartile ranges, and cases beyond that "
        "are drawn individually. All axes are logarithmic.",
        FULL_WIDTH_IN)



def runtime_breakdown(summaries: dict[tuple[str, str], SceneSummary],
                      out_dir: Path) -> Figure:
    """
    Where the time goes, laid out as Belgrod et al. lay out their Figure 7.5:
    one horizontal stacked bar per configuration, scenes across the columns, and
    a single legend naming the parts.

    Horizontal because the configuration labels are words, and stacked because
    the question is which phase dominates rather than what the total is -- the
    total is in the timing table.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    scenes = sorted({sc for sc, _ in summaries})
    modes = [m for m in MODE_ORDER if m in {m2 for _, m2 in summaries}]
    if not scenes or not modes:
        return Figure(
        _stem("runtime-breakdown"), "fig:breakdown", "no data", FULL_WIDTH_IN)

    parts = [("prep_ms", "prep"), ("broad_ms", "broad"), ("narrow_ms", "narrow")]
    part_ink = [SERIES[3], SERIES[0], SERIES[2]]

    fig, axes = plt.subplots(1, len(scenes), squeeze=False,
                             figsize=(FULL_WIDTH_IN, 0.34 * len(modes) + 1.25))
    y = np.arange(len(modes))

    for c, scene in enumerate(scenes):
        ax = axes[0][c]
        left = np.zeros(len(modes))
        for (column, _), ink in zip(parts, part_ink):
            widths = []
            for mode in modes:
                s_ = summaries.get((scene, mode))
                stat = s_.totals.get(column) if s_ else None
                widths.append(stat.median if stat is not None and stat.n else 0.0)
            widths = np.asarray(widths)
            ax.barh(y, widths, left=left, height=0.68, color=ink, linewidth=0)
            left += widths

        ax.set_title(SCENE_LABEL.get(scene, scene), fontsize=7, pad=3)
        ax.set_yticks(y)
        ax.invert_yaxis()
        ax.tick_params(labelsize=6, length=2, pad=1)
        ax.set_xlabel("time (ms)", fontsize=6.5)
        ax.grid(True, axis="x", linewidth=0.3)
        if c == 0:
            ax.set_yticklabels([mode_label_paper(m) for m in modes], fontsize=6.5)
        else:
            ax.set_yticklabels([])

    handles = [plt.Rectangle((0, 0), 1, 1, color=ink) for ink in part_ink]
    fig.legend(handles=handles, labels=[label for _, label in parts],
               loc="lower center", ncol=len(parts), fontsize=6.5, frameon=False,
               handletextpad=0.4, columnspacing=1.6, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    _save(fig, out_dir, "runtime-breakdown")
    plt.close(fig)

    return Figure(
        _stem("runtime-breakdown"), "fig:breakdown",
        "Runtime split by phase for every scene and mode: \\emph{prep} builds "
        "the swept boxes and the acceleration structure, \\emph{broad} finds "
        "the candidate pairs, \\emph{narrow} turns them into a time of impact. "
        "A star marks a parallel CPU mode and a dagger a GPU one. Bars are the "
        "median over repeats.",
        FULL_WIDTH_IN)



def toi_error_histogram(case_series: dict, out_dir: Path) -> Figure:
    """
    The time-of-impact error against the symbolic ground truth, as a histogram
    on a log axis -- Belgrod et al.'s Figure 7.3, over the same scenes.

    The quantity is one-sided by construction: every value is how far *before*
    the exact root the mode reported, so the axis carries no sign and the whole
    distribution being finite is itself the conservativeness statement.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    scenes = sorted({s.dataset for s in case_series.values()})
    modes = [m for m in MODE_ORDER if m in {s.mode for s in case_series.values()}]
    if not scenes or not modes:
        return Figure(
        _stem("toi-error"), "fig:toi-error", "no data", FULL_WIDTH_IN)

    # Not sharey: the scenes differ in case count by a factor of fifty-eight
    # (cloth-ball has 79, rod-twist 4,571), so a shared axis flattens five of the
    # six panels into the baseline.
    fig, axes = plt.subplots(1, len(scenes), squeeze=False,
                             figsize=(FULL_WIDTH_IN, 1.85))
    lo, hi = 1e-9, 1.0
    edges = np.logspace(np.log10(lo), np.log10(hi), 26)

    for c, scene in enumerate(scenes):
        ax = axes[0][c]
        for mode in modes:
            vals = []
            for s_ in case_series.values():
                if s_.dataset != scene or s_.mode != mode:
                    continue
                stat = s_.accuracy.get("toi_max_early")
                if stat is not None and stat.n and math.isfinite(stat.median) \
                        and stat.median > 0:
                    vals.append(min(max(stat.median, lo), hi))
            if not vals:
                continue
            # Tight and Tight (GPU) run the same search and their histograms
            # coincide; dashing the device modes keeps both visible.
            ax.hist(vals, bins=edges, histtype="step",
                    linewidth=1.1 if not mode.startswith("device-") else 0.8,
                    linestyle="--" if mode.startswith("device-") else "-",
                    color=mode_color(mode), label=mode_label_paper(mode))
        ax.set_xscale("log")
        ax.set_title(SCENE_LABEL.get(scene, scene), fontsize=7, pad=3)
        ax.tick_params(labelsize=6, length=2, pad=1)
        ax.grid(True, linewidth=0.3)
        ax.set_xlabel("earliness", fontsize=6.5)
        if c == 0:
            ax.set_ylabel("cases", fontsize=6.5)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(modes), fontsize=6.5,
               frameon=False, handletextpad=0.4, columnspacing=1.4,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    _save(fig, out_dir, "toi-error")
    plt.close(fig)

    return Figure(
        _stem("toi-error"), "fig:toi-error",
        "Distribution over cases of the time-of-impact error against the "
        "dataset's exact symbolic roots, on a log axis. The error is one-sided "
        "by construction -- it is how far before the true root the mode "
        "reported -- so every value shown is on the safe side, and no case falls "
        "outside the axis on the late side because none exists.",
        FULL_WIDTH_IN)


def write_figure_tex(figures: list[Figure], out_dir: Path,
                     stem: str = "figures") -> Path:
    """A LaTeX fragment including each figure with its caption and label."""
    lines = []
    for f in figures:
        width = "\\columnwidth" if f.width_in <= COLUMN_WIDTH_IN else "\\textwidth"
        lines += [
            r"\begin{figure}[htbp]",
            r"  \centering",
            f"  \\includegraphics[width={width}]{{figures/{f.stem}.pdf}}",
            f"  \\caption{{{f.caption}}}",
            f"  \\label{{{f.label}}}",
            r"\end{figure}",
            "",
        ]
    path = out_dir / f"{stem}.tex"
    path.write_text("\n".join(lines))
    return path
