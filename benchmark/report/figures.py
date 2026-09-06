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
from .style import (COLUMN_WIDTH_IN, FULL_WIDTH_IN, MODE_ORDER, SCENE_LABEL,
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


def _save(fig, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf")
    fig.savefig(out_dir / f"{stem}.png")


def phase_breakdown(summaries: dict[tuple[str, str], SceneSummary],
                    out_dir: Path) -> Figure:
    """Broad and narrow phase side by side, per scene and mode, with the range."""
    import matplotlib.pyplot as plt
    import numpy as np

    scenes = sorted({scene for scene, _ in summaries})
    modes = _ordered_modes(summaries)
    fig, ax = plt.subplots(figsize=figsize(FULL_WIDTH_IN, 0.50))

    # Cap the group width so a sweep of one scene does not draw two bars the
    # width of the axes; with several scenes this is the usual 0.8.
    group_width = min(0.8, 0.28 * max(len(modes), 1))
    bar_width = group_width / max(len(modes), 1)
    x = np.arange(len(scenes), dtype=float)

    for i, mode in enumerate(modes):
        offset = -group_width / 2 + bar_width * (i + 0.5)
        broad_med, narrow_med, err_lo, err_hi = [], [], [], []
        for scene in scenes:
            s = summaries.get((scene, mode))
            b = s.totals["broad_ms"] if s else Stat()
            n = s.totals["narrow_ms"] if s else Stat()
            broad_med.append(b.median if b.n else 0.0)
            narrow_med.append(n.median if n else 0.0)
            # The range of the total, which is the bar's full height.
            total = Stat()
            for k in range(min(b.n, n.n)):
                total.add(b.values[k] + n.values[k])
            centre = total.median if total.n else 0.0
            err_lo.append(max(centre - total.lo, 0.0) if total.n else 0.0)
            err_hi.append(max(total.hi - centre, 0.0) if total.n else 0.0)

        colour = mode_color(mode)
        ax.bar(x + offset, broad_med, bar_width * 0.92, color=colour, alpha=0.45,
               label=f"{mode_label(mode)} — broad", linewidth=0)
        ax.bar(x + offset, narrow_med, bar_width * 0.92, bottom=broad_med,
               color=colour, label=f"{mode_label(mode)} — narrow", linewidth=0)
        totals = [b + n for b, n in zip(broad_med, narrow_med)]
        ax.errorbar(x + offset, totals, yerr=[err_lo, err_hi], fmt="none",
                    ecolor="#2A2F35", elinewidth=0.8, capsize=2.0)

    ax.set_xticks(x)
    ax.set_xticklabels([SCENE_LABEL.get(s, s) for s in scenes])
    ax.set_ylabel("time for the whole scene (ms)")
    ax.set_xlim(-0.5, len(scenes) - 0.5)
    # The legend sits above the axes rather than inside them: with the tallest
    # bar near the top of the frame there is no interior corner it does not
    # cover, and covering the data to label it is not a trade worth making.
    ax.legend(ncol=min(len(modes), 3), loc="lower center",
              bbox_to_anchor=(0.5, 1.01), borderaxespad=0.0)
    fig.tight_layout()
    _save(fig, out_dir, "phase-breakdown")
    plt.close(fig)

    return Figure(
        "phase-breakdown", "fig:phase-breakdown",
        "Broad phase (pale) and narrow phase (solid) summed over every case in "
        "each scene, per narrow-phase mode. Bars are the median over independent "
        "repeats; whiskers span the full run-to-run range of the total, so a "
        "difference smaller than a whisker is not a result.",
        FULL_WIDTH_IN)


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
        "narrow-per-case", "fig:narrow-per-case",
        "Narrow-phase time for each individual case against the number of "
        "candidate pairs the broad phase handed it, on log-log axes. Each point "
        "is the median over repeats.",
        COLUMN_WIDTH_IN)


def earliness_distribution(case_series: dict, out_dir: Path) -> Figure:
    """
    How far before the true root each mode reports.

    Drawn on a log axis, always: the quantity spans several decades, and the two
    modes differ by two of them. A previous pair of figures drew the same
    quantity on a log axis in one place and a linear one in another.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    modes = sorted({s.mode for s in case_series.values()})
    ordered = [m for m in MODE_ORDER if m in modes] + \
              sorted(set(modes) - set(MODE_ORDER))
    fig, ax = plt.subplots(figsize=figsize(COLUMN_WIDTH_IN, 0.62))

    drew = False
    for mode in ordered:
        values = []
        for s in case_series.values():
            if s.mode != mode:
                continue
            stat = s.accuracy["toi_med_early"]
            if stat.n and math.isfinite(stat.median) and stat.median > 0:
                values.append(stat.median)
        if not values:
            continue
        drew = True
        values = np.sort(np.asarray(values))
        # An empirical CDF says more than a histogram here: the reader wants
        # "what fraction of cases report within X of the truth", and it does not
        # depend on a bin width nobody chose deliberately.
        ax.step(values, np.arange(1, len(values) + 1) / len(values),
                where="post", color=mode_color(mode), label=mode_label(mode))

    if not drew:
        ax.text(0.5, 0.5, "no case reported a positive earliness",
                ha="center", va="center", transform=ax.transAxes)
    else:
        ax.set_xscale("log")
    ax.set_xlabel("earliness of the reported time of impact")
    ax.set_ylabel("fraction of cases")
    ax.set_ylim(0, 1.02)
    ax.grid(True, which="both", axis="both")
    ax.legend(loc="lower right")
    fig.tight_layout()
    _save(fig, out_dir, "earliness-cdf")
    plt.close(fig)

    return Figure(
        "earliness-cdf", "fig:earliness",
        "Empirical distribution over cases of how far before the true time of "
        "impact each mode reports. Earliness is plotted on a log axis, so a "
        "curve further left reports closer to the true root and a curve "
        "further right reports earlier -- always the safe direction, and "
        "always at the cost of a solver's step size.",
        COLUMN_WIDTH_IN)


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
