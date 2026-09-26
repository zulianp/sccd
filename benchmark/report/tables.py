"""
Tables, emitted twice from one source: LaTeX `booktabs` and Markdown.

The two must not drift, so each table is built once as a `Table` -- a caption, a
label, column specs and rows of already-formatted cells -- and rendered into
either target. Every table names the CSV it came from, so a number in the
article can be traced to the run that produced it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from .data import SceneSummary, Stat, separable
from .style import SCENE_LABEL, mode_label

# The phase tags every caption uses, so "broad" never has to be guessed
# at: it meant the traversal alone in one table and the whole phase in
# another, and a reader had no way to tell which.
TAGS = ("\\emph{BP full} is the whole broad phase, the acceleration structure "
        "plus the queries over it; \\emph{BP prep} is the structure alone and "
        "\\emph{BP queries} the queries alone. The narrow phase is named by the "
        "answer it is asked for: \\emph{NP EToI} returns one earliest time of "
        "impact for the step, so every query prunes against the running minimum, "
        "and \\emph{NP per-pair} returns one per candidate with no shared bound. ")


@dataclass
class Column:
    header: str
    align: str = "r"          # l, c or r
    tex_header: str | None = None

    def tex(self) -> str:
        return self.tex_header if self.tex_header is not None else _tex_escape(self.header)


@dataclass
class Table:
    label: str
    caption: str
    columns: list[Column]
    rows: list[list[str]] = field(default_factory=list)
    source: str = ""
    notes: str = ""

    def add(self, *cells: str) -> None:
        if len(cells) != len(self.columns):
            raise ValueError(
                f"{self.label}: row has {len(cells)} cells, table has "
                f"{len(self.columns)} columns")
        self.rows.append([str(c) for c in cells])


def _tex_escape(text: str) -> str:
    out = str(text)
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
                 ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")):
        out = out.replace(a, b)
    return out


def render_tex(table: Table) -> str:
    spec = "".join(c.align for c in table.columns)
    lines = [
        r"\begin{table}[htbp]",
        r"  \centering",
        f"  \\caption{{{_tex_escape(table.caption)}}}",
        f"  \\label{{{table.label}}}",
        f"  \\begin{{tabular}}{{{spec}}}",
        r"    \toprule",
        "    " + " & ".join(c.tex() for c in table.columns) + r" \\",
        r"    \midrule",
    ]
    for row in table.rows:
        lines.append("    " + " & ".join(_tex_escape(c) for c in row) + r" \\")
    lines += [r"    \bottomrule", r"  \end{tabular}"]
    if table.notes:
        lines.append(f"  \\par\\smallskip\\footnotesize {_tex_escape(table.notes)}")
    if table.source:
        lines.append(f"  \\par\\smallskip\\footnotesize Source: "
                     f"\\texttt{{{_tex_escape(table.source)}}}")
    lines.append(r"\end{table}")
    return "\n".join(lines) + "\n"


def render_markdown(table: Table) -> str:
    widths = [len(c.header) for c in table.columns]
    for row in table.rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def line(cells: list[str]) -> str:
        out = []
        for i, cell in enumerate(cells):
            out.append(cell.rjust(widths[i]) if table.columns[i].align == "r"
                       else cell.ljust(widths[i]))
        return "| " + " | ".join(out) + " |"

    rule = "|" + "|".join(
        ("-" * (widths[i] + 2)) if c.align != "r"
        else ("-" * (widths[i] + 1)) + ":"
        for i, c in enumerate(table.columns)) + "|"

    parts = [line([c.header for c in table.columns]), rule]
    parts += [line(row) for row in table.rows]
    body = "\n".join(parts)
    tail = ""
    if table.notes:
        tail += f"\n\n{table.notes}"
    if table.source:
        tail += f"\n\nSource: `{table.source}`"
    return body + tail + "\n"


def _ms(stat: Stat) -> str:
    """
    Median and slowest repeat, as `median / max`.

    A median alone describes the run you are likely to get and says
    nothing about the run you have to budget for. The two together also
    show the spread directly, which a percentage states but does not let
    the reader check.
    """
    if not stat.n or not math.isfinite(stat.median):
        return "--"
    if stat.n == 1:
        return f"{stat.median:.1f}"
    return f"{stat.median:.1f} / {stat.hi:.1f}"


def dataset_table(summaries: dict[tuple[str, str], SceneSummary],
                  oracle_rows: dict, source: str) -> Table:
    """
    What was measured on, before any result about it.

    A reader's first question is the size and shape of the problems, and it
    should be answerable without reading prose. Cases and candidate pairs come
    from the sweep; queries and roots from the accuracy run, which sees every
    query rather than the subsample a timing run might use.
    """
    table = Table(
        label="tab:dataset",
        caption=("The benchmark problems. \\emph{cases} is simulation steps with "
                 "a runnable query set, \\emph{candidate pairs} the mean number "
                 "the broad phase produces per step, and \\emph{queries} the "
                 "curated per-step query sets that carry exact symbolic roots. Cases and "
                 "queries are counts over the whole scene; candidate pairs is a "
                 "per-step mean."),
        columns=[Column("scene", "l"), Column("cases"),
                 Column("candidate pairs/step", tex_header=r"pairs/step"),
                 Column("queries"), Column("with a root", tex_header=r"w/ root")],
        source=source,
    )
    per_scene: dict[str, SceneSummary] = {}
    for (scene, _), s in summaries.items():
        # any mode will do for case and pair counts; they do not depend on it
        per_scene.setdefault(scene, s)

    # Per mode, then take one mode: every mode sees the same queries, so summing
    # across modes and dividing back by how many there were only invites the
    # count to drift when a scene is swept on one processor and not another.
    by_mode: dict[tuple[str, str], tuple[int, int]] = {}
    for (scene, _, mode), r in (oracle_rows or {}).items():
        if mode == "tight-inclusion":
            continue
        q, g = by_mode.get((scene, mode), (0, 0))
        by_mode[(scene, mode)] = (q + r.queries, g + r.gt_checked)

    roots: dict[str, tuple[int, int]] = {}
    for (scene, _mode), qg in by_mode.items():
        roots[scene] = max(roots.get(scene, (0, 0)), qg)

    for scene in sorted(per_scene):
        s = per_scene[scene]
        pairs = s.queries // s.cases if s.cases else 0
        q, g = roots.get(scene, (0, 0))
        table.add(SCENE_LABEL.get(scene, scene), f"{s.cases:,}", f"{pairs:,}",
                  f"{q:,}" if q else "--",
                  f"{g:,}" if g else "--")
    return table


def throughput_table(summaries: dict[tuple[str, str], SceneSummary], source: str) -> Table:
    """
    Candidate pairs per second, which is what compares across scenes.

    Milliseconds for a whole scene answer "how long did this take"; they cannot
    be compared between a 79-case scene and a 4,571-case one. Throughput can.
    """
    table = Table(
        label="tab:throughput",
        caption=("Broad- and narrow-phase throughput in candidate pairs per "
                 "second, median over repeats. A rate, so it is directly "
                 "comparable between scenes of very different size: a rate, so it is "
                 "neither a whole-scene total nor a per-step figure. " + TAGS),
        columns=[Column("scene", "l"), Column("mode", "l"),
                 Column("broad Mpair/s", tex_header=r"BP full (Mpair/s)"),
                 Column("narrow Mpair/s", tex_header=r"NP EToI (Mpair/s)")],
        source=source,
    )
    for (scene, mode), s in sorted(summaries.items()):
        # Pairs per second of the whole broad phase, structure included.
        prep_s, trav = s.totals["prep_ms"], s.totals["broad_ms"]
        b = Stat()
        for i in range(min(prep_s.n, trav.n)):
            b.add(prep_s.values[i] + trav.values[i])
        n = s.totals["narrow_ms"]
        def rate(stat):
            m = stat.median
            if not stat.n or not math.isfinite(m) or m <= 0 or not s.queries:
                return "--"
            return f"{s.queries / (m / 1000.0) / 1e6:.1f}"
        table.add(SCENE_LABEL.get(scene, scene), mode_label(mode), rate(b), rate(n))
    return table


def timing_table(summaries: dict[tuple[str, str], SceneSummary], source: str) -> Table:
    table = Table(
        label="tab:timing",
        caption=("Wall-clock time per scene and narrow-phase mode, as a "
                 "whole-scene total in milliseconds summed over every case in "
                 "the scene, given as median / slowest over independent "
                 "repeats. " + TAGS +
                 "\\emph{total} is BP full + NP EToI."),
        columns=[
            Column("scene", "l"), Column("mode", "l"), Column("cases"),
            Column("pairs"), Column("rep"),
            Column("broad ms", tex_header=r"BP full (ms)"),
            Column("earliest ms", tex_header=r"NP EToI (ms)"),
            Column("per-pair ms", tex_header=r"NP per-pair (ms)"),
            Column("total ms", tex_header=r"total (ms)"),
        ],
        source=source,
        notes=("Each cell is the median over repeats and the slowest of them. A difference smaller than the gap between the two does not separate "
               "two modes and is not reported as a ratio anywhere in this "
               "document. Which mode is faster depends on the output mode as "
               "well as the scene, so the two are given side by side rather "
               "than one standing for the other."),
    )
    for (scene, mode), s in sorted(summaries.items()):
        prep = s.totals["prep_ms"]
        traversal, narrow = s.totals["broad_ms"], s.totals["narrow_ms"]
        per_pair = s.totals["narrow_ms_s1"]
        # The acceleration structure belongs to the broad phase, so add it in
        # within a repeat and only then across repeats. On rod-twist it is
        # larger than the traversal and the narrow phase together, so a broad
        # phase reported without it understates the pipeline by more than half.
        broad, total = Stat(), Stat()
        for i in range(min(prep.n, traversal.n, narrow.n)):
            broad.add(prep.values[i] + traversal.values[i])
            total.add(prep.values[i] + traversal.values[i] + narrow.values[i])
        table.add(SCENE_LABEL.get(scene, scene), mode_label(mode),
                  f"{s.cases}", f"{s.queries:,}", f"{s.repeats}",
                  _ms(broad), _ms(narrow), _ms(per_pair), _ms(total))
    return table


def conservativeness_table(summaries: dict[tuple[str, str], SceneSummary],
                           source: str) -> Table:
    table = Table(
        label="tab:conservativeness",
        caption=("Conservativeness against the dataset's exact roots, as counts over "
                 "the whole scene summed across its cases. "
                 "\\emph{late} counts queries whose reported time of impact "
                 "falls after the true one; it must be zero, because a late "
                 "time of impact lets a simulation step through the contact, "
                 "which is the failure the search exists to prevent. False "
                 "positives cost work only and are reported for information. "
                 "\\emph{queries} is how many carry ground-truth data at all, "
                 "no-collision cases included; \\emph{toi compared} is how many "
                 "times of impact were actually placed beside an exact root, "
                 "which is what the claim rests on. One row per scene: the two "
                 "processors agree except on false positives, where the host's "
                 "count is given first and the device's in parentheses."),
        columns=[
            Column("scene", "l"),
            Column("queries", tex_header=r"queries"),
            Column("toi compared", tex_header=r"toi compared"),
            Column("late"), Column("false pos."), Column("false neg."),
        ],
        source=source,
        notes=("Measured against the exact roots shipped with the dataset, not "
               "against TightInclusion: TightInclusion's own answer is itself a "
               "lower bound on the truth, so comparing against it over-reports "
               "lateness."),
    )
    # One row per scene. The processors agree on every column but the false
    # positives, and there only on three scenes, so a row each would repeat
    # itself twelve times to show four differing numbers. Where they differ the
    # host's count is given with the device's in parentheses.
    def cell(values: list[str]) -> str:
        first = values[0]
        rest = [v for v in values[1:] if v != first]
        return first if not rest else f"{first} ({', '.join(rest)})"

    by_scene: dict[str, list[tuple[str, SceneSummary]]] = {}
    for (scene, mode), s in sorted(summaries.items()):
        by_scene.setdefault(scene, []).append((mode, s))

    for scene, entries in by_scene.items():
        # Host first, so the parenthesised value is always the device's.
        entries.sort(key=lambda e: e[0].startswith("device-"))
        cols = [[f"{s.gt_queries:,}" for _, s in entries],
                [f"{s.toi_compared:,}" for _, s in entries],
                [f"{s.toi_late}" for _, s in entries],
                [f"{s.fp:,}" for _, s in entries],
                [f"{s.fn}" for _, s in entries]]
        table.add(SCENE_LABEL.get(scene, scene), *(cell(c) for c in cols))
    return table


def accuracy_table(summaries: dict[tuple[str, str], SceneSummary], source: str) -> Table:
    table = Table(
        label="tab:earliness",
        caption=("How far before the true time of impact each mode reports. "
                 "\\emph{median} is the median over cases of the per-case "
                 "median; \\emph{worst case} is the largest earliness reported "
                 "anywhere in the scene. Reporting early is always safe and "
                 "always costs a solver step size, so the worst case is the "
                 "largest step the mode can cost, not a tail to discount."),
        columns=[Column("scene", "l"), Column("mode", "l"),
                 Column("median earliness", tex_header=r"median"),
                 Column("worst case", tex_header=r"worst case")],
        source=source,
    )
    for (scene, mode), s in sorted(summaries.items()):
        if not s.toi_med_early:
            value = "--"
        else:
            ordered = sorted(s.toi_med_early)
            mid = len(ordered) // 2
            med = (ordered[mid] if len(ordered) % 2
                   else 0.5 * (ordered[mid - 1] + ordered[mid]))
            value = f"{med:.2e}"
        worst = f"{max(s.toi_max_early):.2e}" if s.toi_max_early else "--"
        table.add(SCENE_LABEL.get(scene, scene), mode_label(mode), value, worst)
    return table


def comparison_notes(summaries: dict[tuple[str, str], SceneSummary]) -> list[str]:
    """
    Mode-versus-mode statements, each one checked against the noise first.

    This is the discipline docs/BENCHMARKS.md opens with and that none of the
    generated figures honoured: a gap inside the observed run-to-run spread is
    written down as inside noise, not as a ratio.
    """
    notes: list[str] = []
    # Group by processor as well as scene. Ranking every mode of a scene together
    # picks the best and worst overall, which on a scene measured on both
    # processors compares host Relaxed against GPU Tight and reports the sum of
    # two unrelated effects as if it were the mode trade.
    groups: dict[tuple[str, str], dict[str, SceneSummary]] = {}
    for (scene, mode), summary in summaries.items():
        space = "GPU" if mode.startswith("device-") else "CPU"
        groups.setdefault((scene, space), {})[mode] = summary

    # Both output modes, because which narrow-phase mode is faster depends on
    # it. Tight's advantages exist only where a shared running minimum lets a
    # tighter bound prune the queries that follow; with one result per candidate
    # and no shared bound there is nothing for that tightness to buy.
    for (scene, space) in sorted(groups):
        for column, output in (("narrow_ms", "earliest"), ("narrow_ms_s1", "per-pair")):
            _compare(groups[(scene, space)], scene, space, column, output, notes)
    return notes


def _compare(modes, scene, space, column, output, notes) -> None:
        if len(modes) < 2:
            return
        where = f"{SCENE_LABEL.get(scene, scene)} ({space}, {output})"
        ranked = sorted(modes.items(), key=lambda kv: kv[1].totals[column].median)
        (best_mode, best), (worst_mode, worst) = ranked[0], ranked[-1]
        b, w = best.totals[column], worst.totals[column]
        if b.n < 2 or w.n < 2:
            notes.append(
                f"- **{where}**: a single repeat gives no estimate of the noise, "
                f"so {mode_label(best_mode)} and {mode_label(worst_mode)} are not "
                f"separable here.")
            return
        ok, ratio, noise = separable(b, w)
        if ok:
            notes.append(
                f"- **{where}**: {mode_label(best_mode)} is {ratio:.2f}× faster "
                f"than {mode_label(worst_mode)} in the narrow phase "
                f"({b.median:.0f} ms against {w.median:.0f} ms; run-to-run spread "
                f"{noise * 100:.1f}%).")
        else:
            notes.append(
                f"- **{where}**: {mode_label(best_mode)} and "
                f"{mode_label(worst_mode)} are inside noise "
                f"({b.median:.0f} ms against {w.median:.0f} ms, spread "
                f"{noise * 100:.1f}%); this does not separate them.")


def write_tables(tables: list[Table], out_dir: Path, stem: str = "tables") -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    tex_path = out_dir / f"{stem}.tex"
    md_path = out_dir / f"{stem}.md"
    tex_path.write_text("\n".join(render_tex(t) for t in tables))
    md_path.write_text("\n\n".join(
        f"### {t.caption}\n\n{render_markdown(t)}" for t in tables))
    return {"tex": tex_path, "md": md_path}


def broadphase_table(per_strategy: dict[str, dict[tuple[str, str], SceneSummary]],
                     source: str) -> Table:
    """
    The broad-phase strategies over the same scenes.

    They return identical pair sets, so whichever is faster on a given geometry
    is simply the right one. The host fixes that choice from this measurement and
    the device races for it per scene, having no minimum-corner implementation to
    choose between. The preparation is split out from the traversal because that
    is where the strategies differ: the sweep builds its sorted intervals more
    cheaply, the cell list traverses its grid more cheaply, and ordering a cell's
    entries buys traversal at the cost of preparation.

    A strategy missing from a processor drops that processor's rows, so pass only
    strategies measured on both.
    """
    names = sorted(per_strategy)
    table = Table(
        label="tab:broadphase",
        caption=("Broad-phase strategies over the same cases, "
                 "as whole-scene totals in milliseconds, "
                 "median over repeats. " + TAGS +
                 "The two columns per strategy are BP prep, building the sorted "
                 "intervals or the grid, and BP queries over it; BP full is the "
                 "two added, which is what the \\emph{faster} verdict ranks. "
                 "Both strategies report identical "
                 "candidate pairs, so the difference is entirely in how they "
                 "are found."),
        columns=([Column("scene", "l"), Column("mode", "l")]
                 + [Column(f"{n} structure ms", tex_header=f"{n} BP prep") for n in names]
                 + [Column(f"{n} broad ms", tex_header=f"{n} BP queries") for n in names]
                 + [Column("faster", "l")]),
        source=source,
        notes=("`faster` names the winning strategy and by how much on the whole "
               "broad phase. A margin inside the run-to-run spread is a tie."),
    )
    # Both processors. `use_cell2d_` is read by the prep and by both steps on
    # each of them, so SCCD_BROADPHASE selects a real implementation either way
    # and a GPU row compares two of them rather than one against itself.
    keys = sorted({k for s in per_strategy.values() for k in s})
    for scene, mode in keys:
        structure, broad, spreads = {}, {}, {}
        for n in names:
            s = per_strategy[n].get((scene, mode))
            structure[n] = s.totals["prep_ms"] if s else Stat()
            traversal = s.totals["broad_ms"] if s else Stat()
            # The strategies are ranked on the whole broad phase, structure
            # included. Ranking on the traversal alone picks the sweep on most
            # scenes and the cell list on most whole phases, because the sort
            # the sweep saves in traversal it pays for in structure.
            whole = Stat()
            for i in range(min(structure[n].n, traversal.n)):
                whole.add(structure[n].values[i] + traversal.values[i])
            broad[n] = whole
            spreads[n] = whole.spread if whole.n else 0.0
        if not all(broad[n].n for n in names):
            continue
        ranked = sorted(names, key=lambda n: broad[n].median)
        best, worst = ranked[0], ranked[-1]
        gap = broad[worst].median / broad[best].median if broad[best].median else 0.0
        # The same discipline as the mode comparison: a difference smaller than
        # the noise is not a result, and naming a winner there would invent one.
        noise = 1.0 + max(spreads[best], spreads[worst])
        verdict = f"{best} {gap:.2f}x" if gap > noise else "tie"
        table.add(SCENE_LABEL.get(scene, scene), mode_label(mode),
                  *[f"{structure[n].median:.0f}" for n in names],
                  *[f"{broad[n].median:.0f}" for n in names],
                  verdict)
    return table


def broadphase_variant_table(per_strategy: dict[str, dict[tuple[str, str], SceneSummary]],
                             source: str) -> Table:
    """
    The edge-edge queries against each other, on the host.

    A separate question from `broadphase_table`, which asks whether to sort or to
    bin and answers it for both processors. This one takes binning as settled and
    asks which query to run over the cells: every cell a box touches, or the one
    cell its minimum corner is in, with or without each cell ordered on the axis
    the grid does not use. Only the host is asked, because that is where all of
    them exist.

    One column per strategy, the whole broad phase rather than its two halves,
    because the halves trade against each other -- ordering a cell is
    preparation bought back in traversal -- and the sum is what a caller pays.
    """
    names = [n for n in ("sweep", "cell2d", "cell2dmin", "cell2dminsort") if n in per_strategy]
    table = Table(
        label="tab:bpvariant",
        caption=("Edge-edge queries over the same cases on the host, as "
                 "whole-scene totals in milliseconds, median over repeats. " + TAGS +
                 "Every column is BP full, the acceleration structure and the "
                 "queries over it added, because ordering a cell's entries is "
                 "preparation bought back in traversal and only the sum is "
                 "comparable. All of them report identical candidate pairs."),
        columns=([Column("scene", "l")]
                 + [Column(f"{n} BP full ms", tex_header=f"{n} BP full") for n in names]
                 + [Column("fastest", "l")]),
        source=source,
        notes=("`fastest` names the winning strategy and by how much against the "
               "slowest. A margin inside the run-to-run spread is a tie."),
    )
    keys = sorted({k for s in per_strategy.values() for k in s if not k[1].startswith("device-")})
    totals = {n: 0.0 for n in names}
    for scene, mode in keys:
        broad, spreads = {}, {}
        for n in names:
            s_ = per_strategy[n].get((scene, mode))
            structure = s_.totals["prep_ms"] if s_ else Stat()
            traversal = s_.totals["broad_ms"] if s_ else Stat()
            whole = Stat()
            for i in range(min(structure.n, traversal.n)):
                whole.add(structure.values[i] + traversal.values[i])
            broad[n] = whole
            spreads[n] = whole.spread if whole.n else 0.0
        if not all(broad[n].n for n in names):
            continue
        ranked = sorted(names, key=lambda n: broad[n].median)
        best, worst = ranked[0], ranked[-1]
        for n in names:
            totals[n] += broad[n].median
        gap = broad[worst].median / broad[best].median if broad[best].median else 0.0
        noise = 1.0 + max(spreads[best], spreads[worst])
        verdict = f"{best} {gap:.2f}x" if gap > noise else "tie"
        table.add(SCENE_LABEL.get(scene, scene),
                  *[f"{broad[n].median:.0f}" for n in names], verdict)
    if any(totals.values()):
        best = min(names, key=lambda n: totals[n])
        gap = max(totals.values()) / totals[best] if totals[best] else 0.0
        table.add("all scenes", *[f"{totals[n]:.0f}" for n in names], f"{best} {gap:.2f}x")
    return table

def processor_table(summaries: dict[tuple[str, str], SceneSummary],
                    mode: str, source: str) -> Table:
    """
    CPU against GPU for one mode, with the phase that explains the difference.

    A single end-to-end ratio hides the mechanism: the two processors do not win
    the same phase, and on one scene they do not even agree on the sign. Giving
    the per-phase ratios beside the total makes the total readable -- and makes
    it obvious when a scene is the exception.
    """
    gpu = "device-" + mode
    table = Table(
        label="tab:processor",
        caption=("Host against device for the same mode and the same cases. "
                 "Every time is a whole-scene total in milliseconds, summed over "
                 "every case of the scene. " + TAGS +
                 "\\emph{total} is BP full + NP EToI, median "
                 "over repeats. A ratio above one means the GPU is faster."),
        columns=[Column("scene", "l"), Column("CPU ms"), Column("GPU ms"),
                 Column("total", tex_header=r"total$\times$"),
                 Column("broad", tex_header=r"BP full$\times$"),
                 Column("narrow", tex_header=r"NP EToI$\times$")],
        source=source,
        notes=("A ratio is the host median over the device median, so 2.0 means "
               "the device takes half the time. Ratios below 1.0 are the cases "
               "where the host wins and are the ones worth reading."),
    )
    for scene in sorted({s for s, _ in summaries}):
        host, dev = summaries.get((scene, mode)), summaries.get((scene, gpu))
        if not host or not dev:
            continue

        def med(s_, col):
            st = s_.totals.get(col)
            return st.median if st is not None and st.n else float("nan")

        parts = ("prep_ms", "broad_ms", "narrow_ms")
        h_tot = sum(med(host, c) for c in parts)
        d_tot = sum(med(dev, c) for c in parts)
        # Broad phase means the structure and the traversal over it.
        h_b = med(host, "prep_ms") + med(host, "broad_ms")
        d_b = med(dev, "prep_ms") + med(dev, "broad_ms")
        h_n, d_n = med(host, "narrow_ms"), med(dev, "narrow_ms")
        ratio = lambda a, b: f"{a / b:.2f}x" if b else "--"
        table.add(SCENE_LABEL.get(scene, scene),
                  f"{h_tot:,.0f}", f"{d_tot:,.0f}",
                  ratio(h_tot, d_tot), ratio(h_b, d_b), ratio(h_n, d_n))
    return table


def per_frame_table(summaries: dict[tuple[str, str], SceneSummary],
                    source: str) -> Table:
    """
    Mean cost of one simulation step.

    The whole-scene totals answer "what does this dataset cost"; they cannot be
    compared between a 79-case scene and a 4,571-case one, and they are not the
    number a solver author needs. Dividing by the case count gives the per-frame
    figure, which is what a step of that scene costs on this hardware.
    """
    table = Table(
        label="tab:per-frame",
        caption=("Mean time for one simulation step: the scene total, median "
                 "over repeats, divided by the number of frames. A step runs "
                 "both query types, so this is the cost of the vertex-face and "
                 "edge-edge work of that frame together. " + TAGS +
                 "BP full here includes building the swept boxes."),
        columns=[Column("scene", "l"), Column("frames"), Column("mode", "l"),
                 Column("broad ms", tex_header=r"BP full (ms)"), Column("narrow ms", tex_header=r"NP EToI (ms)"),
                 Column("total ms")],
        source=source,
        notes=("A mean rather than a median over steps: the scene total is what "
               "a run costs, and the mean is the only average that divides back "
               "into it."),
    )
    for scene, mode in sorted(summaries):
        s = summaries[(scene, mode)]
        if not s.frames:
            continue

        def per(col):
            stat = s.totals.get(col)
            return (stat.median / s.frames) if stat is not None and stat.n else 0.0

        # The acceleration structure is part of the broad phase, not a phase
        # beside it. The reference study reports one broad-phase time per
        # method for exactly this reason: the thirteen methods it compares
        # build different structures, so only the whole phase is comparable.
        broad = per("prep_ms") + per("broad_ms")
        narrow = per("narrow_ms")
        table.add(SCENE_LABEL.get(scene, scene), f"{s.frames:,}", mode_label(mode),
                  f"{broad:.2f}", f"{narrow:.2f}", f"{broad + narrow:.2f}")
    return table
