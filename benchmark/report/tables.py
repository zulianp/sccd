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
    if not stat.n or not math.isfinite(stat.median):
        return "--"
    if stat.n == 1:
        return f"{stat.median:.1f}"
    return f"{stat.median:.1f} ({stat.spread * 100:.1f} %)"


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
                 "curated per-step query sets that carry exact symbolic roots."),
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
                 "second, median over repeats. Unlike a whole-scene time this is "
                 "comparable between scenes of very different size."),
        columns=[Column("scene", "l"), Column("mode", "l"),
                 Column("broad Mpair/s", tex_header=r"broad (Mpair/s)"),
                 Column("narrow Mpair/s", tex_header=r"narrow (Mpair/s)")],
        source=source,
    )
    for (scene, mode), s in sorted(summaries.items()):
        b, n = s.totals["broad_ms"], s.totals["narrow_ms"]
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
        caption=("Wall-clock time per scene and narrow-phase mode, summed over "
                 "every case in the scene, median over independent repeats with "
                 "the full run-to-run range as a percentage of that median. "
                 "\\emph{prep} is broad-phase preparation, \\emph{broad} the "
                 "overlap query itself. The two narrow-phase columns are the two "
                 "output modes: \\emph{earliest} returns one time of impact for "
                 "the step, so every query prunes against the running minimum, "
                 "while \\emph{per-pair} returns one per candidate with no "
                 "shared bound. \\emph{total} is prep + broad + earliest."),
        columns=[
            Column("scene", "l"), Column("mode", "l"), Column("cases"),
            Column("pairs"), Column("rep"),
            Column("prep ms", tex_header=r"prep (ms)"),
            Column("broad ms", tex_header=r"broad (ms)"),
            Column("earliest ms", tex_header=r"earliest (ms)"),
            Column("per-pair ms", tex_header=r"per-pair (ms)"),
            Column("total ms", tex_header=r"total (ms)"),
        ],
        source=source,
        notes=("A difference smaller than the bracketed spread does not separate "
               "two modes and is not reported as a ratio anywhere in this "
               "document. Which mode is faster depends on the output mode as "
               "well as the scene, so the two are given side by side rather "
               "than one standing for the other."),
    )
    for (scene, mode), s in sorted(summaries.items()):
        prep = s.totals["prep_ms"]
        broad, narrow = s.totals["broad_ms"], s.totals["narrow_ms"]
        per_pair = s.totals["narrow_ms_s1"]
        # Total the phases a caller actually pays for, within a repeat and only
        # then across them. Preparation is part of that: on rod-twist it is
        # larger than the broad and narrow phases together, so a "total" that
        # leaves it out understates the pipeline by more than half.
        total = Stat()
        for i in range(min(prep.n, broad.n, narrow.n)):
            total.add(prep.values[i] + broad.values[i] + narrow.values[i])
        table.add(SCENE_LABEL.get(scene, scene), mode_label(mode),
                  f"{s.cases}", f"{s.queries:,}", f"{s.repeats}",
                  _ms(prep), _ms(broad), _ms(narrow), _ms(per_pair), _ms(total))
    return table


def conservativeness_table(summaries: dict[tuple[str, str], SceneSummary],
                           source: str) -> Table:
    table = Table(
        label="tab:conservativeness",
        caption=("Conservativeness against the dataset's exact roots. "
                 "\\emph{late} counts queries whose reported time of impact "
                 "falls after the true one; it must be zero, because a late "
                 "time of impact lets a simulation step through the contact, "
                 "which is the failure the search exists to prevent. False "
                 "positives cost work only and are reported for information. "
                 "\\emph{queries} is how many carry ground-truth data at all, "
                 "no-collision cases included; \\emph{toi compared} is how many "
                 "times of impact were actually placed beside an exact root, "
                 "which is what the claim rests on."),
        columns=[
            Column("scene", "l"), Column("mode", "l"),
            Column("queries", tex_header=r"queries"),
            Column("toi compared", tex_header=r"toi compared"),
            Column("late"), Column("false pos."), Column("false neg."),
            Column("mesh-path divergence", tex_header=r"mesh div."),
        ],
        source=source,
        notes=("Measured against the exact roots shipped with the dataset, not "
               "against TightInclusion: TightInclusion's own answer is itself a "
               "lower bound on the truth, so comparing against it over-reports "
               "lateness. The last column is not part of the gate. It counts "
               "cases where the earliest-impact answer computed over the *mesh* "
               "is later than the earliest exact root of the *curated queries*, "
               "which are two separately stored geometries: the mesh is read "
               "from PLY, the queries are exact dyadic rationals. It is a "
               "measure of the agreement between those two inputs rather than "
               "of the kernel, and with the mesh stored in double it is zero "
               "everywhere."),
    )
    for (scene, mode), s in sorted(summaries.items()):
        table.add(SCENE_LABEL.get(scene, scene), mode_label(mode),
                  f"{s.gt_queries:,}", f"{s.toi_compared:,}",
                  f"{s.toi_late}", f"{s.fp:,}", f"{s.fn}", f"{s.s0_late}")
    return table


def accuracy_table(summaries: dict[tuple[str, str], SceneSummary], source: str) -> Table:
    table = Table(
        label="tab:earliness",
        caption=("How far before the true time of impact each mode reports, as "
                 "the median over cases of the per-case median earliness. "
                 "Reporting early is always safe and always costs a solver step "
                 "size, so this is the accuracy axis the two modes trade "
                 "against speed."),
        columns=[Column("scene", "l"), Column("mode", "l"),
                 Column("median earliness", tex_header=r"median earliness")],
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
        table.add(SCENE_LABEL.get(scene, scene), mode_label(mode), value)
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
