"""
Generate every benchmark artifact from one sweep CSV.

    python3 -m report <bench.csv> <out-dir> [--scaling <scaling.csv>]

Writes, under <out-dir>:
    figures/*.pdf   figures for LaTeX
    figures/*.png   the same figures for Markdown
    figures.tex     \\includegraphics with a caption and label for each
    tables.tex      booktabs tables
    tables.md       the same tables as Markdown
    summary.md      tables, comparisons and figure references, ready to embed
"""

from __future__ import annotations

import sys
from pathlib import Path

from . import data, figures, style, tables


def main(argv: list[str]) -> int:
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) < 2:
        print(__doc__, file=sys.stderr)
        return 2

    bench_csv = Path(args[0])
    out_dir = Path(args[1])
    if not bench_csv.is_file():
        print(f"error: {bench_csv} does not exist", file=sys.stderr)
        return 2

    figure_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = data.read_rows(bench_csv)
    try:
        data.check_schema(bench_csv, rows)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    scenes = data.by_scene(rows)
    cases = data.by_case(rows)
    source = str(bench_csv)

    style.apply_rcparams()
    drawn = [
        figures.phase_breakdown(scenes, figure_dir),
        figures.narrow_phase_per_case(cases, figure_dir),
        figures.earliness_distribution(cases, figure_dir),
    ]
    figures.write_figure_tex(drawn, out_dir)

    built = [
        tables.timing_table(scenes, source),
        tables.conservativeness_table(scenes, source),
        tables.accuracy_table(scenes, source),
    ]
    tables.write_tables(built, out_dir)

    # A run that reported a late time of impact must be impossible to overlook,
    # so it is stated before anything else and sets the exit status.
    late = sum(s.toi_late + s.s0_late for s in scenes.values())

    lines = ["## Results", ""]
    if late:
        lines += [
            f"> **{late} queries reported a time of impact after the true one.** "
            f"This is a correctness failure, not an accuracy one: a late time of "
            f"impact lets a simulation step through the contact. Every number "
            f"below is suspect until it is fixed.", ""]
    else:
        checked = sum(s.gt_queries for s in scenes.values())
        lines += [
            f"Across {checked:,} queries with an exact root, no mode reported a "
            f"time of impact after the true one and none missed a collision.", ""]

    for table in built:
        lines += [f"### {table.caption}", "", tables.render_markdown(table), ""]

    notes = tables.comparison_notes(scenes)
    if notes:
        lines += ["### Mode comparison", "",
                  "Only differences larger than the run-to-run spread that "
                  "measured them are stated as ratios.", ""] + notes + [""]

    lines += ["### Figures", ""]
    for f in drawn:
        lines += [f"![{f.caption}](figures/{f.stem}.png)", "", f"*{f.caption}*", ""]

    (out_dir / "summary.md").write_text("\n".join(lines))

    print(f"wrote {out_dir}/summary.md, tables.tex, tables.md, figures.tex "
          f"and {len(drawn)} figures")
    if late:
        print(f"FAILED: {late} late times of impact", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
