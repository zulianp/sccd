"""
Generate every benchmark artifact from one sweep CSV.

    python3 -m report <bench.csv> <out-dir> [oracle.csv] [scaling.txt ...]

An oracle.csv is ti_oracle's --csv output; scaling files are
sccd_refine_scaling's stdout. Both are optional -- the timing report is
generated without them, and each adds its own tables or figure when given.

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

from . import data, figures, oracle as oracle_mod, scaling as scaling_mod, style, tables


def main(argv: list[str]) -> int:
    args = [a for a in argv[1:] if not a.startswith("--")]
    if len(args) < 2:
        print(__doc__, file=sys.stderr)
        return 2

    bench_csv = Path(args[0])
    out_dir = Path(args[1])
    extra = [Path(a) for a in args[2:]]
    oracle_csv = next((p for p in extra if p.suffix == ".csv"), None)
    scaling_files = [p for p in extra if p.suffix != ".csv"]
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
    if scaling_files:
        runs = [scaling_mod.parse(p) for p in scaling_files]
        runs = [r for r in runs if r.faces]
        if runs:
            drawn.append(scaling_mod.figure(runs, figure_dir))
    figures.write_figure_tex(drawn, out_dir)

    built = [
        tables.timing_table(scenes, source),
        tables.conservativeness_table(scenes, source),
        tables.accuracy_table(scenes, source),
    ]
    oracle_rows = {}
    if oracle_csv and oracle_csv.is_file():
        oracle_rows = oracle_mod.read(oracle_csv)
        built.append(oracle_mod.gate_table(oracle_rows, str(oracle_csv)))
        built.append(oracle_mod.reference_table(oracle_rows, str(oracle_csv)))
    tables.write_tables(built, out_dir)

    # A run that reported a late time of impact must be impossible to overlook,
    # so it is stated before anything else and sets the exit status.
    late = sum(s.toi_late + s.s0_late for s in scenes.values())
    late += oracle_mod.violations(oracle_rows) if oracle_rows else 0

    lines = ["## Results", ""]
    if late:
        lines += [
            f"> **{late} queries reported a time of impact after the true one.** "
            f"This is a correctness failure, not an accuracy one: a late time of "
            f"impact lets a simulation step through the contact. Every number "
            f"below is suspect until it is fixed.", ""]
    else:
        checked = sum(s.gt_queries for s in scenes.values())
        if oracle_rows:
            checked = max(checked, sum(r.gt_checked for (_, _, m), r in
                                       oracle_rows.items() if m != "tight-inclusion"))
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
