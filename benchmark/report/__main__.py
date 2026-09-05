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

from . import (data, embed, figures, oracle as oracle_mod,
               scaling as scaling_mod, style, tables)


def _repo_relative(path: Path) -> str:
    resolved = Path(path).resolve()
    for parent in [Path.cwd().resolve(), *Path.cwd().resolve().parents]:
        if (parent / ".git").exists():
            try:
                return str(resolved.relative_to(parent))
            except ValueError:
                break
    return str(path)


def _provenance(bench_csv: Path, oracle_csv: Path | None, scenes: dict) -> str:
    """What was run, so a reader can tell what the tables above are of."""
    repeats = sorted({s.repeats for s in scenes.values()})
    cases = sum(s.cases for s in scenes.values()) // max(
        len({m for _, m in scenes}), 1)
    lines = [
        f"- Timings: `{_repo_relative(bench_csv)}`, {cases} cases over "
        f"{len({d for d, _ in scenes})} scenes, "
        f"{repeats[0] if len(repeats) == 1 else f'{min(repeats)}-{max(repeats)}'}"
        f" independent repeats.",
    ]
    if oracle_csv is not None:
        lines.append(f"- Accuracy: `{_repo_relative(oracle_csv)}`, every query of every scene "
                     f"checked against the dataset's exact roots.")
    lines.append("- Regenerate with `python3 -m report <bench.csv> <out> "
                 "<oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to "
                 "assert the document still matches the data.")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    flags = [a for a in argv[1:] if a.startswith("--")]
    args = [a for a in argv[1:] if not a.startswith("--")]
    check_only = "--check" in flags
    embed_into = None
    for flag in flags:
        if flag.startswith("--embed="):
            embed_into = Path(flag.split("=", 1)[1])
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
    # Cite the CSV by its path in the repository, not by wherever it happened to
    # be read from: a "Source:" line naming a scratch directory tells a reader
    # nothing they can act on.
    source = _repo_relative(bench_csv)

    style.apply_rcparams()
    drawn = [
        figures.phase_breakdown(scenes, figure_dir),
        figures.narrow_phase_per_case(cases, figure_dir),
        figures.earliness_distribution(cases, figure_dir),
    ]
    scaling_runs = []
    if scaling_files:
        scaling_runs = [r for r in (scaling_mod.parse(p) for p in scaling_files)
                        if r.faces]
        if scaling_runs:
            drawn.append(scaling_mod.figure(scaling_runs, figure_dir))
    figures.write_figure_tex(drawn, out_dir)

    built = [
        tables.timing_table(scenes, source),
        tables.conservativeness_table(scenes, source),
        tables.accuracy_table(scenes, source),
    ]
    oracle_rows = {}
    if oracle_csv and oracle_csv.is_file():
        oracle_rows = oracle_mod.read(oracle_csv)
        oracle_source = _repo_relative(oracle_csv)
        built.append(oracle_mod.gate_table(oracle_rows, oracle_source))
        built.append(oracle_mod.reference_table(oracle_rows, oracle_source))
    if scaling_runs:
        built.append(scaling_mod.table(
            scaling_runs, ", ".join(_repo_relative(p) for p in scaling_files)))
    tables.write_tables(built, out_dir)

    # A run that reported a late time of impact must be impossible to overlook,
    # so it is stated before anything else and sets the exit status.
    # The gate is the same-geometry comparison only. `s0_late` compares the mesh
    # path against roots belonging to the curated query set, which are different
    # geometries, so it cannot decide conservativeness and does not gate.
    late = sum(s.toi_late for s in scenes.values())
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
        mesh_div = sum(s.s0_late for s in scenes.values())
        if mesh_div:
            lines += [
                f"({mesh_div} cases show the mesh-path earliest-impact answer "
                f"landing after the curated queries' earliest exact root. Those "
                f"are two different geometries -- smesh stores mesh coordinates "
                f"as float32 -- and on every one of them the curated-query "
                f"answer is at or before the root, so it is not a "
                f"conservativeness failure. See the last column of the "
                f"conservativeness table.)", ""]

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

    # Refresh the generated blocks of docs/BENCHMARKS.md, or, with --check,
    # verify that the committed document is what the committed data and this
    # script reproduce. That is the difference between claiming a number is
    # reproducible and demonstrating it.
    if embed_into is not None:
        blocks = {t.label.split(":", 1)[-1]: tables.render_markdown(t) for t in built}
        blocks["comparison"] = "\n".join(notes) if notes else "_No comparison available._"
        blocks["provenance"] = _provenance(bench_csv, oracle_csv, scenes)
        status = embed.apply(embed_into, blocks, check=check_only)
        if status:
            return status

    print(f"wrote {out_dir}/summary.md, tables.tex, tables.md, figures.tex "
          f"and {len(drawn)} figures")
    if late:
        print(f"FAILED: {late} late times of impact", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
