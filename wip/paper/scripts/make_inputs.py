#!/usr/bin/env python3
"""
Regenerate the article's LaTeX tables from the committed benchmark CSVs.

    python3 scripts/make_inputs.py            # write generated/tables/*.tex
    python3 scripts/make_inputs.py --check    # verify the prose against them

The tables are built by `benchmark/report`, the same package that produces
docs/BENCHMARKS.md, through the same call sequence and with the same mode
selection -- so a number in the article and the same number in the documentation
come from one place. Figures are not drawn here: they are committed as PDFs under
figures/, because drawing them needs matplotlib and the article must build
without it. `benchmark/report` imports matplotlib lazily, which is what makes the
table half reachable on a machine that has no plotting stack at all.

Nothing under benchmark/ or docs/ is written.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
REPO = PAPER.parent.parent
sys.path.insert(0, str(REPO / "benchmark"))

RESULTS = REPO / "benchmark" / "results"
# The same CSV `benchmark/scripts/regenerate_reports.sh` builds the documents
# from, so this script and that one cannot disagree about a cell. It holds the
# shipped broad phase and the sweep it is measured against.
BENCH_CSV = REPO / "benchmark" / "assessment" / "broadphase-cell2dmin.csv"
ORACLE_CSV = RESULTS / "oracle-gh200-all.csv"
# Order matters: the report renders one series per file in the order given, and
# these are listed to match docs/BENCHMARKS.md -- host then device -- so the
# article's table and the library's are the same table. Both runs use the
# shipped broad phase, so the one axis the study varies is the processor.
SCALING = [RESULTS / "scaling" / n for n in
           ("host-cell2dmin-mode2.txt", "device-cell2dmin-mode2.txt")]

# The two mode selections the two committed documents are generated with.
# `tight` is the article's subject; `relaxed` supplies the trade-off subsection.
TIGHT = {"tight", "device-tight"}
RELAXED = {"relaxed", "device-relaxed"}


def _build(only_modes: set[str], label_override: dict[str, str], suffix: str):
    """The table half of `python3 -m report`, for one mode selection."""
    from report import data, oracle as oracle_mod, scaling as scaling_mod, style, tables

    style.LABEL_OVERRIDE.update(label_override)

    rows = [r for r in data.read_rows(BENCH_CSV) if r["mode"] in only_modes]
    if not rows:
        raise SystemExit(f"no rows for modes {sorted(only_modes)}")
    data.check_schema(BENCH_CSV, rows)

    # A CSV may hold more than one broad-phase strategy; summarising across them
    # would sum two runs of the same cases into one scene total. The report picks
    # one for the headline tables and compares the rest against it, preferring
    # the cell list because that is what the shipped default probes first.
    strategies = data.broadphases(rows)
    primary = ("cell2d" if "cell2d" in strategies else strategies[0]) if strategies else None
    scenes = data.by_scene(rows, primary)

    source = str(BENCH_CSV.relative_to(REPO))
    oracle_rows = oracle_mod.read(ORACLE_CSV)
    # The reference is kept whatever the mode selection: it is what the article
    # compares against, not one of the subjects.
    keep = set(only_modes) | {"tight-inclusion"}
    oracle_rows = {k: v for k, v in oracle_rows.items() if k[2] in keep}

    built = [
        tables.dataset_table(scenes, oracle_rows, source),
        tables.timing_table(scenes, source),
        tables.throughput_table(scenes, source),
        tables.conservativeness_table(scenes, source),
        tables.accuracy_table(scenes, source),
        tables.per_frame_table(scenes, source),
    ]
    host_modes = sorted({m for _, m in scenes if not m.startswith("device-")})
    if host_modes:
        pick = "tight" if "tight" in host_modes else host_modes[0]
        if any(m == "device-" + pick for _, m in scenes):
            built.append(tables.processor_table(scenes, pick, source))
    if len(strategies) > 1:
        built.append(tables.broadphase_table(
            {n: data.by_scene(rows, n) for n in strategies}, source))

    oracle_source = str(ORACLE_CSV.relative_to(REPO))
    built.append(oracle_mod.gate_table(oracle_rows, oracle_source))
    built.append(oracle_mod.reference_table(oracle_rows, oracle_source))
    built.append(oracle_mod.earliness_table(oracle_rows, oracle_source))

    runs = [r for r in (scaling_mod.parse(p) for p in SCALING) if r.faces]
    if runs and not suffix:
        built.append(scaling_mod.table(
            runs, ", ".join(str(p.relative_to(REPO)) for p in SCALING)))

    # A late time of impact would make every other number meaningless, so it is
    # a hard failure here rather than a line in the output.
    late = sum(s.toi_late for s in scenes.values()) + oracle_mod.violations(oracle_rows)
    if late:
        raise SystemExit(f"FAILED: {late} late times of impact in {sorted(only_modes)}")

    return built, tables


# `benchmark/report.tables.render_tex` now emits the article's LaTeX directly:
# it leaves the caption and the notes as their authors wrote them, gives a
# note's backticks \texttt, and wraps every tabular in \fittable. Nothing is
# rewritten here, so this script and `scripts/make_tables.sh` produce the same
# file for the same table.


def _write(built, tables_mod, suffix: str) -> list[str]:
    out = PAPER / "generated" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for t in built:
        # tab:foo -> tab-foo.tex, so the file name says which \ref resolves to it.
        stem = t.label.replace(":", "-") + suffix
        # Relaxed tables reuse the Tight labels, which would be a duplicate-label
        # warning and a \ref pointing at whichever came last. Suffix both.
        body = tables_mod.render_tex(t)
        if suffix:
            body = body.replace("\\label{" + t.label + "}",
                                "\\label{" + t.label + suffix + "}")
        (out / f"{stem}.tex").write_text(body)
        written.append(stem)
    return written


# Figures quoted in the prose, and the table each must still appear in. A claim
# that drifts from its evidence is the failure this guards against.
PROSE_CLAIMS = [
    ("15.80", "tab-processor"),
    ("2.05", "tab-processor"),
    # The scaling ratios are quoted as divisions of these cells, so the reader
    # can do the arithmetic; checking the cells checks the ratios.
    ("55.7", "tab-scaling"),
    ("75.4", "tab-scaling"),
    ("242.7", "tab-scaling"),
    ("297.5", "tab-scaling"),
    ("441.7", "tab-scaling"),
    ("823.9", "tab-scaling"),
    ("1828.3", "tab-scaling"),
    ("3644.3", "tab-scaling"),
    ("840.2", "tab-scaling"),
    ("1009.4", "tab-scaling"),
    ("2804.1", "tab-scaling"),
    ("818.9", "tab-scaling"),
]

# Totals the prose states that are sums of a generated table's columns rather
# than cells of it, so a grep cannot find them: (column index, expected total).
# The reader can add the column up; this keeps us honest if the data changes.
COLUMN_TOTALS = [
    ("tab-conservativeness", 1, 5_815_032),   # queries carrying ground truth
    ("tab-conservativeness", 2, 5_522_383),   # of those, compared against a root
]


def _check() -> int:
    gen = PAPER / "generated" / "tables"
    if not gen.is_dir():
        print("error: generated/tables does not exist; run without --check first",
              file=sys.stderr)
        return 2
    blob = "\n".join(p.read_text() for p in sorted(gen.glob("*.tex")))
    prose = "\n".join(p.read_text() for p in sorted((PAPER / "sections").glob("*.tex")))

    status = 0
    for needle, table in PROSE_CLAIMS:
        if needle not in prose:
            print(f"note: {needle} no longer appears in the prose", file=sys.stderr)
            continue
        if needle not in (gen / f"{table}.tex").read_text():
            print(f"error: prose claims {needle} but {table}.tex does not contain it",
                  file=sys.stderr)
            status = 1

    # Column totals the prose quotes. The conservativeness table carries one row
    # per scene now -- the processors agree on these columns -- so a column sums
    # to the per-processor figure directly, and the prose may quote that or the
    # doubled one it gets by checking both processors.
    for table, col, expected in COLUMN_TOTALS:
        total = 0
        for line in (gen / f"{table}.tex").read_text().splitlines():
            line = line.strip()
            if not line.endswith(r"\\") or "&" not in line:
                continue
            cells = [c.strip() for c in line[:-2].split("&")]
            if len(cells) <= col or not cells[col].replace(",", "").isdigit():
                continue
            total += int(cells[col].replace(",", ""))
        if total != expected:
            print(f"error: {table} column {col} sums to {total:,}, "
                  f"prose says {expected:,}", file=sys.stderr)
            status = 1
        else:
            # The prose may quote the per-processor figure or the doubled one.
            wanted = {f"{expected:,}", f"{2 * expected:,}"}
            if not any(w.replace(",", "{,}") in prose for w in wanted):
                print(f"note: neither {expected:,} nor {2 * expected:,} appears "
                      f"in the prose", file=sys.stderr)

    # The article's tables and the library's documentation must be the same tables.
    # This is a stronger statement than "both were generated from the same CSV",
    # and it is the one a reader of either document cares about.
    status |= _check_against_docs()

    # Every \ref and \input target must exist.
    for m in re.finditer(r"\\input\{generated/tables/([^}]+)\}", prose):
        if not (gen / f"{m.group(1)}.tex").is_file():
            print(f"error: \\input of a missing table: {m.group(1)}", file=sys.stderr)
            status = 1
    return status


# Generated block in docs/BENCHMARKS.md -> the article's table of the same data.
DOC_TABLES = [
    ("dataset", "tab-dataset"),
    ("conservativeness", "tab-conservativeness"),
    ("earliness-ref", "tab-earliness-ref"),
    ("processor", "tab-processor"),
    ("per-frame", "tab-per-frame"),
    ("broadphase", "tab-broadphase"),
    ("scaling", "tab-scaling"),
]


def _check_against_docs() -> int:
    """Assert every generated table matches the committed documentation row for row."""
    doc_path = REPO / "docs" / "BENCHMARKS.md"
    if not doc_path.is_file():
        print("note: docs/BENCHMARKS.md not found; skipping the cross-check",
              file=sys.stderr)
        return 0
    doc = doc_path.read_text()
    gen = PAPER / "generated" / "tables"

    def markdown_rows(label: str):
        m = re.search(rf"<!-- sccd:begin {label} -->(.*?)<!-- sccd:end {label} -->",
                      doc, re.S)
        if m is None:
            return None
        rows = []
        for line in m.group(1).splitlines():
            line = line.strip()
            # Skip prose, the header rule, and the "Source:" footer.
            if not line.startswith("|") or set(line) <= set("|-: "):
                continue
            rows.append([c.strip() for c in line.strip("|").split("|")])
        return rows[1:]

    def tex_rows(name: str):
        rows = []
        for line in (gen / f"{name}.tex").read_text().splitlines():
            line = line.strip()
            if "&" not in line or not line.endswith(r"\\"):
                continue
            rows.append([c.strip() for c in line[:-2].split("&")])
        return rows[1:]

    status = 0
    for label, name in DOC_TABLES:
        md = markdown_rows(label)
        if md is None:
            print(f"note: no '{label}' block in docs/BENCHMARKS.md", file=sys.stderr)
            continue
        tex = tex_rows(name)
        if md != tex:
            print(f"error: {name}.tex does not match the '{label}' block of "
                  f"docs/BENCHMARKS.md ({len(md)} rows there, {len(tex)} here)",
                  file=sys.stderr)
            for i, (x, y) in enumerate(zip(md, tex)):
                if x != y:
                    print(f"       first difference at row {i}:\n"
                          f"         docs: {x}\n         here: {y}", file=sys.stderr)
                    break
            status = 1
    return status


def main(argv: list[str]) -> int:
    if "--check" in argv[1:]:
        return _check()
    if not BENCH_CSV.is_file():
        print(f"error: {BENCH_CSV} does not exist", file=sys.stderr)
        return 2

    built, tables_mod = _build(TIGHT, {"tight": "CPU", "device-tight": "GPU"}, "")
    written = _write(built, tables_mod, "")
    built_r, _ = _build(RELAXED, {"relaxed": "CPU", "device-relaxed": "GPU"}, "-relaxed")
    written += _write(built_r, tables_mod, "-relaxed")

    print(f"wrote {len(written)} tables into generated/tables:")
    for w in sorted(written):
        print(f"  {w}.tex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
