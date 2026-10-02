#!/bin/bash
# Regenerate the paper's generated tables from the committed benchmark data.
#
#     scripts/make_tables.sh [python]
#
# The report takes flags that change what the tables say, and getting them wrong
# is not a build error -- it produces a table that looks right and is not. The
# two that matter:
#
#   --label   The mode column of this paper reads CPU and GPU. Without the
#             override the report prints the internal mode name instead, so the
#             column comes out "Tight" and "Tight (GPU)". Every table then
#             disagrees with the ones nobody regenerated, and LaTeX has nothing
#             to complain about.
#   <oracle>  The collision counts in tab:dataset come from the oracle CSV.
#             Omit it and those cells become "--", quietly dropping two columns.
#
# So the invocation lives here rather than in a comment, and this script is what
# regenerates the tables.
set -euo pipefail

PAPER="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$PAPER/../.." && pwd)"
PY="${1:-python3}"

BENCH="$REPO/benchmark/assessment/broadphase-cell2dminfv.csv"
ORACLE="$REPO/benchmark/results/oracle-gh200-all.csv"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT

cd "$REPO"

# A table is only as good as the run behind it, and the ways a run can measure
# the wrong thing do not show up in its numbers: a strategy name the build did
# not know and so raced instead, a submit path that lost the thread binding, a
# merge that ate a row per chunk. Those are provenance faults, and this refuses
# to publish a file carrying one.
"$PY" "$REPO/benchmark/scripts/validate_results.py" timings "$BENCH" \
      --expect-bp cell2dminfv || {
    echo "make_tables: refusing to regenerate from $BENCH" >&2
    exit 1
}

"$PY" -m benchmark.report "$BENCH" "$OUT" "$ORACLE" \
      --modes=tight,device-tight \
      --label="tight:CPU,device-tight:GPU"

"$PY" - "$OUT/tables.tex" "$PAPER/generated/tables" "$BENCH" <<'PY'
import re, sys, os
tables_tex, out_dir, bench = sys.argv[1], sys.argv[2], sys.argv[3]
src = open(tables_tex).read()
blocks, pos = [], 0
while True:
    a = src.find('\\begin{table}', pos)
    if a < 0:
        break
    b = src.index('\\end{table}', a) + len('\\end{table}')
    blocks.append(src[a:b]); pos = b

rel = os.path.relpath(bench, os.path.dirname(os.path.dirname(out_dir)) + "/../..")
written, skipped = [], []
for blk in blocks:
    m = re.search(r"label\{tab:([a-z-]+)\}", blk)
    if not m:
        continue
    path = os.path.join(out_dir, f"tab-{m.group(1)}.tex")
    if not os.path.exists(path):
        skipped.append(f"tab-{m.group(1)} (no committed table)")
        continue
    # A cell the report could not fill is "--". Overwriting a table that has the
    # value with one that does not is how content gets lost silently, so it is
    # refused here rather than reviewed later.
    def cells(t):
        body = t[t.index("\\midrule") + 9:t.index("\\bottomrule")]
        return [c.strip() for r in body.strip().split("\\\\") if r.strip()
                for c in r.split("&")]
    old = open(path).read()
    lost = sum(1 for o, n in zip(cells(old), cells(blk)) if n == "--" and o != "--")
    if lost:
        skipped.append(f"tab-{m.group(1)} ({lost} cells would become --)")
        continue
    # Only the Source line, and only when it already names the bench CSV. A
    # blunter rewrite reaches the \texttt in a note -- which made one caption
    # read "<path> names the winning strategy" -- and overwrites the Source of
    # a table whose data comes from the oracle CSV instead.
    # The path comes from the file actually read, not from a literal: a literal
    # here is how every table came to cite a CSV that was not its source after a
    # bulk rename, and `rel` was computed for this and then left unused.
    txt = re.sub(r"(Source: )\\texttt\{(?!benchmark/results/oracle)[^}]*\}",
                 lambda m: m.group(1) + "\\texttt{" + rel + "}", blk)
    open(path, "w").write(txt + "\n")
    written.append(f"tab-{m.group(1)}")

# The mode column is the thing --label controls, so it is checked rather than
# trusted: this is the fault the script exists to prevent.
bad = []
for name in written:
    txt = open(os.path.join(out_dir, name + ".tex")).read()
    if re.search(r"& (Tight|Relaxed)( \(GPU\))? &", txt):
        bad.append(name)
if bad:
    raise SystemExit("mode column is not CPU/GPU in: " + ", ".join(bad))

print("wrote:", ", ".join(written))
if skipped:
    print("skipped:", "; ".join(skipped))
PY
