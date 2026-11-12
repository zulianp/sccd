#!/bin/bash
# Regenerate the benchmark documents and the article's figures from the committed data.
# Run from the repository root with the analysis venv present.
set -euo pipefail
PY="${PY:-./.venv/bin/python}"
export PYTHONPATH=benchmark
# The shipped broad phase on both processors, and nothing else: the paper
# reports no other strategy outside the Scalable CCD comparison.
SC="benchmark/results/scaling/host-cell2dminfv-mode2.txt benchmark/results/scaling/device-cell2dminfv-mode2.txt"
# One CSV behind every table and figure of both documents and the article, so a
# number in one is the number in the others. It holds the shipped broad phase
# and the sweep it is measured against.
BENCH="benchmark/assessment/broadphase-cell2dminfv.csv"

# Same gate as the paper's tables: the documents quote the same measurements, so
# they are refused on the same grounds.
"$PY" benchmark/scripts/validate_results.py timings "$BENCH" --expect-bp cell2dminfv || {
    echo "regenerate_reports: refusing to regenerate from $BENCH" >&2
    exit 1
}

echo "== docs/BENCHMARKS.md =="
$PY -m report "$BENCH" /tmp/report-tight \
    benchmark/results/oracle-gh200-all.csv $SC \
    --modes=tight,device-tight --label="tight:CPU,device-tight:GPU" \
    --embed=docs/BENCHMARKS.md

echo "== docs/BENCHMARKS_RELAXED.md =="
$PY -m report "$BENCH" /tmp/report-relaxed \
    benchmark/results/oracle-gh200-all.csv \
    --modes=relaxed,device-relaxed --label="relaxed:CPU,device-relaxed:GPU" \
    --figure-prefix=relaxed- --embed=docs/BENCHMARKS_RELAXED.md

echo "== figures =="
# The documents carry tables only; every figure lives with the article, in
# zulianp/uvic_papers under sccd/figures/. `python3 -m report` draws its figures
# under its out-dir, and they are copied there from it, or the tables move and the
# plots stay behind. The checkout is PAPER_DIR, or beside this repository; without
# one the figures stay under /tmp and the documents above are unaffected.
PAPER_DIR="${PAPER_DIR:-../uvic_papers/sccd}"
if [ -d "$PAPER_DIR/figures" ]; then
    cp /tmp/report-tight/figures/*.pdf /tmp/report-relaxed/figures/*.pdf "$PAPER_DIR/figures/"
    # The figures the article draws itself, from the profile and competitor data.
    SCCD_ROOT="$PWD" "$PY" "$PAPER_DIR/scripts/make_figures.py"
else
    echo "regenerate_reports: no article checkout at $PAPER_DIR; figures left in /tmp/report-*/figures" >&2
fi

echo "== check both =="
$PY -m report "$BENCH" /tmp/report-tight \
    benchmark/results/oracle-gh200-all.csv $SC \
    --modes=tight,device-tight --label="tight:CPU,device-tight:GPU" \
    --embed=docs/BENCHMARKS.md --check
$PY -m report "$BENCH" /tmp/report-relaxed \
    benchmark/results/oracle-gh200-all.csv \
    --modes=relaxed,device-relaxed --label="relaxed:CPU,device-relaxed:GPU" \
    --figure-prefix=relaxed- --embed=docs/BENCHMARKS_RELAXED.md --check
