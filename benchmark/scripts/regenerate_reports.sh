#!/bin/bash
# Regenerate every document and paper input from the committed data.
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
# `python3 -m report` writes figures under its out-dir; the committed copies in
# docs/figures and the paper have to be refreshed from there, or the tables move
# and the plots stay behind. This is the step that was missing.
cp /tmp/report-tight/figures/*.pdf /tmp/report-tight/figures/*.png docs/figures/
cp /tmp/report-relaxed/figures/*.pdf /tmp/report-relaxed/figures/*.png docs/figures/
# toi-error is not in this list: the article draws its own two-row version,
# with additive CCD beneath ours, from the competitor run.
for f in narrow-per-case reference-speedup refine-scaling results-grid \
         runtime-breakdown; do
    cp docs/figures/$f.pdf wip/paper/figures/$f.pdf
done
# The two the paper draws itself, from the profile data.
( cd wip/paper && ../../.venv/bin/python scripts/make_figures.py )

echo "== check both =="
$PY -m report "$BENCH" /tmp/report-tight \
    benchmark/results/oracle-gh200-all.csv $SC \
    --modes=tight,device-tight --label="tight:CPU,device-tight:GPU" \
    --embed=docs/BENCHMARKS.md --check
$PY -m report "$BENCH" /tmp/report-relaxed \
    benchmark/results/oracle-gh200-all.csv \
    --modes=relaxed,device-relaxed --label="relaxed:CPU,device-relaxed:GPU" \
    --figure-prefix=relaxed- --embed=docs/BENCHMARKS_RELAXED.md --check
