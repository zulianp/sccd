#!/bin/bash
# Regenerate every document and paper input from the committed data.
# Run from the repository root with the analysis venv present.
set -euo pipefail
PY="${PY:-./.venv/bin/python}"
export PYTHONPATH=benchmark
SC="benchmark/results/scaling/host-cell2d-mode2.txt benchmark/results/scaling/host-sweep-mode2.txt benchmark/results/scaling/device-mode2.txt"

echo "== docs/BENCHMARKS.md =="
$PY -m report benchmark/results/sweep-gh200-bp.csv /tmp/report-tight \
    benchmark/results/oracle-gh200-all.csv $SC \
    --modes=tight,device-tight --label="tight:CPU,device-tight:GPU" \
    --embed=docs/BENCHMARKS.md

echo "== docs/BENCHMARKS_RELAXED.md =="
$PY -m report benchmark/results/sweep-gh200-bp.csv /tmp/report-relaxed \
    benchmark/results/oracle-gh200-all.csv \
    --modes=relaxed,device-relaxed --label="relaxed:CPU,device-relaxed:GPU" \
    --figure-prefix=relaxed- --embed=docs/BENCHMARKS_RELAXED.md

echo "== figures =="
# `python3 -m report` writes figures under its out-dir; the committed copies in
# docs/figures and the paper have to be refreshed from there, or the tables move
# and the plots stay behind. This is the step that was missing.
cp /tmp/report-tight/figures/*.pdf /tmp/report-tight/figures/*.png docs/figures/
cp /tmp/report-relaxed/figures/*.pdf /tmp/report-relaxed/figures/*.png docs/figures/
for f in narrow-per-case refine-scaling results-grid runtime-breakdown toi-error; do
    cp docs/figures/$f.pdf wip/paper/figures/$f.pdf
done
# The two the paper draws itself, from the profile data.
( cd wip/paper && ../../.venv/bin/python scripts/make_figures.py )

echo "== check both =="
$PY -m report benchmark/results/sweep-gh200-bp.csv /tmp/report-tight \
    benchmark/results/oracle-gh200-all.csv $SC \
    --modes=tight,device-tight --label="tight:CPU,device-tight:GPU" \
    --embed=docs/BENCHMARKS.md --check
$PY -m report benchmark/results/sweep-gh200-bp.csv /tmp/report-relaxed \
    benchmark/results/oracle-gh200-all.csv \
    --modes=relaxed,device-relaxed --label="relaxed:CPU,device-relaxed:GPU" \
    --figure-prefix=relaxed- --embed=docs/BENCHMARKS_RELAXED.md --check
