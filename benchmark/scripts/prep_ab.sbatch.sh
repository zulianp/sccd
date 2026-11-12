#!/bin/bash
#SBATCH --job-name=sccd-prep-ab
#SBATCH --account=c40
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
# Strong scaling of the host pipeline before and after the preparation work,
# same protocol as the committed strong-*.csv: four cases per scene, Tight, the
# broad-phase tuner left on, two repeats per thread count.
set -uo pipefail
OUT=$SCRATCH/sccd-profile/prep-ab; mkdir -p $OUT
DATA=$SCRATCH/sccd-assess/data
for sc in cloth-funnel armadillo-rollers; do
  for which in base opt; do
    BIN=$SCRATCH/sccd-prep-$which/build-grace/sccd_bench
    f=$OUT/strong-$which-$sc.csv
    : > "$f"
    for t in 1 2 4 8 16 32 64 72; do
      for rep in 1 2; do
        out=$(OMP_NUM_THREADS=$t SCCD_BENCH_MAX_CASES=4 SCCD_NARROWPHASE_MODE=2 \
              timeout 1500 "$BIN" "$DATA" "$sc" 2>/dev/null | grep -v "^dataset,")
        echo "$out" | awk -v t=$t -v r=$rep -F, "NF>5 {print t\",\"r\",\"\$0}" >> "$f"
      done
    done
    echo "$which $sc rows=$(wc -l < "$f")"
  done
done
echo ALLDONE
