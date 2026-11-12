#!/usr/bin/env python3
"""Summarise one comparison run: SCCD against Scalable CCD and Additive CCD.

Timings are per case (one query kind of one simulation step), in milliseconds,
and both the median and the maximum are reported: in continuous collision
detection the worst case is the number a solver has to survive, and a median
hides it by orders of magnitude.
"""
import csv
import statistics as st
import sys
from collections import defaultdict


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fmt(v, w=9, p=3):
    return f"{v:{w}.{p}f}" if v is not None else " " * (w - 1) + "-"


rows = list(csv.DictReader(open(sys.argv[1])))
per_case = [r for r in rows if r.get("type") != "step"]
step_rows = [r for r in rows if r.get("type") == "step"]

# ---- performance ---------------------------------------------------------
print("=" * 96)
print("PERFORMANCE   per case (one query kind of one step), milliseconds")
print("=" * 96)
for scene in sorted({r["dataset"] for r in per_case}):
    print(f"\n{scene}")
    print(f"  {'mode':<22} {'n':>4} {'broad med':>10} {'broad max':>10} "
          f"{'narrow med':>11} {'narrow max':>11} {'total med':>10}")
    acc = defaultdict(lambda: defaultdict(list))
    for r in per_case:
        if r["dataset"] != scene:
            continue
        for col in ("broad_ms", "narrow_ms"):
            v = f(r.get(col))
            if v is not None:
                acc[r["mode"]][col].append(v)
    for mode in sorted(acc):
        b, n = acc[mode]["broad_ms"], acc[mode]["narrow_ms"]
        tot = [x + y for x, y in zip(b, n)] if len(b) == len(n) else []
        print(f"  {mode:<22} {len(b):>4} {fmt(st.median(b), 10)} {fmt(max(b), 10)} "
              f"{fmt(st.median(n), 11)} {fmt(max(n), 11)} "
              f"{fmt(st.median(tot), 10) if tot else '':>10}")

# ---- accuracy, per query -------------------------------------------------
print()
print("=" * 96)
print("ACCURACY   per curated query, against the exact roots")
print("  (Scalable CCD is compared on the earliest time of impact, not per query;")
print("   compare_table.py prints both comparisons.)")
print("=" * 96)
print(f"  {'scene':<20} {'mode':<22} {'queries':>9} {'fp':>7} {'fn':>5} "
      f"{'late':>5} {'max early':>11} {'med early':>11}")
for scene in sorted({r["dataset"] for r in per_case}):
    agg = defaultdict(lambda: defaultdict(float))
    early_max, early_med = defaultdict(list), defaultdict(list)
    for r in per_case:
        if r["dataset"] != scene:
            continue
        m = r["mode"]
        for col in ("queries", "fp", "fn", "toi_late", "toi_n"):
            v = f(r.get(col))
            if v is not None:
                agg[m][col] += v
        if f(r.get("toi_n")):
            for lst, col in ((early_max, "toi_max_early"), (early_med, "toi_med_early")):
                v = f(r.get(col))
                if v is not None:
                    lst[m].append(v)
    for m in sorted(agg):
        if not agg[m]["toi_n"]:
            continue  # no per-query answer to score
        print(f"  {scene:<20} {m:<22} {int(agg[m]['toi_n']):>9} "
              f"{int(agg[m]['fp']):>7} {int(agg[m]['fn']):>5} {int(agg[m]['toi_late']):>5} "
              f"{fmt(max(early_max[m]) if early_max[m] else None, 11, 3+3)} "
              f"{fmt(st.median(early_med[m]) if early_med[m] else None, 11, 3+3)}")

# ---- accuracy, per step --------------------------------------------------
if step_rows:
    print()
    print("=" * 96)
    print("ACCURACY   earliest time of impact per simulation step")
    print("  The unit Scalable CCD's API answers in. `late` counts steps whose")
    print("  reported time of impact is after the earliest exact root, which is")
    print("  the failure that lets a solver step through a contact.")
    print("=" * 96)
    print(f"  {'scene':<20} {'mode':<22} {'steps':>6} {'late':>5} {'worst margin':>14}")
    agg = defaultdict(lambda: [0, 0, None])
    for r in step_rows:
        k = (r["dataset"], r["mode"])
        agg[k][0] += 1
        agg[k][1] += int(f(r.get("s0_late")) or 0)
        got, ref = f(r.get("s0_toi")), f(r.get("gt_earliest"))
        if got is not None and ref is not None and ref < 1:
            d = got - ref
            if agg[k][2] is None or d > agg[k][2]:
                agg[k][2] = d
    for (scene, m), (n, late, worst) in sorted(agg.items()):
        print(f"  {scene:<20} {m:<22} {n:>6} {late:>5} {fmt(worst, 14, 6)}")
