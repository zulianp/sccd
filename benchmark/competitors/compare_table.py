#!/usr/bin/env python3
"""The two competitor comparisons, one table each.

Usage: compare_table.py <comparison.csv> [<more.csv> ...]

Each competitor is compared on the question it answers.

* **Earliest time of impact: SCCD against Scalable CCD.** Both run their
  earliest-time-of-impact path per case -- one query kind of one simulation step,
  from a bound of 1, over identical broad-phase candidates. `narrow ms` is SCCD's
  `narrow_ms` (`ToiOutput::Earliest`) and Scalable CCD's default build. The error
  is per step and signed: the minimum of the step's cases minus the earliest
  exact root, scored for every repeat. Negative is conservative; `err worst` is
  the largest signed error and `late` counts steps where it is positive, per
  pass.

* **Per collision pair: SCCD against Additive CCD.** `narrow ms` is the per-pair
  path over the broad-phase candidates: SCCD's `narrow_ms_s1`
  (`ToiOutput::PerPair`) and ACCD's `narrow_ms` over SCCD's host candidates.
  `fp`, `fn` and the time-of-impact error score every curated query at its own
  coordinates, where the exact roots were computed. The error is signed, reported
  time of impact minus exact root, so negative is conservative: `err med` is the
  median, `worst late` the largest positive one (`0` means no query was late) and
  `worst early` the most negative, the largest gap on the safe side. `fp`, `fn` and `late` are counts
  per pass.

Timings are per case in milliseconds, median / maximum. `ns/q` is narrow-phase
time over the candidates it was handed, summed over the run.
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


def step_of(case):
    d = "".join(c for c in case if c.isdigit())
    return int(d) if d else -1


EARLIEST = [
    ("SCCD host Relaxed", "relaxed"),
    ("SCCD host Tight", "tight"),
    ("SCCD device Relaxed", "device-relaxed"),
    ("SCCD device Tight", "device-tight"),
    ("Scalable CCD device", "scalable-ccd-device"),
    ("Scalable CCD host", "scalable-ccd-host"),
]

PER_PAIR = [
    # (name, mode, narrow-phase column)
    ("SCCD host Relaxed", "relaxed", "narrow_ms_s1"),
    ("SCCD host Tight", "tight", "narrow_ms_s1"),
    ("SCCD device Relaxed", "device-relaxed", "narrow_ms_s1"),
    ("SCCD device Tight", "device-tight", "narrow_ms_s1"),
    ("ACCD host", "accd", "narrow_ms"),
]


def g(v, w, p=3):
    return f"{v:{w}.{p}g}" if v is not None else " " * (w - 1) + "-"


def passes_of(rows):
    """Repeats in a mode's rows: how many times each case appears."""
    cases = {r["case"] for r in rows}
    return max(1, len(rows) // max(1, len(cases)))


def timing(rows, narrow_col):
    pr = [f(r["prep_ms"]) or 0.0 for r in rows]
    b = [f(r["broad_ms"]) or 0.0 for r in rows]
    n = [f(r[narrow_col]) or 0.0 for r in rows]
    tot = [x + y + z for x, y, z in zip(pr, b, n)]
    nq = sum(f(r["queries"]) or 0 for r in rows)
    nsq = (sum(n) / nq * 1e6) if nq and sum(n) > 0 else None
    return (f"{st.median(pr):>8.3f} | {st.median(b):>7.3f} /{max(b):>8.2f} | "
            f"{st.median(n):>7.3f} /{max(n):>8.2f} | {st.median(tot):>8.3f} | {g(nsq, 8, 4)}")


def table_header(cols):
    line = "| " + " | ".join(cols) + " |"
    print(line)
    print("|" + "|".join("-" * (len(c) + 2) for c in cols) + "|")


def earliest_table(scene, rows):
    # Earliest time of impact per step per repeat. A row's repeat is how many
    # times its (mode, case) has been seen, so a library whose answer varies
    # between runs is scored on every run rather than its best one.
    step_ref = {}
    got = defaultdict(dict)  # mode -> (step, repeat) -> earliest
    seen = defaultdict(int)
    for r in rows:
        s = step_of(r["case"])
        if s < 0:
            continue
        ref = f(r.get("gt_earliest"))
        if ref is not None and ref < 1:
            step_ref[s] = min(step_ref.get(s, 1.0), ref)
        k = (r["mode"], r["case"])
        rep = seen[k]
        seen[k] += 1
        if r["mode"] == "scalable-ccd-host":
            continue  # no narrow phase
        t = f(r.get("s0_toi"))
        if t is None:
            continue
        key = (s, rep)
        prev = got[r["mode"]].get(key)
        got[r["mode"]][key] = t if prev is None else min(prev, t)

    print(f"\n### {scene} -- earliest time of impact\n")
    table_header([f"{'library':<21}", f"{'prep ms':>8}", f"{'broad ms':>16}", f"{'narrow ms':>16}",
                  f"{'total':>8}", f"{'ns/q':>8}", f"{'err med':>11}", f"{'err worst':>13}",
                  f"{'late':>4}"])
    for name, mode in EARLIEST:
        mr = [r for r in rows if r["mode"] == mode]
        if not mr:
            continue
        errs, late = [], 0
        for (s, _), t in got.get(mode, {}).items():
            if s in step_ref:
                errs.append(t - step_ref[s])
                late += t > step_ref[s]
        reps = max(1, len({rep for (_, rep) in got.get(mode, {})}))
        print(f"| {name:<21} | {timing(mr, 'narrow_ms')} | "
              f"{g(st.median(errs) if errs else None, 11)} | {g(max(errs) if errs else None, 13)} | "
              f"{round(late / reps) if errs else '-':>4} |")
    print(f"\n  steps with a known contact: {len(step_ref)}")


def per_pair_table(scene, rows):
    print(f"\n### {scene} -- per collision pair\n")
    table_header([f"{'library':<21}", f"{'prep ms':>8}", f"{'broad ms':>16}", f"{'narrow ms':>16}",
                  f"{'total':>8}", f"{'ns/q':>8}", f"{'fp':>5}", f"{'fn':>5}", f"{'late':>4}",
                  f"{'err med':>10}", f"{'worst late':>10}", f"{'worst early':>11}"])
    queries = None
    for name, mode, col in PER_PAIR:
        mr = [r for r in rows if r["mode"] == mode]
        if not mr:
            continue
        p = passes_of(mr)
        fp = sum(f(r["fp"]) or 0 for r in mr) / p
        fn = sum(f(r["fn"]) or 0 for r in mr) / p
        late = sum(f(r["toi_late"]) or 0 for r in mr) / p
        if queries is None:
            queries = sum((f(r["toi_n"]) or 0) + (f(r["fn"]) or 0) for r in mr) / p
        scored = [r for r in mr if f(r["toi_n"])]
        # Signed: the harness records how far each answer fell before its root
        # (early) and after it (late), so the error is one negated and the other
        # as it stands. Negative is conservative.
        med = -st.median([f(r["toi_med_early"]) for r in scored]) if scored else None
        worst_late = max([f(r["toi_max_late"]) for r in scored], default=None)
        mx = -max([f(r["toi_max_early"]) for r in scored], default=0.0) if scored else None
        print(f"| {name:<21} | {timing(mr, col)} | {round(fp):>5} | {round(fn):>5} | {round(late):>4} | "
              f"{g(med, 10)} | {g(worst_late, 10)} | {g(mx, 11)} |")
    if queries is not None:
        print(f"\n  curated queries with a contact per pass: {round(queries)}")


def main(paths):
    rows = []
    for p in paths:
        rows += [r for r in csv.DictReader(open(p)) if r.get("type") in ("vf", "ee")]
    for scene in sorted({r["dataset"] for r in rows}):
        mine = [r for r in rows if r["dataset"] == scene]
        earliest_table(scene, mine)
        per_pair_table(scene, mine)
    print("\n  Time-of-impact error is signed: reported minus exact root, so negative is")
    print("  conservative and a positive value is a contact reported after the true one.")
    print("  Counts are per pass over the case list. The mesh path and the exact roots")
    print("  describe the same geometry here, because smesh is built with double geom_t;")
    print("  a float32 mesh answers about different coordinates and makes correct kernels")
    print("  look late.")


if __name__ == "__main__":
    main(sys.argv[1:])
