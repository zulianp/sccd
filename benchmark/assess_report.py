#!/usr/bin/env python3
"""Turn assessment.csv into the tables wip/ASSESSMENT.md cites.

Every keep-or-demote call in the cleanup has to point at a row of the CSV, so
this deliberately does very little: it takes the median over repeats, puts
variants of one component side by side, and flags the two things that would make
a comparison meaningless.

The median rather than the mean because the failure mode here is a single slow
run -- another job landing on the node, a page fault storm on first touch -- and
one such run moves a three-sample mean by more than most of the differences
being judged.

Two checks run over every comparison, and both exist because this project got
caught by them:

  * **Pair counts must agree across the variants of a component.** A broad phase
    that finds a different pair set is not faster, it is wrong, and a table of
    milliseconds alone cannot tell those two apart.
  * **A difference smaller than the spread between repeats is not a result.**
    Several conclusions in this work were retracted after being read off
    differences that the noise covered, so a comparison whose gap is inside the
    observed spread is printed as "inside noise" rather than as a ratio.

Usage:
    python3 benchmark/assess_report.py benchmark/assessment/assessment.csv
"""

import csv
import sys
from collections import defaultdict
from statistics import median


def load(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                r["ms"] = float(r["ms"])
                r["pairs"] = int(r["pairs"])
            except (ValueError, KeyError):
                continue
            rows.append(r)
    return rows


def group(rows):
    """(hardware, scene, level, component, phase) -> variant -> [rows]."""
    out = defaultdict(lambda: defaultdict(list))
    for r in rows:
        key = (r["hardware"], r["scene"], r["level"], r["component"], r["phase"])
        out[key][r["variant"]].append(r)
    return out


def spread(values):
    """Observed run-to-run spread, as a fraction of the median."""
    if len(values) < 2:
        return None
    m = median(values)
    if m <= 0:
        return None
    return (max(values) - min(values)) / m


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1

    rows = load(argv[1])
    if not rows:
        print(f"no usable rows in {argv[1]}", file=sys.stderr)
        return 1

    grouped = group(rows)
    problems = []

    print("# Assessment\n")
    print(f"Source: `{argv[1]}` — {len(rows)} rows, "
          f"{len({r['repeat'] for r in rows})} repeats.\n")

    for key in sorted(grouped):
        hardware, scene, level, component, phase = key
        variants = grouped[key]
        if len(variants) < 2:
            continue

        lvl = "" if level == "-" else f" level {level}"
        print(f"## {component} — {scene}{lvl}, {phase}, {hardware}\n")
        print("| variant | median ms | spread | pairs | notes |")
        print("|---|---:|---:|---:|---|")

        medians = {}
        pair_counts = set()
        spreads = []
        for variant in sorted(variants):
            vrows = variants[variant]
            times = [r["ms"] for r in vrows]
            med = median(times)
            medians[variant] = med
            sp = spread(times)
            if sp is not None:
                spreads.append(sp)
            pairs = vrows[0]["pairs"]
            pair_counts.add(pairs)
            sp_txt = "-" if sp is None else f"{sp * 100:.0f}%"
            print(f"| {variant} | {med:.2f} | {sp_txt} | {pairs} | "
                  f"{vrows[0].get('notes', '')} |")
        print()

        # A component whose variants disagree on the pair set is not being
        # compared on equal terms; that is a correctness report, not a timing one.
        if len(pair_counts) > 1:
            msg = (f"{component}/{scene}{lvl}/{phase}/{hardware}: variants "
                   f"disagree on the pair count {sorted(pair_counts)}")
            problems.append(msg)
            print(f"> **Pair sets differ across variants** ({sorted(pair_counts)}). "
                  f"These variants are not computing the same thing, so the times "
                  f"below are not comparable and one of them is wrong.\n")

        best = min(medians, key=medians.get)
        worst = max(medians, key=medians.get)
        if medians[best] > 0:
            ratio = medians[worst] / medians[best]
            gap = (medians[worst] - medians[best]) / medians[worst]
            if not spreads:
                # One repeat gives no estimate of the noise, and this harness is
                # noisy enough that a single pair of numbers separates nothing.
                print(f"**{best} vs {worst}: single repeat, not separable.** "
                      f"{medians[best]:.2f} vs {medians[worst]:.2f} ms with no "
                      f"spread to compare against; re-run with ASSESS_REPEATS>1 "
                      f"before reading anything into this.\n")
                continue
            noise = max(spreads)
            if gap <= noise:
                print(f"**{best} vs {worst}: inside noise.** The gap is "
                      f"{gap * 100:.0f}% and the observed run-to-run spread is "
                      f"{noise * 100:.0f}%, so this does not separate them.\n")
            else:
                print(f"**{best} is {ratio:.2f}x faster than {worst}** "
                      f"({medians[best]:.2f} vs {medians[worst]:.2f} ms; "
                      f"spread {noise * 100:.0f}%).\n")

    # Pair counts must also agree ACROSS hardware for the same scene and phase.
    # Grouping above puts host and device in different buckets, so a device broad
    # phase that finds fewer pairs than the host would never be compared with it
    # -- and that is the direction that matters, because a pair the device does
    # not report is a collision the narrow phase never gets the chance to find.
    cross = defaultdict(lambda: defaultdict(set))
    for r in rows:
        cross[(r["scene"], r["level"], r["phase"])][r["hardware"]].add(r["pairs"])
    for (scene, level, phase), by_hw in sorted(cross.items()):
        counts = {hw: sorted(v) for hw, v in by_hw.items() if v}
        if len(counts) < 2:
            continue
        flat = {c for v in counts.values() for c in v}
        if len(flat) > 1:
            problems.append(
                f"PAIR COUNT DIFFERS ACROSS HARDWARE: {scene} level {level} "
                f"{phase}: " + ", ".join(f"{hw}={v}" for hw, v in sorted(counts.items()))
                + " -- a pair one side does not report is a collision its narrow "
                  "phase never sees, so this needs explaining before any timing "
                  "on these rows is used")

    # Missed collisions are illegal; report them regardless of any timing.
    for r in rows:
        notes = r.get("notes", "")
        if "fn=" in notes:
            fn = notes.split("fn=")[1].split(";")[0]
            if fn.isdigit() and int(fn) > 0:
                problems.append(
                    f"MISSED COLLISIONS: {r['hardware']}/{r['scene']}/"
                    f"{r['component']}/{r['variant']} reports fn={fn}")

    if problems:
        print("## Problems\n")
        for p in sorted(set(problems)):
            print(f"- {p}")
        print()
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
