#!/usr/bin/env python3
"""The two competitor tables, from the comparison CSV.

Each competitor answers a different question, so it gets its own table rather
than a shared one:

* Scalable CCD returns one earliest time of impact per simulation step, so both
  libraries run their earliest-time-of-impact path per case from a fresh bound
  and are scored on that number.
* Additive CCD answers per pair, so both narrow phases are timed over the same
  broad-phase candidates and scored on the curated queries at the coordinates
  their exact roots were computed for.

Usage: make_competitors.py [--check]
"""
import csv
import gzip
import pathlib
import statistics as st
import sys
from collections import defaultdict

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
DATA = ROOT / "benchmark" / "competitors" / "results"
OUT = HERE.parent / "generated" / "tables"

SCENES = ["armadillo-rollers", "cloth-ball", "cloth-funnel"]
PRETTY = {"armadillo-rollers": "armadillo-rollers", "cloth-ball": "cloth-ball",
          "cloth-funnel": "cloth-funnel"}


def rows():
    """Every per-case row of the newest comparison CSV, gzipped or not."""
    candidates = sorted(DATA.glob("compare-gh200-full-*.csv.gz")) + \
        sorted(DATA.glob("compare-gh200-full-*.csv"))
    if not candidates:
        sys.exit(f"error: no comparison CSV under {DATA}")
    path = candidates[-1]
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        data = [r for r in csv.DictReader(fh) if r.get("type") in ("vf", "ee")]
    return path, data


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def step_of(case):
    d = "".join(c for c in case if c.isdigit())
    return int(d) if d else -1


def passes(rs):
    """Repeats in a mode's rows: how many times each case appears."""
    return max(1, len(rs) // max(1, len({r["case"] for r in rs})))


def fmt(v, digits=3):
    if v is None:
        return "--"
    if v == 0:
        return "0"
    return f"{v:.{digits}g}"


def earliest_table(data):
    """Per step: the earliest time of impact each library reports."""
    lines = []
    for scene in SCENES:
        mine = [r for r in data if r["dataset"] == scene]
        # The earliest exact root of each step, and each library's answer for it,
        # kept per repeat: a library whose answer varies between runs is scored
        # on every run rather than on its best one.
        ref, got, seen = {}, defaultdict(dict), defaultdict(int)
        for r in mine:
            s = step_of(r["case"])
            if s < 0:
                continue
            g = num(r.get("gt_earliest"))
            if g is not None and g < 1:
                ref[s] = min(ref.get(s, 1.0), g)
            k = (r["mode"], r["case"])
            rep = seen[k]
            seen[k] += 1
            if r["mode"] == "scalable-ccd-host":
                continue
            t = num(r.get("s0_toi"))
            if t is None:
                continue
            prev = got[r["mode"]].get((s, rep))
            got[r["mode"]][(s, rep)] = t if prev is None else min(prev, t)

        for label, mode in (("SCCD (CPU)", "tight"), ("SCCD (GPU)", "device-tight"),
                            ("Scalable CCD (GPU)", "scalable-ccd-device")):
            mr = [r for r in mine if r["mode"] == mode]
            if not mr:
                continue
            broad = st.median([(num(r["prep_ms"]) or 0) + (num(r["broad_ms"]) or 0) for r in mr])
            narrow = st.median([num(r["narrow_ms"]) or 0 for r in mr])
            cand = sum(num(r["queries"]) or 0 for r in mr)
            nsq = sum(num(r["narrow_ms"]) or 0 for r in mr) / cand * 1e6 if cand else None
            errs = [t - ref[s] for (s, _), t in got[mode].items() if s in ref]
            late = sum(1 for e in errs if e > 0)
            reps = max(1, len({rep for (_, rep) in got[mode]}))
            lines.append(
                f"    {PRETTY[scene] if label.startswith('SCCD (CPU)') else ''} & {label} & "
                f"{broad:.1f} & {narrow:.1f} & {broad + narrow:.1f} & {fmt(nsq, 3)} & "
                f"{fmt(st.median(errs) if errs else None, 3)} & {fmt(max(errs) if errs else None, 3)} & "
                f"{round(late / reps)} \\\\")
        lines.append("    \\midrule")
    if lines and lines[-1].strip() == "\\midrule":
        lines.pop()

    return f"""\\begin{{table}}[htbp]
  \\centering
  \\caption{{Earliest time of impact per simulation step, against Scalable
    CCD~\\citep{{belgrod2025toi}}. Both libraries run their
    earliest-time-of-impact path once per case from a bound of $1$, over
    identical broad-phase candidates, and a step's answer is the minimum over
    its cases. \\emph{{broad}} is the acceleration structure and the traversal
    together and \\emph{{narrow}} the root finding, both medians per case in
    milliseconds; \\emph{{ns/cand.}} is narrow-phase time over the candidates it
    was handed. The error is signed, reported minus exact root, so negative is
    conservative; \\emph{{late}} counts steps where it is positive, per pass over
    the case list.}}
  \\label{{tab:competitor-earliest}}
  \\fittable{{%
\\begin{{tabular}}{{llrrrrrrr}}
    \\toprule
    scene & library & broad (ms) & narrow (ms) & total (ms) & ns/cand. & err.\\ med. & err.\\ worst & late \\\\
    \\midrule
{chr(10).join(lines)}
    \\bottomrule
  \\end{{tabular}}}}
  \\par\\smallskip\\footnotesize Scalable CCD's narrow phase is CUDA only, so it has
  no CPU row; its host broad phase is discussed in the text.
  \\par\\smallskip\\footnotesize Source: \\texttt{{benchmark/competitors/results/}}
\\end{{table}}
"""


def pair_table(data):
    """Per curated query: what each library reports for the same pair."""
    lines = []
    for scene in SCENES:
        mine = [r for r in data if r["dataset"] == scene]
        for label, mode, col in (("SCCD (CPU)", "tight", "narrow_ms_s1"),
                                 ("SCCD (GPU)", "device-tight", "narrow_ms_s1"),
                                 ("Additive CCD (CPU)", "accd", "narrow_ms")):
            mr = [r for r in mine if r["mode"] == mode]
            if not mr:
                continue
            p = passes(mr)
            narrow = st.median([num(r[col]) or 0 for r in mr])
            cand = sum(num(r["queries"]) or 0 for r in mr)
            nsq = sum(num(r[col]) or 0 for r in mr) / cand * 1e6 if cand else None
            fp = sum(num(r["fp"]) or 0 for r in mr) / p
            fn = sum(num(r["fn"]) or 0 for r in mr) / p
            late = sum(num(r["toi_late"]) or 0 for r in mr) / p
            scored = [r for r in mr if num(r["toi_n"])]
            med = -st.median([num(r["toi_med_early"]) for r in scored]) if scored else None
            worst = -max([num(r["toi_max_early"]) for r in scored], default=0.0) if scored else None
            lines.append(
                f"    {PRETTY[scene] if label.startswith('SCCD (CPU)') else ''} & {label} & "
                f"{narrow:.2f} & {fmt(nsq, 3)} & {round(fp)} & {round(fn)} & {round(late)} & "
                f"{fmt(med, 3)} & {fmt(worst, 3)} \\\\")
        lines.append("    \\midrule")
    if lines and lines[-1].strip() == "\\midrule":
        lines.pop()

    return f"""\\begin{{table}}[htbp]
  \\centering
  \\caption{{Per collision pair, against additive CCD~\\citep{{li2021codim}} as the
    IPC toolkit implements it. Both narrow phases are timed over the same
    broad-phase candidates -- additive CCD answers one pair at a time and carries
    no parallelism of its own, so it is driven by the parallel loop the toolkit's
    own stepsize query uses, on the same $72$ threads -- and both are scored on
    every curated query at the coordinates its exact root was computed for.
    \\emph{{f.p.}}, \\emph{{missed}} and \\emph{{late}} are counts per pass over the
    case list. The error is signed, so negative is conservative.}}
  \\label{{tab:competitor-pair}}
  \\fittable{{%
\\begin{{tabular}}{{llrrrrrrr}}
    \\toprule
    scene & library & narrow (ms) & ns/pair & f.p. & missed & late & err.\\ med. & err.\\ worst \\\\
    \\midrule
{chr(10).join(lines)}
    \\bottomrule
  \\end{{tabular}}}}
  \\par\\smallskip\\footnotesize Source: \\texttt{{benchmark/competitors/results/}}
\\end{{table}}
"""


def main():
    path, data = rows()
    OUT.mkdir(parents=True, exist_ok=True)
    wrote = []
    for name, text in (("tab-competitor-earliest", earliest_table(data)),
                       ("tab-competitor-pair", pair_table(data))):
        (OUT / f"{name}.tex").write_text(text)
        wrote.append(name)
    print(f"wrote {', '.join(wrote)} from {path.name}")


if __name__ == "__main__":
    main()
