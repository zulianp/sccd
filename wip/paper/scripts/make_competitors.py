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

Both tables report earliness, the exact root minus the time of impact reported
for it, so a positive value is conservative and a negative one is an answer
after the root.

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

SCENES = ["armadillo-rollers", "cloth-ball", "cloth-funnel", "n-body-simulation",
          "puffer-ball", "rod-twist"]
# The name the paper's other tables use for each scene.
PRETTY = {"armadillo-rollers": "armadillo-rollers", "cloth-ball": "cloth-ball",
          "cloth-funnel": "cloth-funnel", "n-body-simulation": "n-body",
          "puffer-ball": "puffer-ball", "rod-twist": "rod-twist"}


def rows():
    """Every per-case row of the newest comparison CSV, gzipped or not."""
    # Newest by name, whichever form it is in. Sorting the two globs separately
    # and taking the last would prefer an ungzipped leftover over a newer
    # archive, which is how a stale file got read once.
    candidates = sorted(DATA.glob("compare-gh200-full-*.csv.gz")) + \
        sorted(DATA.glob("compare-gh200-full-*.csv"))
    if not candidates:
        sys.exit(f"error: no comparison CSV under {DATA}")
    path = max(candidates, key=lambda p: p.name.replace(".csv.gz", "").replace(".csv", ""))
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        data = [r for r in csv.DictReader(fh) if r.get("type") in ("vf", "ee")]

    # The comparison measures our device over both broad-phase strategies, so
    # that they and Scalable CCD's sit in one allocation. Only one of them
    # belongs in these tables, or a per-case median would be taken over two
    # different runs of the same case and describe neither. The cell list is the
    # faster of the two on every scene of the device, so it is the one a caller
    # gets and the one the comparison reports.
    kept = [r for r in data
            if not r["mode"].startswith("device-") or r.get("broadphase") == "cell2d"]
    dropped = len(data) - len(kept)
    if dropped:
        print(f"  ({dropped} device rows of the other strategy left out)")
    return path, kept


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


def whole_scene(rs, cols):
    """Whole-scene total in ms: per case the median over repeats, summed.

    A median over repeats first keeps one slow case on one pass from moving the
    scene's total, which is the convention the paper's other timing tables use.
    """
    per_case = defaultdict(list)
    for r in rs:
        per_case[r["case"]].append(sum(num(r[c]) or 0.0 for c in cols))
    return sum(st.median(v) for v in per_case.values()), len(per_case)


def fmt(v, digits=3):
    if v is None:
        return "--"
    if v == 0:
        return "0"
    return f"{v:.{digits}g}"


def ms(v):
    """A millisecond total, grouped in thousands."""
    return f"{v:,.0f}" if v >= 1000 else f"{v:.1f}"


def ratio(ours, theirs):
    """How many times further from the root the competitor's median sits.

    Left blank where there is no such number: where the competitor answers after
    the root, which is a failure rather than a looser answer, and where both
    sides are exactly $0$ because every step of the scene begins in contact.

    Nothing but a number or a blank may be returned. A word here would be read
    as a property of the row it is printed on, which is ours, and the one word
    worth saying about a competitor -- that it is late -- would then libel the
    library that never is. It is said on the competitor's own row instead.
    """
    if ours is None or theirs is None or theirs < 0 or ours == 0 or theirs == 0:
        return "--"
    v = theirs / ours
    return f"{v:,.0f}$\\times$" if v >= 100 else f"{v:.1f}$\\times$"


def earliest_table(data):
    """Per step: the earliest time of impact each library reports."""
    lines, touching = [], []
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

        # A step whose earliest exact root is 0 begins with its primitives
        # already touching, so every correct answer for it is 0 and it carries
        # no information about tightness. Say how many of each scene's steps
        # those are, so the column is read for what it measures.
        zero = sum(1 for v in ref.values() if v == 0.0)
        if zero:
            touching.append(f"{PRETTY[scene]} {zero} of {len(ref)}")

        ROWS = (("SCCD (CPU)", "tight"), ("SCCD (GPU)", "device-tight"),
                ("Scalable CCD (GPU)", "scalable-ccd-device"))
        stat = {}
        for _, mode in ROWS:
            mr = [r for r in mine if r["mode"] == mode]
            if not mr:
                continue
            total, ncases = whole_scene(mr, ("prep_ms", "broad_ms", "narrow_ms"))
            # Earliness: the exact root minus the answer, so positive is
            # conservative and a negative value is an answer after the root.
            early = [ref[s] - t for (s, _), t in got[mode].items() if s in ref]
            reps = max(1, len({rep for (_, rep) in got[mode]}))
            stat[mode] = {
                "total": total, "avg": total / ncases,
                "med": st.median(early) if early else None,
                "max": max(early) if early else None,
                "late": round(sum(1 for e in early if e < 0) / reps),
            }

        them = stat.get("scalable-ccd-device")
        for label, mode in ROWS:
            s = stat.get(mode)
            if s is None:
                continue
            if them is None or mode == "scalable-ccd-device":
                # "late" belongs here, on the row of the library it describes,
                # and nowhere else.
                speedup = "1.0$\\times$ (ref)"
                tighter = ("late" if s["med"] is not None and s["med"] < 0
                           else "1.0$\\times$ (ref)")
            else:
                speedup = f"{them['total'] / s['total']:.2f}$\\times$"
                tighter = ratio(s["med"], them["med"])
            lines.append(
                f"    {PRETTY[scene] if label.startswith('SCCD (CPU)') else ''} & {label} & "
                f"{fmt(s['med'], 3)} & {fmt(s['max'], 3)} & {s['late']} & "
                f"{ms(s['total'])} & {s['avg']:.2f} & {speedup} & {tighter} \\\\")
        lines.append("    \\midrule")
    if lines and lines[-1].strip() == "\\midrule":
        lines.pop()

    return f"""\\begin{{table}}[htbp]
  \\centering
  \\caption{{Earliest time of impact per simulation step, against Scalable
    CCD~\\citep{{belgrod2025toi}}. Both libraries run their
    earliest-time-of-impact path once per case from a bound of $1$, over
    identical broad-phase candidates, and a step's answer is the minimum over
    its cases. \\emph{{earliness}} is the step's earliest exact root minus the
    answer reported for it, so a positive value is conservative and a negative
    one is an answer after the root; \\emph{{late}} counts the steps where it is
    negative, per pass over the case list. \\emph{{total}} is the whole scene,
    prep and broad phase and narrow phase together, summed over every case with
    the median over repeats taken first; \\emph{{avg}} divides it by the cases of
    the scene. \\emph{{speedup}} is Scalable CCD's total over ours, so above one
    is our lead, and \\emph{{tighter}} is its median earliness over ours, so
    above one is how many times further from the root its median answer sits.
    Both are quoted on our rows against the reference row, which reads
    \\emph{{late}} on the scene where Scalable CCD's own median answer falls
    after the root; \\emph{{--}} marks a tightness ratio there is no number for.
    Our rows never read \\emph{{late}}: no SCCD answer in this table or the next
    falls after a root.}}
  \\label{{tab:competitor-earliest}}
  \\fittable{{%
\\begin{{tabular}}{{llrrrrrrr}}
    \\toprule
    scene & library & earl.\\ med. & earl.\\ max & late & total (ms) & avg (ms) & speedup & tighter \\\\
    \\midrule
{chr(10).join(lines)}
    \\bottomrule
  \\end{{tabular}}}}
  \\par\\smallskip\\footnotesize Scalable CCD's narrow phase is CUDA only, so it has
  no CPU row; its host broad phase is discussed in the text.
  \\par\\smallskip\\footnotesize Some steps begin with their primitives already
  touching, so their earliest exact root is $0$ and every correct answer for them
  is $0$: {', '.join(touching)}. Those steps hold the earliness columns at $0$
  without saying anything about tightness.
  \\par\\smallskip\\footnotesize Source: \\texttt{{benchmark/competitors/results/}}
\\end{{table}}
"""


def pair_table(data):
    """Per curated query: what each library reports for the same pair."""
    lines, late_total = [], 0
    for scene in SCENES:
        mine = [r for r in data if r["dataset"] == scene]
        ROWS = (("SCCD (CPU)", "tight", "narrow_ms_s1"),
                ("SCCD (GPU)", "device-tight", "narrow_ms_s1"),
                ("Additive CCD (CPU)", "accd", "narrow_ms"))
        stat = {}
        for _, mode, col in ROWS:
            mr = [r for r in mine if r["mode"] == mode]
            if not mr:
                continue
            p = passes(mr)
            total, _ = whole_scene(mr, (col,))
            cand = sum(num(r["queries"]) or 0 for r in mr) / p
            late_total += sum(num(r["toi_late"]) or 0 for r in mr) / p
            scored = [r for r in mr if num(r["toi_n"])]
            stat[mode] = {
                "total": total,
                "nsq": total / cand * 1e6 if cand else None,
                "fp": round(sum(num(r["fp"]) or 0 for r in mr) / p),
                "fn": round(sum(num(r["fn"]) or 0 for r in mr) / p),
                "med": st.median([num(r["toi_med_early"]) for r in scored]) if scored else None,
                "max": max([num(r["toi_max_early"]) for r in scored], default=None) if scored else None,
            }

        them = stat.get("accd")
        for label, mode, _ in ROWS:
            s = stat.get(mode)
            if s is None:
                continue
            if them is None or mode == "accd":
                slowdown = "1.0$\\times$ (ref)"
                tighter = ("late" if s["med"] is not None and s["med"] < 0
                           else "1.0$\\times$ (ref)")
            else:
                slowdown = f"{s['total'] / them['total']:.2f}$\\times$"
                tighter = ratio(s["med"], them["med"])
            lines.append(
                f"    {PRETTY[scene] if label.startswith('SCCD (CPU)') else ''} & {label} & "
                f"{s['fp']} & {s['fn']} & {fmt(s['med'], 3)} & {fmt(s['max'], 3)} & "
                f"{ms(s['total'])} & {fmt(s['nsq'], 3)} & {slowdown} & {tighter} \\\\")
        lines.append("    \\midrule")
    if lines and lines[-1].strip() == "\\midrule":
        lines.pop()

    late_note = ("no query of either library is answered after its root."
                 if late_total == 0 else
                 f"{round(late_total)} queries are answered after their root.")
    return f"""\\begin{{table}}[htbp]
  \\centering
  \\caption{{Per collision pair, against additive CCD~\\citep{{li2021codim}} as the
    IPC toolkit implements it. Both narrow phases are timed over the same
    broad-phase candidates -- additive CCD answers one pair at a time and carries
    no parallelism of its own, so it is driven by the parallel loop the toolkit's
    own stepsize query uses, on the same $72$ threads -- and both are scored on
    every curated query at the coordinates its exact root was computed for.
    \\emph{{f.p.}} and \\emph{{missed}} are counts per pass over the case list.
    \\emph{{earliness}} is the query's exact root minus the time of impact
    reported for it, so a positive value is conservative, and {late_note}
    \\emph{{total}} is the narrow phase over the whole scene, summed over every case with the median
    over repeats taken first, and \\emph{{avg}} divides it by the candidates it
    was handed. \\emph{{slowdown}} is our total over additive CCD's, so above one
    is what the tighter answer costs, and \\emph{{tighter}} is its median
    earliness over ours, so above one is how many times further from the root its
    median answer sits.}}
  \\label{{tab:competitor-pair}}
  \\fittable{{%
\\begin{{tabular}}{{llrrrrrrrr}}
    \\toprule
    scene & library & f.p. & missed & earl.\\ med. & earl.\\ max & total (ms) & avg (ns/pair) & slowdown & tighter \\\\
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
