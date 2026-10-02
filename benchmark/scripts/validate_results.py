#!/usr/bin/env python3
"""Refuse a result file whose provenance is not what was asked for.

Row counts and the conservativeness gates cannot catch a run that measured the
wrong thing: every answer in it is correct and only the conditions were wrong.
Four such runs were published or nearly published in one week -- a strategy name
the build did not know and so raced instead, a submit path that lost the thread
binding, another that lost the CSV header so the merge ate a row per chunk, and a
scaling study against a synthesised mesh instead of the dataset's.

This checks the three things that distinguish those from a good run.

    validate_results.py timings NEW.csv [--expect-bp NAME] [--control OLD.csv]
    validate_results.py scaling NEW.txt [--expect-bp NAME] [--expect-mesh SUBSTR]
                                        [--expect-topology TRI3] [--reference OLD.txt]

Exits non-zero on any failure, so it can gate a regeneration.
"""
import argparse
import collections
import csv
import gzip
import statistics as st
import sys

# The query types a scene's cases divide into. A run missing one of them whole
# is the failure mode that hid behind a zero in a published table.
TYPES = ("ee", "vf")

# The modes SCCD's own pipeline reports. A comparison CSV also carries a mode per
# competitor library, and those have their own case coverage and their own misses;
# counting them here would both misfire and attribute another library's failure to
# SCCD, which its own rows must never carry.
SCCD_MODES = ("tight", "relaxed", "device-tight", "device-relaxed")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load(path):
    with _open(path) as fh:
        return list(csv.DictReader(fh))


def fail(msgs):
    for m in msgs:
        print(f"FAIL: {m}")
    return 1


def check_timings(args):
    rows = load(args.csv)
    bad = []
    print(f"{len(rows)} rows from {args.csv}")

    # 1. Provenance: the strategy the file claims, against the one asked for.
    bps = collections.Counter(r.get("broadphase", "") for r in rows)
    print(f"broadphase column: {dict(bps)}")
    if "auto" in bps:
        bad.append(f"{bps['auto']} rows recorded broadphase=auto, which means the "
                   "name asked for was not one the build knew and it raced instead")
    if args.expect_bp and args.expect_bp not in bps:
        bad.append(f"no row records broadphase={args.expect_bp}; the run did not "
                   f"measure what was asked for (saw {sorted(bps)})")

    # 2. Completeness by category. A mode or a query type short of the others is
    #    a dropped slice, not a property of that mode.
    ours = [r for r in rows if r["mode"] in SCCD_MODES]
    other = len(rows) - len(ours)
    if other:
        print(f"({other} rows belong to another library and are not checked here)")
    per_mode = collections.Counter(r["mode"] for r in ours)
    print(f"rows per mode: {dict(per_mode)}")
    if len(set(per_mode.values())) > 1:
        bad.append(f"modes hold different row counts {dict(per_mode)}; a slice was "
                   "dropped, and the merge eating one row per chunk looks exactly "
                   "like this")
    for scene in sorted({r["dataset"] for r in ours}):
        for mode in sorted({r["mode"] for r in ours}):
            seen = {t for r in ours
                    if r["dataset"] == scene and r["mode"] == mode for t in [r["type"]]}
            missing = set(TYPES) - seen
            if missing and seen:
                bad.append(f"{scene} {mode} has no {'/'.join(sorted(missing))} rows "
                           "at all; a whole query type is absent")

    # 3. Conservativeness, which is a hard invariant rather than a statistic.
    fatal = [r for r in ours if any(int(r.get(k, 0) or 0)
                                   for k in ("fn", "toi_late", "broad_fn", "s0_late"))]
    print(f"missed collisions / late times of impact: {len(fatal)}")
    if fatal:
        bad.append(f"{len(fatal)} rows report a missed pair or a late time of impact")

    # 4. An untouched control. The sweep broad phase shares the file and no
    #    cell-list change can affect it, so if its times moved the run measured a
    #    different machine and nothing in it is comparable to the old numbers.
    if args.control:
        old = load(args.control)

        def totals(rs):
            runs = collections.defaultdict(list)
            for r in rs:
                if r.get("broadphase") != "sweep":
                    continue
                runs[(r["dataset"], r["mode"], r["case"], r["type"])].append(
                    float(r["prep_ms"]) + float(r["broad_ms"]))
            out = collections.defaultdict(float)
            for (d, m, _c, t), v in runs.items():
                out[(d, m, t)] += st.median(v)
            return out

        a, b = totals(old), totals(rows)
        shared = sorted(set(a) & set(b))
        if not shared:
            bad.append("no sweep control series shared with the reference file")
        worst, worst_key = 0.0, None
        for k in shared:
            if a[k] <= 0:
                continue
            dev = abs(b[k] / a[k] - 1.0)
            if dev > worst:
                worst, worst_key = dev, k
        print(f"sweep control: {len(shared)} series, largest deviation "
              f"{worst:.0%} ({worst_key})")
        if worst > args.tolerance:
            bad.append(f"the sweep control moved {worst:.0%} at {worst_key}, over the "
                       f"{args.tolerance:.0%} tolerance; this run measured a "
                       "different machine than the reference")

    return fail(bad) if bad else 0


def check_scaling(args):
    with open(args.txt) as fh:
        lines = fh.read().splitlines()
    if not lines or not lines[0].startswith("#"):
        return fail([f"{args.txt} has no provenance header"])
    head = lines[0]
    print(f"header: {head[:160]}")
    bad = []

    # The scaling driver synthesises a mesh of its own when none is passed and
    # says so only here, in a line that is easy not to read.
    for want, label in ((args.expect_bp, "broadphase"),
                        (args.expect_topology, "base_topology")):
        if want and f"{label}={want}" not in head:
            bad.append(f"header does not say {label}={want}")
    if args.expect_mesh and args.expect_mesh not in head:
        bad.append(f"header does not name the mesh {args.expect_mesh!r}; "
                   "`t0=generated t1=synthesized` means the driver invented one")
    if "t0=generated" in head or "t1=synthesized" in head:
        bad.append("the driver synthesised its own mesh, so the levels do not "
                   "describe the dataset")

    data = [l.split() for l in lines[2:] if l.strip()]
    print(f"levels: {[(d[0], d[1]) for d in data]}")
    if args.reference:
        with open(args.reference) as fh:
            ref = [l.split() for l in fh.read().splitlines()[2:] if l.strip()]
        got = {d[0]: d[1] for d in data}
        exp = {d[0]: d[1] for d in ref}
        if set(exp) - set(got):
            bad.append(f"levels missing against the reference: "
                       f"{sorted(set(exp) - set(got))}")
        for lvl, faces in exp.items():
            if lvl in got and got[lvl] != faces:
                bad.append(f"level {lvl} has {got[lvl]} faces where the reference "
                           f"has {faces}; this is different geometry")
    return fail(bad) if bad else 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="kind", required=True)

    t = sub.add_parser("timings", help="a per-case benchmark CSV")
    t.add_argument("csv")
    t.add_argument("--expect-bp", help="the broad phase the run was asked for")
    t.add_argument("--control", help="the previous CSV, for the sweep control")
    t.add_argument("--tolerance", type=float, default=0.15)

    s = sub.add_parser("scaling", help="a refine_scaling table")
    s.add_argument("txt")
    s.add_argument("--expect-bp")
    s.add_argument("--expect-mesh")
    s.add_argument("--expect-topology")
    s.add_argument("--reference", help="the previous table, for levels and faces")

    args = p.parse_args()
    rc = check_timings(args) if args.kind == "timings" else check_scaling(args)
    print("PASS" if rc == 0 else "REFUSED")
    return rc


if __name__ == "__main__":
    sys.exit(main())
