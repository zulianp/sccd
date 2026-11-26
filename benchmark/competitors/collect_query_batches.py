#!/usr/bin/env python3
"""Collect the narrow-phase batch timings on the query-only sets into one CSV.

    collect_query_batches.py <narrow-repeats-dir> > query-batches.csv

<narrow-repeats-dir> holds one directory per set, each with the logs of three
repeats run as separate processes (narrow-repeats.sbatch on Fir):

  sccd-<r>.log                ti_oracle --bench 1 --bench-device-only
  scalable-earliest-<r>.log   scalable_ccd_bench --queries, default build
  scalable-perquery-<r>.log   scalable_ccd_bench --queries, SCALABLE_CCD_TOI_PER_QUERY

Every timing is one batch per query type: all of a set's vertex-face queries in
one call, then all of its edge-edge queries, after one untimed warm-up. Each row
is one repeat; the median and maximum over repeats are left to the reader.

Columns: set, type (vf|ee), library, answer (earliest|per-pair), repeat, ms,
calls, failed_call (the call size Scalable CCD failed to allocate for, if any),
earliest_toi, collisions (the per-query build's returned collisions).
"""
import csv
import re
import sys
from pathlib import Path

FIELDS = ["set", "type", "library", "answer", "repeat", "ms", "calls", "failed_call",
          "earliest_toi", "collisions"]

SCCD_HEAD = re.compile(r"^(\S+) / (VF|EE)\s+throughput:")
SCCD_ROW = re.compile(r"^(device-relaxed|device-tight)\s+([\d.]+)\s+([\d.]+)\s+[\d.]+x\s+earliest toi (\S+)")
SC_FAIL = re.compile(r"^(\S+) / (VF|EE) batch: their narrow phase failed .* on a call of (\d+) queries")
SC_OK = re.compile(r"^(\S+) / (VF|EE) batch: (\d+) queries in (\d+) call\(s\) of at most (\d+)"
                   r"(?: \(after failing on calls of (\d+)\))?, (earliest|per-query) build, "
                   r"(?:earliest toi (\S+)|(\d+) collisions returned), best ([\d.]+) ms")


def sccd_rows(path, repeat):
    kind = None
    for line in path.read_text().splitlines():
        m = SCCD_HEAD.match(line)
        if m:
            kind = m.group(2).lower()
            continue
        m = SCCD_ROW.match(line)
        if m and kind:
            mode, per_pair, earliest, toi = m.groups()
            lib = "sccd-" + mode
            yield dict(type=kind, library=lib, answer="per-pair", repeat=repeat, ms=per_pair)
            yield dict(type=kind, library=lib, answer="earliest", repeat=repeat, ms=earliest,
                       earliest_toi=toi)


def scalable_rows(path, repeat):
    failed = {}
    for line in path.read_text().splitlines():
        m = SC_FAIL.match(line)
        if m:
            kind = m.group(2).lower()
            failed[kind] = max(failed.get(kind, 0), int(m.group(3)))
            continue
        m = SC_OK.match(line)
        if m:
            kind = m.group(2).lower()
            build = m.group(7)
            yield dict(type=kind, library="scalable-ccd",
                       answer="earliest" if build == "earliest" else "per-pair",
                       repeat=repeat, ms=m.group(10), calls=m.group(4),
                       failed_call=m.group(6) or "", earliest_toi=m.group(8) or "",
                       collisions=m.group(9) or "")


def main(argv):
    if len(argv) != 2:
        sys.exit(__doc__)
    root = Path(argv[1])
    out = csv.DictWriter(sys.stdout, fieldnames=FIELDS)
    out.writeheader()
    n = 0
    for set_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        for log in sorted(set_dir.glob("*.log")):
            m = re.match(r"(sccd|scalable-earliest|scalable-perquery)-(\d+)\.log$", log.name)
            if not m:
                continue
            gen = sccd_rows if m.group(1) == "sccd" else scalable_rows
            for row in gen(log, int(m.group(2))):
                out.writerow({"set": set_dir.name, **row})
                n += 1
    print(f"{n} rows", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
