#!/usr/bin/env python3
"""
Build the query-difficulty table from the committed box counts.

    python3 scripts/make_difficulty.py

Reads benchmark/results/boxes/ and writes generated/tables/tab-difficulty.tex,
the per-query box-count distribution on cloth-funnel for both processors. Like
the other tables, nothing here is typed by hand.

Where the input comes from. Both narrow phases carry a per-query box counter
behind -DSCCD_NP_COUNT_BOXES, and each call prints its distribution as a
twenty-four-bucket log2 histogram -- the same bucketing code on both processors,
so the two columns measure the same thing. The counters put a global atomic on
the device's hot path, so an instrumented build's timings mean nothing and only
its counts are kept. On one GH200 module:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSCCD_ENABLE_CUDA=ON \\
          -DSCCD_ENABLE_SMESH=ON -DSCCD_NP_COUNT_BOXES=ON
    cmake --build build --target sccd_bench
    for space in host device; do
      SCCD_BENCH_EXECUTION_SPACE=$space SCCD_NARROWPHASE_MODE=2 \\
        build/sccd_bench "$DATA" cloth-funnel 2> $space.err
    done

and the two stderr streams become the two CSVs here, one row per narrow-phase
call: `cloth-funnel-boxes.csv` carries each call's histogram and
`cloth-funnel-boxes-totals.csv` its totals.

Two filters are applied, and both matter.

The benchmark runs the first case of a scene once untimed to pay for allocation
and first touch, so the first call of each query type appears twice. The second
occurrence is the timed one; the first is dropped here, which leaves exactly the
scene's 577 cases.

Only the earliest-time-of-impact calls are tabulated. The per-pair calls are a
different measurement -- no bound is shared, so nothing prunes across queries --
and the same run's per-pair rows also include the small curated accuracy sets,
which are different query data. Both are summarised in the prose instead, from
the totals file.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
BOXES = PAPER.parent.parent / "benchmark" / "results" / "boxes"
SCENE = "cloth-funnel"

ETOI = 0  # ToiOutput::Earliest
NBUCKET = 24

# A bucket b holds the queries that examined between 2^b and 2^(b+1)-1 boxes,
# except bucket 0, which holds those that examined one box or none. So the row
# "at least 2^b" is the sum of buckets b upward.
ROWS = [(0, "$0$--$1$"),
        (3, r"$\ge 8$"),
        (6, r"$\ge 64$"),
        (10, r"$\ge 1{,}024$"),
        (14, r"$\ge 16{,}384$"),
        (20, r"$\ge 1{,}048{,}576$")]

PROCS = ["host", "device"]


def load() -> tuple[dict, dict, dict]:
    """Per-processor bucket totals, worst query and case count, EToI calls only."""
    buckets = {p: [0] * NBUCKET for p in PROCS}
    worst = {p: 0 for p in PROCS}
    queries = {p: 0 for p in PROCS}
    cases = {p: 0 for p in PROCS}
    seen: set[tuple[str, str]] = set()
    with open(BOXES / f"{SCENE}-boxes.csv", newline="") as f:
        for row in csv.DictReader(f):
            if int(row["toi_output"]) != ETOI:
                continue
            p, qt = row["processor"], row["query_type"]
            # The untimed warm-up: the first call of each query type.
            if (p, qt) not in seen:
                seen.add((p, qt))
                continue
            for b in range(NBUCKET):
                buckets[p][b] += int(row[f"b{b}"])
            worst[p] = max(worst[p], int(row["worst"]))
            queries[p] += int(row["queries"])
            cases[p] += 1
    return buckets, worst, (queries, cases)


def totals() -> dict:
    """Boxes and bound-discarded boxes per processor, EToI calls only, warm-up kept
    out the same way."""
    out = defaultdict(lambda: [0, 0, 0])
    seen: set[tuple[str, str]] = set()
    with open(BOXES / f"{SCENE}-boxes-totals.csv", newline="") as f:
        for row in csv.DictReader(f):
            if int(row["toi_output"]) != ETOI:
                continue
            p, qt = row["processor"], row["query_type"]
            if (p, qt) not in seen:
                seen.add((p, qt))
                continue
            acc = out[p]
            acc[0] += int(row["queries"])
            acc[1] += int(row["boxes"])
            acc[2] += int(row["bound_discarded"])
    return out


def group(n: int) -> str:
    return f"{n:,}".replace(",", "{,}")


def main() -> int:
    buckets, worst, (queries, cases) = load()
    agg = totals()
    assert queries["host"] == queries["device"], "the two processors ran different sets"
    assert cases["host"] == cases["device"]

    body = []
    for b, label in ROWS:
        counts = [buckets[p][0] if b == 0 else sum(buckets[p][b:]) for p in PROCS]
        body.append(f"    {label:<24}" + "".join(f" & ${group(c)}$" for c in counts) + r" \\")
    body.append(r"    \midrule")
    body.append(r"    worst single query      "
                + "".join(f" & ${group(worst[p])}$" for p in PROCS) + r" \\")
    body.append(r"    \midrule")
    body.append(r"    boxes in total          "
                + "".join(f" & ${group(agg[p][1])}$" for p in PROCS) + r" \\")
    body.append(r"    mean per query          "
                + "".join(f" & ${agg[p][1] / agg[p][0]:.2f}$" for p in PROCS) + r" \\")

    discarded = agg["device"][2]
    share = 100.0 * discarded / agg["device"][1]

    out = PAPER / "generated" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    (out / "tab-difficulty.tex").write_text(rf"""\begin{{table}}[htbp]
  \centering
  \caption{{Boxes examined per query on {SCENE}, over every one of the
    ${group(queries['host'])}$ candidate pairs the broad phase produced in the
    scene's ${group(cases['host'])}$ steps. Both processors are handed the same
    candidates and asked for one earliest time of impact per step, so both
    prune against a shared running minimum. Each narrow phase counts its own
    boxes per query and buckets them by $\lfloor\log_{{2}}\rfloor$ with the same
    code, so the two columns measure the same quantity; the device's count also
    includes boxes it evaluated and then discarded because the bound had already
    passed them, which the host discards before counting, and which is
    ${share:.2f}\%$ of the device total here. Timings are not reported from this
    run: the counters serialise the device kernel.}}
  \label{{tab:difficulty}}
  \begin{{tabular}}{{lrr}}
    \toprule
    boxes examined & host queries & device queries \\
    \midrule
{chr(10).join(body)}
    \bottomrule
  \end{{tabular}}
  \par\smallskip\footnotesize Source: \texttt{{benchmark/results/boxes/}}
\end{{table}}
""")
    print(f"wrote tab-difficulty.tex ({cases['host']} cases, "
          f"{queries['host']:,} queries per processor)")
    for p in PROCS:
        q, bx, bk = agg[p]
        print(f"  {p:<7} boxes {bx:,} ({bx / q:.2f} per query), "
              f"bound-discarded {bk:,}, worst {worst[p]:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
