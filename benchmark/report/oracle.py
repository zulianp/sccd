"""
The accuracy side: ti_oracle's per-mode comparison against the exact roots.

Two references appear here and they are not interchangeable, which is the whole
reason this is a separate table from the timings.

`gt_*` compares against the **dataset's exact roots**, computed symbolically.
That is the conservativeness gate: `gt_missed` is a collision the mode did not
report and `gt_late` is a time of impact after the true one. Both must be zero.

The unprefixed columns compare against **TightInclusion**, whose own answer is
itself a conservative lower bound on the truth rather than the truth. A mode that
reports earlier than TightInclusion is counted "late" or "false negative" there
while being perfectly safe, so those columns over-report and are shown for
context, not as a verdict.

TightInclusion also appears as a row of its own, which is what makes it a timed
reference and not only a correctness oracle.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .data import Stat
from .style import SCENE_LABEL, mode_label, normalise_mode
from .tables import Column, Table


@dataclass
class OracleRow:
    dataset: str
    phase: str
    mode: str
    queries: int = 0
    hits: int = 0
    gt_checked: int = 0
    gt_missed: int = 0
    gt_late: int = 0
    gt_worst_overshoot: float = 0.0
    relerr_median: float = math.nan
    abserr_max: float = math.nan
    ms: Stat = field(default_factory=Stat)


def _f(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def _repeat_count(values: list[int]) -> int:
    """
    How many times the chunk sequence repeats in `values`.

    A scene large enough to be split runs as R identical passes over C chunks,
    and the merged CSV holds all R*C rows under the same (dataset, phase, mode)
    key. Nothing in the row says which chunk or which repeat it is, so the
    period has to be recovered from the sequence: the largest R dividing the
    length for which the values are R copies of the leading C = len/R block.

    Counting a scene's queries has to know this. Summing every row multiplies
    by the repeat count; taking the last row keeps only the final chunk, which
    is what the dataset census used to do -- it reported rod-twist, the one
    scene the oracle splits, as 61,920 queries instead of 549,208.
    """
    n = len(values)
    for r in range(n, 0, -1):
        if n % r:
            continue
        c = n // r
        if all(values[i] == values[i % c] for i in range(n)):
            return r
    return 1


def read(csv_path: Path) -> dict[tuple[str, str, str], OracleRow]:
    """
    Rows keyed by (dataset, phase, mode).

    Per-repeat quantities (`ms`) accumulate into a Stat; per-query totals are
    summed over the chunks of a single pass, so a scene split across chunks is
    counted once and in full.
    """
    out: dict[tuple[str, str, str], OracleRow] = {}
    per_key: dict[tuple[str, str, str], list[dict]] = {}
    with Path(csv_path).open(newline="") as f:
        for raw in csv.DictReader(f):
            if not raw.get("dataset"):
                continue
            mode = normalise_mode(raw.get("mode", ""))
            key = (raw["dataset"], raw.get("phase", ""), mode)
            row = out.get(key)
            if row is None:
                row = OracleRow(raw["dataset"], raw.get("phase", ""), mode)
                out[key] = row
                per_key[key] = []
            per_key[key].append(raw)
            overshoot = _f(raw.get("gt_worst_overshoot"))
            if math.isfinite(overshoot):
                row.gt_worst_overshoot = max(row.gt_worst_overshoot, overshoot)
            row.relerr_median = _f(raw.get("relerr_median"))
            row.abserr_max = _f(raw.get("abserr_max"))
            row.ms.add(_f(raw.get("ms")))

    for key, raws in per_key.items():
        row = out[key]
        counts = [int(_f(r.get("queries")) or 0) for r in raws]
        repeats = _repeat_count(counts)
        per_pass = len(raws) // repeats
        # One pass covers the scene exactly once.
        for r in raws[:per_pass]:
            row.queries += int(_f(r.get("queries")) or 0)
            row.hits += int(_f(r.get("hits")) or 0)
            row.gt_checked += int(_f(r.get("gt_checked")) or 0)
        # Failures are summed over every repeat instead: a violation that shows
        # up in one pass out of three is still a violation, and averaging or
        # sampling it away is the one mistake this table must not make.
        for r in raws:
            row.gt_missed += int(_f(r.get("gt_missed")) or 0)
            row.gt_late += int(_f(r.get("gt_late")) or 0)
    return out


def _ms(stat: Stat) -> str:
    if not stat.n or not math.isfinite(stat.median):
        return "--"
    if stat.n == 1:
        return f"{stat.median:.0f}"
    # median / slowest, as everywhere else a repeat is reported
    return f"{stat.median:.0f} / {stat.hi:.0f}"


def gate_table(rows: dict[tuple[str, str, str], OracleRow], source: str) -> Table:
    table = Table(
        label="tab:gate",
        caption=("The conservativeness gate: every mode's answer checked against "
                 "the dataset's exact symbolic roots. \\emph{missed} counts "
                 "collisions not reported and \\emph{late} counts times of "
                 "impact after the true one. Both must be zero; a late time of "
                 "impact lets a simulation step through the contact."),
        columns=[Column("scene", "l"), Column("phase", "l"), Column("mode", "l"),
                 Column("queries checked"), Column("missed"), Column("late")],
        source=source,
        notes=("TightInclusion is excluded from this table: it is the reference "
               "the queries were selected against, not a subject of it."),
    )
    for (scene, phase, mode), r in sorted(rows.items()):
        if mode == "tight-inclusion" or r.gt_checked == 0:
            continue
        table.add(SCENE_LABEL.get(scene, scene), phase, mode_label(mode),
                  f"{r.gt_checked:,}", f"{r.gt_missed}", f"{r.gt_late}")
    return table


def reference_table(rows: dict[tuple[str, str, str], OracleRow], source: str) -> Table:
    table = Table(
        label="tab:reference",
        caption=("SCCD against TightInclusion over the same queries: hits "
                 "reported, and time taken. TightInclusion is the reference "
                 "implementation of a certified conservative narrow phase, so "
                 "matching its hit count is the strongest agreement available; "
                 "reporting more hits is a false positive, which costs work but "
                 "is never unsafe."),
        columns=[Column("scene", "l"), Column("phase", "l"), Column("mode", "l"),
                 Column("queries"), Column("hits"),
                 Column("time ms", tex_header="time (ms)"),
                 Column("vs. TI", tex_header=r"vs.\ TI")],
        source=source,
    )
    for (scene, phase, mode), r in sorted(rows.items()):
        ti = rows.get((scene, phase, "tight-inclusion"))
        speedup = "--"
        if ti and ti.ms.n and r.ms.n and r.ms.median > 0 and mode != "tight-inclusion":
            speedup = f"{ti.ms.median / r.ms.median:.1f}×"
        elif mode == "tight-inclusion":
            speedup = "1.0× (ref)"
        table.add(SCENE_LABEL.get(scene, scene), phase, mode_label(mode),
                  f"{r.queries:,}", f"{r.hits:,}", _ms(r.ms), speedup)
    return table


def violations(rows: dict[tuple[str, str, str], OracleRow]) -> int:
    return sum(r.gt_missed + r.gt_late
               for (_, _, mode), r in rows.items() if mode != "tight-inclusion")
