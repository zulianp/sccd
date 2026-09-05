"""
Read a sweep CSV and aggregate it without throwing the spread away.

Repeats arrive as separate rows sharing (dataset, mode, case, type), because the
driver measures each case once per process. That is what makes run-to-run
variation visible at all, so it is carried through aggregation rather than
collapsed into a mean: every quantity comes out as a `Stat` with a median and a
spread, and the report refuses to call a difference a result when it is smaller
than the spread it was measured against.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .style import normalise_mode

# Columns the driver emits as timings, in the order a breakdown should show them.
TIMING_COLUMNS = ["prep_ms", "broad_ms", "narrow_ms", "query_narrow_ms", "narrow_ms_s1"]

# Accuracy and conservativeness columns. `toi_late` and `s0_late` are the gate:
# a non-zero value in either is a time of impact reported after the exact root.
ACCURACY_COLUMNS = [
    "fp", "fn", "broad_fp", "broad_fn",
    "toi_n", "toi_late", "toi_max_late", "toi_max_early", "toi_med_early",
    "s0_late", "s0_margin", "s0_toi", "gt_earliest", "root_n", "s1_min",
]


@dataclass
class Stat:
    """One quantity measured several times."""
    values: list[float] = field(default_factory=list)

    def add(self, value: float) -> None:
        if value is not None and math.isfinite(value):
            self.values.append(value)

    @property
    def n(self) -> int:
        return len(self.values)

    @property
    def median(self) -> float:
        if not self.values:
            return math.nan
        ordered = sorted(self.values)
        mid = len(ordered) // 2
        if len(ordered) % 2:
            return ordered[mid]
        return 0.5 * (ordered[mid - 1] + ordered[mid])

    @property
    def mean(self) -> float:
        return sum(self.values) / len(self.values) if self.values else math.nan

    @property
    def std(self) -> float:
        if len(self.values) < 2:
            return 0.0
        m = self.mean
        return math.sqrt(sum((v - m) ** 2 for v in self.values) / (len(self.values) - 1))

    @property
    def lo(self) -> float:
        return min(self.values) if self.values else math.nan

    @property
    def hi(self) -> float:
        return max(self.values) if self.values else math.nan

    @property
    def spread(self) -> float:
        """Full range as a fraction of the median -- the harness's noise floor."""
        m = self.median
        if not self.values or not math.isfinite(m) or m <= 0:
            return 0.0
        return (self.hi - self.lo) / m

    @property
    def total(self) -> float:
        return sum(self.values)


def _as_float(text: str) -> float | None:
    if text is None or text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


# The columns the driver emits. A CSV missing any of them was written by a
# script that had its own idea of the schema, which is exactly how twelve
# columns -- every accuracy column among them -- went unnamed and were filed
# under the None key where nothing ever read them.
EXPECTED_COLUMNS = [
    "dataset", "mode", "case", "type", "queries",
    "prep_ms", "broad_ms", "narrow_ms", "query_narrow_ms",
    "fp", "fn", "broad_fp", "broad_fn", "narrow_ms_s1",
    "toi_n", "toi_late", "toi_max_late", "toi_max_early", "toi_med_early",
    "s0_late", "s0_margin", "s0_toi", "gt_earliest", "root_n", "s1_min",
]

# Columns that must actually carry values, not merely be named. These are the
# accuracy and conservativeness data; a sweep where they are all empty looks
# perfectly healthy and proves nothing.
REQUIRED_NON_EMPTY = ["narrow_ms", "toi_late", "s0_late", "root_n"]


def check_schema(csv_path: Path, rows: list[dict]) -> None:
    """Raise if the CSV is not the schema the driver emits, or is empty where it matters."""
    if not rows:
        raise ValueError(f"{csv_path}: no data rows")

    present = set(rows[0].keys())
    missing = [c for c in EXPECTED_COLUMNS if c not in present]
    if missing:
        raise ValueError(
            f"{csv_path}: missing {len(missing)} column(s): {', '.join(missing)}. "
            f"The header does not match the driver's; regenerate the CSV with a "
            f"header from `sccd_bench --header`.")

    # A row shorter than the header leaves its trailing keys as None, which is
    # not an empty cell -- it means rows from a different schema were
    # concatenated in. That is not hypothetical: a merge step once swept
    # ti_oracle's 20-column output into the 25-column timing CSV, and its
    # timings landed silently in the `toi_late` column, where they read as
    # thousands of late times of impact.
    for i, row in enumerate(rows):
        missing_cells = [c for c in EXPECTED_COLUMNS if row.get(c) is None]
        if missing_cells:
            raise ValueError(
                f"{csv_path}: row {i + 2} has fewer fields than the header "
                f"({len(missing_cells)} column(s) unset, first "
                f"'{missing_cells[0]}'). Rows from another schema look to have "
                f"been concatenated in.")
        if None in row:
            raise ValueError(
                f"{csv_path}: row {i + 2} has more fields than the header names.")

    for column in REQUIRED_NON_EMPTY:
        if not any((row.get(column) or "").strip() for row in rows):
            raise ValueError(
                f"{csv_path}: column '{column}' is empty on every row. The sweep "
                f"produced no {column} data, so nothing computed from it means "
                f"anything.")


def read_rows(csv_path: Path) -> list[dict]:
    """Rows of a sweep CSV, with the mode name mapped onto its current spelling."""
    rows: list[dict] = []
    with Path(csv_path).open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path}: no header")
        if None in reader.fieldnames:
            raise ValueError(
                f"{csv_path}: the header names fewer columns than the rows carry, "
                f"so some columns would be dropped. Regenerate it with "
                f"`sccd_bench --header`.")
        for raw in reader:
            if not raw.get("dataset"):
                continue
            row = dict(raw)
            row["mode"] = normalise_mode(raw.get("mode", ""))
            rows.append(row)
    return rows


@dataclass
class CaseSeries:
    """One case measured under one mode, over however many repeats there were."""
    dataset: str
    mode: str
    case: str
    kind: str
    queries: int = 0
    timings: dict[str, Stat] = field(default_factory=lambda: defaultdict(Stat))
    accuracy: dict[str, Stat] = field(default_factory=lambda: defaultdict(Stat))

    @property
    def repeats(self) -> int:
        for stat in self.timings.values():
            if stat.n:
                return stat.n
        return 0


def by_case(rows: list[dict]) -> dict[tuple[str, str, str], CaseSeries]:
    out: dict[tuple[str, str, str], CaseSeries] = {}
    for row in rows:
        key = (row["dataset"], row["mode"], row["case"])
        series = out.get(key)
        if series is None:
            series = CaseSeries(row["dataset"], row["mode"], row["case"],
                                row.get("type", ""))
            out[key] = series
        queries = _as_float(row.get("queries"))
        if queries is not None:
            series.queries = int(queries)
        for column in TIMING_COLUMNS:
            series.timings[column].add(_as_float(row.get(column)))
        for column in ACCURACY_COLUMNS:
            series.accuracy[column].add(_as_float(row.get(column)))
    return out


@dataclass
class SceneSummary:
    """A whole scene under one mode: totals over cases, spread over repeats."""
    dataset: str
    mode: str
    cases: int = 0
    queries: int = 0
    repeats: int = 0
    # Per-repeat scene totals, so the spread reported is the spread of the number
    # actually quoted -- summing medians would report a total nobody measured.
    totals: dict[str, Stat] = field(default_factory=lambda: defaultdict(Stat))
    fp: int = 0
    fn: int = 0
    broad_fp: int = 0
    broad_fn: int = 0
    toi_late: int = 0
    s0_late: int = 0
    # Two different denominators, and conflating them overstates the evidence.
    # `gt_queries` is how many queries carry ground-truth data at all, NaN
    # (no collision) included. `toi_compared` is how many times of impact were
    # actually put side by side with an exact root, which is the number the
    # conservativeness claim rests on.
    gt_queries: int = 0
    toi_compared: int = 0
    toi_med_early: list[float] = field(default_factory=list)


def by_scene(rows: list[dict]) -> dict[tuple[str, str], SceneSummary]:
    """
    Scene-level summary per mode.

    Timings are totalled within a repeat first and only then compared across
    repeats. Adding up per-case medians would produce a scene total that no
    single run ever produced, and would understate the spread besides.
    """
    # (dataset, mode) -> repeat index -> column -> running total
    per_repeat: dict[tuple[str, str], list[dict[str, float]]] = defaultdict(list)
    seen_case_repeat: dict[tuple[str, str, str], int] = defaultdict(int)

    summaries: dict[tuple[str, str], SceneSummary] = {}
    cases_seen: dict[tuple[str, str], set] = defaultdict(set)
    queries_seen: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)

    for row in rows:
        key = (row["dataset"], row["mode"])
        case = row["case"]
        repeat = seen_case_repeat[(row["dataset"], row["mode"], case)]
        seen_case_repeat[(row["dataset"], row["mode"], case)] += 1

        bucket = per_repeat[key]
        while len(bucket) <= repeat:
            bucket.append(defaultdict(float))
        for column in TIMING_COLUMNS:
            value = _as_float(row.get(column))
            if value is not None and math.isfinite(value):
                bucket[repeat][column] += value

        summary = summaries.get(key)
        if summary is None:
            summary = SceneSummary(row["dataset"], row["mode"])
            summaries[key] = summary
        cases_seen[key].add(case)

        queries = _as_float(row.get("queries"))
        if queries is not None:
            queries_seen[key][case] = int(queries)

        # Counts are properties of the run, not timings: take them once per case
        # (they are identical across repeats) except the violation counters,
        # which are summed over every row so a single bad repeat cannot hide.
        if repeat == 0:
            for name in ("fp", "fn", "broad_fp", "broad_fn"):
                value = _as_float(row.get(name))
                if value is not None:
                    setattr(summary, name, getattr(summary, name) + int(value))
            root_n = _as_float(row.get("root_n"))
            if root_n is not None and root_n > 0:
                summary.gt_queries += int(root_n)
            toi_n = _as_float(row.get("toi_n"))
            if toi_n is not None and toi_n > 0:
                summary.toi_compared += int(toi_n)
            med_early = _as_float(row.get("toi_med_early"))
            if med_early is not None and med_early > 0:
                summary.toi_med_early.append(med_early)
        for name in ("toi_late", "s0_late"):
            value = _as_float(row.get(name))
            if value is not None:
                setattr(summary, name, getattr(summary, name) + int(value))

    for key, summary in summaries.items():
        summary.cases = len(cases_seen[key])
        summary.queries = sum(queries_seen[key].values())
        bucket = per_repeat[key]
        summary.repeats = len(bucket)
        for repeat_totals in bucket:
            for column in TIMING_COLUMNS:
                summary.totals[column].add(repeat_totals.get(column, 0.0))
    return summaries


def separable(best: Stat, worst: Stat) -> tuple[bool, float, float]:
    """
    Is the gap between two measurements larger than the noise that measured them?

    Returns (separable, ratio, noise). A gap inside the observed run-to-run
    spread is reported as inside noise rather than as a ratio -- the discipline
    docs/BENCHMARKS.md states in its first paragraph and that none of the figures
    used to honour.
    """
    b, w = best.median, worst.median
    if not (math.isfinite(b) and math.isfinite(w)) or b <= 0 or w <= 0:
        return (False, math.nan, math.nan)
    ratio = w / b
    gap = (w - b) / w
    noise = max(best.spread, worst.spread)
    return (gap > noise, ratio, noise)
