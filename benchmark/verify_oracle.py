"""
Check that a prepared dataset's ground truth is complete and self-consistent.

The benchmark's accuracy numbers are only as good as the oracle behind them, and
an incomplete oracle does not announce itself: a root that failed to convert
lands as NaN, and NaN is indistinguishable downstream from "this query has no
collision". cloth-funnel sat at 47% converted for exactly that reason, and every
accuracy figure published for it was computed against the other 53% without a
word of warning.

So the two sources are cross-checked against each other. The dataset ships both a
boolean per query (`mma_bool`) and a time of impact (`roots`), and they must
agree query for query:

    mma_bool[i] == isfinite(toi[i])

A `true` with a NaN is a root that never converted; a `false` with a finite time
is a root that should not exist. Either is a defect in the prepared data, and
both exit non-zero -- a benchmark run on top of one is not worth publishing.

Stdlib only, deliberately: this is a gate, and a gate that cannot run because
numpy is missing is a gate that gets skipped.
"""

from __future__ import annotations

import array
import math
import sys
from pathlib import Path


def count_rows(path: Path) -> int:
    """Non-empty line count, read in blocks -- the query CSVs reach a gigabyte."""
    rows = 0
    trailing_newline = True
    with path.open("rb") as f:
        while True:
            block = f.read(1 << 22)
            if not block:
                break
            rows += block.count(b"\n")
            trailing_newline = block.endswith(b"\n")
    if not trailing_newline:
        rows += 1  # a final line without its newline is still a line
    return rows


class SceneReport:
    def __init__(self, name: str) -> None:
        self.name = name
        self.cases = 0
        self.queries = 0
        self.covered = 0
        self.mma_true = 0
        self.missing_toi: list[str] = []
        self.missing_mma: list[str] = []
        self.length_mismatch: list[str] = []
        self.unconverted = 0        # mma_bool true, toi NaN
        self.spurious = 0           # mma_bool false, toi finite
        self.unconverted_cases: list[str] = []
        self.spurious_cases: list[str] = []

    @property
    def gaps(self) -> int:
        return (self.unconverted + self.spurious + len(self.missing_toi)
                + len(self.missing_mma) + len(self.length_mismatch))


def verify_scene(data_dir: Path, scene: str, check_csv: bool) -> SceneReport | None:
    base = data_dir / scene
    queries_dir = base / "queries"
    if not queries_dir.is_dir():
        return None

    report = SceneReport(scene)
    for query_csv in sorted(queries_dir.glob("*.csv")):
        key = query_csv.stem
        toi_path = base / "roots" / key / "toi.float64"
        mma_path = base / "mma_bool" / key / "mma_bool.uint8"

        has_archive = ((base / "roots" / f"{key}_roots.tar.gz").exists()
                       or (base / "roots" / f"{key}_roots.tar").exists())
        if not toi_path.exists():
            # Only a missing conversion is a defect. A scene that ships no roots
            # archive for a case has nothing to convert and nothing to report.
            if has_archive:
                report.missing_toi.append(key)
            continue

        report.cases += 1
        toi = array.array("d")
        with toi_path.open("rb") as f:
            toi.frombytes(f.read())
        report.queries += len(toi)
        report.covered += sum(1 for t in toi if math.isfinite(t))

        if check_csv:
            rows = count_rows(query_csv)
            if rows != 8 * len(toi):
                report.length_mismatch.append(
                    f"{key}: {rows} csv rows is {rows / 8:g} queries, "
                    f"toi has {len(toi)}")

        if not mma_path.exists():
            report.missing_mma.append(key)
            continue

        mma = array.array("B")
        with mma_path.open("rb") as f:
            mma.frombytes(f.read())
        if len(mma) != len(toi):
            report.length_mismatch.append(
                f"{key}: mma_bool has {len(mma)} entries, toi has {len(toi)}")
            continue

        report.mma_true += sum(1 for b in mma if b)
        unconverted = sum(1 for b, t in zip(mma, toi) if b and not math.isfinite(t))
        spurious = sum(1 for b, t in zip(mma, toi) if not b and math.isfinite(t))
        report.unconverted += unconverted
        report.spurious += spurious
        if unconverted:
            report.unconverted_cases.append(f"{key} ({unconverted})")
        if spurious:
            report.spurious_cases.append(f"{key} ({spurious})")

    return report


def main(argv: list[str]) -> int:
    args = argv[1:]
    check_csv = "--check-csv" in args
    args = [a for a in args if a != "--check-csv"]
    if len(args) < 2:
        print(f"usage: {argv[0]} [--check-csv] <data-dir> <scene> [scene ...]",
              file=sys.stderr)
        return 2

    data_dir = Path(args[0])
    scenes = args[1:]

    header = (f"{'scene':<22}{'cases':>7}{'queries':>10}{'roots':>10}"
              f"{'mma_true':>10}{'gaps':>7}")
    print(header)
    print("-" * len(header))

    failed = False
    seen = False
    for scene in scenes:
        report = verify_scene(data_dir, scene, check_csv)
        if report is None:
            print(f"{scene:<22}{'-- not present --':>44}")
            continue
        seen = True
        print(f"{report.name:<22}{report.cases:>7}{report.queries:>10}"
              f"{report.covered:>10}{report.mma_true:>10}{report.gaps:>7}")

        if report.gaps:
            failed = True
            if report.unconverted:
                print(f"    {report.unconverted} queries: mma_bool says a collision "
                      f"exists but the time of impact is NaN (unconverted root)")
                print(f"      {', '.join(report.unconverted_cases[:8])}"
                      + (" ..." if len(report.unconverted_cases) > 8 else ""))
            if report.spurious:
                print(f"    {report.spurious} queries: a finite time of impact "
                      f"where mma_bool says there is no collision")
                print(f"      {', '.join(report.spurious_cases[:8])}"
                      + (" ..." if len(report.spurious_cases) > 8 else ""))
            for label, items in (("no toi.float64", report.missing_toi),
                                 ("no mma_bool.uint8", report.missing_mma),
                                 ("length mismatch", report.length_mismatch)):
                if items:
                    print(f"    {len(items)} cases with {label}: "
                          f"{', '.join(items[:8])}"
                          + (" ..." if len(items) > 8 else ""))

    if not seen:
        print("\nno scene was present; nothing was verified", file=sys.stderr)
        return 2
    if failed:
        print("\nFAILED: the ground truth is incomplete or inconsistent. Accuracy\n"
              "numbers computed against it would silently score unconverted roots\n"
              "as 'no collision'. Re-run the conversion before benchmarking.",
              file=sys.stderr)
        return 1

    print("\nOK: every query has a time of impact exactly where mma_bool says a "
          "collision exists.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
