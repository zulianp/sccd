"""
Pack the query CSVs into a binary form the benchmark can read without parsing.

Each query CSV row is three rationals written as six integers, and the driver
turns each row into three coordinates. Text parsing dominates the cost of
starting a case: cloth-ball's 5.3 M rows take about 1.5 s to parse on a machine
where `long double` is `double`, and the driver re-reads them for every mode of
every repeat. On AArch64 Linux -- the GH200 the GPU sweep runs on -- `long
double` was software-emulated 128-bit quad, which made it far worse. Since the
sweep is chunked into thirty-minute jobs, that time comes straight out of the
budget for measurement.

The packed form is the divided coordinates as float64, in CSV row order, with no
header: the row count follows from the file size, the same convention
`toi.float64` uses. The division is done in double, which is exactly what the
driver now does, so a run over the cache and a run over the text agree bit for
bit rather than approximately.

Conversion is skipped when the output is newer than its input, so this resumes
after an interruption without redoing finished work.
"""

from __future__ import annotations

import array
import sys
from pathlib import Path


def stale(output: Path, source: Path) -> bool:
    if not output.exists():
        return True
    return source.stat().st_mtime > output.stat().st_mtime


def pack_csv(source: Path, output: Path) -> int:
    coords = array.array("d")
    with source.open("rb") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(b",")
            if len(parts) != 6:
                raise ValueError(f"{source}: expected 6 fields, got {len(parts)}")
            values = [int(p) for p in parts]
            coords.append(values[0] / values[1])
            coords.append(values[2] / values[3])
            coords.append(values[4] / values[5])

    rows = len(coords) // 3
    if rows % 8 != 0:
        raise ValueError(f"{source}: expected a multiple of 8 rows, got {rows}")

    output.parent.mkdir(parents=True, exist_ok=True)
    # Write to a neighbouring temporary and rename, so an interrupted run leaves
    # either the previous cache or none -- never a half-written one that would
    # then look current to the driver.
    tmp = output.with_suffix(output.suffix + ".partial")
    with tmp.open("wb") as f:
        coords.tofile(f)
    tmp.replace(output)
    return rows // 8


def convert_scene(data_dir: Path, scene: str, verbose: bool = True) -> tuple[int, int]:
    queries_dir = data_dir / scene / "queries"
    if not queries_dir.is_dir():
        return (0, 0)

    out_dir = data_dir / scene / "queries_raw"
    converted = skipped = 0
    for source in sorted(queries_dir.glob("*.csv")):
        output = out_dir / f"{source.stem}.f64"
        if not stale(output, source):
            skipped += 1
            continue
        pack_csv(source, output)
        converted += 1

    if verbose:
        print(f"  {scene}: packed {converted} query files, {skipped} already current")
    return (converted, skipped)


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(f"usage: {argv[0]} <data-dir> <scene> [scene ...]", file=sys.stderr)
        return 2

    data_dir = Path(argv[1])
    if not data_dir.is_dir():
        print(f"error: {data_dir} is not a directory", file=sys.stderr)
        return 2

    for scene in argv[2:]:
        convert_scene(data_dir, scene)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
