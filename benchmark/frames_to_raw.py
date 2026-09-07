"""
Convert a scene's PLY frames into the raw-array directories smesh reads.

The driver does this lazily, the first time it needs a frame, by shelling out to
smesh's `db_to_raw`. That has two problems for a prepared dataset: `db_to_raw`
is a Python script needing `meshio` and is not on PATH inside a scheduler job --
a scene whose frames were not already converted then yields a chunk with a
header and no rows -- and it writes only `x.float32`, so a smesh built with
`SMESH_GEOM_TYPE=float64` finds no coordinates it can read.

Doing it here instead settles both. `ply_to_smesh.py` needs nothing outside the
standard library, and it writes whichever precision the smesh you will run
against was built with. Verified against `db_to_raw`: the float32 output is
bit-identical, and the float64 output is the PLY's own doubles exactly for the
scenes that ship `property double`.

**One precision per directory.** Writing `x.float64` beside `x.float32` does not
give a directory that serves both builds -- it gives one that serves neither
reliably. Measured: a float32 build over three armadillo-rollers cases takes one
second against a float32-only tree and does not finish in three hundred seconds
when the float64 files are also present. smesh does not select the file matching
its `geom_t`, so the coordinates it ends up reading are not the ones intended.
The other precision's files are therefore removed from a directory as it is
written, and switching a prepared dataset between builds means reconverting.

Frames are independent, so this runs a process per core, and a frame whose
output is newer than its PLY is skipped -- an interrupted run resumes.
"""

from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ply_to_smesh import read_ply, write_smesh_folder  # noqa: E402

DEFAULT_PRECISION = "float32"
OTHER = {"float32": "float64", "float64": "float32"}


def is_current(out_dir: Path, ply: Path, precision: str) -> bool:
    """
    Present, no older than the PLY, and not contaminated by the other precision.

    The last clause is not fussiness: a directory holding both precisions makes
    the reader take the wrong coordinates, so one left over from an earlier
    conversion has to be treated as stale even though the wanted files are there.
    """
    if any((out_dir / f"{n}.{OTHER[precision]}").exists() for n in "xyz"):
        return False
    needed = [out_dir / f"{n}.{precision}" for n in "xyz"]
    needed += [out_dir / f"i{k}.int32" for k in range(3)]
    if not all(f.exists() for f in needed):
        return False
    ply_mtime = ply.stat().st_mtime
    return all(f.stat().st_mtime >= ply_mtime for f in needed)


def frame_plys(frames: Path) -> list[Path]:
    """
    The PLY files, without the archive's leavings.

    Some of the dataset tarballs were made on macOS and carry an AppleDouble
    resource fork beside each file -- `._balls16_0.ply` next to
    `balls16_0.ply`. They match `*.ply` and are not PLYs, so a converter that
    globs naively tries to parse one and fails.
    """
    return [p for p in sorted(frames.glob("*.ply")) if not p.name.startswith("._")]


def convert_one(job: tuple[str, str, str]) -> tuple[str, str]:
    """Returns (name, error): an empty error means it converted."""
    ply, out_dir, precision = job
    try:
        _, coords, faces = read_ply(ply)
        write_smesh_folder(out_dir, coords, faces, (precision,))
        # Leave no coordinates of the other precision behind: a mixed directory
        # is read wrongly rather than flexibly.
        for n in "xyz":
            stale = Path(out_dir) / f"{n}.{OTHER[precision]}"
            if stale.exists():
                stale.unlink()
        return (Path(ply).name, "")
    except Exception as exc:  # noqa: BLE001 - one bad frame must not stop the rest
        # Report and carry on. A pool that dies on the first unreadable file
        # leaves every later scene unconverted, and the failure surfaces as a
        # scene that mysteriously has no frames rather than as a bad file.
        return (Path(ply).name, f"{type(exc).__name__}: {exc}")


def convert_scene(data_dir: Path, scene: str, workers: int,
                  precision: str = DEFAULT_PRECISION) -> int:
    frames = data_dir / scene / "frames"
    if not frames.is_dir():
        return 0
    out_root = data_dir / scene / "frames_raw"

    plys = frame_plys(frames)
    jobs = []
    for ply in plys:
        out_dir = out_root / ply.stem
        if is_current(out_dir, ply, precision):
            continue
        jobs.append((str(ply), str(out_dir), precision))

    total = len(plys)
    if not jobs:
        print(f"  {scene}: {total} frames already current ({precision})", flush=True)
        return 0

    print(f"  {scene}: converting {len(jobs)} of {total} frames to {precision} "
          f"on {workers} process(es)", flush=True)

    failures: list[tuple[str, str]] = []
    if workers <= 1:
        results = (convert_one(job) for job in jobs)
    else:
        pool = ProcessPoolExecutor(max_workers=workers)
        results = pool.map(convert_one, jobs)
    for name, error in results:
        if error:
            failures.append((name, error))
    if workers > 1:
        pool.shutdown()

    print(f"  {scene}: {len(jobs) - len(failures)} frames converted"
          + (f", {len(failures)} failed" if failures else ""), flush=True)
    for name, error in failures[:5]:
        print(f"    {name}: {error}", flush=True)
    return len(failures)


def main(argv: list[str]) -> int:
    args = [a for a in argv[1:] if not a.startswith("--")]
    workers = 0
    precision = DEFAULT_PRECISION
    for flag in argv[1:]:
        if flag.startswith("--jobs="):
            workers = int(flag.split("=", 1)[1])
        elif flag.startswith("--precision="):
            precision = flag.split("=", 1)[1]
            if precision not in OTHER:
                print(f"error: precision must be float32 or float64",
                      file=sys.stderr)
                return 2
    if len(args) < 2:
        print(f"usage: {argv[0]} [--jobs=N] [--precision=float32|float64] "
              f"<data-dir> <scene> [scene ...]", file=sys.stderr)
        return 2
    if workers <= 0:
        workers = os.cpu_count() or 1

    data_dir = Path(args[0])
    failures = 0
    for scene in args[1:]:
        failures += convert_scene(data_dir, scene, workers, precision)
    if failures:
        print(f"{failures} frame(s) could not be converted", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
