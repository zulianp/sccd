"""
Convert a scene's PLY frames into the raw-array directories smesh reads.

The driver does this lazily, the first time it needs a frame, by shelling out to
smesh's `db_to_raw`. That has two problems for a prepared dataset: `db_to_raw`
is a Python script needing `meshio` and is not on PATH inside a scheduler job --
a scene whose frames were not already converted then yields a chunk with a
header and no rows -- and it writes only `x.float32`, so a smesh built with
`SMESH_GEOM_TYPE=float64` finds no coordinates it can read.

Doing it here instead settles both. `ply_to_smesh.py` needs nothing outside the
standard library, and it writes both precisions: smesh reads whichever matches
the `geom_t` it was built with, so one prepared dataset serves a float32 and a
float64 build without reconverting. Verified against `db_to_raw`: the float32
output is bit-identical, and the float64 output is the PLY's own doubles exactly
for the scenes that ship `property double`.

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

PRECISIONS = ("float64", "float32")


def is_current(out_dir: Path, ply: Path) -> bool:
    """Every file the reader needs is present and no older than the PLY."""
    needed = [out_dir / f"{n}.{p}" for n in "xyz" for p in PRECISIONS]
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


def convert_one(job: tuple[str, str]) -> tuple[str, str]:
    """Returns (name, error): an empty error means it converted."""
    ply, out_dir = job
    try:
        _, coords, faces = read_ply(ply)
        write_smesh_folder(out_dir, coords, faces, PRECISIONS)
        return (Path(ply).name, "")
    except Exception as exc:  # noqa: BLE001 - one bad frame must not stop the rest
        # Report and carry on. A pool that dies on the first unreadable file
        # leaves every later scene unconverted, and the failure surfaces as a
        # scene that mysteriously has no frames rather than as a bad file.
        return (Path(ply).name, f"{type(exc).__name__}: {exc}")


def convert_scene(data_dir: Path, scene: str, workers: int) -> int:
    frames = data_dir / scene / "frames"
    if not frames.is_dir():
        return 0
    out_root = data_dir / scene / "frames_raw"

    plys = frame_plys(frames)
    jobs = []
    for ply in plys:
        out_dir = out_root / ply.stem
        if is_current(out_dir, ply):
            continue
        jobs.append((str(ply), str(out_dir)))

    total = len(plys)
    if not jobs:
        print(f"  {scene}: {total} frames already current", flush=True)
        return 0

    print(f"  {scene}: converting {len(jobs)} of {total} frames "
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
    for flag in argv[1:]:
        if flag.startswith("--jobs="):
            workers = int(flag.split("=", 1)[1])
    if len(args) < 2:
        print(f"usage: {argv[0]} [--jobs=N] <data-dir> <scene> [scene ...]",
              file=sys.stderr)
        return 2
    if workers <= 0:
        workers = os.cpu_count() or 1

    data_dir = Path(args[0])
    failures = 0
    for scene in args[1:]:
        failures += convert_scene(data_dir, scene, workers)
    if failures:
        print(f"{failures} frame(s) could not be converted", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
