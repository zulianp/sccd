from __future__ import annotations

import array
import csv
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


def query_count(path: Path) -> int:
    with path.open(newline="") as f:
        rows = sum(1 for row in csv.reader(f) if row)
    if rows % 8 != 0:
        raise ValueError(f"{path}: expected a multiple of 8 rows, got {rows}")
    return rows // 8


def roots_key(path: Path) -> str:
    name = path.name
    for suffix in ("_roots.tar.gz", "_roots.tar"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"{path}: unsupported roots archive suffix")


def stale(output: Path, inputs: tuple[Path, ...]) -> bool:
    if not output.exists():
        return True
    output_mtime = output.stat().st_mtime
    return any(input_path.stat().st_mtime > output_mtime for input_path in inputs)


def convert_one(job: tuple[str, str, str, str]) -> tuple[str, int]:
    """
    Convert one archive. Runs in a worker process, so it takes and returns only
    picklable values and imports the reader itself.
    """
    archive, query_path, out_path, python_dir = job
    if python_dir not in sys.path:
        sys.path.insert(0, python_dir)
    import sccd_read_roots

    n_queries = query_count(Path(query_path))
    toi = array.array("d", [math.nan]) * n_queries
    found = 0
    for query_id, root in sccd_read_roots.read_wxf_roots(archive).items():
        if 0 <= query_id < n_queries:
            toi[query_id] = root["t"]
            found += 1

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Write through a temporary and rename: an interrupted run must leave either
    # the previous conversion or none, never a short file that later looks
    # current to `stale()` and is silently treated as complete.
    tmp = out.with_suffix(out.suffix + ".partial")
    with tmp.open("wb") as f:
        toi.tofile(f)
    tmp.replace(out)
    return (Path(archive).name, found)


def pending_jobs(data_dir: Path, dataset: str, python_dir: str) -> list[tuple[str, str, str, str]]:
    base = data_dir / dataset
    roots_dir = base / "roots"
    queries_dir = base / "queries"
    if not roots_dir.is_dir() or not queries_dir.is_dir():
        return []

    archives = sorted(roots_dir.glob("*_roots.tar.gz")) + sorted(roots_dir.glob("*_roots.tar"))
    jobs = []
    for archive in archives:
        key = roots_key(archive)
        query_path = queries_dir / f"{key}.csv"
        if not query_path.exists():
            continue
        out_path = roots_dir / key / "toi.float64"
        if not stale(out_path, (archive, query_path)):
            continue
        jobs.append((str(archive), str(query_path), str(out_path), python_dir))
    return jobs


def convert_dataset(data_dir: Path, dataset: str, python_dir: str,
                    workers: int) -> None:
    """
    Convert every stale archive of a dataset.

    Archives are independent, and the work is symbolic evaluation that pins one
    core, so this runs a process per core. It matters: puffer-ball is 712 MB of
    root archives and takes over three hours in one process. Already-converted
    archives are skipped before the pool is built, so an interrupted run resumes
    and costs only what was in flight.
    """
    jobs = pending_jobs(data_dir, dataset, python_dir)
    if not jobs:
        return

    print(f"  {dataset}: converting {len(jobs)} root archive(s) "
          f"on {workers} process(es)", flush=True)

    if workers <= 1:
        for i, job in enumerate(jobs, 1):
            name, found = convert_one(job)
            print(f"    [{i}/{len(jobs)}] {name}: {found} root(s)", flush=True)
        return

    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, (name, found) in enumerate(pool.map(convert_one, jobs), 1):
            print(f"    [{i}/{len(jobs)}] {name}: {found} root(s)", flush=True)


def main(argv: list[str]) -> int:
    args = [a for a in argv[1:] if not a.startswith("--")]
    workers = 0
    for flag in argv[1:]:
        if flag.startswith("--jobs="):
            workers = int(flag.split("=", 1)[1])

    if len(args) < 2:
        print(f"usage: {argv[0]} [--jobs=N] <data-dir> <python-dir> [dataset ...]",
              file=sys.stderr)
        return 1

    data_dir = Path(args[0])
    python_dir = str(Path(args[1]).resolve())
    datasets = args[2:]
    sys.path.insert(0, python_dir)

    if workers <= 0:
        workers = os.cpu_count() or 1

    # sccd_read_roots pulls in sympy, and is only needed for an archive that is
    # actually stale. Check for it once, and only once there is work to do:
    # importing it eagerly made this script fail on machines where every root
    # file was already converted and there was nothing to convert.
    if any(pending_jobs(data_dir, dataset, python_dir) for dataset in datasets):
        try:
            import sccd_read_roots  # noqa: F401
        except ImportError as exc:
            print(f"error: converting roots needs the 'sccd_read_roots' module and "
                  f"its sympy dependency ({exc}). Install sympy, or leave the "
                  f"already-converted toi.float64 files in place so no conversion "
                  f"is required.", file=sys.stderr)
            return 1

    for dataset in datasets:
        convert_dataset(data_dir, dataset, python_dir, workers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
