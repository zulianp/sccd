from __future__ import annotations

import array
import csv
import math
import sys
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


def convert_dataset(data_dir: Path, dataset: str, read_wxf_module: object) -> None:
    base = data_dir / dataset
    roots_dir = base / "roots"
    queries_dir = base / "queries"
    if not roots_dir.is_dir() or not queries_dir.is_dir():
        return

    archives = sorted(roots_dir.glob("*_roots.tar.gz")) + sorted(roots_dir.glob("*_roots.tar"))
    for archive in archives:
        key = roots_key(archive)
        query_path = queries_dir / f"{key}.csv"
        if not query_path.exists():
            continue

        out_dir = roots_dir / key
        out_path = out_dir / "toi.float64"
        if not stale(out_path, (archive, query_path)):
            continue

        n_queries = query_count(query_path)
        toi = array.array("d", [math.nan]) * n_queries
        for query_id, root in read_wxf_module.read_wxf_roots(str(archive)).items():
            if 0 <= query_id < n_queries:
                toi[query_id] = root["t"]

        out_dir.mkdir(parents=True, exist_ok=True)
        with out_path.open("wb") as f:
            toi.tofile(f)


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(f"usage: {argv[0]} <data-dir> <python-dir> [dataset ...]", file=sys.stderr)
        return 1

    data_dir = Path(argv[1])
    python_dir = argv[2]
    datasets = argv[3:]
    sys.path.insert(0, python_dir)

    import read_wxf

    for dataset in datasets:
        convert_dataset(data_dir, dataset, read_wxf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
