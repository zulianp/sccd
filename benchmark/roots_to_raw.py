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

    # sccd_read_roots pulls in sympy, which is only needed to convert an archive that is
    # actually stale. Importing it eagerly made this script fail on machines
    # where every root file was already converted and there was no work to do.
    class LazyRootReader:
        module = None

        def read_wxf_roots(self, path):
            if LazyRootReader.module is None:
                try:
                    import sccd_read_roots
                except ImportError as exc:
                    raise SystemExit(
                        f"error: converting {path} needs the 'sccd_read_roots' module and its sympy "
                        f"dependency ({exc}). Install sympy, or leave the already-converted "
                        f"toi.float64 files in place so no conversion is required."
                    ) from exc
                LazyRootReader.module = sccd_read_roots
            return LazyRootReader.module.read_wxf_roots(path)

    reader = LazyRootReader()
    for dataset in datasets:
        convert_dataset(data_dir, dataset, reader)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
