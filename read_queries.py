from __future__ import annotations

from typing import Dict, List

import numpy as np


def read_queries(path: str) -> Dict[str, List[List[float]]]:
    """
    Read a query CSV as documented in data/README.md and return a structure of
    arrays. For vertex-face queries the returned entries are, in order:
    v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1.
    Edge-edge files use: ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1,
    eb0_t1, eb1_t1.
    Each entry stores three lists: x, y, z for every query (SoA).
    """
    data = np.loadtxt(path, delimiter=",", dtype=np.float64, ndmin=2)

    if data.shape[1] != 6:
        raise ValueError("Expected 6 columns per row (numerator/denominator pairs for xyz).")
    if data.shape[0] % 8 != 0:
        raise ValueError("Row count must be a multiple of 8 (8 rows per query).")

    queries = data.reshape(-1, 8, 3, 2)
    coords = queries[..., 0] / queries[..., 1]

    lower_path = path.lower()
    if "vf" in lower_path:
        names = [
            "v_t0",
            "f0_t0",
            "f1_t0",
            "f2_t0",
            "v_t1",
            "f0_t1",
            "f1_t1",
            "f2_t1",
        ]
    else:
        names = [
            "ea0_t0",
            "ea1_t0",
            "eb0_t0",
            "eb1_t0",
            "ea0_t1",
            "ea1_t1",
            "eb0_t1",
            "eb1_t1",
        ]

    soa: Dict[str, List[List[float]]] = {}
    for idx, name in enumerate(names):
        soa[name] = [
            coords[:, idx, 0].tolist(),  # x
            coords[:, idx, 1].tolist(),  # y
            coords[:, idx, 2].tolist(),  # z
        ]
    return soa


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python read_queries.py <path_to_query_file>")
        sys.exit(1)

    queries = read_queries(sys.argv[1])
    for key, (xs, ys, zs) in queries.items():
        print(f"{key}:")
        print(f"  x: {xs}")
        print(f"  y: {ys}")
        print(f"  z: {zs}")
