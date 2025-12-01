from __future__ import annotations

from typing import List

import numpy as np


def read_queries(path: str) -> List[List[float]]:
    """
    Read a query CSV as documented in data/README.md and return a structure of
    arrays: one list per coordinate per vertex in file order (8 vertices).
    """
    data = np.loadtxt(path, delimiter=",", dtype=np.float64, ndmin=2)

    if data.shape[1] != 6:
        raise ValueError("Expected 6 columns per row (numerator/denominator pairs for xyz).")
    if data.shape[0] % 8 != 0:
        raise ValueError("Row count must be a multiple of 8 (8 rows per query).")

    queries = data.reshape(-1, 8, 3, 2)
    coords = queries[..., 0] / queries[..., 1]

    soa: List[List[float]] = []
    for vertex in range(8):
        soa.append(coords[:, vertex, 0].tolist())
        soa.append(coords[:, vertex, 1].tolist())
        soa.append(coords[:, vertex, 2].tolist())
    return soa


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python read_queries.py <path_to_query_file>")
        sys.exit(1)

    queries = read_queries(sys.argv[1])
    print(queries)
