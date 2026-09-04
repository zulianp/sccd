#!/usr/bin/env python3
"""Cross-check the C ABI against the datasets' exact roots.

    python3 sccd_crosscheck.py <dataset-dir> <vf|ee>

Reads three files per frame from a benchmark dataset -- the query CSV, the
Mathematica root archive, and the hit/miss booleans -- runs every query through
the ctypes binding, and compares. It writes `<dataset>_<vf|ee>_table.csv`, which
sccd_plot_toi_error.py and sccd_plot_fp_fn.py render.

This is not a unit test and is not registered with ctest: it needs the datasets,
which are gigabytes and not in the repository. sccd_binding_test.py is the test.

What counts as a failure here is the conservativeness invariant, and only that:

  a missed collision, or a time of impact later than the exact root, fails and
  the script exits non-zero;

  a false positive, or a time of impact earlier than the exact root, is reported
  and is not a failure -- reporting a contact that is not there, or one before it
  happens, costs work and never safety.

Requires sympy, for reading the root archives. See README.md.
"""

import sccd_read_queries
import sccd_read_roots
from sccd_reference import vf_F_3d, ee_F_3d, find_root_dfs_3D
import sccd_read_mma
import sccd

import time

class CCDTable:
    def __init__(self, dataset_name: str):
        self.table = []
        self.table.append(["frame", "query_id", "hit", "expected_hit", "toi", "expected_toi", "time"])

    def add_row(self, frame: int, query_id: int, hit: bool, expected_hit: bool, toi: float, expected_toi: float, tts: float):
        self.table.append([frame, query_id, hit, expected_hit, toi, expected_toi, tts])

    def save(self, filename: str):
        with open(filename, "w") as f:
            import csv
            writer = csv.writer(f, delimiter=",")
            writer.writerows(self.table)

class Stats:
    def __init__(self):
        self.num_collision = 0
        self.num_nocollision = 0
        self.time_collision = 0
        self.time_nocollision = 0

    def print(self):
        print(f'YES: {self.num_collision}, {self.time_collision} [s], {self.num_collision/(self.time_collision + 1e-16)} QPS')
        print(f'NO:  {self.num_nocollision}, {self.time_nocollision} [s], {self.num_nocollision/(self.time_collision + 1e-16)} QPS')

if __name__ == "__main__":
    import os
    import glob
    import re as _re
    import sys

    if len(sys.argv) != 3:
        print(f"usage {sys.argv[0]} <data-set> <ee|vf>")
        exit(1) 

    base_folder = sys.argv[1]
    prefix = sys.argv[2]

    if prefix not in ("ee", "vf"):
        print(f"Unsupported prefix {prefix}")
        exit(1)

    # Collect triples for vertex-face (vf) datasets
    query_paths = sorted(glob.glob(os.path.join(base_folder, "queries", f"*{prefix}.csv")))
    roots_paths = sorted(glob.glob(os.path.join(base_folder, "roots", f"*{prefix}_roots.tar.gz")))
    mma_paths   = sorted(glob.glob(os.path.join(base_folder, "mma_bool", f"*{prefix}_mma_bool.json")))

    # Build index by dataset key (e.g., "74vf")
    def vf_key_from_path(p: str) -> str:
        name = os.path.basename(p)
        if prefix == "vf":
            m = _re.match(r"(\d+vf)", name)
        elif prefix == "ee":
            m = _re.match(r"(\d+ee)", name)
        return m.group(1) if m else ""

    queries_by_key = {vf_key_from_path(p): p for p in query_paths}
    roots_by_key   = {vf_key_from_path(p): p for p in roots_paths}
    mma_by_key     = {vf_key_from_path(p): p for p in mma_paths}

    all_keys = sorted(set(queries_by_key) & set(roots_by_key) & set(mma_by_key))
    if not all_keys:
        print(f"No {prefix} dataset triples found.")
        raise SystemExit(0)


    tol_t = 1e-2
    tol_uv = 1e-2

    total_cases = 0
    mismatches = 0
    false_positives = 0
    false_negatives = 0
    # Late is the failure the whole design exists to prevent, and it is invisible
    # to an unsigned error: a query can agree on hit-versus-miss and still report
    # a time of impact after the true one. Count it separately from `mismatches`,
    # which lumps early and late together.
    late = 0
    worst_late = 0.0
    
    funs = {
        "vf" : (sccd.find_root_vf, vf_F_3d, 40, 1e-14),
        "ee" : (sccd.find_root_ee, ee_F_3d, 40, 1e-14),
        # "vf" : (sccd.find_root_tight_inclusion_vf, vf_F_3d, 10000, 1e-8),
        # "ee" : (sccd.find_root_tight_inclusion_ee, ee_F_3d, 10000, 1e-8),
    }

    find_root, F_eval, max_iter, tol = funs[prefix]

    table = CCDTable(os.path.basename(base_folder))

    for key in all_keys:
        query_file = queries_by_key[key]
        root_file  = roots_by_key[key]
        mma_file   = mma_by_key[key]
        
        print(f"Dataset {key}:")
        
        query_data = sccd_read_queries.read_queries(query_file)
        root_map   = sccd_read_roots.read_wxf_roots(root_file)  # Dict[int] -> {t,a,b}
        mma_bool   = sccd_read_mma.read_mma_bool(mma_file)    # List[bool]

        n = len(query_data["s0"][0])
        if len(mma_bool) != n:
            print(f"  Warning: mma_bool length {len(mma_bool)} != queries {n}")

        total_time = 0
        min_toi = 10000
        min_toi_expected = 10000

        stats = Stats()

        print("Started!")
        for i in range(n):
            # print(f"{i}/{n})")
            s0 = [ query_data["s0"][0][i], query_data["s0"][1][i], query_data["s0"][2][i]]
            s1 = [ query_data["s1"][0][i], query_data["s1"][1][i], query_data["s1"][2][i]]
            s2 = [ query_data["s2"][0][i], query_data["s2"][1][i], query_data["s2"][2][i]]
            s3 = [ query_data["s3"][0][i], query_data["s3"][1][i], query_data["s3"][2][i]]
            e0 = [ query_data["e0"][0][i], query_data["e0"][1][i], query_data["e0"][2][i]]
            e1 = [ query_data["e1"][0][i], query_data["e1"][1][i], query_data["e1"][2][i]]
            e2 = [ query_data["e2"][0][i], query_data["e2"][1][i], query_data["e2"][2][i]]
            e3 = [ query_data["e3"][0][i], query_data["e3"][1][i], query_data["e3"][2][i]]

            time_start = time.perf_counter()
            ret = find_root(max_iter, tol, s0, s1, s2, s3, e0, e1, e2, e3)
            time_end = time.perf_counter()
            total_time += time_end - time_start
            expected_hit = bool(mma_bool[i]) 

            if ret[0]:
                stats.time_collision += time_end - time_start
                stats.num_collision += 1
            else:
                stats.time_nocollision += time_end - time_start
                stats.num_nocollision += 1

            if ret[0]:
                min_toi = min(min_toi, ret[1])

            total_cases += 1
            if ret[0] != expected_hit:
                print(f'  {key}:{i}/{n}) hit_mismatch: got {ret[0]} expected {expected_hit}')
                false_positives += 1

                if ret[0]:
                    print("-"*80)
                    print(f'{key}:{i}/{n}) false positive: ret={ret[1:]}')
                    Fx, Fy, Fz = F_eval(s0, s1, s2, s3, e0, e1, e2, e3, ret[1], ret[2], ret[3])
                    print(f'  ({Fx}, {Fy}, {Fz}) residual')
                    print("-"*80)

                if expected_hit:
                    gt = root_map[i]
                    eFx, eFy, eFz = F_eval(s0, s1, s2, s3, e0, e1, e2, e3, gt["t"], gt["a"], gt["b"])
                    print(f'{key}:{i}/{n}) false negative: ret={ret[1:]}, gt=({gt["t"]}, {gt["a"]}, {gt["b"]}), F=({eFx}, {eFy}, {eFz})')
                    false_negatives += 1
                    assert False
                # continue

            expected_toi = 10000
            if expected_hit:
                if i in root_map:
                    gt = root_map[i]
                    t_diff = abs(ret[1] - gt["t"])
                    a_diff = abs(ret[2] - gt["a"])
                    b_diff = abs(ret[3] - gt["b"])
                    expected_toi = gt["t"]
                    min_toi_expected = min(min_toi_expected, expected_toi)
                    if ret[1] > gt["t"]:
                        late += 1
                        worst_late = max(worst_late, ret[1] - gt["t"])
                        print(f'  {key}:{i}/{n}) LATE: toi={ret[1]} after the exact '
                              f'root {gt["t"]} by {ret[1] - gt["t"]:.3e}')
                    if t_diff > tol_t or a_diff > tol_uv or b_diff > tol_uv:
                        print("-"*80)
                        print(f'  {key}:{i}/{n}) root_mismatch: ret={ret[1:]}, gt=({gt["t"]}, {gt["a"]}, {gt["b"]})')
                        print(f'             diffs: t={t_diff} u={a_diff} v={b_diff}')

                        eFx, eFy, eFz = F_eval(s0, s1, s2, s3, e0, e1, e2, e3, gt["t"], gt["a"], gt["b"])
                        Fx, Fy, Fz = F_eval(s0, s1, s2, s3, e0, e1, e2, e3, ret[1], ret[2], ret[3])
                        print(f'  ({Fx}, {Fy}, {Fz}) vs expected ({eFx}, {eFy}, {eFz})')
                        print("-"*80)

                        assert ret[1] <= gt["t"]
                        mismatches += 1
                else:
                    false_positives += 1

            table.add_row(os.path.basename(query_file), i, ret[0], expected_hit, ret[1], expected_toi, time_end - time_start)

        
        print(f"  Done {key}, {mismatches} mismatches {false_positives} false positives, {false_negatives} false negatives. Time: {total_time} [s] Queries/s: {n / (total_time)}");
        print(f"  Min TOI {min_toi} (Conservative? {min_toi <= min_toi_expected} Expected: {min_toi_expected})")
        stats.print()
    print(f"Summary: {total_cases} cases, {mismatches} mismatches, {false_positives} false positives, {false_negatives} false negatives.")

    table.save(f"{os.path.basename(base_folder)}_{prefix}_table.csv")

    if late or false_negatives:
        print(f"FAIL: {false_negatives} missed collisions, {late} times of impact "
              f"after the exact root (worst by {worst_late:.3e}). "
              f"Both break the conservativeness invariant.")
        raise SystemExit(1)
    print(f"OK: no missed collisions and no late times of impact over "
          f"{total_cases} cases. {false_positives} false positives, which are allowed.")
    raise SystemExit(0)

