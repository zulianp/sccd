# Read query and root files and test the CCD algorithm
# use read_queries.py and read_wxf.py to read the files
# use ccd3D.py to test the CCD algorithm
# use the results to test the CCD algorithm

import read_queries
import read_wxf
from ccd3D import find_root_vf
from numeric_roots import vf_F_3d, find_root_dfs_3D
import read_mma
import sccd_py

import time

if __name__ == "__main__":
    import os
    import glob
    import re as _re
    import sys

    if len(sys.argv) != 2:
        print(f"usage {sys.argv[0]} <data-set>")
        exit(1) 

    base_folder = sys.argv[1]

    # Collect triples for vertex-face (vf) datasets
    query_paths = sorted(glob.glob(os.path.join(base_folder, "queries", "*vf.csv")))
    roots_paths = sorted(glob.glob(os.path.join(base_folder, "roots", "*vf_roots.tar.gz")))
    mma_paths   = sorted(glob.glob(os.path.join(base_folder, "mma_bool", "*vf_mma_bool.json")))

    # Build index by dataset key (e.g., "74vf")
    def vf_key_from_path(p: str) -> str:
        name = os.path.basename(p)
        m = _re.match(r"(\d+vf)", name)
        return m.group(1) if m else ""

    queries_by_key = {vf_key_from_path(p): p for p in query_paths}
    roots_by_key   = {vf_key_from_path(p): p for p in roots_paths}
    mma_by_key     = {vf_key_from_path(p): p for p in mma_paths}

    all_keys = sorted(set(queries_by_key) & set(roots_by_key) & set(mma_by_key))
    if not all_keys:
        print("No vf dataset triples found.")
        raise SystemExit(0)

    tol_t = 1e-2
    tol_uv = 1e-2
    total_cases = 0
    mismatches = 0
    false_positives = 0
    false_negatives = 0

    for key in all_keys:
        query_file = queries_by_key[key]
        root_file  = roots_by_key[key]
        mma_file   = mma_by_key[key]
        
        print(f"Dataset {key}:")
        
        query_data = read_queries.read_queries(query_file)
        root_map   = read_wxf.read_wxf_roots(root_file)  # Dict[int] -> {t,a,b}
        mma_bool   = read_mma.read_mma_bool(mma_file)    # List[bool]

        n = len(query_data["v_t0"][0])
        if len(mma_bool) != n:
            print(f"  Warning: mma_bool length {len(mma_bool)} != queries {n}")

        total_time = 0
        min_toi = 10000
        min_toi_expected = 10000
        for i in range(n):
            sv_3d = [ query_data["v_t0"][0][i], query_data["v_t0"][1][i], query_data["v_t0"][2][i]]
            s1_3d = [ query_data["f0_t0"][0][i], query_data["f0_t0"][1][i], query_data["f0_t0"][2][i]]
            s2_3d = [ query_data["f1_t0"][0][i], query_data["f1_t0"][1][i], query_data["f1_t0"][2][i]]
            s3_3d = [ query_data["f2_t0"][0][i], query_data["f2_t0"][1][i], query_data["f2_t0"][2][i]]
            ev_3d = [ query_data["v_t1"][0][i], query_data["v_t1"][1][i], query_data["v_t1"][2][i]]
            e1_3d = [ query_data["f0_t1"][0][i], query_data["f0_t1"][1][i], query_data["f0_t1"][2][i]]
            e2_3d = [ query_data["f1_t1"][0][i], query_data["f1_t1"][1][i], query_data["f1_t1"][2][i]]
            e3_3d = [ query_data["f2_t1"][0][i], query_data["f2_t1"][1][i], query_data["f2_t1"][2][i]]

            time_start = time.perf_counter()
            ret = sccd_py.find_root_vf_d(1000, 1e-12, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d)
            # ret = sccd_py.find_root_bisection_vf_d(10000, 1e-10, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d)
            time_end = time.perf_counter()
            total_time += time_end - time_start
            expected_hit = bool(mma_bool[i]) 

            if ret[0]:
                min_toi = min(min_toi, ret[1])

            total_cases += 1
            if ret[0] != expected_hit:
                print(f'  {key}:{i}) hit_mismatch: got {ret[0]} expected {expected_hit}')
                false_positives += 1

                if ret[0]:
                    print("-"*80)
                    print(f'{key}:{i}) false positive: ret={ret[1:]}')
                    Fx, Fy, Fz = vf_F_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, ret[1], ret[2], ret[3])
                    print(f'  ({Fx}, {Fy}, {Fz}) residual')
                    print("-"*80)

                if expected_hit:
                    gt = root_map[i]
                    eFx, eFy, eFz = vf_F_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, gt["t"], gt["a"], gt["b"])
                    print(f'{key}:{i}) false negative: ret={ret[1:]}, gt=({gt["t"]}, {gt["a"]}, {gt["b"]}), F=({eFx}, {eFy}, {eFz})')
                    false_negatives += 1
                    assert False
                continue


            if expected_hit:
                # Compare with root map if available (keys are per-query indices from _q<number>_)
                if i in root_map:
                    gt = root_map[i]
                    t_diff = abs(ret[1] - gt["t"])
                    a_diff = abs(ret[2] - gt["a"])
                    b_diff = abs(ret[3] - gt["b"])
                    min_toi_expected = min(min_toi_expected, gt["t"])
                    if t_diff > tol_t or a_diff > tol_uv or b_diff > tol_uv:
                        print("-"*80)
                        print(f'  {key}:{i}) root_mismatch: ret={ret[1:]}, gt=({gt["t"]}, {gt["a"]}, {gt["b"]})')
                        print(f'             diffs: t={t_diff} a={a_diff} b={b_diff}')

                        eFx, eFy, eFz = vf_F_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, gt["t"], gt["a"], gt["b"])
                        Fx, Fy, Fz = vf_F_3d(sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d, ret[1], ret[2], ret[3])
                        print(f'  ({Fx}, {Fy}, {Fz}) vs expected ({eFx}, {eFy}, {eFz})')
                        print("-"*80)
                        mismatches += 1
                else:
                    # Root not present; report but don't fail loudly
                    false_positives += 1

        
        print(f"  Done {key}, {mismatches} mismatches {false_positives} false positives, {false_negatives} false negatives. Time: {total_time} [s] Queries/s: {n / (total_time)}, Min TOI {min_toi} (Expected: {min_toi_expected})")
    print(f"Summary: {total_cases} cases, {mismatches} mismatches, {false_positives} false positives, {false_negatives} false negatives.")

