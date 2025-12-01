# Read query and root files and test the CCD algorithm
# use read_queries.py and read_wxf.py to read the files
# use ccd3D.py to test the CCD algorithm
# use the results to test the CCD algorithm

import read_queries
import read_wxf
from ccd3D import find_root_vf
import read_mma


if __name__ == "__main__":
    # query_file = "data/armadillo-rollers/queries/5vf.csv"
    # root_file  = "data/armadillo-rollers/roots/5vf_roots.tar.gz"
    # mma_bool_file = "data/armadillo-rollers/mma_bool/5vf_mma_bool.json"

    query_file = "data/armadillo-rollers/queries/74vf.csv"
    root_file  = "data/armadillo-rollers/roots/74vf_roots.tar.gz"
    mma_bool_file = "data/armadillo-rollers/mma_bool/74vf_mma_bool.json"
    query_data = read_queries.read_queries(query_file)
    root_data = read_wxf.read_wxf_roots(root_file)
    mma_bool_data = read_mma.read_mma_bool(mma_bool_file)


    for i in range(len(query_data["v_t0"][0])):
        sv_3d = [ query_data["v_t0"][0][i], query_data["v_t0"][1][i], query_data["v_t0"][2][i]]
        s1_3d = [ query_data["f0_t0"][0][i], query_data["f0_t0"][1][i], query_data["f0_t0"][2][i]]
        s2_3d = [ query_data["f1_t0"][0][i], query_data["f1_t0"][1][i], query_data["f1_t0"][2][i]]
        s3_3d = [ query_data["f2_t0"][0][i], query_data["f2_t0"][1][i], query_data["f2_t0"][2][i]]
        ev_3d = [ query_data["v_t1"][0][i], query_data["v_t1"][1][i], query_data["v_t1"][2][i]]
        e1_3d = [ query_data["f0_t1"][0][i], query_data["f0_t1"][1][i], query_data["f0_t1"][2][i]]
        e2_3d = [ query_data["f1_t1"][0][i], query_data["f1_t1"][1][i], query_data["f1_t1"][2][i]]
        e3_3d = [ query_data["f2_t1"][0][i], query_data["f2_t1"][1][i], query_data["f2_t1"][2][i]]
        found = find_root_vf(1000, 1e-6, sv_3d, s1_3d, s2_3d, s3_3d, ev_3d, e1_3d, e2_3d, e3_3d)

        if found[0] == mma_bool_data[i]:

            if found[0]:
                t_diff = abs(found[1] - root_data[i]["t"])
                a_diff = abs(found[2] - root_data[i]["a"])
                b_diff = abs(found[3] - root_data[i]["b"])
                if t_diff > 1e-6 or a_diff > 1e-6 or b_diff > 1e-6:
                    print(f'{i}) {found[1:]} != ({root_data[i]["t"]}, {root_data[i]["a"]}, {root_data[i]["b"]})')
                    print(f't_diff: {t_diff}, a_diff: {a_diff}, b_diff: {b_diff}')
        else:
            print(f'{i}) {found} != {mma_bool_data[i]}')

