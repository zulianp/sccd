#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == "__main__":
    import sys
    import os
    if len(sys.argv) != 3:
        print("Usage: python plot_ccd.py <name> <table.csv>")
        sys.exit(1)


    name = sys.argv[1]
    figure_path = os.path.join("figures", name)
    if not os.path.exists("./figures"):
        os.mkdir("./figures")

    
        
    table_file = sys.argv[2]
    # read the table.csv file
    table = pd.read_csv(table_file)
    # plot the table

    expected_hits = table["expected_hit"]
    hits = table["hit"]

    toi = table["toi"][expected_hits].to_numpy()
    expected_toi = table["expected_toi"][expected_hits].to_numpy()
    time = table["time"]

    diff = expected_toi - toi
    normalization = expected_toi + 1e-12
    diff /= normalization

    print(f'Mean:   {np.mean(diff)}')
    print(f'Std:    {np.std(diff)}')
    print(f'Min:    {np.min(diff)}')
    print(f'Max:    {np.max(diff)}')
    print(f'Median: {np.median(diff)}')

    plt.hist(diff, bins=100, log=True)
    plt.ylabel("Count")
    plt.xlabel("Relative Difference Between Expected TOI and TOI")
    plt.savefig(f"{figure_path}_diff_histogram.pdf")
    plt.close()
    # plt.show()

    ind = np.argsort(expected_toi)
    # plt.semilogy(diff)

    
    
    plt.plot(expected_toi[ind], toi[ind]/normalization[ind], '.', label="TOI", alpha=0.2)
    plt.plot([np.min(expected_toi), np.max(expected_toi)], [1, 1], '-', color="red", linewidth=1.0, label="Reference")

    # plt.plot(expected_toi[ind], expected_toi[ind]/normalization[ind], '.', label="Expected TOI", alpha=0.2)

    plt.legend()
    plt.yscale("log")
    plt.ylabel("Relative TOI (1 is reference)")
    plt.xlabel("Time of Impact (Expected TOI)")
    plt.title("Time of Impact (TOI) Difference")
    plt.savefig(f"{figure_path}_diff.png")
    # plt.show()

    hit_diff = expected_hits != hits

    print(f"Number of hit differences: {np.count_nonzero(hit_diff)}")
    print(f"Number of False Positives: {np.count_nonzero(hits[hit_diff])}/{np.count_nonzero(~expected_hits)}")
    print(f"Number of False Negatives: {np.count_nonzero(~hits[hit_diff])}/{np.count_nonzero(expected_hits)}")
    print(f"Queries per second {len(time)}/{time.sum():.5f} = {len(time)/time.sum():.5f} [qxs]")