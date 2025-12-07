#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python plot_ccd.py <name> <table.csv>")
        sys.exit(1)
    name = sys.argv[1]
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


    print(f'Mean:   {np.mean(diff)}')
    print(f'Std:    {np.std(diff)}')
    print(f'Min:    {np.min(diff)}')
    print(f'Max:    {np.max(diff)}')
    print(f'Median: {np.median(diff)}')

    # TODO create a histogram of the difference semilogy
    plt.hist(diff, bins=100, log=True)
    plt.ylabel("Frequency")
    plt.xlabel("Difference Between Expected TOI and TOI")
    plt.title("Time of Impact (TOI) Difference Histogram")
    plt.savefig(f"{name}_diff_histogram.pdf")
    # plt.show()

    ind = np.argsort(expected_toi)
    # plt.semilogy(diff)
    
    plt.plot(toi[ind], label="TOI")
    plt.plot(expected_toi[ind], label="Expected TOI")
    plt.plot(diff[ind], label="Difference")
    plt.legend()
    plt.yscale("log")
    plt.ylabel("Difference Between Expected TOI and TOI")
    plt.xlabel("Query (sorted by Expected TOI)")
    plt.title("Time of Impact (TOI) Difference")
    plt.savefig(f"{name}_diff.pdf")
    # plt.show()

    hit_diff = expected_hits != hits

    print(f"Number of hit differences: {np.count_nonzero(hit_diff)}")
    print(f"Number of False Positives: {np.count_nonzero(hits[hit_diff])}")
    print(f"Number of False Negatives: {np.count_nonzero(~hits[hit_diff])}")
    print(f"Queries per second {len(time)}/{time.sum():.5f} = {len(time)/time.sum():.5f} [qxs]")