#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plot_ccd.py <table.csv>")
        sys.exit(1)
    table_file = sys.argv[1]
    # read the table.csv file
    table = pd.read_csv(table_file)
    # plot the table

    expected_hits = table["expected_hit"]
    hits = table["hit"]
    toi = table["toi"]
    expected_toi = table["expected_toi"]
    time = table["time"]

    diff = np.abs(toi[expected_hits] - expected_toi[expected_hits])
    hit_diff = expected_hits != hits

    print(f"Number of hit differences: {np.count_nonzero(hit_diff)}")
    print(f"Number of False Positives: {np.count_nonzero(hits[hit_diff])}")
    print(f"Number of False Negatives: {np.count_nonzero(~hits[hit_diff])}")
    print(f"Queries per second {len(time)}/{time.sum():.5f} = {len(time)/time.sum():.5f} [qxs]")

    print(f'Mean:   {np.mean(diff)}')
    print(f'Std:    {np.std(diff)}')
    print(f'Min:    {np.min(diff)}')
    print(f'Max:    {np.max(diff)}')
    print(f'Median: {np.median(diff)}')

    # TODO create a histogram of the difference semilogy
    plt.hist(diff, bins=100, log=True)
    plt.ylabel("Frequency")
    plt.xlabel("Difference Between TOI and Expected TOI")
    plt.title("Time of Impact (TOI) Difference Histogram")
    plt.savefig("diff_histogram.pdf")
    plt.show()

    # plt.plot(diff)
    plt.semilogy(diff)
    plt.ylabel("Difference Between TOI and Expected TOI")
    plt.xlabel("Query")
    plt.title("Time of Impact (TOI) Difference")
    plt.savefig("diff.pdf")
    plt.show()