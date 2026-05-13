#!/usr/bin/env python3
"""Bar charts of total false positives and false negatives per simulation scene."""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def totals_from_csv(path):
    t = pd.read_csv(path)
    h = t["hit"].to_numpy()
    e = t["expected_hit"].to_numpy()
    fp = int(np.count_nonzero(h & ~e))
    fn = int(np.count_nonzero(~h & e))
    return fp, fn


def main():
    argv = sys.argv[1:]
    if len(argv) < 3 or len(argv) % 2 != 1:
        print(
            "Usage: python plot_ccd_fp_fn_bars.py <output_stem> "
            "<scene_label1> <table1.csv> [ <scene_label2> <table2.csv> ... ]",
            file=sys.stderr,
        )
        sys.exit(1)

    out_stem = argv[0]
    pairs = list(zip(argv[1::2], argv[2::2]))
    labels = [lbl for lbl, _ in pairs]
    fp_totals = []
    fn_totals = []
    for _, path in pairs:
        fp, fn = totals_from_csv(path)
        fp_totals.append(fp)
        fn_totals.append(fn)

    x = np.arange(len(labels))
    fig, (ax_fp, ax_fn) = plt.subplots(2, 1, figsize=(max(6.0, 2.0 * len(labels)), 8), sharex=True)

    ax_fp.bar(x, fp_totals, color="tab:orange", edgecolor="black", linewidth=0.5)
    ax_fp.set_ylabel("Total false positives")
    ax_fp.set_title("False positives per scene")

    ax_fn.bar(x, fn_totals, color="tab:blue", edgecolor="black", linewidth=0.5)
    ax_fn.set_ylabel("Total false negatives")
    ax_fn.set_title("False negatives per scene")

    ax_fn.set_xticks(x)
    ax_fn.set_xticklabels(labels)
    plt.setp(ax_fn.xaxis.get_majorticklabels(), rotation=15, ha="right")

    fig.suptitle("Classification errors (full result tables)", y=1.02)
    fig.tight_layout()
    fig.savefig(f"{out_stem}.pdf", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
