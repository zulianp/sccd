#!/usr/bin/env python3

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


def mean(values):
    return sum(values) / len(values) if values else 0.0


def std(values):
    if len(values) < 2:
        return 0.0
    avg = mean(values)
    return math.sqrt(sum((value - avg) * (value - avg) for value in values) / (len(values) - 1))


def stats(values):
    return (mean(values), std(values), min(values) if values else 0.0, max(values) if values else 0.0)


def relative_to_or_self(path, base):
    try:
        return path.relative_to(base)
    except ValueError:
        return path


def tex_path(path):
    return r"\detokenize{" + str(path) + "}"


def read_rows(bench_csv):
    rows = []
    with bench_csv.open(newline="") as f:
        for row in csv.DictReader(f):
            if not row.get("dataset"):
                continue
            converted = {
                "dataset": row["dataset"],
                "case": row["case"],
                "type": row["type"],
                "queries": float(row["queries"]),
                "broad_ms": float(row["broad_ms"]),
                "narrow_ms": float(row["narrow_ms"]),
                "fp": float(row["fp"]),
                "fn": float(row["fn"]),
                "broad_fp": float(row["broad_fp"]),
                "broad_fn": float(row["broad_fn"]),
            }
            match = re.match(r"^(\d+)", converted["case"])
            converted["step"] = int(match.group(1)) if match else len(rows)
            rows.append(converted)
    return rows


def pair_rows(rows):
    paired = defaultdict(lambda: {
        "queries": 0.0,
        "broad_ms": 0.0,
        "narrow_ms": 0.0,
        "fp": 0.0,
        "fn": 0.0,
        "broad_fp": 0.0,
        "broad_fn": 0.0,
    })
    for row in rows:
        group = paired[(row["dataset"], row["step"])]
        for key in ("queries", "broad_ms", "narrow_ms", "fp", "fn", "broad_fp", "broad_fn"):
            group[key] += row[key]

    by_dataset = defaultdict(list)
    for (dataset, step), values in paired.items():
        by_dataset[dataset].append({"dataset": dataset, "step": step, **values})
    return by_dataset


def write_aggregate_csv(agg_csv, by_dataset):
    fields = [
        "dataset",
        "cases",
        "queries_mean",
        "queries_std",
        "queries_min",
        "queries_max",
        "broad_ms_mean",
        "broad_ms_std",
        "broad_ms_min",
        "broad_ms_max",
        "narrow_ms_mean",
        "narrow_ms_std",
        "narrow_ms_min",
        "narrow_ms_max",
        "fp_mean",
        "fp_std",
        "fp_min",
        "fp_max",
        "fn_mean",
        "fn_std",
        "fn_min",
        "fn_max",
        "broad_fp_mean",
        "broad_fp_std",
        "broad_fp_min",
        "broad_fp_max",
        "broad_fn_mean",
        "broad_fn_std",
        "broad_fn_min",
        "broad_fn_max",
    ]

    agg_rows = []
    for dataset in sorted(by_dataset):
        values = by_dataset[dataset]
        out = {"dataset": dataset, "cases": len(values)}
        for key in ("queries", "broad_ms", "narrow_ms", "fp", "fn", "broad_fp", "broad_fn"):
            out[f"{key}_mean"], out[f"{key}_std"], out[f"{key}_min"], out[f"{key}_max"] = stats(
                [value[key] for value in values]
            )
        agg_rows.append(out)

    with agg_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in agg_rows:
            writer.writerow({
                key: (f"{row[key]:.9g}" if isinstance(row.get(key), float) else row.get(key, ""))
                for key in fields
            })


def tex_escape(text):
    return str(text).replace("\\", r"\textbackslash{}").replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")


def write_pgf_plot(figure_dir, dataset, values):
    values = sorted(values, key=lambda row: row["step"])
    if not values:
        return None
    broad = " ".join(f"({row['step']},{row['broad_ms']:.12g})" for row in values)
    narrow = " ".join(f"({row['step']},{row['narrow_ms']:.12g})" for row in values)
    path = figure_dir / f"{dataset}_average_timings.tex"
    with path.open("w") as f:
        f.write(r"\begin{tikzpicture}" "\n")
        f.write(r"\begin{axis}[" "\n")
        f.write(r"    width=\linewidth,height=6cm," "\n")
        f.write(f"    title={{{tex_escape(dataset)}}},\n")
        f.write(r"    xlabel={Case}," "\n")
        f.write(r"    ylabel={EE + VF time [ms]}," "\n")
        f.write(r"    grid=both," "\n")
        f.write(r"    legend pos=north west" "\n")
        f.write(r"]" "\n")
        f.write(r"\addplot+[mark=none,thick] coordinates {" + broad + "};\n")
        f.write(r"\addlegendentry{Broad phase}" "\n")
        f.write(r"\addplot+[mark=none,thick] coordinates {" + narrow + "};\n")
        f.write(r"\addlegendentry{Narrow phase}" "\n")
        f.write(r"\end{axis}" "\n")
        f.write(r"\end{tikzpicture}" "\n")
    return path


def write_figures(figure_dir, by_dataset):
    figures = []
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for dataset in sorted(by_dataset):
            values = sorted(by_dataset[dataset], key=lambda row: row["step"])
            if not values:
                continue
            steps = [row["step"] for row in values]
            broad = [row["broad_ms"] for row in values]
            narrow = [row["narrow_ms"] for row in values]
            fig, ax = plt.subplots(figsize=(7.0, 3.8), constrained_layout=True)
            ax.plot(steps, broad, label=f"Broad phase mean {mean(broad):.3g} ms", linewidth=1.5)
            ax.plot(steps, narrow, label=f"Narrow phase mean {mean(narrow):.3g} ms", linewidth=1.5)
            ax.set_title(dataset)
            ax.set_xlabel("Case")
            ax.set_ylabel("EE + VF time [ms]")
            ax.grid(True, alpha=0.25)
            ax.legend()
            pdf = figure_dir / f"{dataset}_average_timings.pdf"
            png = figure_dir / f"{dataset}_average_timings.png"
            fig.savefig(pdf)
            fig.savefig(png, dpi=180)
            plt.close(fig)
            figures.append(pdf)
    except Exception as exc:
        print(f"warning: matplotlib timing plots were not generated; using PGFPlots fallback: {exc}", file=sys.stderr)
        for dataset in sorted(by_dataset):
            figure = write_pgf_plot(figure_dir, dataset, by_dataset[dataset])
            if figure is not None:
                figures.append(figure)
    return figures


def write_report(report_tex, agg_csv, figures):
    rel_agg = relative_to_or_self(agg_csv, report_tex.parent)
    with report_tex.open("w") as f:
        f.write(r"""\documentclass{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{pgfplots}
\usepackage{pgfplotstable}
\pgfplotsset{compat=1.18}
\begin{document}
\section*{SCCD Benchmark Summary}
\pgfplotstabletypeset[
    col sep=comma,
    columns={dataset,cases,broad_ms_mean,narrow_ms_mean,fp_mean,fn_mean,broad_fp_mean,broad_fn_mean},
    columns/dataset/.style={string type,column name=Dataset},
    columns/cases/.style={column name=Cases},
    columns/broad_ms_mean/.style={fixed,precision=3,column name={Broad ms}},
    columns/narrow_ms_mean/.style={fixed,precision=3,column name={Narrow ms}},
    columns/fp_mean/.style={fixed,precision=3,column name={NP FP}},
    columns/fn_mean/.style={fixed,precision=3,column name={NP FN}},
    columns/broad_fp_mean/.style={fixed,precision=3,column name={BP FP}},
    columns/broad_fn_mean/.style={fixed,precision=3,column name={BP FN}},
    every head row/.style={before row=\toprule,after row=\midrule},
    every last row/.style={after row=\bottomrule}
]{""")
        f.write(tex_path(rel_agg))
        f.write(r"""}

""")
        for figure in figures:
            rel = relative_to_or_self(figure, report_tex.parent)
            f.write(r"\begin{figure}[htbp]" "\n")
            f.write(r"\centering" "\n")
            if figure.suffix == ".tex":
                f.write(r"\input{" + tex_path(rel) + "}\n")
            else:
                f.write(r"\includegraphics[width=\linewidth]{" + tex_path(rel) + "}\n")
            f.write(r"\end{figure}" "\n\n")
        f.write(r"\end{document}" "\n")


def main(argv):
    if len(argv) != 5:
        print("usage: bench_postprocess.py <bench.csv> <aggregate.csv> <figure-dir> <report.tex>", file=sys.stderr)
        return 2

    bench_csv = Path(argv[1])
    agg_csv = Path(argv[2])
    figure_dir = Path(argv[3])
    report_tex = Path(argv[4])
    agg_csv.parent.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_tex.parent.mkdir(parents=True, exist_ok=True)

    by_dataset = pair_rows(read_rows(bench_csv))
    write_aggregate_csv(agg_csv, by_dataset)
    figures = write_figures(figure_dir, by_dataset)
    write_report(report_tex, agg_csv, figures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
