#!/usr/bin/env python3

import csv
import math
import re
import sys
from array import array
from collections import defaultdict
from pathlib import Path


def mean(values):
    return sum(values) / len(values) if values else 0.0


def std(values):
    if len(values) < 2:
        return 0.0
    avg = mean(values)
    return math.sqrt(
        sum((value - avg) * (value - avg) for value in values) / (len(values) - 1)
    )


def stats(values):
    return (
        mean(values),
        std(values),
        min(values) if values else 0.0,
        max(values) if values else 0.0,
    )


def relative_to_or_self(path, base):
    try:
        return path.relative_to(base)
    except ValueError:
        return path


def tex_path(path):
    return r"\detokenize{" + str(path) + "}"


def read_float64_raw(path):
    values = array("d")
    with path.open("rb") as f:
        values.fromfile(f, path.stat().st_size // values.itemsize)
    return values


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
    grouped = defaultdict(dict)
    for row in rows:
        key = (row["dataset"], row["step"])
        if row["type"] in grouped[key]:
            print(
                f"warning: duplicate {row['type']} row for {row['dataset']} step {row['step']}; summing duplicates",
                file=sys.stderr,
            )
            prev = grouped[key][row["type"]]
            for value_key in (
                "queries",
                "broad_ms",
                "narrow_ms",
                "fp",
                "fn",
                "broad_fp",
                "broad_fn",
            ):
                prev[value_key] += row[value_key]
        else:
            grouped[key][row["type"]] = dict(row)

    paired = []
    incomplete = []
    for (dataset, step), by_type in sorted(grouped.items()):
        if "ee" not in by_type or "vf" not in by_type:
            incomplete.append((dataset, step, ",".join(sorted(by_type))))
            continue

        pair = {
            "dataset": dataset,
            "step": step,
            "ee_case": by_type["ee"]["case"],
            "vf_case": by_type["vf"]["case"],
        }
        for key in (
            "queries",
            "broad_ms",
            "narrow_ms",
            "fp",
            "fn",
            "broad_fp",
            "broad_fn",
        ):
            pair[f"ee_{key}"] = by_type["ee"][key]
            pair[f"vf_{key}"] = by_type["vf"][key]
            pair[key] = by_type["ee"][key] + by_type["vf"][key]
        paired.append(pair)

    by_dataset = defaultdict(list)
    for row in paired:
        by_dataset[row["dataset"]].append(row)

    if incomplete:
        examples = ", ".join(
            f"{dataset}/{step} ({types})" for dataset, step, types in incomplete[:5]
        )
        suffix = "" if len(incomplete) <= 5 else f", ... +{len(incomplete) - 5} more"
        print(
            f"warning: skipped {len(incomplete)} incomplete EE/VF pairs: {examples}{suffix}",
            file=sys.stderr,
        )

    return by_dataset, paired


def write_paired_csv(paired_csv, paired_rows):
    fields = ["dataset", "step", "ee_case", "vf_case"]
    for prefix in ("", "ee_", "vf_"):
        fields.extend(
            f"{prefix}{key}"
            for key in (
                "queries",
                "broad_ms",
                "narrow_ms",
                "fp",
                "fn",
                "broad_fp",
                "broad_fn",
            )
        )
    with paired_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in paired_rows:
            writer.writerow(
                {
                    key: (
                        f"{row[key]:.9g}"
                        if isinstance(row.get(key), float)
                        else row.get(key, "")
                    )
                    for key in fields
                }
            )


def write_aggregate_csv(agg_csv, by_dataset):
    fields = [
        "dataset",
        "cases",
        "queries_mean",
        "queries_sum",
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
        "fp_sum",
        "fn_sum",
        "broad_fp_sum",
        "broad_fn_sum",
    ]

    agg_rows = []
    for dataset in sorted(by_dataset):
        values = by_dataset[dataset]
        out = {"dataset": dataset, "cases": len(values)}
        for key in ("queries", "broad_ms", "narrow_ms"):
            values_for_key = [value[key] for value in values]
            (
                out[f"{key}_mean"],
                out[f"{key}_std"],
                out[f"{key}_min"],
                out[f"{key}_max"],
            ) = stats(values_for_key)
            if key == "queries":
                out[f"{key}_sum"] = sum(values_for_key)
        for key in ("fp", "fn", "broad_fp", "broad_fn"):
            out[f"{key}_sum"] = sum(value[key] for value in values)
        agg_rows.append(out)

    with agg_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in agg_rows:
            writer.writerow(
                {
                    key: (
                        f"{row[key]:.9g}"
                        if isinstance(row.get(key), float)
                        else row.get(key, "")
                    )
                    for key in fields
                }
            )


def tex_escape(text):
    return (
        str(text)
        .replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("&", r"\&")
        .replace("%", r"\%")
    )


def write_pgf_plot(figure_dir, dataset, values):
    values = sorted(values, key=lambda row: row["step"])
    if not values:
        return None
    broad = "\n".join(
        f"        ({row['step']},{row['broad_ms']:.12g})" for row in values
    )
    narrow = "\n".join(
        f"        ({row['step']},{row['narrow_ms']:.12g})" for row in values
    )
    total = "\n".join(
        f"        ({row['step']},{(row['broad_ms'] + row['narrow_ms']):.12g})"
        for row in values
    )
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
        f.write(
            r"\addplot+[mark=none,thick] coordinates {"
            "\n" + broad + "\n"
            r"    };"
            "\n"
        )
        f.write(r"\addlegendentry{Broad phase}" "\n")
        f.write(
            r"\addplot+[mark=none,thick] coordinates {"
            "\n" + narrow + "\n"
            r"    };"
            "\n"
        )
        f.write(r"\addlegendentry{Narrow phase}" "\n")
        f.write(
            r"\addplot+[green!60!black,mark=none,thick] coordinates {"
            "\n" + total + "\n"
            r"    };"
            "\n"
        )
        f.write(r"\addlegendentry{Total}" "\n")
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
            total = [row["broad_ms"] + row["narrow_ms"] for row in values]
            fig, ax = plt.subplots(figsize=(7.0, 3.8), constrained_layout=True)
            ax.plot(
                steps,
                broad,
                label=f"Broad phase mean {mean(broad):.3g} ms",
                linewidth=1.5,
            )
            ax.plot(
                steps,
                narrow,
                label=f"Narrow phase mean {mean(narrow):.3g} ms",
                linewidth=1.5,
            )
            ax.plot(
                steps,
                total,
                label=f"Total mean {mean(total):.3g} ms",
                color="green",
                linewidth=1.5,
            )
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
        print(
            f"warning: matplotlib timing plots were not generated; using PGFPlots fallback: {exc}",
            file=sys.stderr,
        )
        for dataset in sorted(by_dataset):
            figure = write_pgf_plot(figure_dir, dataset, by_dataset[dataset])
            if figure is not None:
                figures.append(figure)
    return figures


def collect_toi_points(data_dir, rows):
    points_by_dataset = defaultdict(list)
    missing = 0
    mismatched = 0
    for row in rows:
        roots_dir = data_dir / row["dataset"] / "roots" / row["case"]
        reference_path = roots_dir / "toi.float64"
        sccd_path = roots_dir / "sccd_toi.float64"
        if not reference_path.exists() or not sccd_path.exists():
            missing += 1
            continue

        reference = read_float64_raw(reference_path)
        sccd = read_float64_raw(sccd_path)
        if len(reference) != len(sccd):
            mismatched += 1
            continue

        points = points_by_dataset[row["dataset"]]
        for reference_toi, sccd_toi in zip(reference, sccd):
            if not math.isfinite(reference_toi) or not math.isfinite(sccd_toi):
                continue
            denominator = reference_toi + 1.0e-12
            points.append((reference_toi, sccd_toi / denominator))

    if missing:
        print(
            f"warning: skipped {missing} cases without TOI raw files", file=sys.stderr
        )
    if mismatched:
        print(
            f"warning: skipped {mismatched} cases with mismatched TOI raw sizes",
            file=sys.stderr,
        )
    return points_by_dataset


def write_toi_pgf_plot(figure_dir, dataset, points):
    if not points:
        return None
    csv_path = figure_dir / f"{dataset}_toi_ratio.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("reference_toi", "toi_over_reference"))
        writer.writerows(
            (f"{reference_toi:.12g}", f"{ratio:.12g}")
            for reference_toi, ratio in points
        )

    tex_path_out = figure_dir / f"{dataset}_toi_ratio.tex"
    with tex_path_out.open("w") as f:
        f.write(r"\begin{tikzpicture}" "\n")
        f.write(r"\begin{axis}[" "\n")
        f.write(r"    width=\linewidth,height=6cm," "\n")
        f.write(f"    title={{{tex_escape(dataset)}}},\n")
        f.write(r"    xlabel={Reference time of impact}," "\n")
        f.write(r"    ylabel={Time of impact / reference time of impact}," "\n")
        f.write(r"    grid=both," "\n")
        f.write(r"    legend pos=north west" "\n")
        f.write(r"]" "\n")
        f.write(
            r"\addplot+[only marks,mark=*,mark size=0.35pt,opacity=0.25] "
            r"table[x=reference_toi,y=toi_over_reference,col sep=comma] {"
            + tex_path(relative_to_or_self(csv_path, tex_path_out.parent))
            + r"};"
            "\n"
        )
        f.write(r"\addlegendentry{SCCD}" "\n")
        f.write(r"\addplot+[red,mark=none,thick] coordinates {(0,1) (1,1)};" "\n")
        f.write(r"\addlegendentry{Reference}" "\n")
        f.write(r"\end{axis}" "\n")
        f.write(r"\end{tikzpicture}" "\n")
    return tex_path_out


def write_toi_figures(data_dir, figure_dir, rows):
    points_by_dataset = collect_toi_points(data_dir, rows)
    figures = []
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for dataset in sorted(points_by_dataset):
            points = points_by_dataset[dataset]
            if not points:
                continue
            x = [point[0] for point in points]
            y = [point[1] for point in points]
            fig, ax = plt.subplots(figsize=(7.0, 3.8), constrained_layout=True)
            ax.plot(x, y, ".", alpha=0.25, markersize=1.5, label="SCCD")
            ax.plot(
                [0.0, 1.0],
                [1.0, 1.0],
                "-",
                color="red",
                linewidth=1.0,
                label="Reference",
            )
            ax.set_title(dataset)
            ax.set_xlabel("Reference time of impact")
            ax.set_ylabel("Time of impact / reference time of impact")
            ax.grid(True, alpha=0.25)
            ax.legend()
            pdf = figure_dir / f"{dataset}_toi_ratio.pdf"
            png = figure_dir / f"{dataset}_toi_ratio.png"
            fig.savefig(pdf)
            fig.savefig(png, dpi=180)
            plt.close(fig)
            figures.append(png)
    except Exception as exc:
        print(
            f"warning: matplotlib TOI plots were not generated; using PGFPlots fallback: {exc}",
            file=sys.stderr,
        )
        for dataset in sorted(points_by_dataset):
            figure = write_toi_pgf_plot(figure_dir, dataset, points_by_dataset[dataset])
            if figure is not None:
                figures.append(figure)
    return figures


def write_report(report_tex, agg_csv, paired_csv, figures):
    rel_agg = relative_to_or_self(agg_csv, report_tex.parent)
    rel_paired = relative_to_or_self(paired_csv, report_tex.parent)
    with report_tex.open("w") as f:
        f.write(
            r"""\documentclass{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{pgfplots}
\usepackage{pgfplotstable}
\pgfplotsset{compat=1.18}
\pgfkeys{/pgf/number format/1000 sep={}}
\begin{document}
\section*{SCCD Benchmark Summary}
\noindent Paired data source: \texttt{"""
        )
        f.write(tex_escape(rel_paired))
        f.write(
            r"""}. Timing columns aggregate complete-case EE+VF sums.

{\footnotesize
\pgfplotstabletypeset[
    col sep=comma,
    columns={dataset,cases,broad_ms_mean,broad_ms_min,broad_ms_max,broad_fp_sum,broad_fn_sum},
    columns/dataset/.style={string type,column name=Dataset},
    columns/cases/.style={column name=Cases},
    columns/broad_ms_mean/.style={fixed,precision=3,column name={BP mean}},
    columns/broad_ms_min/.style={fixed,precision=3,column name={BP min}},
    columns/broad_ms_max/.style={fixed,precision=3,column name={BP max}},
    columns/broad_fp_sum/.style={fixed,precision=0,column name={BP FP}},
    columns/broad_fn_sum/.style={fixed,precision=0,column name={BP FN}},
    every head row/.style={before row=\toprule,after row=\midrule},
    every last row/.style={after row=\bottomrule}
]{"""
        )
        f.write(tex_path(rel_agg))
        f.write(
            r"""}
}

{\footnotesize
\pgfplotstabletypeset[
    col sep=comma,
    columns={dataset,cases,queries_sum,narrow_ms_mean,narrow_ms_min,narrow_ms_max,fp_sum,fn_sum},
    columns/dataset/.style={string type,column name=Dataset},
    columns/cases/.style={column name=Cases},
    columns/queries_sum/.style={fixed,precision=0,column name={Queries}},
    columns/narrow_ms_mean/.style={fixed,precision=3,column name={NP mean}},
    columns/narrow_ms_min/.style={fixed,precision=3,column name={NP min}},
    columns/narrow_ms_max/.style={fixed,precision=3,column name={NP max}},
    columns/fp_sum/.style={fixed,precision=0,column name={NP FP}},
    columns/fn_sum/.style={fixed,precision=0,column name={NP FN}},
    every head row/.style={before row=\toprule,after row=\midrule},
    every last row/.style={after row=\bottomrule}
]{"""
        )
        f.write(tex_path(rel_agg))
        f.write(
            r"""}
}

"""
        )
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
    if len(argv) not in (5, 6):
        print(
            "usage: bench_postprocess.py <bench.csv> <aggregate.csv> <figure-dir> <report.tex> [data-dir]",
            file=sys.stderr,
        )
        return 2

    bench_csv = Path(argv[1])
    agg_csv = Path(argv[2])
    paired_csv = agg_csv.with_name(agg_csv.name.replace("_aggregate", "_paired"))
    if paired_csv == agg_csv:
        paired_csv = agg_csv.with_name(f"{agg_csv.stem}_paired{agg_csv.suffix}")
    figure_dir = Path(argv[3])
    report_tex = Path(argv[4])
    data_dir = Path(argv[5]) if len(argv) == 6 else bench_csv.parent.parent / "data"
    agg_csv.parent.mkdir(parents=True, exist_ok=True)
    paired_csv.parent.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    report_tex.parent.mkdir(parents=True, exist_ok=True)

    rows = read_rows(bench_csv)
    by_dataset, paired_rows = pair_rows(rows)
    write_paired_csv(paired_csv, paired_rows)
    write_aggregate_csv(agg_csv, by_dataset)
    figures = write_figures(figure_dir, by_dataset)
    if data_dir.exists():
        figures.extend(write_toi_figures(data_dir, figure_dir, rows))
    else:
        print(
            f"warning: skipped TOI plots because data directory does not exist: {data_dir}",
            file=sys.stderr,
        )
    write_report(report_tex, agg_csv, paired_csv, figures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
