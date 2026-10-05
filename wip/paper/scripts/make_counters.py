#!/usr/bin/env python3
"""
Build the hardware-counter table from the committed profiling output.

    python3 scripts/make_counters.py

Reads benchmark/results/profile/ -- Nsight Compute CSV for the Hopper
narrow-phase kernel and `perf stat` output for Grace -- and writes
generated/tables/tab-counters.tex. Like the benchmark tables, nothing is typed
by hand: the numbers in the article are the numbers the profilers wrote.

Scenes appear as rows rather than columns, because there are six of them and a
transposed table would not fit the page. A scene with no file yet is simply
absent, so this runs against a partial sweep.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
PROF = PAPER.parent.parent / "benchmark" / "results" / "profile"

# Benchmark order, and the short labels the other tables use.
SCENES = [
    ("armadillo-rollers", "armadillo-rollers"),
    ("cloth-ball", "cloth-ball"),
    ("cloth-funnel", "cloth-funnel"),
    ("n-body-simulation", "n-body"),
    ("puffer-ball", "puffer-ball"),
    ("rod-twist", "rod-twist"),
]


def ncu(tag: str, scene: str, mode: str = "2") -> dict[str, float]:
    """Mean of each metric over the profiled launches of the search kernel."""
    path = PROF / f"hopper-{tag}-{scene}-m{mode}.csv"
    if not path.is_file():
        return {}
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith('"ID"'))
    except StopIteration:
        return {}
    acc: dict[str, list[float]] = {}
    for row in csv.DictReader(lines[start:]):
        kernel = re.sub(r"<.*", "", (row.get("Kernel Name") or "").replace("void ", ""))
        # The search kernel only. The launch also contains a trivial
        # initialisation kernel whose counters say nothing about the search.
        if "dfs" not in kernel:
            continue
        try:
            acc.setdefault(row["Metric Name"], []).append(
                float((row.get("Metric Value") or "").replace(",", "")))
        except (ValueError, KeyError):
            continue
    return {k: sum(v) / len(v) for k, v in acc.items() if v}


def perf(kind: str, scene: str) -> dict[str, float]:
    """The counters and the derived rates `perf stat` prints beside them."""
    path = PROF / f"grace-{kind}-{scene}.txt"
    if not path.is_file():
        return {}
    out: dict[str, float] = {}
    for line in path.read_text(errors="replace").splitlines():
        m = re.match(r"\s*([\d,]+)\s+([\w\-.]+)", line)
        if m:
            out[m.group(2)] = float(m.group(1).replace(",", ""))
        m = re.search(r"#\s+([\d.]+)\s+insn per cycle", line)
        if m:
            out["ipc"] = float(m.group(1))
        m = re.search(r"#\s+([\d.]+)%\s+of all L1-dcache accesses", line)
        if m:
            out["l1_miss_pct"] = float(m.group(1))
    return out


def fmt(v, spec: str) -> str:
    return "--" if v is None else format(v, spec)


def device_table() -> str:
    rows = []
    for key, label in SCENES:
        m, o, s = ncu("mem", key), ncu("occ", key), ncu("stall", key)
        if not (m or o or s):
            continue
        dram = None
        if m:
            dram = (m.get("dram__bytes_read.sum", 0.0)
                    + m.get("dram__bytes_write.sum", 0.0)) / 1e6
        rows.append(" & ".join([
            label,
            fmt(dram, ".2f"),
            fmt(m.get("dram__throughput.avg.pct_of_peak_sustained_elapsed"), ".3f"),
            fmt(o.get("sm__warps_active.avg.pct_of_peak_sustained_active"), ".1f"),
            fmt(o.get("smsp__thread_inst_executed_per_inst_executed.ratio"), ".1f"),
            fmt(o.get("smsp__inst_executed.avg.per_cycle_active"), ".2f"),
            fmt(o.get("launch__waves_per_multiprocessor"), ".2f"),
            fmt(s.get("smsp__average_warps_issue_stalled_wait_per_issue_active.ratio"), ".2f"),
            fmt(s.get("smsp__average_warps_issue_stalled_barrier_per_issue_active.ratio"), ".2f"),
            fmt(s.get("smsp__average_warps_issue_stalled_long_scoreboard_per_issue_active.ratio"), ".2f"),
        ]) + r" \\")
    body = "\n".join("    " + r for r in rows)
    return rf"""\begin{{table}}[!tb]
  \centering
  \caption{{Nsight Compute counters for the device narrow-phase kernel,
    at the shipped parameters. \emph{{DRAM}} is read and write together for
    one launch and the share of sustained peak bandwidth it represents;
    \emph{{lanes}} is how many of a warp's $32$ are doing useful work when an
    instruction issues; \emph{{issue}} is instructions per scheduler cycle out of a
    possible one; \emph{{waves}} is how many times the grid covers the device. The
    last three columns are warps stalled per issue-active cycle, by reason. A
    profiler serialises launches, so the ratios here are comparable with
    \cref{{sec:results}} and the absolute times are not.}}
  \label{{tab:counters}}
  \fittable{{%
\begin{{tabular}}{{lrrrrrrrrr}}
    \toprule
    & \multicolumn{{2}}{{c}}{{DRAM}} & \multicolumn{{4}}{{c}}{{utilisation}}
      & \multicolumn{{3}}{{c}}{{warps stalled on}} \\
    \cmidrule(lr){{2-3}}\cmidrule(lr){{4-7}}\cmidrule(lr){{8-10}}
    scene & MB & \% peak & occ.\ \% & lanes & issue & waves
          & depend. & barrier & memory \\
    \midrule
{body}
    \bottomrule
  \end{{tabular}}}}
  \par\smallskip\footnotesize Source: \texttt{{benchmark/results/profile/}}
\end{{table}}
"""


def host_table() -> str:
    rows = []
    for key, label in SCENES:
        p, n = perf("pipeline", key), perf("narrow", key)
        if not (p or n):
            continue
        rows.append(" & ".join([
            label,
            fmt(p.get("ipc"), ".2f"),
            fmt(p.get("l1_miss_pct"), ".2f"),
            fmt(n.get("ipc"), ".2f"),
            fmt(n.get("l1_miss_pct"), ".2f"),
        ]) + r" \\")
    if not rows:
        return ""
    body = "\n".join("    " + r for r in rows)
    return rf"""\begin{{table}}[!tb]
  \centering
  \caption{{\texttt{{perf stat}} counters on Grace. \emph{{pipeline}} is the whole
    host pipeline from a build with the device path compiled out, so it is host
    work throughout; \emph{{narrow phase}} is the accuracy oracle, which runs the
    narrow phase over the curated query sets and no broad phase at all. Both are
    whole-process counters.}}
  \label{{tab:counters-host}}
  \begin{{tabular}}{{lrrrr}}
    \toprule
    & \multicolumn{{2}}{{c}}{{pipeline}} & \multicolumn{{2}}{{c}}{{narrow phase}} \\
    \cmidrule(lr){{2-3}}\cmidrule(lr){{4-5}}
    scene & inst/cycle & L1 miss \% & inst/cycle & L1 miss \% \\
    \midrule
{body}
    \bottomrule
  \end{{tabular}}
  \par\smallskip\footnotesize Source: \texttt{{benchmark/results/profile/}}
\end{{table}}
"""


def pipe_table() -> str:
    """Which functional unit is busy, and how busy the SM is overall.

    The pipe percentages are of peak while the SM is *active*; the throughput
    figure is of peak over *elapsed* time. The gap between them is the machine
    standing idle, not the pipes running slowly, which is the point the table
    is here to make.
    """
    rows = []
    for key, label in SCENES:
        m = ncu("pipe", key)
        if not m:
            continue
        add = m.get("smsp__sass_thread_inst_executed_op_dadd_pred_on.sum", 0.0)
        mul = m.get("smsp__sass_thread_inst_executed_op_dmul_pred_on.sum", 0.0)
        fma = m.get("smsp__sass_thread_inst_executed_op_dfma_pred_on.sum", 0.0)
        integer = m.get("smsp__sass_thread_inst_executed_op_integer_pred_on.sum", 0.0)
        fp64 = add + mul + fma
        rows.append(" & ".join([
            label,
            fmt(m.get("sm__inst_executed_pipe_fp64.avg.pct_of_peak_sustained_active"), ".1f"),
            fmt(m.get("sm__inst_executed_pipe_alu.avg.pct_of_peak_sustained_active"), ".1f"),
            fmt(m.get("sm__inst_executed_pipe_lsu.avg.pct_of_peak_sustained_active"), ".1f"),
            fmt(m.get("sm__throughput.avg.pct_of_peak_sustained_elapsed"), ".2f"),
            fmt(integer / fp64 if fp64 else None, ".2f"),
        ]) + r" \\")
    if not rows:
        return ""
    body = "\n".join("    " + r for r in rows)
    return rf"""\begin{{table}}[!tb]
  \centering
  \caption{{Functional-unit use in the device search kernel. The first three
    columns are each pipe's share of its peak issue rate \emph{{while the SM is
    active}}; \emph{{SM}} is the whole multiprocessor's share of peak over
    \emph{{elapsed}} time, so the gap between them is the machine standing idle
    rather than the pipes running slowly. \emph{{int/fp64}} is integer
    thread-instructions per double-precision one. No flop rate is reported: the
    throughput figure for this pipeline is candidate pairs resolved per second
    (\cref{{tab:throughput}}), and an implementation that raised its arithmetic
    rate while returning the same answers would be a worse one.}}
  \label{{tab:pipes}}
  \begin{{tabular}}{{lrrrrr}}
    \toprule
    & \multicolumn{{3}}{{c}}{{pipe, \% of peak while active}}
      & \multicolumn{{2}}{{c}}{{}} \\
    \cmidrule(lr){{2-4}}
    scene & FP64 & ALU & LSU & SM \% & int/fp64 \\
    \midrule
{body}
    \bottomrule
  \end{{tabular}}
  \par\smallskip\footnotesize Source: \texttt{{benchmark/results/profile/}}
\end{{table}}
"""


def main() -> int:
    if not PROF.is_dir():
        print(f"error: {PROF} does not exist", file=sys.stderr)
        return 2
    out = PAPER / "generated" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    (out / "tab-counters.tex").write_text(device_table())
    (out / "tab-pipes.tex").write_text(pipe_table())
    host = host_table()
    (out / "tab-counters-host.tex").write_text(host)
    have = sum(1 for k, _ in SCENES if ncu("mem", k))
    print(f"wrote tab-counters.tex ({have} of {len(SCENES)} scenes) "
          f"and tab-counters-host.tex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
