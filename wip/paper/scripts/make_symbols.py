#!/usr/bin/env python3
"""
Fold the Grace sampling profiles into a per-object breakdown.

    python3 scripts/make_symbols.py

`perf stat` gives whole-process counters, which cannot say where they were
spent. `perf record` can. An earlier version of this script classified by
symbol, which left a large unattributed remainder: a sampling profile of an
optimised build has many local symbols with no name, and guessing which library
they belong to from the address is exactly the kind of inference this project
does not make. Sorting by shared object instead settles it, so that is what the
table reports, with the symbol-level profile used only to name the largest
entries inside our own binary.

Writes generated/tables/tab-symbols.tex.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
PROF = PAPER.parent.parent / "benchmark" / "results" / "profile"

SCENES = [
    ("armadillo-rollers", "armadillo-rollers"),
    ("cloth-ball", "cloth-ball"),
    ("cloth-funnel", "cloth-funnel"),
    ("n-body-simulation", "n-body"),
    ("puffer-ball", "puffer-ball"),
    ("rod-twist", "rod-twist"),
]

# Shared object -> the column it belongs in. SCCD is header-heavy and inlines
# into the driver, so `sccd_bench` and `libsccd` are the same code to a reader.
COLUMNS = [
    ("SCCD", (r"^sccd_bench$", r"^libsccd")),
    ("OpenMP", (r"^libgomp",)),
    ("libm", (r"^libm\.",)),
    ("libc", (r"^libc\.", r"^ld-linux")),
    ("kernel, other", (r"^\[unknown\]$", r"^\[kernel", r".*")),
]


def dso(scene: str) -> dict[str, float] | None:
    path = PROF / f"grace-dso-{scene}.txt"
    if not path.is_file():
        return None
    out = {name: 0.0 for name, _ in COLUMNS}
    seen = False
    for line in path.read_text(errors="replace").splitlines():
        m = re.match(r"\s*([\d.]+)%\s+(\S+)", line)
        if not m or line.lstrip().startswith("#"):
            continue
        pct, obj = float(m.group(1)), m.group(2)
        seen = True
        for name, pats in COLUMNS:
            if any(re.match(p, obj) for p in pats):
                out[name] += pct
                break
    return out if seen else None


def main() -> int:
    rows = []
    for key, label in SCENES:
        d = dso(key)
        if d is None:
            continue
        rows.append("    " + " & ".join(
            [label] + [f"{d[name]:.1f}" for name, _ in COLUMNS]) + r" \\")
    if not rows:
        print("no DSO profiles found", file=sys.stderr)
        return 1
    body = "\n".join(rows)
    head = " & ".join(name for name, _ in COLUMNS)
    tex = rf"""\begin{{table}}[!tb]
  \centering
  \caption{{Where the host spends its cycles, from a sampling profile of the
    CUDA-free build attributed by shared object. \emph{{SCCD}} is our own code,
    which inlines into the driver; \emph{{OpenMP}} is the runtime behind the
    parallel loops; \emph{{libm}} is almost entirely the call that raises $\delta$
    to the third power in \cref{{eq:err}}. The profile covers the whole process,
    so mesh loading and allocation are in it as well as the pipeline, which is
    most of what the last column holds.}}
  \label{{tab:symbols}}
  \begin{{tabular}}{{lrrrrr}}
    \toprule
    scene & {head} \\
    \midrule
{body}
    \bottomrule
  \end{{tabular}}
  \par\smallskip\footnotesize Percentages of sampled cycles.
  Source: \texttt{{benchmark/results/profile/}}
\end{{table}}
"""
    out = PAPER / "generated" / "tables" / "tab-symbols.tex"
    out.write_text(tex)
    print(f"wrote {out.name} ({len(rows)} scenes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
