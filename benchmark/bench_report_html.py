#!/usr/bin/env python3
"""Render the SCCD benchmark results as a self-contained HTML report.

Reads whichever inputs are present and renders what it can:

  --aggregate   bench_aggregate.csv   pipeline timing, one row per dataset+mode
  --toi-error   bench_toi_error.csv   time-of-impact error distribution
  --oracle      oracle.csv            per-mode comparison against TightInclusion

The output has no external assets beyond a web font, so it can be mailed,
committed, or opened straight off a cluster filesystem. Charts are emitted as
inline SVG rather than through a plotting library: matplotlib is frequently
absent on the machines this runs on, and the fallback used to silently degrade.
"""

import argparse
import csv
import datetime as _dt
import html
import math
import subprocess
from pathlib import Path

# --- design tokens ---------------------------------------------------------
# Series colours are validated for colour-vision deficiency separation and for
# contrast against each mode's surface; dark is a separate set of steps, not an
# inversion of light.
SERIES_LIGHT = ["#0E8AA0", "#C07C10", "#8A3C74"]
SERIES_DARK = ["#0096AE", "#B8862F", "#9C5A88"]
REFERENCE_INK = "#7C8894"

MODE_ORDER = ["scalar", "fast-vector", "vector", "conservative", "ti-vec", "ti-compat", "ti-reference"]
# The oracle and the benchmark driver name the same kernels differently.
MODE_ALIASES = {"vector": "fast-vector", "ti-vec": "conservative"}
MODE_BLURB = {
    "scalar": "Scalar reference search.",
    "fast-vector": "Lane-packed vectorized kernel, vertex-face only.",
    "conservative": "Lane-packed kernel reproducing TightInclusion exactly.",
    "ti-compat": "Fast kernel corrected by TightInclusion. An oracle, not a code path.",
    "ti-reference": "TightInclusion itself, the reference implementation.",
}
CONSERVATIVE_MODES = {"conservative", "ti-compat", "ti-reference"}


def canonical_mode(name):
    return MODE_ALIASES.get(name, name)


def mode_sort_key(name):
    return MODE_ORDER.index(name) if name in MODE_ORDER else len(MODE_ORDER)


def read_csv(path):
    if not path or not Path(path).is_file():
        return []
    with open(path, newline="") as f:
        return [row for row in csv.DictReader(f) if any(v for v in row.values())]


def as_float(row, key, default=0.0):
    try:
        return float(row.get(key, default) or default)
    except (TypeError, ValueError):
        return default


def as_int(row, key, default=0):
    return int(as_float(row, key, default))


def fmt_int(value):
    return f"{int(value):,}"


def fmt_ms(value):
    if value >= 10000:
        return f"{value / 1000:.1f} s"
    if value >= 100:
        return f"{value:.0f} ms"
    return f"{value:.1f} ms"


def fmt_sci(value):
    if value == 0:
        return "0"
    if value < 1e-3:
        exponent = int(math.floor(math.log10(abs(value))))
        return f"{value / (10 ** exponent):.1f}e{exponent}"
    return f"{value:.3g}"


# --- chart primitives ------------------------------------------------------
# Marks are thin, ends are rounded and anchored to the baseline, adjacent fills
# keep a 2px surface gap, and grid lines stay recessive.

def svg_grouped_bars(groups, series, values, *, unit, log=False, height=250, label_fmt=fmt_ms):
    """Grouped bar chart.

    groups  outer category labels (one cluster each)
    series  series names, in fixed order -- colour follows the series, never rank
    values  {(group, series): value}
    """
    if not groups or not series:
        return '<p class="empty">No data for this chart.</p>'

    pad_l, pad_r, pad_t, pad_b = 58, 14, 14, 46
    width = max(360, 96 * len(groups) + pad_l + pad_r)
    plot_w = width - pad_l - pad_r
    plot_h = height - pad_t - pad_b

    finite = [v for v in values.values() if v is not None and math.isfinite(v) and v > 0]
    vmax = max(finite) if finite else 1.0
    vmin = min(finite) if finite else 1.0

    if log:
        lo = math.log10(max(vmin, 1e-12)) - 0.15
        hi = math.log10(vmax) + 0.15

        def scale(v):
            if v is None or v <= 0:
                return 0.0
            return (math.log10(v) - lo) / (hi - lo) * plot_h
        ticks = [10 ** e for e in range(int(math.floor(lo)), int(math.ceil(hi)) + 1)]
    else:
        hi = vmax * 1.12 or 1.0

        def scale(v):
            return 0.0 if not v else (v / hi) * plot_h
        ticks = [hi * f for f in (0, 0.25, 0.5, 0.75, 1.0)]

    out = [f'<svg viewBox="0 0 {width} {height}" role="img" class="chart" '
           f'preserveAspectRatio="xMinYMid meet">']

    for t in ticks:
        y = pad_t + plot_h - scale(t)
        if not (pad_t - 1 <= y <= pad_t + plot_h + 1):
            continue
        out.append(f'<line class="grid" x1="{pad_l}" x2="{width - pad_r}" y1="{y:.1f}" y2="{y:.1f}"/>')
        out.append(f'<text class="tick" x="{pad_l - 8}" y="{y + 3.5:.1f}" text-anchor="end">'
                   f'{html.escape(label_fmt(t) if not log else fmt_sci(t))}</text>')

    group_w = plot_w / len(groups)
    bar_gap = 2.0  # surface gap between adjacent fills
    bar_w = max(6.0, (group_w * 0.72 - bar_gap * (len(series) - 1)) / len(series))

    for gi, group in enumerate(groups):
        cx = pad_l + group_w * (gi + 0.5)
        span = bar_w * len(series) + bar_gap * (len(series) - 1)
        x0 = cx - span / 2
        for si, name in enumerate(series):
            v = values.get((group, name))
            x = x0 + si * (bar_w + bar_gap)
            h = scale(v)
            y = pad_t + plot_h - h
            label = "no data" if v is None else f"{label_fmt(v)}"
            out.append(
                f'<rect class="bar s{si}" x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" '
                f'height="{max(h, 0.8):.1f}" rx="2" '
                f'data-tip="{html.escape(group)} &middot; {html.escape(name)}: {html.escape(label)}">'
                f'<title>{html.escape(group)} - {html.escape(name)}: {html.escape(label)}</title></rect>')
        out.append(f'<text class="axis" x="{cx:.1f}" y="{height - pad_b + 18:.1f}" '
                   f'text-anchor="middle">{html.escape(group)}</text>')

    out.append(f'<line class="axis-line" x1="{pad_l}" x2="{width - pad_r}" '
               f'y1="{pad_t + plot_h}" y2="{pad_t + plot_h}"/>')
    out.append(f'<text class="unit" x="{pad_l - 8}" y="{pad_t - 2}" text-anchor="end">'
               f'{html.escape(unit)}</text>')
    out.append("</svg>")
    return "".join(out)


def legend(series):
    items = "".join(
        f'<span class="key"><i class="sw s{i}"></i>{html.escape(name)}</span>'
        for i, name in enumerate(series))
    return f'<div class="legend">{items}</div>'


def data_table(headers, rows, *, caption):
    head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>" for row in rows)
    return (f'<details class="tableview"><summary>{html.escape(caption)}</summary>'
            f'<div class="scroll"><table><thead><tr>{head}</tr></thead>'
            f"<tbody>{body}</tbody></table></div></details>")


# --- page ------------------------------------------------------------------

STYLE = """
:root{
  color-scheme: light dark;
  --surface:#FAFBFC; --raised:#FFFFFF; --sunken:#F1F4F6;
  --ink:#141A20; --ink-2:#4E5A66; --ink-3:#7C8894;
  --rule:#E1E7EC; --rule-strong:#CBD4DC;
  --accent:#0E8AA0; --accent-soft:#E2F2F5;
  --good:#1F7A4C; --good-soft:#E4F2EA;
  --bad:#B03A32; --bad-soft:#FBE9E7;
  --warn:#9A6B12; --warn-soft:#FBF0DC;
  --s0:#0E8AA0; --s1:#C07C10; --s2:#8A3C74; --s3:#7C8894;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --surface:#101418; --raised:#171D23; --sunken:#0B0F13;
    --ink:#E6EDF3; --ink-2:#A3B3C0; --ink-3:#71818F;
    --rule:#232C34; --rule-strong:#33404A;
    --accent:#0096AE; --accent-soft:#0C2A31;
    --good:#4FB07A; --good-soft:#12261C;
    --bad:#E07069; --bad-soft:#2B1512;
    --warn:#D3A03F; --warn-soft:#2A2113;
    --s0:#0096AE; --s1:#B8862F; --s2:#9C5A88; --s3:#71818F;
  }
}
:root[data-theme="dark"]{
  --surface:#101418; --raised:#171D23; --sunken:#0B0F13;
  --ink:#E6EDF3; --ink-2:#A3B3C0; --ink-3:#71818F;
  --rule:#232C34; --rule-strong:#33404A;
  --accent:#0096AE; --accent-soft:#0C2A31;
  --good:#4FB07A; --good-soft:#12261C;
  --bad:#E07069; --bad-soft:#2B1512;
  --warn:#D3A03F; --warn-soft:#2A2113;
  --s0:#0096AE; --s1:#B8862F; --s2:#9C5A88; --s3:#71818F;
}

*{box-sizing:border-box}
body{
  margin:0; background:var(--surface); color:var(--ink);
  font-family:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
  font-size:15px; line-height:1.62;
  -webkit-font-smoothing:antialiased;
}
.wrap{max-width:60rem;margin:0 auto;padding:3.5rem 1.5rem 5rem}
h1,h2,h3{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;text-wrap:balance;margin:0}
h1{font-size:2.3rem;line-height:1.18;letter-spacing:-.015em}
h2{font-size:1.3rem;margin-bottom:.4rem}
h3{font-size:1rem;margin-bottom:.2rem}
p{margin:.5rem 0;max-width:65ch;color:var(--ink-2)}
.eyebrow{
  font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:.7rem;letter-spacing:.14em;text-transform:uppercase;
  color:var(--ink-3);margin:0 0 .6rem
}
header.masthead{display:flex;flex-direction:column;gap:.4rem;padding-bottom:1.75rem;
  border-bottom:1px solid var(--rule)}
.meta{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.75rem;color:var(--ink-3);
  display:flex;flex-wrap:wrap;gap:.35rem 1.25rem;margin-top:.6rem}
section{margin-top:3rem}

/* verdict: the summary before the detail */
.verdict{display:flex;gap:1rem;align-items:flex-start;margin-top:2rem;padding:1.15rem 1.3rem;
  border-radius:10px;border:1px solid var(--rule);background:var(--raised);
  border-left:3px solid var(--accent)}
.verdict.pass{border-left-color:var(--good);background:var(--good-soft)}
.verdict.fail{border-left-color:var(--bad);background:var(--bad-soft)}
.verdict h2{font-size:1.05rem}
.verdict p{margin:.25rem 0 0;color:var(--ink-2)}

.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(9.5rem,1fr));gap:.75rem;margin-top:1.5rem}
.tile{background:var(--raised);border:1px solid var(--rule);border-radius:9px;padding:.85rem .95rem}
.tile .n{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:1.5rem;font-weight:500;
  letter-spacing:-.02em;font-variant-numeric:tabular-nums;display:block;line-height:1.25}
.tile .k{font-size:.72rem;color:var(--ink-3);letter-spacing:.04em;text-transform:uppercase;
  font-family:"IBM Plex Mono",ui-monospace,monospace}
.tile.is-bad .n{color:var(--bad)} .tile.is-good .n{color:var(--good)}

.card{background:var(--raised);border:1px solid var(--rule);border-radius:10px;padding:1.25rem;margin-top:1rem}
.scroll{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:.85rem}
th,td{text-align:right;padding:.42rem .6rem;border-bottom:1px solid var(--rule);white-space:nowrap}
th:first-child,td:first-child{text-align:left}
th{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.68rem;letter-spacing:.07em;
  text-transform:uppercase;color:var(--ink-3);font-weight:500;border-bottom:1px solid var(--rule-strong)}
td{font-family:"IBM Plex Mono",ui-monospace,monospace;font-variant-numeric:tabular-nums;color:var(--ink-2)}
td.name{font-family:"IBM Plex Sans",sans-serif;color:var(--ink)}
tbody tr:hover{background:var(--sunken)}

.pill{display:inline-flex;align-items:center;gap:.35rem;padding:.1rem .5rem;border-radius:999px;
  font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.7rem;letter-spacing:.03em;
  border:1px solid transparent}
.pill.ok{background:var(--good-soft);color:var(--good);border-color:var(--good)}
.pill.no{background:var(--bad-soft);color:var(--bad);border-color:var(--bad)}
.pill.ref{background:var(--sunken);color:var(--ink-3);border-color:var(--rule-strong)}
.pill::before{content:"";width:.45rem;height:.45rem;border-radius:50%;background:currentColor}

.chart{width:100%;height:auto;display:block;overflow:visible}
.chart .grid{stroke:var(--rule);stroke-width:1}
.chart .axis-line{stroke:var(--rule-strong);stroke-width:1}
.chart .tick,.chart .axis,.chart .unit{font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:10px;fill:var(--ink-3)}
.chart .axis{fill:var(--ink-2);font-size:10.5px}
.chart .bar{transition:opacity .12s ease}
.chart:hover .bar{opacity:.45}
.chart .bar:hover{opacity:1}
.chart .bar.s0{fill:var(--s0)} .chart .bar.s1{fill:var(--s1)}
.chart .bar.s2{fill:var(--s2)} .chart .bar.s3{fill:var(--s3)}
.legend{display:flex;flex-wrap:wrap;gap:.35rem 1rem;margin:.7rem 0 0;
  font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.74rem;color:var(--ink-2)}
.key{display:inline-flex;align-items:center;gap:.4rem}
.sw{width:.7rem;height:.7rem;border-radius:2px;display:inline-block}
.sw.s0{background:var(--s0)} .sw.s1{background:var(--s1)}
.sw.s2{background:var(--s2)} .sw.s3{background:var(--s3)}

.tableview{margin-top:.9rem;border-top:1px solid var(--rule);padding-top:.5rem}
.tableview summary{cursor:pointer;font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:.74rem;color:var(--ink-3);letter-spacing:.04em}
.tableview summary:hover{color:var(--ink-2)}
.tableview[open] summary{margin-bottom:.6rem}
.empty{color:var(--ink-3);font-style:italic}
.modes{display:grid;gap:.6rem;margin-top:1rem}
.mode{display:grid;grid-template-columns:auto 1fr auto;gap:.85rem;align-items:baseline;
  padding:.7rem .9rem;border:1px solid var(--rule);border-radius:8px;background:var(--raised)}
.mode code{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.78rem;color:var(--accent)}
.mode p{margin:0;font-size:.86rem}
footer{margin-top:3.5rem;padding-top:1.25rem;border-top:1px solid var(--rule);
  font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.73rem;color:var(--ink-3)}
footer p{color:inherit;font-family:inherit;font-size:inherit}
#tip{position:fixed;z-index:20;pointer-events:none;opacity:0;transition:opacity .1s ease;
  background:var(--ink);color:var(--surface);padding:.3rem .55rem;border-radius:5px;
  font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.72rem;white-space:nowrap;
  box-shadow:0 4px 14px rgba(0,0,0,.18)}
a{color:var(--accent)}
:focus-visible{outline:2px solid var(--accent);outline-offset:2px;border-radius:3px}
@media (prefers-reduced-motion:reduce){*{transition:none!important}}
@media (max-width:640px){.wrap{padding:2.25rem 1rem 3rem}h1{font-size:1.75rem}}
"""

SCRIPT = """
(function(){
  var tip=document.getElementById('tip');
  document.addEventListener('mouseover',function(e){
    var t=e.target.closest('[data-tip]'); if(!t) return;
    tip.innerHTML=t.getAttribute('data-tip'); tip.style.opacity='1';
  });
  document.addEventListener('mousemove',function(e){
    if(tip.style.opacity!=='1') return;
    var x=e.clientX+14, y=e.clientY-30;
    if(x+tip.offsetWidth>window.innerWidth-8) x=e.clientX-tip.offsetWidth-14;
    tip.style.left=x+'px'; tip.style.top=Math.max(8,y)+'px';
  });
  document.addEventListener('mouseout',function(e){
    if(e.target.closest('[data-tip]')) tip.style.opacity='0';
  });
})();
"""


def join_names(names):
    """Join a list into prose: a, b and c."""
    names = list(names)
    if len(names) <= 1:
        return "".join(names)
    return ", ".join(names[:-1]) + " and " + names[-1]


def git_revision():
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return "unknown"


def build_report(oracle_rows, agg_rows, toi_rows, title):
    # --- accuracy, from the oracle -----------------------------------------
    modes = sorted({canonical_mode(r["mode"]) for r in oracle_rows}, key=mode_sort_key)
    compare_modes = [m for m in modes if m != "ti-reference"]
    groups = []
    for r in oracle_rows:
        g = f"{r['dataset']} {r['phase']}"
        if g not in groups:
            groups.append(g)

    by_key = {}
    for r in oracle_rows:
        by_key[(f"{r['dataset']} {r['phase']}", canonical_mode(r["mode"]))] = r

    total_queries = sum(as_int(r, "queries") for r in oracle_rows
                        if canonical_mode(r["mode"]) == (compare_modes[0] if compare_modes else ""))
    # The invariant is stated against the true time of impact. The datasets ship
    # exact roots, so use those; TightInclusion's own answer is a conservative
    # lower bound on the truth and flags safe results as late.
    have_truth = any(as_int(r, "gt_checked") for r in oracle_rows)
    violations = {m: 0 for m in compare_modes}
    missed = {m: 0 for m in compare_modes}
    late_vs_ti = {m: 0 for m in compare_modes}
    for r in oracle_rows:
        m = canonical_mode(r["mode"])
        if m not in violations:
            continue
        late_vs_ti[m] += as_int(r, "late")
        if have_truth:
            violations[m] += as_int(r, "gt_missed") + as_int(r, "gt_late")
            missed[m] += as_int(r, "gt_missed")
        else:
            violations[m] += as_int(r, "false_negative") + as_int(r, "late")
            missed[m] += as_int(r, "false_negative")

    safe = [m for m in compare_modes if violations.get(m, 0) == 0]
    unsafe = [m for m in compare_modes if violations.get(m, 0) > 0]

    parts = []

    # --- masthead ----------------------------------------------------------
    parts.append('<div class="wrap"><header class="masthead">')
    parts.append('<p class="eyebrow">SCCD &middot; narrow phase</p>')
    parts.append(f"<h1>{html.escape(title)}</h1>")
    reference = ("the exact roots shipped with each dataset" if have_truth
                 else "TightInclusion")
    parts.append("<p>Every narrow-phase mode measured query by query against "
                 f"{reference}. A mode is usable only if it never misses a collision and never "
                 "reports a time of impact later than the true one; speed is the second "
                 "question, not the first.</p>")
    parts.append('<div class="meta">')
    parts.append(f"<span>generated {_dt.datetime.now().strftime('%Y-%m-%d %H:%M')}</span>")
    parts.append(f"<span>revision {html.escape(git_revision())}</span>")
    if total_queries:
        parts.append(f"<span>{fmt_int(total_queries)} queries per mode</span>")
    parts.append("</div></header>")

    # --- verdict -----------------------------------------------------------
    if oracle_rows:
        ok = not unsafe
        cls = "pass" if ok else "fail"
        if ok:
            head = "Every mode measured is conservative."
            body = ("No missed collisions and no time of impact later than the true one, "
                    f"across {fmt_int(total_queries)} queries per mode.")
        else:
            head = (f"{join_names(unsafe)} " +
                    ("is" if len(unsafe) == 1 else "are") + " not conservative.")
            worst = max(unsafe, key=lambda m: violations[m])
            body = (f"{fmt_int(violations[worst])} of {fmt_int(total_queries)} queries break the "
                    f"invariant in <strong>{html.escape(worst)}</strong>"
                    + (f", including {fmt_int(missed[worst])} missed collisions" if missed[worst] else "")
                    + ". "
                    + (f"{join_names(safe)} passed cleanly." if safe else ""))
        parts.append(f'<div class="verdict {cls}"><div><h2>{head}</h2><p>{body}</p></div></div>')

        parts.append('<div class="tiles">')
        parts.append(f'<div class="tile"><span class="n">{fmt_int(total_queries)}</span>'
                     f'<span class="k">queries checked</span></div>')
        for m in compare_modes:
            v = violations.get(m, 0)
            state = "is-good" if v == 0 else "is-bad"
            parts.append(f'<div class="tile {state}"><span class="n">{fmt_int(v)}</span>'
                         f'<span class="k">{html.escape(m)} violations</span></div>')
        parts.append("</div>")

    # --- accuracy table ----------------------------------------------------
    if oracle_rows:
        parts.append('<section><h2>Conservativeness</h2>')
        parts.append("<p>A violation is a missed collision, or a time of impact after the true one. "
                     "False positives are listed but are not failures: reporting a contact that is "
                     "not there, or one earlier than it happens, costs work and never safety.</p>")
        if have_truth:
            parts.append("<p><em>Late vs TI</em> counts results falling between TightInclusion's "
                         "answer and the true root. TightInclusion is itself conservative, so those "
                         "are inside the safe band and are shown only for context.</p>")
        rows = []
        for m in compare_modes:
            v = violations.get(m, 0)
            pill = ('<span class="pill ok">passes</span>' if v == 0
                    else '<span class="pill no">fails</span>')
            fp = sum(as_int(by_key[(g, m)], "false_positive") for g in groups if (g, m) in by_key)
            rel = [as_float(by_key[(g, m)], "relerr_median") for g in groups if (g, m) in by_key]
            rel_med = sorted(rel)[len(rel) // 2] if rel else 0.0
            rows.append([f'<span class="name">{html.escape(m)}</span>', pill,
                         fmt_int(missed.get(m, 0)),
                         fmt_int(v - missed.get(m, 0)),
                         fmt_int(late_vs_ti.get(m, 0)),
                         fmt_int(fp), fmt_sci(rel_med)])
        head = "".join(f"<th>{h}</th>" for h in
                       ["mode", "status", "missed", "late vs truth", "late vs TI",
                        "false positives", "median TOI error"])
        body = "".join("<tr><td class='name'>" + r[0] + "</td><td>" + r[1] + "</td>" +
                       "".join(f"<td>{c}</td>" for c in r[2:]) + "</tr>" for r in rows)
        parts.append(f'<div class="card"><div class="scroll"><table><thead><tr>{head}</tr></thead>'
                     f"<tbody>{body}</tbody></table></div></div></section>")

    # --- violations chart --------------------------------------------------
    if oracle_rows and len(compare_modes) > 1:
        vals = {}
        for g in groups:
            for m in compare_modes:
                r = by_key.get((g, m))
                if r is None:
                    vals[(g, m)] = None
                elif have_truth:
                    vals[(g, m)] = as_int(r, "gt_missed") + as_int(r, "gt_late")
                else:
                    vals[(g, m)] = as_int(r, "false_negative") + as_int(r, "late")
        parts.append('<section><h2>Where the violations are</h2>')
        parts.append("<p>Queries per dataset and phase whose result is unsafe. Zero is the only "
                     "acceptable height.</p>")
        parts.append('<div class="card">')
        parts.append(svg_grouped_bars(groups, compare_modes, vals, unit="queries",
                                      label_fmt=lambda v: fmt_int(round(v))))
        parts.append(legend(compare_modes))
        parts.append(data_table(["dataset / phase"] + compare_modes,
                                [[html.escape(g)] + [fmt_int(vals.get((g, m)) or 0) for m in compare_modes]
                                 for g in groups],
                                caption="View as table"))
        parts.append("</div></section>")

    # --- timing chart ------------------------------------------------------
    if oracle_rows and any(as_float(r, "ms") for r in oracle_rows):
        tvals = {(f"{r['dataset']} {r['phase']}", canonical_mode(r["mode"])): as_float(r, "ms")
                 for r in oracle_rows}
        tmodes = [m for m in modes]
        parts.append('<section><h2>Cost of each mode</h2>')
        parts.append("<p>Wall time for the same query set, including TightInclusion itself. "
                     "Log scale: the modes differ by more than an order of magnitude.</p>")
        parts.append('<div class="card">')
        parts.append(svg_grouped_bars(groups, tmodes, tvals, unit="ms (log)", log=True, height=270))
        parts.append(legend(tmodes))
        parts.append(data_table(["dataset / phase"] + tmodes,
                                [[html.escape(g)] + [fmt_ms(tvals.get((g, m)) or 0) for m in tmodes]
                                 for g in groups],
                                caption="View as table"))
        parts.append("</div></section>")

    # --- pipeline timing from the benchmark driver -------------------------
    if agg_rows:
        agg_modes = sorted({r.get("mode") or "scalar" for r in agg_rows}, key=mode_sort_key)
        agg_groups = []
        for r in agg_rows:
            g = r.get("dataset_name") or r["dataset"]
            if g not in agg_groups:
                agg_groups.append(g)
        vals = {}
        for r in agg_rows:
            g = r.get("dataset_name") or r["dataset"]
            vals[(g, r.get("mode") or "scalar")] = as_float(r, "narrow_ms_median")
        parts.append('<section><h2>Narrow phase in the full pipeline</h2>')
        parts.append("<p>Median narrow-phase time per step from the benchmark driver, over the "
                     "broadphase output rather than a fixed query list.</p>")
        parts.append('<div class="card">')
        parts.append(svg_grouped_bars(agg_groups, agg_modes, vals, unit="ms / step"))
        parts.append(legend(agg_modes))
        parts.append("</div></section>")

    # --- mode reference ----------------------------------------------------
    parts.append('<section><h2>The modes</h2>')
    parts.append('<div class="modes">')
    for i, m in enumerate(modes or list(MODE_BLURB)):
        pill = ('<span class="pill ref">reference</span>' if m == "ti-reference"
                else '<span class="pill ok">conservative</span>' if m in CONSERVATIVE_MODES
                else '<span class="pill no">not conservative</span>')
        env = {"scalar": "0", "fast-vector": "1", "conservative": "2", "ti-compat": "3"}.get(m)
        code = f"<code>SCCD_NARROWPHASE_MODE={env}</code>" if env else "<code>&mdash;</code>"
        parts.append(f'<div class="mode">{code}<p>{html.escape(MODE_BLURB.get(m, ""))}</p>{pill}</div>')
    parts.append("</div></section>")

    parts.append("<footer><p>Generated by <code>benchmark/bench_report_html.py</code>. "
                 "Accuracy from <code>ti_oracle</code>, timing from <code>sccd_bench</code>. "
                 "Series colours are checked for colour-vision separation and surface contrast; "
                 "every chart has a table view.</p></footer>")
    parts.append("</div><div id='tip'></div>")
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aggregate")
    ap.add_argument("--toi-error")
    ap.add_argument("--oracle")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="Narrow Phase Conformance")
    ap.add_argument("--artifact", action="store_true",
                    help="emit only the page content (title, style, body) for hosts that "
                         "supply their own document skeleton")
    args = ap.parse_args()

    oracle_rows = read_csv(args.oracle)
    agg_rows = read_csv(args.aggregate)
    toi_rows = read_csv(args.toi_error)

    if not (oracle_rows or agg_rows):
        print("bench_report_html: no input CSVs found; nothing to report")
        return 1

    body = build_report(oracle_rows, agg_rows, toi_rows, args.title)
    fonts = ("<link rel='preconnect' href='https://fonts.googleapis.com'>"
             "<link rel='preconnect' href='https://fonts.gstatic.com' crossorigin>"
             "<link rel='stylesheet' href='https://fonts.googleapis.com/css2?"
             "family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&"
             "family=IBM+Plex+Serif:wght@500;600&display=swap'>")

    if args.artifact:
        page = (f"<title>{html.escape(args.title)}</title>{fonts}"
                f"<style>{STYLE}</style>{body}<script>{SCRIPT}</script>")
    else:
        page = ("<!doctype html><html lang='en'><head><meta charset='utf-8'>"
                "<meta name='viewport' content='width=device-width,initial-scale=1'>"
                f"<title>{html.escape(args.title)}</title>{fonts}"
                f"<style>{STYLE}</style></head><body>{body}"
                f"<script>{SCRIPT}</script></body></html>")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    print(f"bench_report_html: wrote {out} ({len(page) / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
