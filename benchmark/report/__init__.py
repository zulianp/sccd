"""
One reporting module for the SCCD benchmark.

Three generators came before this and shared no code: bench_postprocess.py
(matplotlib with a PGFPlots fallback), bench_report_html.py (stdlib, inline SVG)
and assess_report.py (stdlib, Markdown). The pieces worth keeping were in the
wrong places -- the only authoritative map from every historical mode name onto
the current ones lived in the HTML generator, where the matplotlib code could not
import it, and the discipline that refuses to report a ratio smaller than the
measured noise lived in the assessment reporter, where the figures could not.

Everything is here now, and one run over one CSV emits every artifact: PDF
figures for LaTeX, PNGs of the same figures for Markdown, booktabs tables, and
the Markdown that goes into docs/BENCHMARKS.md.
"""

from . import data, figures, style, tables

__all__ = ["data", "figures", "style", "tables"]
