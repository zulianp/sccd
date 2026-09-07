"""
Design tokens and matplotlib defaults, in one place.

The three report generators that came before this shared no code, so a figure's
appearance was whatever the function drawing it happened to do: three aspect
ratios across one figure set, no font control anywhere, and `tab10` sampled
continuously over [0, 1] -- which for three series picks entries 0, 5 and 9 and
throws away the separation the palette was designed to have. At the size a
7.0 x 3.8 inch figure is reproduced in a journal column, the default tick labels
land near 4 pt.

Everything here is set once and applied to every figure.
"""

from __future__ import annotations

# --- page geometry ---------------------------------------------------------
# SIAM's article class sets a 4.8 inch text column in the single-column layouts
# these figures target. Drawing at the width it will be printed at is the whole
# point: a figure scaled down in LaTeX takes its fonts down with it.
COLUMN_WIDTH_IN = 4.8
FULL_WIDTH_IN = 6.5
GOLDEN = 0.618

def figsize(width_in: float = COLUMN_WIDTH_IN, ratio: float = GOLDEN) -> tuple[float, float]:
    return (width_in, width_in * ratio)


# --- colour ----------------------------------------------------------------
# Checked for colour-vision-deficiency separation and for contrast on white.
# Used as a discrete sequence -- index it, never sample a colormap over [0, 1].
SERIES = ["#0E8AA0", "#C07C10", "#8A3C74", "#3E6E4B", "#8C5A3C"]
REFERENCE_INK = "#7C8894"
GRID_INK = "#D6DBE0"
RULE_INK = "#333A41"

# Mode -> colour, fixed so a mode keeps its colour across every figure.
MODE_ORDER = ["relaxed", "tight", "device-relaxed", "device-tight", "tight-inclusion"]
MODE_COLOR = {
    "relaxed": SERIES[0],
    "tight": SERIES[1],
    "device-relaxed": SERIES[2],
    "device-tight": SERIES[3],
    "tight-inclusion": REFERENCE_INK,
}
MODE_MARKER = {
    "relaxed": "o",
    "tight": "s",
    "device-relaxed": "^",
    "device-tight": "D",
    "tight-inclusion": "x",
}

# Committed CSVs keep whatever mode name was current when they were written, so
# a recorded run keeps saying what was run and the mapping lives here instead.
# "vector" and "ti-compat" were pinned to narrow-phase modes 1 and 3, which no
# longer exist -- setting either warns and runs Relaxed -- so rows carrying those
# names are a second measurement of Relaxed and are labelled as such rather than
# folded into it.
MODE_ALIASES = {
    "scalar": "relaxed",
    "fast": "relaxed",
    "conservative": "tight",
    "ti-vec": "tight",
    "vector": "relaxed (retired mode 1)",
    "fast-vector": "relaxed (retired mode 1)",
    "ti-compat": "relaxed (retired mode 3)",
    "device": "device-relaxed",
    "device-fast": "device-relaxed",
    "device-ti": "device-tight",
    "ti-reference": "tight-inclusion",
}

MODE_LABEL = {
    "relaxed": "Relaxed",
    "tight": "Tight",
    "device-relaxed": "Relaxed (GPU)",
    "device-tight": "Tight (GPU)",
    "tight-inclusion": "TightInclusion",
}

SCENE_LABEL = {
    "armadillo-rollers": "armadillo-rollers",
    "cloth-ball": "cloth-ball",
    "cloth-funnel": "cloth-funnel",
    "puffer-ball": "puffer-ball",
    "n-body-simulation": "n-body",
    "rod-twist": "rod-twist",
}


def normalise_mode(name: str) -> str:
    key = (name or "").strip().lower()
    return MODE_ALIASES.get(key, key)


def mode_color(mode: str) -> str:
    return MODE_COLOR.get(mode, SERIES[hash(mode) % len(SERIES)])


def mode_label(mode: str) -> str:
    return MODE_LABEL.get(mode, mode)


def apply_rcparams() -> None:
    """Set every figure's typography and furniture. Call once, before plotting."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import rcParams

    rcParams.update({
        # Type. A serif face to sit with the body text of the article, and
        # mathtext from the same family so an axis label with a symbol in it
        # does not change typeface halfway through.
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,

        # Furniture. The identification of a figure belongs in its \caption, so
        # axes carry no title; a title drawn into the image duplicates the
        # caption and cannot be referenced.
        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": GRID_INK,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "axes.edgecolor": RULE_INK,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,

        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.columnspacing": 1.2,

        "lines.linewidth": 1.4,
        "lines.markersize": 4,

        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
