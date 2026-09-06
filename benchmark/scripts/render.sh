#!/usr/bin/env bash
set -euo pipefail

# Render a Markdown document to a styled, self-contained HTML page.
#
#   benchmark/scripts/render.sh                     # docs/BENCHMARKS.md
#   benchmark/scripts/render.sh docs/API.md         # any document
#   benchmark/scripts/render.sh --link docs/X.md    # link the stylesheet
#
# The CSS is inlined by default, so the result is one file that can be mailed
# or opened from anywhere without carrying assets alongside. `--link` keeps a
# <link> to docs/assets/style.css instead, which is what you want while editing
# the stylesheet itself.
#
# The generated document is table-heavy, and two things the Markdown converter
# cannot do are applied afterwards: every table is wrapped so a wide one scrolls
# in its own box rather than pushing the page sideways, and the zero in a
# failure column is marked so it reads as the good case it is.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CSS="${ROOT_DIR}/docs/assets/style.css"

INLINE=1
INPUTS=()
for arg in "$@"; do
    case "${arg}" in
        --link) INLINE=0 ;;
        --inline) INLINE=1 ;;
        -h|--help) sed -n '3,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*) printf 'error: unknown option %s\n' "${arg}" >&2; exit 2 ;;
        *) INPUTS+=("${arg}") ;;
    esac
done
[[ "${#INPUTS[@]}" -eq 0 ]] && INPUTS=("${ROOT_DIR}/docs/BENCHMARKS.md")

# python-markdown, from wherever this machine keeps it. Not `markdown_py` off
# the PATH: on this workstation that resolves to a Python 2.7 build that cannot
# run the extensions.
pick_python() {
    if [[ -n "${PYTHON:-}" ]]; then printf '%s\n' "${PYTHON}"; return 0; fi
    local candidate
    for candidate in "${ROOT_DIR}/data/venv/bin/python" python3 python3.13 python3.12 python3.11; do
        command -v "${candidate}" >/dev/null 2>&1 || [[ -x "${candidate}" ]] || continue
        if "${candidate}" -c 'import markdown' >/dev/null 2>&1; then
            printf '%s\n' "${candidate}"; return 0
        fi
    done
    return 1
}

if ! PY="$(pick_python)"; then
    echo "error: no Python with the 'markdown' package. Install it with:" >&2
    echo "         python3 -m pip install markdown pygments" >&2
    exit 1
fi

[[ -f "${CSS}" ]] || { printf 'error: missing %s\n' "${CSS}" >&2; exit 1; }

for input in "${INPUTS[@]}"; do
    [[ -f "${input}" ]] || { printf 'error: %s does not exist\n' "${input}" >&2; exit 1; }
    output="${input%.md}.html"
    CSS_PATH="${CSS}" INLINE_CSS="${INLINE}" "${PY}" - "${input}" "${output}" <<'PY'
import html
import os
import re
import sys

import markdown

src, dst = sys.argv[1], sys.argv[2]
text = open(src, encoding="utf-8").read()

# The generated blocks are delimited by HTML comments; they are instructions to
# the report generator, not content, so they do not belong in the page.
text = re.sub(r"^<!--\s*sccd:(begin|end)\s+[\w-]+\s*-->\n?", "", text, flags=re.M)

title = "Document"
for line in text.splitlines():
    if line.startswith("# "):
        title = line[2:].strip()
        break

body = markdown.markdown(
    text,
    extensions=["tables", "fenced_code", "codehilite", "toc", "attr_list"],
    extension_configs={"codehilite": {"guess_lang": False, "noclasses": True}},
)

# A wide table has to scroll inside its own box, or it widens the page and every
# paragraph on it goes with it. python-markdown emits a bare <table>.
body = body.replace("<table>", '<div class="table-scroll"><table>')
body = body.replace("</table>", "</table></div>")

# In a column counting failures, 0 is the result worth seeing. Mark the zeros in
# the columns whose header says so, and anything non-zero there as its opposite.
FAIL_COLUMNS = ("late", "missed", "false neg.", "fn")


def mark_failure_cells(match: "re.Match[str]") -> str:
    table = match.group(0)
    header = re.search(r"<thead>.*?</thead>", table, re.S)
    if not header:
        return table
    names = [re.sub(r"<[^>]+>", "", c).strip().lower()
             for c in re.findall(r"<th[^>]*>(.*?)</th>", header.group(0), re.S)]
    watch = {i for i, n in enumerate(names) if n in FAIL_COLUMNS}
    if not watch:
        return table

    body_only = re.search(r"<tbody>.*?</tbody>", table, re.S)
    if not body_only:
        return table

    def do_row(row: "re.Match[str]") -> str:
        cells = re.findall(r"<td[^>]*>.*?</td>", row.group(0), re.S)
        out = []
        for i, cell in enumerate(cells):
            value = re.sub(r"<[^>]+>", "", cell).strip()
            if i in watch and value not in ("", "--"):
                cls = "zero" if value in ("0", "0.0") else "nonzero-bad"
                cell = cell.replace("<td", f'<td class="{cls}"', 1)
            out.append(cell)
        return "<tr>" + "".join(out) + "</tr>"

    marked = re.sub(r"<tr>(?:(?!</tr>).)*</tr>", do_row, body_only.group(0), flags=re.S)
    return table.replace(body_only.group(0), marked)


# Over the whole table, not the body: the column names are in <thead>, which
# is outside <tbody>, so matching on the body alone never found them and the
# marking silently did nothing.
body = re.sub(r"<table>.*?</table>", mark_failure_cells, body, flags=re.S)

# Figure captions: the paragraph straight after a paragraph holding only an
# image. Marked here rather than selected in CSS, because CSS cannot express
# "the paragraph after an image-only paragraph" and the obvious near-miss --
# `p > em:only-child` -- also matches ordinary mid-sentence emphasis, since
# :only-child counts element siblings and ignores the text around them.
body = re.sub(
    r"(<p><img[^>]*/?></p>)\s*<p>((?:(?!</p>).)*)</p>",
    lambda m: f'{m.group(1)}<p class="caption">{m.group(2)}</p>',
    body, flags=re.S)

# A paragraph that is nothing but emphasis is a caption too, where a document
# writes them that way.
body = re.sub(r"<p><em>((?:(?!</em>).)*)</em></p>",
              r'<p class="caption">\1</p>', body, flags=re.S)

css_path = os.environ["CSS_PATH"]
if os.environ.get("INLINE_CSS") == "1":
    style = "<style>\n" + open(css_path, encoding="utf-8").read() + "\n</style>"
else:
    rel = os.path.relpath(css_path, os.path.dirname(os.path.abspath(dst)))
    style = f'<link rel="stylesheet" href="{rel}">'

page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
{style}
</head>
<body>
{body}
</body>
</html>
"""
open(dst, "w", encoding="utf-8").write(page)
print(f"  {src} -> {dst}  ({len(page) / 1024:.0f} KB)")
PY
done
