"""
Refresh the generated blocks of a Markdown document in place.

docs/BENCHMARKS.md is part prose and part measurement, and the two want opposite
treatment: the prose is written and should not be regenerated, the numbers are
measured and should never be typed. So the document carries named blocks

    <!-- sccd:begin timing -->
    ...anything here is replaced...
    <!-- sccd:end timing -->

and this module rewrites the contents of each block from the sweep, leaving
everything outside them untouched.

That makes "regenerate it and diff" a real check rather than an aspiration: if
the committed CSV and the committed script do not reproduce the committed
document, the diff is non-empty and the claim that the numbers are reproducible
is false. `--check` does exactly that and returns non-zero, which is the form a
CI gate wants.
"""

from __future__ import annotations

import re
from pathlib import Path

_BLOCK = re.compile(
    r"(?P<begin><!--\s*sccd:begin\s+(?P<name>[A-Za-z0-9_-]+)\s*-->)"
    r"(?P<body>.*?)"
    r"(?P<end><!--\s*sccd:end\s+(?P=name)\s*-->)",
    re.DOTALL)


def block_names(text: str) -> list[str]:
    return [m.group("name") for m in _BLOCK.finditer(text)]


def render(text: str, blocks: dict[str, str]) -> tuple[str, list[str]]:
    """
    Replace each named block's body. Returns the new text and the names that the
    document asked for but the caller did not supply -- a block left stale is
    worse than one that is obviously missing, so the caller is told.
    """
    unknown: list[str] = []

    def substitute(match: re.Match[str]) -> str:
        name = match.group("name")
        if name not in blocks:
            unknown.append(name)
            return match.group(0)
        body = blocks[name].strip("\n")
        return f"{match.group('begin')}\n\n{body}\n\n{match.group('end')}"

    return _BLOCK.sub(substitute, text), unknown


def apply(path: Path, blocks: dict[str, str], check: bool = False) -> int:
    """
    Rewrite `path`'s generated blocks. With `check`, write nothing and return 1
    if the file is not already what regeneration would produce.
    """
    path = Path(path)
    original = path.read_text()
    updated, unknown = render(original, blocks)

    if unknown:
        print(f"warning: {path} has blocks nothing generated: "
              f"{', '.join(sorted(set(unknown)))}")

    unused = sorted(set(blocks) - set(block_names(original)))
    if unused:
        print(f"warning: generated blocks the document does not embed: "
              f"{', '.join(unused)}")

    if check:
        if updated != original:
            print(f"FAILED: {path} is not what the committed CSV and this script "
                  f"produce. Regenerate it, or explain the difference.")
            return 1
        print(f"OK: {path} matches what the committed data reproduces.")
        return 0

    if updated != original:
        path.write_text(updated)
        print(f"updated {path}")
    else:
        print(f"{path} was already current")
    return 0
