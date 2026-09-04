from __future__ import annotations

import argparse
from pathlib import Path


_ASCII_REPLACEMENT_TABLE = bytes(range(128)) + (b" " * 128)


def replace_nonascii_header(path: str | Path) -> int:
    path = Path(path)
    header = bytearray()

    with path.open("r+b") as file:
        for line in file:
            header.extend(line)
            if line.strip() == b"end_header":
                break
        else:
            raise ValueError(f"{path}: missing PLY end_header")

        filtered = header.translate(_ASCII_REPLACEMENT_TABLE)
        if filtered != header:
            file.seek(0)
            file.write(filtered)

    return sum(byte >= 128 for byte in header)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replace non-ASCII bytes in PLY headers in place."
    )
    parser.add_argument("paths", nargs="+", help="PLY files to update")
    args = parser.parse_args()

    for path in args.paths:
        replace_nonascii_header(path)


if __name__ == "__main__":
    main()
