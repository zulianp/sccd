"""
Convert a query archive of the 2020 CCD benchmark (NYU handle 2451/61518) into a
scene directory that ti_oracle reads.

    python3 ccd_queries_to_scene.py <archive.tar.gz> <out-scene-dir>

The archives are the Handcrafted dataset (2451/61519), the Simulation dataset
(2451/61520, e.g. cow-heads) and their Rounded versions (2451/63808). Each holds
`[prefix/]<sub-scene>/{vertex-face,edge-edge}/<kind>-NNNN.csv`, where a query is
eight rows -- the four points at t=0, then at t=1 -- of

    x_num,x_den,y_num,y_den,z_num,z_den,truth

with `truth` the exact answer, the same on all eight rows. Their README says it
is 0 or 1, but the files also hold -1, 2 and 3 (cow-heads, the Rounded
handcrafted set). The authors' own reader (CCD-Wrapper, read_rational_csv)
prints an error for those and stores the value in a std::vector<bool>, so any
non-zero value is a collision there, and it is one here: for checking missed
collisions that is the strict reading. The raw value is kept beside the
answer, and the conversion reports how many queries carry each value.

There are no times of impact, so these sets can check missed collisions and
false positives but not lateness.

The scene written is the layout of the 2451/74508 datasets:

    queries/<sub-scene>-NNNN{vf,ee}.csv      the first six columns
    mma_bool/<same key>/mma_bool.uint8       one byte per query: 1 if truth != 0
    mma_bool/<same key>/truth_code.int8      one byte per query: truth as written

The archive is read in place; nothing is extracted. Any query whose rows disagree
on the truth stops the conversion, since there is then no answer to check against.

It also counts the coordinates a double cannot hold exactly (a numerator above
2^53, or a denominator that is not a power of two). ti_oracle divides numerator by
denominator in double, so for such a set a "missed collision" may be a collision
of the rounded geometry, not of the one the truth describes. The Rounded datasets
exist for that reason; their counts should be zero.
"""
from __future__ import annotations

import sys
import tarfile
from collections import Counter
from pathlib import Path

KINDS = {"vertex-face": "vf", "edge-edge": "ee"}


def exactly_double(num: int, den: int) -> bool:
    return abs(num) <= (1 << 53) and den > 0 and (den & (den - 1)) == 0


def convert(archive: Path, out: Path) -> int:
    n_files = n_queries = n_true = inexact = 0
    codes = Counter()
    seen = set()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar:
            parts = member.name.split("/")
            if not member.isfile() or not parts[-1].endswith(".csv") or parts[-1].startswith("._"):
                continue
            if len(parts) < 3 or parts[-2] not in KINDS:
                print(f"error: unexpected path {member.name}", file=sys.stderr)
                return 1
            sub, kind, name = parts[-3], KINDS[parts[-2]], parts[-1]
            number = name.rsplit("-", 1)[-1][: -len(".csv")]
            key = f"{sub}-{number}{kind}"
            if key in seen:
                print(f"error: two files map to {key}", file=sys.stderr)
                return 1
            seen.add(key)

            rows = [r for r in tar.extractfile(member).read().decode().splitlines() if r.strip()]
            if not rows or len(rows) % 8:
                print(f"error: {member.name} has {len(rows)} rows, not a multiple of 8", file=sys.stderr)
                return 1
            coords, truth, raw = [], bytearray(), bytearray()
            for q in range(len(rows) // 8):
                answers = set()
                for row in rows[8 * q : 8 * q + 8]:
                    cols = row.split(",")
                    if len(cols) != 7:
                        print(f"error: {member.name}: expected 7 columns, got {len(cols)}", file=sys.stderr)
                        return 1
                    nums = [int(c) for c in cols[:6]]
                    inexact += sum(not exactly_double(nums[i], nums[i + 1]) for i in (0, 2, 4))
                    coords.append(",".join(cols[:6]))
                    answers.add(int(cols[6]))
                if len(answers) != 1:
                    print(f"error: {member.name} query {q}: truth column {sorted(answers)}", file=sys.stderr)
                    return 1
                code = answers.pop()
                codes[code] += 1
                raw.append(code & 0xFF)
                truth.append(1 if code != 0 else 0)

            (out / "queries").mkdir(parents=True, exist_ok=True)
            (out / "mma_bool" / key).mkdir(parents=True, exist_ok=True)
            (out / "queries" / f"{key}.csv").write_text("\n".join(coords) + "\n")
            (out / "mma_bool" / key / "mma_bool.uint8").write_bytes(bytes(truth))
            (out / "mma_bool" / key / "truth_code.int8").write_bytes(bytes(raw))
            n_files += 1
            n_queries += len(truth)
            n_true += sum(truth)

    print(f"{archive.name} -> {out}: {n_files} files, {n_queries} queries, "
          f"{n_true} colliding; {inexact} coordinate(s) not exactly representable in double; "
          f"queries per truth value: {dict(sorted(codes.items()))}")
    return 0 if n_files else 1


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[3], file=sys.stderr)
        return 2
    return convert(Path(argv[1]), Path(argv[2]))


if __name__ == "__main__":
    sys.exit(main(sys.argv))
