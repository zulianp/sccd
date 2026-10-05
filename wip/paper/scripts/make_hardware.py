#!/usr/bin/env python3
"""
Build the hardware table from what was queried on the compute node.

    python3 scripts/make_hardware.py

Reads benchmark/results/hardware-gh200.txt -- lscpu, numactl, sysfs and
cudaGetDeviceProperties as they answered on an Alps compute node -- and writes
generated/tables/tab-hardware.tex. It was the last table in the article typed by
hand, which matters more here than the count of rows suggests: the peak
arithmetic figures are derived from the clocks and unit counts in the same
table, and the article reads a ratio off them, so a stale clock would propagate
into a claim about which processor should be ahead.

Three things in the table are not queried and are named here as constants: the
width of a Neoverse V2 SVE2 pipe, the number of double-precision units on a
Hopper multiprocessor, and the memory technologies. They are architectural
facts, cited in the article to the vendor documents, and no interface on the
node reports them.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
SRC = PAPER.parent.parent / "benchmark" / "results" / "hardware-gh200.txt"

# Architecture, not measurement. A Neoverse V2 core issues four 128-bit SVE2
# operations per cycle and a Hopper multiprocessor has 64 FP64 units; both do a
# fused multiply-add, so each counts two flops.
SVE_PIPES = 4
FP64_UNITS_PER_SM = 64
FLOPS_PER_FMA = 2
WARP = 32
HOST_MEMORY_TECH = "LPDDR5X"
DEVICE_MEMORY_TECH = "HBM3"


def parse(text: str) -> tuple[dict[str, str], dict[str, str], dict[str, float]]:
    """The file's three blocks: the host block, the device block, the triad lines.

    Each block is `key<two or more spaces>value`, with the key optionally ending
    in a colon and optionally prefixed by `---` for the lines the collecting
    script added to what lscpu printed. The host and the device both have a key
    called `L2 cache`, which is why the two blocks stay separate.
    """
    host: dict[str, str] = {}
    device: dict[str, str] = {}
    triad: dict[str, float] = {}
    block = None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("##"):
            low = line.lower()
            block = host if "host" in low else device if "device" in low else None
            continue
        m = re.match(r"^(Grace|Hopper)\s+STREAM triad\s+([\d.]+) GB/s", line)
        if m:
            triad[m.group(1)] = float(m.group(2))
            continue
        if block is None or not line or line.startswith("#"):
            continue
        line = re.sub(r"^---\s*", "", line)
        m = re.match(r"^(.+?):\s+(.+)$", line) or re.match(r"^(.+?)\s{2,}(.+)$", line)
        if m:
            block[m.group(1).strip()] = m.group(2).strip()
    return host, device, triad


def number(s: str) -> float:
    return float(re.sub(r"[^\d.]", "", s))


def cache_each(spec: str) -> float:
    """`18 MiB (288 instances)` as MiB per instance."""
    m = re.match(r"([\d.]+)\s*MiB\s*\((\d+) instances\)", spec)
    if not m:
        raise ValueError(f"cannot read a per-instance cache size from {spec!r}")
    return float(m.group(1)) / int(m.group(2))


def mib(x: float) -> str:
    return f"${x:g}$ MiB" if x >= 1 else f"${x * 1024:g}$ KiB"


def main() -> int:
    if not SRC.is_file():
        print(f"error: {SRC} does not exist", file=sys.stderr)
        return 2
    host, device, triad = parse(SRC.read_text())

    cores = int(number(host["Core(s) per socket"]))
    host_ghz = number(host["CPU max MHz"]) / 1000.0
    sve_bits = int(number(host["SVE bytes"])) * 8
    l1_each = cache_each(host["L1d cache"])
    l2_each = cache_each(host["L2 cache"])
    l3_each = cache_each(host["L3 cache"])
    host_mem_gib = number(host["numa0 MB"]) / 1024.0
    host_tflops = (cores * host_ghz * 1e9 * SVE_PIPES
                   * (sve_bits // 64) * FLOPS_PER_FMA) / 1e12

    sms = int(number(device["multiprocessors"]))
    per_sm = int(number(device["max threads / SM"]))
    dev_ghz = number(device["clock rate (SM)"])
    dev_tflops = sms * FP64_UNITS_PER_SM * FLOPS_PER_FMA * dev_ghz * 1e9 / 1e12

    rows = [
        ("architecture", "Arm Neoverse V2",
         f"NVIDIA Hopper, sm\\_{device['compute capability'].replace('.', '')}"),
        ("parallel units", f"${cores}$ cores", f"${sms}$ multiprocessors"),
        ("maximum clock", f"${host_ghz:.2f}$ GHz", f"${dev_ghz:.2f}$ GHz"),
        ("vector width", f"${sve_bits}$-bit SVE2", f"${WARP}$-thread warp"),
        ("threads in flight", f"${cores}$", f"${sms * per_sm:,}$".replace(",", "{,}")),
        ("L1 / shared memory", f"{mib(l1_each)} per core",
         f"${number(device['shared mem / SM']):g}$ KiB per multiprocessor"),
        ("L2 cache", f"{mib(l2_each)} per core", f"${number(device['L2 cache']):g}$ MiB"),
        ("L3 cache", f"{mib(l3_each)}", "--"),
        ("memory", f"${host_mem_gib:.0f}$ GiB {HOST_MEMORY_TECH}",
         f"${number(device['global memory']):.0f}$ GiB {DEVICE_MEMORY_TECH}"),
        ("peak \\textsc{fp64}", f"${host_tflops:.1f}$ TFLOP/s",
         f"${dev_tflops:.1f}$ TFLOP/s"),
        ("triad bandwidth", f"${triad['Grace']:.0f}$ GB/s",
         f"${triad['Hopper']:,.0f}$ GB/s".replace(",", "{,}")),
    ]

    body = "\n".join(f"    {a} & {b} & {c} \\\\" for a, b, c in rows)
    out = PAPER / "generated" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    (out / "tab-hardware.tex").write_text(rf"""\begin{{table}}[htbp]
  \centering
  \caption{{The two processors of one \textsc{{GH200}} module, as reported by the
    operating system and the CUDA device properties on the compute nodes used.
    Peak \textsc{{fp64}} is derived from those clocks and the per-unit counts of
    the two architectures~\citep{{arm2022neoversev2,nvidia2022hopper}}: ${SVE_PIPES}$
    ${sve_bits}$-bit SVE2 operations per cycle on a core and ${FP64_UNITS_PER_SM}$
    double-precision units on a multiprocessor, each a fused multiply-add. The
    bandwidth row is the best of ten \textsc{{stream}} triad passes.}}
  \label{{tab:hardware}}
  \fittable{{%
  \begin{{tabular}}{{lll}}
    \toprule
     & Grace (host) & Hopper (device) \\
    \midrule
{body}
    \bottomrule
  \end{{tabular}}}}
  \par\smallskip\footnotesize Source:
  \texttt{{benchmark/results/hardware-gh200.txt}}
\end{{table}}
""")
    print(f"wrote tab-hardware.tex (host {host_tflops:.1f} TFLOP/s, "
          f"device {dev_tflops:.1f} TFLOP/s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
