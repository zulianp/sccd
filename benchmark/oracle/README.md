# Narrow-phase accuracy oracle

`ti_oracle` compares each narrow-phase mode against TightInclusion, query by
query, on the CCD benchmark sets under `data/`. It is the acceptance gate for
the TI-accurate vectorized kernel: no kernel change lands without a run here.

## The invariant it enforces

SCCD must be **conservative**:

* a collision must never be missed (a false negative is a correctness failure);
* a reported time of impact must never be **later** than the true one, because a
  solver would then step straight through the contact;
* a false positive, or a time of impact that is too **early**, is acceptable --
  it only costs work.

The benchmark datasets ship exact roots (`data/<scene>/roots/<key>/toi.float64`,
NaN where there is no collision), so the invariant is checked directly:
`gt_missed` counts queries with a real root that a mode reported as no
collision, and `gt_late` counts times of impact after the true one. Those two
are the gate.

**Do not gate on TightInclusion.** Its answer is a *conservative lower bound* on
the true time of impact, not the truth, so a result landing between its answer
and the real root is inside the safe band and is not a violation. The `late`
column measures exactly that and over-reports; it is kept for context only. The
same caveat applies to `false_negative`: TightInclusion reports hits on queries
the exact roots say have no collision at all.

`ti_oracle` exits non-zero on `gt_missed` or `gt_late`, falling back to the
TightInclusion comparison only for datasets that ship no roots. Pass
`--no-strict` to report without failing.

## Build and run

```sh
cmake -S . -B build -DSCCD_ENABLE_TIGHT_INCLUSION=ON
cmake --build build -j --target ti_oracle
./build/ti_oracle data/cloth-funnel --csv benchmark/oracle/cloth-funnel.csv
```

Exits non-zero on `gt_missed` or `gt_late`. Pass `--no-strict` to report without
failing.

Recorded runs, throughput figures and the reasoning behind the gate are in
[`../../wip/ORACLE_NOTES.md`](../../wip/ORACLE_NOTES.md).
