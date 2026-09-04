# Three edge-edge queries the device narrow phase cannot finish cheaply

`worst-query-w{1,2,3}ee.csv`, in the datasets' own query format: four points at
the start of the step, four at the end, one row each, `numerator,denominator` per
axis. They came out of cloth-funnel's broad-phase candidates via
`SCCD_NP_WORST_CSV`, which dumps the geometry of the costliest queries by
per-query box count.

One query per call, `max_toi = 1`, `toi_stride = 1` — no shared bound on either
side, nothing else in flight:

| | host `Tight` | device `Tight` | | device reports |
|---|---:|---:|---:|---|
| `w1ee` | **55** | **53,015,756** | 964,000× | no collision |
| `w2ee` | 543 | 47,578,924 | 87,600× | no collision |
| `w3ee` | 645 | 2,466,808 | 3,800× | `toi = 0.163146973` |

`w1ee` has **no collision**. The host proves that in 55 boxes. The device reaches
the same answer after 53 million.

Both sides are conservative and both agree on the answer, so this is not a
correctness bug — it is the whole of the device narrow phase's remaining cost,
reduced to one query that runs in seconds.

## Running them

```sh
mkdir -p /tmp/worst/queries
cp benchmark/assessment/worst-query-w1ee.csv /tmp/worst/queries/w1ee.csv
sccd_np_trace /tmp/worst ee --file w1 --mode 2 --top 3            # host
sccd_np_trace /tmp/worst ee --file w1 --mode 2 --device --no-isolated --top 0
```

Build with `-DSCCD_NP_COUNT_BOXES` or the counts read zero.

## Where to look first

The host terminates a branch when the domain widths fall inside the per-axis
domain tolerances. 53 million boxes is roughly a binary tree to depth 26, which
is what a search that *never* satisfies that condition looks like — it runs until
the geometry stops it rather than until the tolerance does. So the first thing to
compare is the per-axis tolerance each side computes for this query:
`compute_edge_edge_tolerance` on the host against `load_query_and_tol`'s `atol` on
the device. They are supposed to be transcriptions of each other.
