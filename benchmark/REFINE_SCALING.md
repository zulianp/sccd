# Refinement scaling study

`sccd_refine_scaling` takes a closed surface, refines it repeatedly, and runs one
collision step per level. It answers a question neither `ti_oracle` nor
`bench.exe` can: how cost responds to surface element count. Both of those run
fixed query sets.

```sh
cmake -S . -B build -DSCCD_ENABLE_SMESH=ON -Dsmesh_DIR=<prefix>/lib/cmake/smesh
cmake --build build -j --target sccd_refine_scaling
./build/sccd_refine_scaling 4                    # generated cube, synthetic motion
./build/sccd_refine_scaling 4 <mesh_t0> <mesh_t1>   # a real pair of frames
```

Env: `SCCD_NARROWPHASE_MODE`, `SCCD_MAX_DEPTH`, `SCCD_TOL`, `SCCD_SCALE`,
`SCCD_BENCH_EXECUTION_SPACE=device`.

A surface mesh is refined directly; a volume mesh is refined and then skinned.
smesh's `refine` rejects `TRISHELL3`, the type `skin` produces, which is why the
volume is refined before skinning rather than after.

## What it measures, and why the columns are what they are

`vf_pairs` and `ee_pairs` are printed alongside the element count because they
need not grow at the same rate, and it is the pair count that drives the narrow
phase. `ns_per_pair` is printed so a change in total time can be attributed to
the pair count or to the per-pair cost. `toi` is printed because it is the check
that every level is solving the same problem -- see below.

## Result on real data: cloth-ball

`benchmark/ply_to_smesh.py` converts the benchmark frames, which are PLY, into the
raw-array folders smesh reads. It handles ascii, little- and big-endian binary,
because the scenes are not consistent -- cloth-funnel is little-endian, cloth-ball
is big-endian.

```sh
python3 benchmark/ply_to_smesh.py data/cloth-ball/frames/cloth_ball40.ply /tmp/cb-40
python3 benchmark/ply_to_smesh.py data/cloth-ball/frames/cloth_ball41.ply /tmp/cb-41
./build/sccd_refine_scaling 2 /tmp/cb-40 /tmp/cb-41
```

Frames 40 and 41, host, mode 0:

| level | faces | nodes | vf_pairs | ee_pairs | broad_ms | narrow_ms | ns/pair |
|---|---|---|---|---|---|---|---|
| 0 | 92,230 | 46,598 | 12,104 | 56,881 | 50 | 2.0 | 29 |
| 1 | 368,920 | 185,423 | 122,626 | 505,058 | 418 | 14.6 | 23 |
| 2 | 1,475,680 | 739,763 | 2,080,634 | 6,481,554 | 4,353 | 175 | 21 |

Three things, and the first two are the ones that matter for where to spend
effort.

**The broad phase dominates, by 25x.** At 1.5M triangles it is 4,353 ms against
the narrow phase's 175 ms. Every optimisation in this repository so far has gone
into the narrow phase; on a refined real scene that is the 4% part. Broken down
and diagnosed in `BROADPHASE.md`: the sweep scans ~2,100 candidates per pair it
finds, and there is an unused grid implementation already in the tree.

**Candidate pairs grow superlinearly in elements.** Faces go up 4x per level by
construction, but vertex-face pairs go up 10x then 17x, and edge-edge 9x then
13x. Refining a cloth that is already in near-contact packs many more elements
into the same proximity, so pair count scales closer to `n^1.7` than to `n`. A
narrow phase that is linear *in pairs* is therefore still superlinear in
elements, and quoting "linear" without saying in what is misleading.

**Per-pair cost is flat, and slightly improving with size**: 29, 23, 21 ns. So
the narrow phase itself has no superlinear term -- consistent with the synthetic
run below, and now confirmed on a real scene.

### Mode 2 is not the expensive one here

Same frames, same levels, mode 2 (TightInclusion-exact) against mode 0:

| level | mode 0 narrow_ms | mode 2 narrow_ms |
|---|---|---|
| 0 | 2.0 | 2.2 |
| 1 | 14.6 | 13.7 |
| 2 | 175.1 | 149.5 |

The conservative kernel is level with the scalar reference and slightly ahead at
the largest size. Worth stating plainly because the armadillo-rollers edge-edge
case in `benchmark/oracle/README.md` shows mode 2 losing badly, and it would be
easy to generalise from that: on this scene it costs nothing.

## Synthetic geometry: the narrow phase is linear in candidate pairs

Generated cube, host, mode 0, on an M-series laptop. Kept because it reaches
far larger pair counts than the real scene does:

| level | faces | vf_pairs | ee_pairs | broad_ms | narrow_ms | ns/pair |
|---|---|---|---|---|---|---|
| 0 | 48 | 1,104 | 2,200 | 2.4 | 0.3 | 83 |
| 1 | 192 | 18,240 | 39,808 | 2.3 | 1.8 | 31 |
| 2 | 768 | 294,144 | 656,776 | 3.9 | 23.3 | 25 |
| 3 | 3,072 | 4,715,520 | 10,589,560 | 27.4 | 394.6 | 26 |
| 4 | 12,288 | 75,485,184 | 169,759,960 | 302.6 | 5,259 | 21 |
| 5 | 49,152 | 1,207,910,400 | 2,717,471,128 | 8,234 | 112,472 | 29 |

**Per-pair cost is flat at 21-31 ns across five orders of magnitude in pair
count**, after the first level where the mesh is too small to amortize anything.
So the narrow phase is linear in candidate pairs and there is no hidden
superlinear term in it. Mode 1 measures the same shape (20-27 ns/pair).

## The synthetic motion degenerates the broad phase -- read the caveat

The pair counts above are `n^2/2`: at level 5, 49,152 faces give 1.2e9
vertex-face pairs. That is not a property of the broad phase, it is the motion.

A convex surface cannot self-collide under any gentle deformation, so making a
generated cube collide at all requires a violent motion, and any violent motion
gives every element a swept box spanning the model, so everything overlaps
everything. **Total cost therefore scales as `n^2` here, and a real simulation
step would not.** Use the two-mesh form for a representative workload.

Three motions were tried and the `toi` column caught two of them:

- **Pushing the halves through each other**: toi came out exactly 1, 1/2, 1/4,
  1/8, 1/16 over five levels. Not distant geometry meeting -- adjacent elements
  either side of a discontinuity interpenetrating at once, so contact time falls
  with element size and every level measures a different problem.
- **A large rigid rotation**: never collided (toi stayed 1), and pair counts grew
  15x per level against 4x for the elements.
- **Point reflection through the centre** (the default, `SCCD_SCALE=-0.5`): toi
  is 0.666667 at every level, which is what a resolution-independent contact time
  looks like. This is the one to keep.

## Mode 2 does not complete on this input

The TightInclusion-exact mode does not finish level 2 -- 950k pairs -- in 900 s,
where mode 0 takes 23 ms. Do not read that as a mode-2 scaling result. A cube
reflected through its own centre is a degenerate configuration, with faces landing
exactly coplanar, which is the classic worst case for an interval root finder:
boxes stay ambiguous and the search drives to its depth cap. It is a real input a
solver could encounter, and worth keeping in mind as a robustness case, but it
says nothing about mode 2's cost on ordinary geometry. Compare against
`ti_oracle --bench`, where mode 2 runs within 2x of mode 0 on real query sets.
