# Additive CCD

`accd_bench.exe.cpp` runs `ipc::AdditiveCCD` from the IPC toolkit over the
benchmark, the same cases and CSV columns as `sccd_bench`, for the per-pair
comparison. Cost is measured over SCCD's sweep broad-phase candidates, so both
narrow phases time the same list on the same cores -- the loop over candidates is
the `tbb::parallel_for` the toolkit's own stepsize query uses, because the method
itself is a per-pair kernel; accuracy is scored on every curated query at
the coordinates its exact root was computed for, as `sccd_bench` scores SCCD.

`conservative_rescaling` stays at the toolkit's default of `0.9`, the value its
header gives for comparison with the method's paper. `min_distance` is `0`.

The method stops short of contact by design, so `toi_late` must be zero and the
binary treats anything else as a harness error. `toi_med_early` and
`toi_max_early` measure how far short it stops.

SCCD's broad phase reports indices into its own edge list, which is not the
benchmark's ordering. The harness maps them through the same edge-id table
`bench.exe.cpp` builds; passing unmapped indices hands the method edges that do
not exist, and it responds with degenerate-input warnings and non-convergence.
