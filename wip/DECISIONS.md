# Decision record

Claims this project made, measured, and in several cases withdrew. Kept because a
retraction is only useful if it is as findable as the claim was, and because
every one of these was believed at the time on evidence that looked sufficient.

## Withdrawn or corrected

- [Withdrawn: "mode 2 is about 100× slower on armadillo edge-edge"](ASSESSMENT.md)
- [Retracted: "the device narrow phase loses on every scene"](ASSESSMENT.md)
- [Fixed on the way past: an unsound rejection in the device's mode-0 kernel](ASSESSMENT.md)
- [Withdrawn: "mode 2's earliest impact is late on armadillo-rollers"](ASSESSMENT.md)
- [Retracted: every `fn=0` in this document before the matrix above](ASSESSMENT.md)
- [Resolved: the sweep dropped touching pairs](ASSESSMENT.md)
- Reversed: demoting `external/json` to `spikes/` — see below

The full argument and the numbers behind each sit in
[`ASSESSMENT.md`](ASSESSMENT.md).

## What kept going wrong

Five distinct failures, all of the same shape: a measurement that answered a
different question than the one being asked.

1. **Comparing across code paths.** "The device narrow phase loses by 87x" raced
   the host's *fastest* kernel against the device's *slowest*, because
   `SCCD_NARROWPHASE_MODE=2` names a different kernel on each side.
2. **Comparing across populations.** "The device does 520x the host's work" was
   an average over two paths differing by 15x; split by `toi_stride` it is 94x on
   the interface anyone calls and 1397x on a path handling 274 queries.
3. **Comparing against a reference that does not apply.** "Mode 2 reports late on
   armadillo-rollers" compared a mesh-path result against roots computed for the
   query geometry. smesh stores coordinates as `float`, so those are different
   numbers, and only armadillo's coordinates need more than 24 mantissa bits.
4. **Comparing against nothing at all.** Every `fn=0` in this branch before the
   full matrix was vacuous: that data tree ships no exact roots, so `expected`
   was false for every query and a false negative was impossible by construction.

Two independent instruments agreeing is not enough — in (2) they agreed because
both were making the same mistake. The check that has actually caught things is
asking what the number would look like if the code were wrong.

## 5. Measuring the wrong object entirely

Added after the fact, because it nearly cost a history rewrite across 12
branches.

"`git clone` of this repository transfers 4.6 GiB" was recorded as the
top-priority item, with a filter-repo plan and a force-push sequence. It came
from `git count-objects -vH` on a local working clone and from
`git clone --local .`, which hardlinks the whole local object store whether or not
it is reachable. A real clone from the URL is **24 MB**; the remote history holds
2 `data/` objects against 24,963 locally.

The 4.56 GiB was local, and nearly all of it was pinned by one stale
`refs/codex/...` ref left by tooling. Dropping it and running `git gc` gave
24.03 MiB with every branch intact.

Same shape as the other four: a number that answered a different question than
the one being asked. It survived longer than the others because it was never
challenged — it was measured once, written down with a table of largest blobs
that made it look thoroughly established, and cited unchallenged for several
turns.

## 6. Checking one caller and calling it none

`external/json` was demoted to `spikes/` in `e313109` on the stated grounds that
it was "never `add_subdirectory`'d, so it was not in the build at all". That
sentence is true and the conclusion drawn from it is not. The main
`CMakeLists.txt` does not reference it, which is what I checked; `benchmark/bench.sh`
configures it as a **standalone CMake project**, which I did not.

Its two programs turn the datasets' `boxes/*.json` and `mma_bool/*.json` into the
raw arrays `bench.exe.cpp` reads — `boxes/<key>/c0.int32`, `c1.int32` and
`mma_bool/<key>/mma_bool.uint8`. Without them `sccd_bench` has no expected pair
sets and no hit/miss booleans, so the whole accuracy half of the benchmark is
gone. And `bench.sh` did not degrade gracefully: it configures that path
unconditionally, so from `e313109` until now the harness failed on its second
step with "The source directory .../external/json does not exist".

It is now `benchmark/json/`, beside its only caller. It stays a separate CMake
project rather than joining the main one, because it fetches simdjson at
configure time and the shipped library must build with nothing fetched.

The shape is the same as the other five, one level up: the *check* answered a
different question than the one being asked. "Is it in the build?" is not "does
anything use it?", and a grep of `CMakeLists.txt` cannot tell them apart.

## A note on the C ABI

Both are the kind of thing that only shows up once something is actually tested,
which this ABI was not until `src/tests/c_abi_test.exe.cpp` was written.

**It shipped the losing splitter.** Every C entry point defaulted
`SCCD_ADAPTIVE_SPLIT` to `0`, so the installed ABI ran uniform interval splitting
while the C++ path ran the adaptive splitter. The assessment measured uniform
behind adaptive on every real scene, and it has since been demoted to a spike.
The C entry points now use the adaptive splitter, like everything else.

**The bisecting variants could miss a collision.** At the subdivision cap, the
search discarded the box instead of accepting it — an unsound rejection, and the
one way this algorithm can lose a root. It was reachable: a vertex dropping
straight through a stationary triangle was missed by
`sccd_find_root_bisection_vf_d` while `sccd_find_root_bisection_vf_f` found it,
because in single precision the box met a tolerance condition before the cap and
in double it did not. Exhaustion now accepts at the box's `t` lower bound, as
every other termination path in that loop already did.
