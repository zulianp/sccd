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
- Withdrawn: "the quad device kernel's local stack costs it 255 registers and a
  spill" — see below

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
`CMakeLists.txt` does not reference it, which is what I checked; `benchmark/scripts/bench.sh`
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

## 7. Reading a symptom off the wrong line of the compiler's output

`narrow_phase_vq_kernel` reported, at `-Xptxas -v`, an 8128-byte stack frame,
112 bytes of spill stores, 216 of spill loads and 255 registers — the only
kernel in the build that spilled. Its `Domain stack[140]` was 7840 of those
8128 bytes, so the diagnosis wrote itself: the thread-local stack is the
problem, move it to a block-shared pool with a global overflow queue the way
the triangle kernel does. That was the plan, and it was approved.

It was wrong on every count, and one compile settled it. Rebuilding at stack
caps from 4 to 513 entries:

    cap    frame      registers   spill
      4      448 B    254           0/0
      8      736 B    255       112/216
     16     1184 B    255       112/216
     32     2080 B    255       112/216
    140     8128 B    255       112/216
    513    29024 B    255       112/216

The registers and the spill do not move with the array at all. They are the
inclusion function's working set — thirty coordinates, two `Frame`s, eight
corner triples — and no stack change touches them. The frame is the array, and
the frame is not the spill; they had been read as one number.

Nor does the frame cost time. Only the entries a search actually reaches are
ever written; the rest is address space. On GH200, 400k queries at depth 69:
cap 140 runs 7.589 ms, cap 257 runs 7.441 ms, cap 513 runs 7.355 ms — the
largest is the fastest, which is to say it is noise.

**The near-miss.** A first timing run appeared to confirm the restructure
handsomely: cap 140 at 7.613 ms against cap 32 at 3.681 ms, a 2.07x win for the
smaller stack. It measured the wrong thing. `kMaxDeviceDepth` is *derived* from
the cap, so shrinking the array shortens the search: cap 32 means depth 15, not
depth 69. Holding depth fixed at 15, cap 140 and cap 32 run 3.719 ms and
3.680 ms to a bit-identical answer. This is failure mode 5 again, measuring the
wrong object entirely, and it would have "validated" a rewrite that could not
have delivered what it promised.

What the exercise did produce is the opposite change from the one planned.
Since depth is what costs and headroom for it is free, the stack is now sized
*from* a depth (`SCCD_VQ_MAX_DEPTH`, 128) rather than being an entry count that
a depth is derived from. The old 140 entries capped the device at 69, which
happens to equal the host's default `SCCD_MAX_DEPTH` — so the two agreed by
coincidence, and raising `SCCD_MAX_DEPTH` gave a host that searched deeper and
a device that quietly did not. The clamp is safe, since exhaustion accepts at
the box's `t` lower bound, but it is an accuracy divergence and it is now
reported on stderr instead of being silent.

The rule to carry forward: **when a plan rests on a number, re-derive the
number before executing the plan, and check that the knob you are varying moves
only the thing you think it moves.**

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
