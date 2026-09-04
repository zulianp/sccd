// The SIMD AABB kernels against the scalar predicate they are an optimisation of.
//
// vaabb.hpp carries four implementations of the same test -- AVX-512, AVX2, NEON
// and a scalar fallback -- selected by preprocessor. Which one a build gets
// depends on SCCD_ENABLE_NATIVE_ARCH, so turning that option off silently swaps
// the kernel underneath the broad phase, and nothing checked that the answers
// match.
//
// They must agree exactly. `disjoint` is the definition; everything else is a
// faster way of computing it, and a lane that disagrees is a pair the broad phase
// either invents or, far worse, drops.
//
// The touching case gets its own attention. The predicate rejects only on a
// strict `amin > bmax`, so boxes that meet exactly are overlapping -- and a sweep
// whose candidate window disagreed with that on the sort axis silently dropped
// real pairs, which is a bug this project shipped. Any kernel that got the
// boundary wrong would do the same, so it is exercised deliberately rather than
// left to chance in random data.

#include "sccd_aabb.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace {

    int failures = 0;

    template <typename T>
    struct Boxes {
        std::vector<T> lo[3], hi[3];
        ptrdiff_t n = 0;
    };

    // A mixture on purpose: coordinates drawn from a small integer grid so exact
    // ties are common, plus zero-extent boxes. Continuous random coordinates
    // essentially never produce a tie, which is exactly the case that matters.
    template <typename T>
    Boxes<T> make_boxes(std::mt19937& rng, const ptrdiff_t n, const int grid) {
        std::uniform_int_distribution<int> pos(0, grid);
        std::uniform_int_distribution<int> ext(0, 2);
        Boxes<T> b;
        b.n = n;
        for (int d = 0; d < 3; ++d) {
            b.lo[d].resize(n);
            b.hi[d].resize(n);
        }
        for (ptrdiff_t i = 0; i < n; ++i) {
            for (int d = 0; d < 3; ++d) {
                const T l = (T)pos(rng);
                b.lo[d][i] = l;
                b.hi[d][i] = l + (T)ext(rng);   // extent 0 is deliberate
            }
        }
        return b;
    }

    template <typename T>
    int check(const char* name, const ptrdiff_t n, const int grid) {
        std::mt19937 rng(2026);
        Boxes<T> a = make_boxes<T>(rng, n, grid);
        Boxes<T> b = make_boxes<T>(rng, n, grid);

        int mismatches = 0;
        long long overlaps = 0, touching = 0;

        // vaabb_overlap_one_to_many_bits works a chunk at a time and returns a
        // bitmask of the lanes that overlap, so it is compared lane by lane
        // against the scalar predicate on the same inputs.
        for (ptrdiff_t i = 0; i < a.n; ++i) {
            for (ptrdiff_t start = 0; start < b.n; start += SCCD_AABB_DISJOINT_CHUNK_SIZE) {
                const int len = (int)std::min<ptrdiff_t>(SCCD_AABB_DISJOINT_CHUNK_SIZE, b.n - start);

                T* bp[6] = {b.lo[0].data(), b.lo[1].data(), b.lo[2].data(),
                            b.hi[0].data(), b.hi[1].data(), b.hi[2].data()};

                const std::uint32_t bits = sccd::vaabb_overlap_one_to_many_bits<T>(
                    a.lo[0][i], a.lo[1][i], a.lo[2][i], a.hi[0][i], a.hi[1][i], a.hi[2][i],
                    bp[0] + start, bp[1] + start, bp[2] + start,
                    bp[3] + start, bp[4] + start, bp[5] + start, len);

                for (int lane = 0; lane < len; ++lane) {
                    const ptrdiff_t j = start + lane;
                    const bool vector_says = ((bits >> lane) & 1u) != 0;
                    const bool scalar_says = sccd::disjoint<T>(
                        a.lo[0][i], a.lo[1][i], a.lo[2][i], a.hi[0][i], a.hi[1][i], a.hi[2][i],
                        b.lo[0][j], b.lo[1][j], b.lo[2][j],
                        b.hi[0][j], b.hi[1][j], b.hi[2][j]) == 0;

                    if (scalar_says) ++overlaps;
                    for (int d = 0; d < 3; ++d) {
                        if (a.lo[d][i] == b.hi[d][j] || b.lo[d][j] == a.hi[d][i]) { ++touching; break; }
                    }

                    if (vector_says != scalar_says) {
                        if (mismatches < 5) {
                            std::printf("    lane disagrees: box %ld vs %ld -- vector says %s, "
                                        "scalar says %s\n",
                                        (long)i, (long)j,
                                        vector_says ? "overlap" : "disjoint",
                                        scalar_says ? "overlap" : "disjoint");
                        }
                        ++mismatches;
                    }
                }
            }
        }

        std::printf("  %-30s n=%-6ld overlaps=%-9lld touching=%-8lld %s\n",
                    name, (long)n, overlaps, touching,
                    mismatches ? "MISMATCH" : "ok");
        if (mismatches) {
            std::printf("    %d lanes disagree with the scalar predicate\n", mismatches);
            ++failures;
        }
        if (touching == 0) {
            std::printf("    WARNING: no touching boxes in this case, so the boundary "
                        "was not exercised\n");
            ++failures;
        }
        return mismatches;
    }

    // vaabb_disjoint_one_to_many writes a mask array rather than returning bits.
    template <typename T>
    int check_mask(const char* name, const ptrdiff_t n, const int grid) {
        std::mt19937 rng(99);
        Boxes<T> a = make_boxes<T>(rng, n, grid);
        Boxes<T> b = make_boxes<T>(rng, n, grid);

        // It writes exactly SCCD_AABB_DISJOINT_CHUNK_SIZE lanes and takes no count, so
        // the second list is walked a whole chunk at a time; n is chosen to be a
        // multiple of the chunk.
        std::vector<std::uint32_t> mask(SCCD_AABB_DISJOINT_CHUNK_SIZE);
        int mismatches = 0;

        for (ptrdiff_t i = 0; i < a.n; ++i) {
          for (ptrdiff_t start = 0; start + SCCD_AABB_DISJOINT_CHUNK_SIZE <= b.n;
               start += SCCD_AABB_DISJOINT_CHUNK_SIZE) {
            sccd::vaabb_disjoint_one_to_many<T>(
                a.lo[0][i], a.lo[1][i], a.lo[2][i], a.hi[0][i], a.hi[1][i], a.hi[2][i],
                b.lo[0].data() + start, b.lo[1].data() + start, b.lo[2].data() + start,
                b.hi[0].data() + start, b.hi[1].data() + start, b.hi[2].data() + start,
                mask.data());

            for (ptrdiff_t lane = 0; lane < SCCD_AABB_DISJOINT_CHUNK_SIZE; ++lane) {
                const ptrdiff_t j = start + lane;
                const bool vector_disjoint = mask[lane] != 0;
                const bool scalar_disjoint = sccd::disjoint<T>(
                    a.lo[0][i], a.lo[1][i], a.lo[2][i], a.hi[0][i], a.hi[1][i], a.hi[2][i],
                    b.lo[0][j], b.lo[1][j], b.lo[2][j],
                    b.hi[0][j], b.hi[1][j], b.hi[2][j]) != 0;
                if (vector_disjoint != scalar_disjoint) {
                    if (mismatches < 5) {
                        std::printf("    mask disagrees: box %ld vs %ld\n", (long)i, (long)j);
                    }
                    ++mismatches;
                }
            }
          }
        }

        std::printf("  %-30s n=%-6ld %s\n", name, (long)n, mismatches ? "MISMATCH" : "ok");
        if (mismatches) ++failures;
        return mismatches;
    }

    // Boxes that meet exactly, built by hand so the boundary is certain to be hit
    // on every axis rather than left to the generator.
    template <typename T>
    int check_touching(const char* name) {
        const T alo[3] = {0, 0, 0};
        const T ahi[3] = {1, 1, 1};

        struct Case { T blo[3], bhi[3]; bool overlap; const char* what; };
        const Case cases[] = {
            {{1, 0, 0}, {2, 1, 1}, true,  "face-to-face on x"},
            {{0, 1, 0}, {1, 2, 1}, true,  "face-to-face on y"},
            {{0, 0, 1}, {1, 1, 2}, true,  "face-to-face on z"},
            {{1, 1, 1}, {2, 2, 2}, true,  "corner-to-corner"},
            {{-1, 0, 0}, {0, 1, 1}, true, "face-to-face on -x"},
            {{1, 1, 1}, {1, 1, 1}, true,  "zero-extent box on the corner"},
            {{2, 0, 0}, {3, 1, 1}, false, "clear of it on x"},
        };

        int bad = 0;
        for (const Case& c : cases) {
            const bool scalar_says = sccd::disjoint<T>(alo[0], alo[1], alo[2], ahi[0], ahi[1], ahi[2],
                                                       c.blo[0], c.blo[1], c.blo[2],
                                                       c.bhi[0], c.bhi[1], c.bhi[2]) == 0;
            const T* bl[3] = {&c.blo[0], &c.blo[1], &c.blo[2]};
            const T* bh[3] = {&c.bhi[0], &c.bhi[1], &c.bhi[2]};
            const std::uint32_t bits = sccd::vaabb_overlap_one_to_many_bits<T>(
                alo[0], alo[1], alo[2], ahi[0], ahi[1], ahi[2],
                bl[0], bl[1], bl[2], bh[0], bh[1], bh[2], 1);
            const bool vector_says = (bits & 1u) != 0;

            const bool ok = (scalar_says == c.overlap) && (vector_says == c.overlap);
            std::printf("  %-30s %-30s expected %-8s scalar %-8s vector %-8s %s\n",
                        name, c.what,
                        c.overlap ? "overlap" : "disjoint",
                        scalar_says ? "overlap" : "disjoint",
                        vector_says ? "overlap" : "disjoint",
                        ok ? "ok" : "MISMATCH");
            if (!ok) { ++failures; ++bad; }
        }
        return bad;
    }

}  // namespace

int main() {
    std::printf("touching boxes must count as overlapping:\n");
    check_touching<double>("double");
    check_touching<float>("float");

    std::printf("vaabb_overlap_one_to_many_bits vs the scalar predicate:\n");
    check<double>("double, dense grid", 400, 6);
    check<double>("double, sparse grid", 400, 40);
    check<float>("float, dense grid", 400, 6);
    check<float>("float, sparse grid", 400, 40);

    std::printf("vaabb_disjoint_one_to_many vs the scalar predicate:\n");
    check_mask<double>("double", 320, 6);
    check_mask<float>("float", 320, 6);

    std::printf("%s\n", failures ? "FAIL"
                                 : "OK: every SIMD kernel agrees with the scalar predicate");
    return failures ? 1 : 0;
}
