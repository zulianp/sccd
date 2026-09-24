#ifndef SCCD_BROADPHASE_SWEEP_HPP
#define SCCD_BROADPHASE_SWEEP_HPP

#include <algorithm>
#include <cfloat>
#include <memory>
#include <cstdint>
#include <cstring>
#include <vector>

#include "sccd_parallel.hpp"
#include "sccd_aabb.hpp"

namespace sccd {

    /**
     * \brief Choose the axis (0=x,1=y,2=z) with largest variance of AABB centers.
     * \param n Number of AABBs.
     * \param aabb SoA arrays of size 6: minx,miny,minz,maxx,maxy,maxz; each of
     * length n. \return Axis index in {0,1,2}.
     */
    /**
     * \brief Accumulate per-axis center means and variances in two sweeps.
     *
     * Both sweeps visit all three axes per element, so the AABB arrays are
     * streamed twice instead of six times, and both are deterministic tiled
     * reductions: the tiling is a function of \p n alone, so the variances --
     * and therefore the axis chosen from them -- do not move with the worker
     * count. They are not bit-identical to a single serial accumulation, which
     * only matters when two axes tie to the last bit, and either is a valid
     * choice then.
     */
    template <typename T>
    static void center_variance(const ptrdiff_t n, T **const SCCD_RESTRICT aabb, T (&var)[3]) {
        for (int d = 0; d < 3; d++) {
            var[d] = 0;
        }
        if (n <= 0) {
            return;
        }

        // Three per-axis sums, the accumulator both passes reduce into. Local
        // because nothing outside this function has any use for it.
        struct Sum3 {
            T v[3];
        };

        const auto join = [](const Sum3 a, const Sum3 b) {
            Sum3 r;
            for (int d = 0; d < 3; d++) r.v[d] = a.v[d] + b.v[d];
            return r;
        };

        const Sum3 total = sccd::parallel_tiled_reduce<Sum3>(
            0,
            n,
            [&](const ptrdiff_t lo, const ptrdiff_t hi) {
                Sum3 acc = {{0, 0, 0}};
                for (ptrdiff_t i = lo; i < hi; i++) {
                    for (int d = 0; d < 3; d++) {
                        acc.v[d] += (aabb[d + 3][i] + aabb[d][i]) / 2;
                    }
                }
                return acc;
            },
            join);

        T mean[3];
        for (int d = 0; d < 3; d++) {
            mean[d] = total.v[d] / (T)n;
        }

        const Sum3 sq = sccd::parallel_tiled_reduce<Sum3>(
            0,
            n,
            [&](const ptrdiff_t lo, const ptrdiff_t hi) {
                Sum3 acc = {{0, 0, 0}};
                for (ptrdiff_t i = lo; i < hi; i++) {
                    for (int d = 0; d < 3; d++) {
                        const T c = (aabb[d + 3][i] + aabb[d][i]) / 2;
                        acc.v[d] += (c - mean[d]) * (c - mean[d]);
                    }
                }
                return acc;
            },
            join);

        for (int d = 0; d < 3; d++) {
            var[d] = sq.v[d];
        }
    }

    template <typename T>
    static int choose_axis(const ptrdiff_t n, T **const SCCD_RESTRICT aabb) {
        T var[3];
        center_variance(n, aabb, var);

        int fargmax = 0;
        T fmax = var[0];

        for (int d = 1; d < 3; d++) {
            if (fmax < var[d]) {
                fmax = var[d];
                fargmax = d;
            }
        }

        return fargmax;
    }

    template <typename T>
    static void largest_variance_axes_sort(const ptrdiff_t n, T **const SCCD_RESTRICT aabb, int *axes) {
        T var[3];
        center_variance(n, aabb, var);

        T vars[3] = {var[0], var[1], var[2]};
        // printf("vars: %g, %g, %g\n", vars[0], vars[1], vars[2]);

        axes[0] = 0;
        axes[1] = 1;
        axes[2] = 2;
        std::sort(axes, axes + 3, [&](const int a, const int b) { return vars[a] > vars[b]; });

        // printf("axes: %d, %d, %d\n", axes[0], axes[1], axes[2]);
    }

    namespace detail {

        /**
         * \brief Map a float to an unsigned integer of the same width, order preserving.
         *
         * Sorting the images by magnitude sorts the originals by value, which is
         * what lets a radix pass replace a comparison. For a non-negative float
         * the IEEE bit pattern already orders correctly once the sign bit is
         * set; for a negative one every bit is flipped, which reverses the order
         * within the negatives and places them below the positives.
         *
         * A NaN key has no position in any ordering and none is defined here.
         * Swept boxes are built from mesh coordinates by min and max, so a NaN
         * key means a NaN coordinate reached the broad phase.
         */
        static inline uint32_t radix_key(const float v) {
            uint32_t b;
            std::memcpy(&b, &v, sizeof(b));
            return (b & 0x80000000u) ? ~b : (b | 0x80000000u);
        }

        static inline uint64_t radix_key(const double v) {
            uint64_t b;
            std::memcpy(&b, &v, sizeof(b));
            return (b & 0x8000000000000000ull) ? ~b : (b | 0x8000000000000000ull);
        }

        /**
         * \brief Least-significant-digit radix sort of (key, index) pairs.
         *
         * One pass per byte of the key, each pass a privatised histogram over
         * $256$ buckets, a serial offset pass over the histogram, and a scatter.
         * The histogram is over buckets and not over keys, so the private copies
         * cost workers times $256$ entries however long the input is, which is
         * what a counting sort over the keys themselves could not afford.
         *
         * The sort is stable, and the offsets are laid out bucket-major and
         * chunk-minor, so equal keys keep the order they arrived in. The pairs
         * arrive in index order, so stability supplies the tie-break on the
         * index that the comparison form spells out, and the result is the same
         * total order.
         *
         * \param tmp Scratch of \p n pairs. The number of passes is even for
         *        every key width used here, so the sorted result ends in \p a.
         */
        template <typename Pair, typename KeyOf>
        static void radix_sort_pairs(const ptrdiff_t n,
                                     Pair* const SCCD_RESTRICT a,
                                     Pair* const SCCD_RESTRICT tmp,
                                     KeyOf key_of) {
            using U = decltype(key_of(a[0]));
            static const int RADIX = 256;
            static const int PASSES = (int)sizeof(U);

            const int nchunks = sccd::max<int>(1, sccd::min<int>(sccd::max_concurrency() * 2, (int)(n / 4096 + 1)));
            const ptrdiff_t chunk = (n + nchunks - 1) / nchunks;
            std::vector<ptrdiff_t> hist((size_t)nchunks * RADIX);

            Pair* src = a;
            Pair* dst = tmp;

            for (int pass = 0; pass < PASSES; ++pass) {
                const int shift = pass * 8;
                std::fill(hist.begin(), hist.end(), (ptrdiff_t)0);

                sccd::parallel_for_chunks(0, nchunks, [&](const ptrdiff_t c) {
                    ptrdiff_t* const row = hist.data() + c * RADIX;
                    const ptrdiff_t begin = c * chunk;
                    const ptrdiff_t end = sccd::min<ptrdiff_t>(begin + chunk, n);
                    for (ptrdiff_t i = begin; i < end; ++i) {
                        row[(key_of(src[i]) >> shift) & (RADIX - 1)] += 1;
                    }
                });

                // Bucket-major, chunk-minor: a chunk's slot inside a bucket
                // follows the chunks before it, which is what keeps the sort
                // stable.
                ptrdiff_t running = 0;
                for (int b = 0; b < RADIX; ++b) {
                    for (int c = 0; c < nchunks; ++c) {
                        ptrdiff_t& slot = hist[(size_t)c * RADIX + b];
                        const ptrdiff_t take = slot;
                        slot = running;
                        running += take;
                    }
                }

                sccd::parallel_for_chunks(0, nchunks, [&](const ptrdiff_t c) {
                    ptrdiff_t* const row = hist.data() + c * RADIX;
                    const ptrdiff_t begin = c * chunk;
                    const ptrdiff_t end = sccd::min<ptrdiff_t>(begin + chunk, n);
                    for (ptrdiff_t i = begin; i < end; ++i) {
                        dst[row[(key_of(src[i]) >> shift) & (RADIX - 1)]++] = src[i];
                    }
                });

                Pair* const swap = src;
                src = dst;
                dst = swap;
            }
        }

    }  // namespace detail

    /**
     * \brief Reusable working storage for sort_along_axis.
     *
     * The two arrays are sized by the list being sorted, so they depend on the
     * mesh and not on the step. A caller that sorts the same mesh every step
     * should hold one of these and pass it in; otherwise every step allocates
     * and frees 2n pairs once per list, which is per-mesh work charged to a
     * per-step measurement.
     *
     * Passing nothing keeps the allocating behaviour, which is what a one-shot
     * caller wants.
     */
    template <typename T, typename I>
    struct SortScratch {
        struct KeyIndex {
            T key;
            I idx;
        };

        std::unique_ptr<KeyIndex[]> keys;
        std::unique_ptr<KeyIndex[]> tmp;
        ptrdiff_t capacity = 0;

        /** \brief Grow to hold \p n pairs. Never shrinks, so a steady state is free. */
        void reserve(const ptrdiff_t n) {
            if (n <= capacity) return;
            // new[] on a trivial type leaves the storage uninitialised, where a
            // std::vector would value-initialise every entry in a serial pass
            // and then have it overwritten.
            keys.reset(new KeyIndex[(size_t)n]);
            tmp.reset(new KeyIndex[(size_t)n]);
            capacity = n;
        }
    };

    /**
     * \brief Sort AABBs along \p sort_axis and permute all six SoA arrays
     * coherently. \param n Number of AABBs. \param sort_axis Axis to sort by
     * (0=x,1=y,2=z). \param arrays SoA arrays [6][n]:
     * minx,miny,minz,maxx,maxy,maxz. \param idx Output permutation array of size n
     * (initialized to 0..n-1 then sorted). \param scratch Scratch buffer of length
     * n used to permute each component.
     *
     * Arrays must be valid and sufficiently sized. Uses parallel sort.
     */
    template <typename T, typename I>
    static void sort_along_axis(const ptrdiff_t n,
                                const int sort_axis,
                                T **const SCCD_RESTRICT arrays,
                                I *const SCCD_RESTRICT idx,
                                T *const SCCD_RESTRICT scratch,
                                SortScratch<T, I> *const sort_scratch = nullptr) {
        // Sort (key, index) pairs rather than indices with an indirect
        // comparator: every comparison in the old form was a random load out of
        // arrays[sort_axis], which is the dominant cost of the sort.
        using Scratch = SortScratch<T, I>;
        using KeyIndex = typename Scratch::KeyIndex;

        Scratch owned;
        Scratch &ks = sort_scratch ? *sort_scratch : owned;
        ks.reserve(n);
        KeyIndex *const keys = ks.keys.get();
        const T *const SCCD_RESTRICT x = arrays[sort_axis];

        sccd::parallel_for_br(0, n, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t i = rbegin; i < rend; i++) {
                keys[i].key = x[i];
                keys[i].idx = I(i);
            }
        });

        // A radix sort, because this is the dominant cost of the sweep's
        // preparation and a comparison sort does not scale through it: the
        // merge sort behind parallel_sort stops improving at about 2.4x on this
        // input, where the passes below are linear and independent.
        // A fixed size, not a per-worker one: this chooses between two sorting
        // algorithms rather than deciding how wide to go, and the radix sort
        // beats the comparison sort on a single worker too.
        static const ptrdiff_t RADIX_CUTOFF = 4096;
        if (n >= RADIX_CUTOFF) {
            detail::radix_sort_pairs(
                n, keys, ks.tmp.get(), [](const KeyIndex &k) { return detail::radix_key(k.key); });
        } else {
            std::sort(keys, keys + n, [](const KeyIndex &l, const KeyIndex &r) {
                if (l.key < r.key) return true;
                if (r.key < l.key) return false;
                return l.idx < r.idx;
            });
        }

        sccd::parallel_for_br(0, n, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t i = rbegin; i < rend; i++) {
                idx[i] = keys[i].idx;
            }
        });

        for (int d = 0; d < 6; d++) {
            memcpy(scratch, arrays[d], sizeof(T) * n);
            T *const SCCD_RESTRICT dst = arrays[d];
            sccd::parallel_for_br(0, n, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                for (ptrdiff_t i = rbegin; i < rend; i++) {
                    dst[i] = scratch[idx[i]];
                }
            });
        }
    }


    /**
     * \namespace detail
     * \brief Internal helpers for vectorized AABB disjoint tests and overlap
     * filtering.
     */
    namespace detail {

        /**
         * \brief Load the \p nxe vertex indices of element \p elem_idx.
         * \tparam nxe Number of indices per element.
         * \param elements SoA arrays for element indices.
         * \param elem_idx Logical element index.
         * \param element_stride Stride between consecutive elements in the arrays.
         * \param out Output array of size nxe.
         */
        template <int nxe, typename I>
        static inline void load_ev(I **const SCCD_RESTRICT elements,
                                   const I elem_idx,
                                   const ptrdiff_t element_stride,
                                   I (&out)[nxe]) {
            for (int v = 0; v < nxe; ++v) {
                out[v] = elements[v][elem_idx * element_stride];
            }
        }

        /**
         * \brief Test if two index tuples share any vertex id.
         * \tparam n1 Length of first tuple.
         * \tparam n2 Length of second tuple.
         * \param a First tuple of indices.
         * \param b Second tuple of indices.
         * \return True if any index matches.
         */
        template <int n1, int n2, typename I>
        static inline bool shares_vertex(const I (&a)[n1], const I (&b)[n2]) {
            for (int i = 0; i < n1; ++i) {
                for (int j = 0; j < n2; ++j) {
                    if (a[i] == b[j]) {
                        return true;
                    }
                }
            }
            return false;
        }

        /**
         * \brief Advance begin and compute end of candidate window along a sorted
         * axis.
         * \param fimin Min coordinate of first AABB on the sort axis.
         * \param fimax Max coordinate of first AABB on the sort axis.
         * \param second_xmax Sorted xmax array of the second list along the sort
         * axis.
         * \param second_xmin Sorted xmin array of the second list along the sort
         * axis.
         * \param second_count Number of AABBs in the second list.
         * \param begin In/out: advanced to first index where second_xmax[begin] >
         * fimin. \param end Out: first index where second_xmin[end] > fimax.
         */
        template <typename T>
        static inline void compute_candidate_window_progressive(const T fimin,
                                                                const T fimax,
                                                                const T *const SCCD_RESTRICT second_xmax,
                                                                const T *const SCCD_RESTRICT second_xmin,
                                                                const ptrdiff_t second_count,
                                                                ptrdiff_t &begin,
                                                                ptrdiff_t &end) {
            for (; begin < second_count; ++begin) {
                // Inclusive: the overlap predicate counts touching boxes as
                // overlapping (it rejects only on a strict `amin > bmax`), so the
                // window that feeds it has to keep a box whose xmax equals fimin.
                // With a strict `<` here the two disagreed, and the sweep dropped
                // real pairs -- a false negative, which the conservativeness
                // invariant does not allow. Boxes are sorted by xmin, so this can
                // only bite when xmin == xmax == fimin, i.e. a zero-extent box
                // sitting exactly at another box's lower bound on the sort axis.
                // That is what an axis-aligned face produces, and it cost 20 of
                // 2220 edge-edge pairs on a refined cube.
                if (fimin <= second_xmax[begin]) {
                    break;
                }
            }
            end = begin;
            for (; end < second_count; ++end) {
                if (fimax < second_xmin[end]) {
                    break;
                }
            }
        }


        /**
         * \brief Packed overlap mask of AABB \p fi against a block of the second list.
         * \return Bit \p i set when second_aabbs[start+i] is not axis-disjoint from A.
         */
        template <typename T>
        static inline uint32_t overlap_bits_for_block(T **const SCCD_RESTRICT second_aabbs,
                                                      const ptrdiff_t start,
                                                      const int chunk_len,
                                                      const T aminx,
                                                      const T aminy,
                                                      const T aminz,
                                                      const T amaxx,
                                                      const T amaxy,
                                                      const T amaxz) {
            return vaabb_overlap_one_to_many_bits<T>(aminx,
                                                     aminy,
                                                     aminz,
                                                     amaxx,
                                                     amaxy,
                                                     amaxz,
                                                     second_aabbs[0] + start,
                                                     second_aabbs[1] + start,
                                                     second_aabbs[2] + start,
                                                     second_aabbs[3] + start,
                                                     second_aabbs[4] + start,
                                                     second_aabbs[5] + start,
                                                     chunk_len);
        }

        /**
         * \brief Clear the bits of \p bits whose pair shares a vertex (two lists).
         *
         * Only visits lanes that survived the AABB filter, via a bit-scan, so the
         * indirect element loads are never issued for rejected candidates.
         */
        template <int F, int S, typename I>
        static inline uint32_t mask_out_shared_two_lists_bits(uint32_t bits,
                                                              const ptrdiff_t noffset,
                                                              const I (&ev)[F],
                                                              const I *const SCCD_RESTRICT second_idx,
                                                              I **const SCCD_RESTRICT second_elements,
                                                              const ptrdiff_t second_element_stride) {
            uint32_t remaining = bits;
            while (remaining) {
                const int lane = sccd::ctz32(remaining);
                remaining &= remaining - 1;

                const I jidx = second_idx[noffset + lane];
                int match = 0;
                if constexpr (S > 1) {
                    I sev[S];
                    for (int v = 0; v < S; ++v) {
                        sev[v] = second_elements[v][jidx * second_element_stride];
                    }
                    for (int a = 0; a < F; ++a) {
                        for (int b = 0; b < S; ++b) {
                            match |= (ev[a] == sev[b]);
                        }
                    }
                } else {
                    for (int a = 0; a < F; ++a) {
                        match |= (ev[a] == jidx);
                    }
                }
                if (match) {
                    bits &= ~(uint32_t(1) << lane);
                }
            }
            return bits;
        }

        /**
         * \brief Clear the bits of \p bits whose pair shares a vertex (self path).
         */
        template <int N, typename I>
        static inline uint32_t mask_out_shared_self_bits(uint32_t bits,
                                                         const ptrdiff_t noffset,
                                                         const I (&ev)[N],
                                                         const I *const SCCD_RESTRICT idx,
                                                         I **const SCCD_RESTRICT elements,
                                                         const ptrdiff_t element_stride) {
            uint32_t remaining = bits;
            while (remaining) {
                const int lane = sccd::ctz32(remaining);
                remaining &= remaining - 1;

                const I jidx = idx[noffset + lane];
                I sev[N];
                load_ev<N>(elements, jidx, element_stride, sev);
                if (shares_vertex<N, N>(ev, sev)) {
                    bits &= ~(uint32_t(1) << lane);
                }
            }
            return bits;
        }

        // -----------------------------

        /**
         * \brief Scalar reference: count candidate overlaps in [begin,end) for two
         * lists. \return Number of non-disjoint, non-shared-vertex candidates.
         */

        /**
         * \brief Scalar reference: collect candidate overlaps in [begin,end) for
         * two lists. \return Number of pairs written to \p first_out and \p
         * second_out.
         */

        // -----------------------------

    }  // namespace detail

    template <typename T>
    void cummax(const ptrdiff_t n, const T *const SCCD_RESTRICT in, T *const SCCD_RESTRICT out) {
        T acc = in[0];
        out[0] = acc;

        for (ptrdiff_t i = 1; i < n; ++i) {
            acc = sccd::max(acc, in[i]);
            out[i] = acc;
        }
    }

    /**
     * \brief Count candidate overlaps between two sorted AABB lists.
     *
     * Excludes pairs that are axis-disjoint or share a vertex. Both AABB lists
     * must be sorted by \p sort_axis and their index arrays aligned accordingly.
     *
     * \tparam first_nxe Number of vertices per element in the first list.
     * \tparam second_nxe Number of vertices per element in the second list.
     * \param sort_axis Sort axis (0=x,1=y,2=z).
     * \param first_count Number of AABBs in the first list.
     * \param first_aabbs SoA arrays [6][first_count].
     * \param first_idx Mapping from sorted position to element id for the first
     * list. \param first_element_stride Stride between elements in the first element
     * arrays. \param first_elements SoA element-vertex arrays for the first list.
     * \param second_count Number of AABBs in the second list.
     * \param second_aabbs SoA arrays [6][second_count].
     * \param second_idx Mapping from sorted position to element id for the second
     * list. \param second_element_stride Stride between elements in the second element
     * arrays. \param second_elements SoA element-vertex arrays for the second list.
     * \param ccdptr Prefix sum array of size first_count+1. On return:
     *               ccdptr[i+1]-ccdptr[i] = candidates for first i, and
     *               ccdptr[first_count] = total candidates.
     * \return True if any candidates exist.
     */
    template <int first_nxe, int second_nxe, typename T, typename I>
    bool count_overlaps(const int sort_axis,
                        const ptrdiff_t first_count,
                        T **const SCCD_RESTRICT first_aabbs,
                        I *const SCCD_RESTRICT first_idx,
                        const ptrdiff_t first_element_stride,
                        I **const SCCD_RESTRICT first_elements,
                        const ptrdiff_t second_count,
                        T **const SCCD_RESTRICT second_aabbs,
                        I *const SCCD_RESTRICT second_idx,
                        const ptrdiff_t second_element_stride,
                        I **const SCCD_RESTRICT second_elements,
                        ptrdiff_t *const SCCD_RESTRICT ccdptr,
                        const T *const SCCD_RESTRICT second_xmax_running) {
        const T *const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
        const T *const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
        const T *const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
        const T *const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

        ccdptr[0] = 0;

        if (first_xmax[first_count - 1] < second_xmin[0]) return false;
        if (second_xmax[second_count - 1] < first_xmin[0]) return false;

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            ptrdiff_t ni = std::lower_bound(second_xmax_running, second_xmax_running + second_count, first_xmin[rbegin]) - second_xmax_running;

            for (ptrdiff_t fi = rbegin; fi < rend; fi++) {
                const T fimin = first_xmin[fi];
                const T fimax = first_xmax[fi];
                const I first_idxi = first_idx[fi];

                I ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_element_stride];
                }

                ptrdiff_t end = ni;
                detail::compute_candidate_window_progressive(
                    fimin, fimax, second_xmax, second_xmin, second_count, ni, end);

                if (ni >= end) {
                    ccdptr[fi + 1] = 0;
                    continue;
                }

                ptrdiff_t count = 0;

                for (ptrdiff_t noffset = ni; noffset < end;) {
                    const int chunk_len = (int)sccd::min((ptrdiff_t)SCCD_AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    uint32_t bits = detail::overlap_bits_for_block(second_aabbs,
                                                                        noffset,
                                                                        chunk_len,
                                                                        first_aabbs[0][fi],
                                                                        first_aabbs[1][fi],
                                                                        first_aabbs[2][fi],
                                                                        first_aabbs[3][fi],
                                                                        first_aabbs[4][fi],
                                                                        first_aabbs[5][fi]);

                    bits = detail::mask_out_shared_two_lists_bits<first_nxe, second_nxe>(
                        bits, noffset, ev, second_idx, second_elements, second_element_stride);

                    count += sccd::popcount32(bits);
                    noffset += chunk_len;
                }

                ccdptr[fi + 1] = count;
            }
        });

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + first_count + 1);
        return ccdptr[first_count] > 0;
    }

    /**
     * \brief Collect candidate overlaps between two sorted AABB lists.
     *
     * Uses the prefix offsets computed by count_overlaps to write pairs.
     *
     * \param sort_axis Sort axis (0=x,1=y,2=z).
     * \param first_count Number of AABBs in the first list.
     * \param first_aabbs SoA arrays [6][first_count].
     * \param first_idx Mapping from sorted position to element id for the first
     * list. \param first_element_stride Stride between elements in the first element
     * arrays. \param first_elements SoA element-vertex arrays for the first list.
     * \param second_count Number of AABBs in the second list.
     * \param second_aabbs SoA arrays [6][second_count].
     * \param second_idx Mapping from sorted position to element id for the second
     * list. \param second_element_stride Stride between elements in the second element
     * arrays. \param second_elements SoA element-vertex arrays for the second list.
     * \param ccdptr Prefix offsets from the count pass (size first_count+1).
     * \param first_out Output array (size ccdptr[first_count]) for first indices.
     * \param second_out Output array (size ccdptr[first_count]) for second indices.
     */
    template <int first_nxe, int second_nxe, typename T, typename I>
    void collect_overlaps(const int sort_axis,
                          const ptrdiff_t first_count,
                          T **const SCCD_RESTRICT first_aabbs,
                          I *const SCCD_RESTRICT first_idx,
                          const ptrdiff_t first_element_stride,
                          I **SCCD_RESTRICT const first_elements,
                          const ptrdiff_t second_count,
                          T **const SCCD_RESTRICT second_aabbs,
                          I *const SCCD_RESTRICT second_idx,
                          const ptrdiff_t second_element_stride,
                          I **SCCD_RESTRICT const second_elements,
                          const ptrdiff_t *const SCCD_RESTRICT ccdptr,
                          const T *const SCCD_RESTRICT second_xmax_running,
                          I *SCCD_RESTRICT first_out,
                          I *SCCD_RESTRICT second_out) {
        const T *const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
        const T *const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
        const T *const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
        const T *const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

        if (first_xmax[first_count - 1] < second_xmin[0]) return;
        if (second_xmax[second_count - 1] < first_xmin[0]) return;

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            ptrdiff_t ni = std::lower_bound(second_xmax_running, second_xmax_running + second_count, first_xmin[rbegin]) - second_xmax_running;

            for (ptrdiff_t fi = rbegin; fi < rend; fi++) {
                const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
                if (expected_count == 0) {
                    continue;
                }

                const T fimin = first_xmin[fi];
                const T fimax = first_xmax[fi];
                const I first_idxi = first_idx[fi];

                I *SCCD_RESTRICT const first_local_elements = &first_out[ccdptr[fi]];
                I *SCCD_RESTRICT const second_local_elements = &second_out[ccdptr[fi]];

                I ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_element_stride];
                }

                ptrdiff_t end = ni;
                detail::compute_candidate_window_progressive(
                    fimin, fimax, second_xmax, second_xmin, second_count, ni, end);

                if (ni >= end) {
                    continue;
                }

                ptrdiff_t count = 0;

                for (ptrdiff_t noffset = ni; noffset < end;) {
                    const int chunk_len = (int)sccd::min((ptrdiff_t)SCCD_AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    uint32_t bits = detail::overlap_bits_for_block(second_aabbs,
                                                                        noffset,
                                                                        chunk_len,
                                                                        first_aabbs[0][fi],
                                                                        first_aabbs[1][fi],
                                                                        first_aabbs[2][fi],
                                                                        first_aabbs[3][fi],
                                                                        first_aabbs[4][fi],
                                                                        first_aabbs[5][fi]);

                    bits = detail::mask_out_shared_two_lists_bits<first_nxe, second_nxe>(
                        bits, noffset, ev, second_idx, second_elements, second_element_stride);

                    while (bits) {
                        const int lane = sccd::ctz32(bits);
                        bits &= bits - 1;
                        first_local_elements[count] = first_idxi;
                        second_local_elements[count] = second_idx[noffset + lane];
                        count += 1;
                    }

                    noffset += chunk_len;
                }

                assert(expected_count == count);
            }
        });
    }

    // --------------------------------------

    /**
     * \brief Count candidate self-overlaps (upper triangle) within one sorted AABB
     * list.
     *
     * Excludes pairs that are axis-disjoint or share a vertex. The AABBs must be
     * sorted by \p sort_axis and \p idx maps sorted positions to element ids.
     *
     * \tparam nxe Number of vertices per element.
     * \param sort_axis Sort axis (0=x,1=y,2=z).
     * \param element_count Number of AABBs/elements.
     * \param aabbs SoA arrays [6][element_count].
     * \param idx Mapping from sorted position to element id.
     * \param element_stride Stride between elements in the vertex arrays.
     * \param elements SoA element-vertex arrays.
     * \param ccdptr Prefix sum array size element_count+1; filled as in the
     * two-lists case. \return True if any candidates exist.
     */
    template <int nxe, typename T, typename I>
    bool count_self_overlaps(const int sort_axis,
                             const ptrdiff_t element_count,
                             T **const SCCD_RESTRICT aabbs,
                             I *const SCCD_RESTRICT idx,
                             const ptrdiff_t element_stride,
                             I **const SCCD_RESTRICT elements,
                             ptrdiff_t *const SCCD_RESTRICT ccdptr) {
        const T *const SCCD_RESTRICT xmin = aabbs[sort_axis];
        const T *const SCCD_RESTRICT xmax = aabbs[3 + sort_axis];

        ccdptr[0] = 0;

        sccd::parallel_for_br(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; fi++) {
                const T fimin = xmin[fi];
                const T fimax = xmax[fi];
                const I idxi = idx[fi];

                assert(idxi >= 0);
                assert(idxi < element_count);

                I ev[nxe];
                for (int v = 0; v < nxe; v++) {
                    ev[v] = elements[v][idxi * element_stride];
                }

                ptrdiff_t noffset = fi + 1;
                ptrdiff_t end = noffset;
                detail::compute_candidate_window_progressive(
                    fimin, fimax, xmax, xmin, element_count, noffset, end);

                if (noffset >= end) {
                    ccdptr[fi + 1] = 0;
                    continue;
                }

                ptrdiff_t count = 0;

                for (; noffset < end;) {
                    const int chunk_len = (int)sccd::min((ptrdiff_t)SCCD_AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    uint32_t bits = detail::overlap_bits_for_block(aabbs,
                                                                        noffset,
                                                                        chunk_len,
                                                                        aabbs[0][fi],
                                                                        aabbs[1][fi],
                                                                        aabbs[2][fi],
                                                                        aabbs[3][fi],
                                                                        aabbs[4][fi],
                                                                        aabbs[5][fi]);

                    bits = detail::mask_out_shared_self_bits<nxe>(bits, noffset, ev, idx, elements, element_stride);

                    count += sccd::popcount32(bits);
                    noffset += chunk_len;
                }

                ccdptr[fi + 1] = count;
            }
        });

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + element_count + 1);
        return ccdptr[element_count] > 0;
    }

    /**
     * \brief Collect candidate self-overlap pairs using prefix offsets.
     *
     * Writes pairs (min(idxi,jidx), max(idxi,jidx)) for j>i into output arrays.
     *
     * \param sort_axis Sort axis (0=x,1=y,2=z).
     * \param element_count Number of elements/AABBs.
     * \param aabbs SoA arrays [6][element_count].
     * \param idx Mapping from sorted position to element id.
     * \param element_stride Stride between elements in the vertex arrays.
     * \param elements SoA element-vertex arrays.
     * \param ccdptr Prefix offsets from the self count pass (size element_count+1).
     * \param first_out Output array of first indices.
     * \param second_out Output array of second indices.
     */
    template <int nxe, typename T, typename I>
    void collect_self_overlaps(const int sort_axis,
                               const ptrdiff_t element_count,
                               T **const SCCD_RESTRICT aabbs,
                               I *const SCCD_RESTRICT idx,
                               const ptrdiff_t element_stride,
                               I **const elements,
                               const ptrdiff_t *const SCCD_RESTRICT ccdptr,
                               I *SCCD_RESTRICT first_out,
                               I *SCCD_RESTRICT second_out) {
        const T *const SCCD_RESTRICT xmin = aabbs[sort_axis];
        const T *const SCCD_RESTRICT xmax = aabbs[3 + sort_axis];

        sccd::parallel_for_br(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; fi++) {
                const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
                if (!expected_count) continue;

                const T fimin = xmin[fi];
                const T fimax = xmax[fi];
                const I idxi = idx[fi];

                I *SCCD_RESTRICT const first_local_elements = &first_out[ccdptr[fi]];

                I *SCCD_RESTRICT const second_local_elements = &second_out[ccdptr[fi]];

                I ev[nxe];
                for (int v = 0; v < nxe; v++) {
                    ev[v] = elements[v][idxi * element_stride];
                }

                ptrdiff_t noffset = fi + 1;
                ptrdiff_t end = noffset;
                detail::compute_candidate_window_progressive(
                    fimin, fimax, xmax, xmin, element_count, noffset, end);

                if (noffset >= end) {
                    continue;
                }

                ptrdiff_t count = 0;
                for (; noffset < end;) {
                    const int chunk_len = (int)sccd::min((ptrdiff_t)SCCD_AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    uint32_t bits = detail::overlap_bits_for_block(aabbs,
                                                                        noffset,
                                                                        chunk_len,
                                                                        aabbs[0][fi],
                                                                        aabbs[1][fi],
                                                                        aabbs[2][fi],
                                                                        aabbs[3][fi],
                                                                        aabbs[4][fi],
                                                                        aabbs[5][fi]);

                    bits = detail::mask_out_shared_self_bits<nxe>(bits, noffset, ev, idx, elements, element_stride);

                    while (bits) {
                        const int lane = sccd::ctz32(bits);
                        bits &= bits - 1;
                        const I jidx = idx[noffset + lane];
                        first_local_elements[count] = sccd::min(idxi, jidx);
                        second_local_elements[count] = sccd::max(idxi, jidx);
                        count += 1;
                    }

                    noffset += chunk_len;
                }

                assert(expected_count == count);
            }
        });
    }

}  // namespace sccd
#endif
