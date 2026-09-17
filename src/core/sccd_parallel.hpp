#ifndef SCCD_PARALLEL_HPP
#define SCCD_PARALLEL_HPP

// sccd_base.hpp pulls in the generated sccd_config.hpp, which is what defines
// SCCD_ENABLE_TBB/SCCD_ENABLE_OPENMP. It must come before the checks below,
// otherwise the TBB backend silently compiles out.
#include "sccd_base.hpp"

#ifdef SCCD_ENABLE_TBB
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_scan.h>
#include <tbb/parallel_sort.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <vector>

#include "sccd_math.hpp"

namespace sccd {
    /**
     * \brief Shortest range worth opening a parallel region for.
     *
     * Forking and joining costs a fixed amount, measured in tens of
     * microseconds, while a pass over the range costs time proportional to its
     * length. Below the crossover a parallel loop is slower than a serial one,
     * and the broad phase runs several short passes per step on a small mesh,
     * so the guard is worth having.
     *
     * A fixed length rather than a length per worker, which is the form that
     * looks more principled and measures worse. Sizing the team from the work
     * instead -- one worker per few hundred elements -- was four times slower
     * than this on armadillo-rollers at seventy-two threads, because the
     * preparation passes are tens of thousands of elements and the narrowed
     * team left most of the machine idle on every one of them.
     */
    static const ptrdiff_t SCCD_MIN_PARALLEL_N = 4096;

    /**
     * \brief How many workers a parallel_for here can expect to get.
     *
     * Used to size work partitions, not to schedule anything: a partition that
     * ignores the worker count either leaves most of them idle or cuts the work
     * so finely that the bookkeeping costs more than the work.
     */
    static inline int max_concurrency() {
#if defined(SCCD_ENABLE_TBB)
        return (int)tbb::this_task_arena::max_concurrency();
#elif defined(_OPENMP)
        return omp_get_max_threads();
#else
        return 1;
#endif
    }

    template <typename F>
    void parallel_for_br(const ptrdiff_t start, const ptrdiff_t end, F fun) {
#ifdef SCCD_ENABLE_TBB
        tbb::parallel_for(tbb::blocked_range<ptrdiff_t>(start, end),
                          [&](const tbb::blocked_range<ptrdiff_t>& r) { fun(r.begin(), r.end()); });
#else
        static const ptrdiff_t TILE_SIZE = 128;
#pragma omp parallel for if (end - start >= SCCD_MIN_PARALLEL_N)
        for (ptrdiff_t i = start; i < end; i += TILE_SIZE) {
            ptrdiff_t iend = min(i + TILE_SIZE, end);
            fun(i, iend);
        }
#endif
    }

    /**
     * \brief Parallel loop over a short range, one index per worker.
     *
     * parallel_for_br cuts its range into fixed tiles of 128, so a range shorter
     * than that runs entirely on one worker. That is the right default for a
     * loop over elements and the wrong one for a loop over work chunks the
     * caller has already sized, where the range length is the worker count and
     * every index is a large piece of work. Use this for the second kind.
     */
    template <typename F>
    void parallel_for_chunks(const ptrdiff_t start, const ptrdiff_t end, F fun) {
#ifdef SCCD_ENABLE_TBB
        tbb::parallel_for(
            tbb::blocked_range<ptrdiff_t>(start, end, 1),
            [&](const tbb::blocked_range<ptrdiff_t>& r) {
                for (ptrdiff_t i = r.begin(); i < r.end(); ++i) fun(i);
            },
            tbb::simple_partitioner());
#elif defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1)
        for (ptrdiff_t i = start; i < end; ++i) {
            fun(i);
        }
#else
        for (ptrdiff_t i = start; i < end; ++i) {
            fun(i);
        }
#endif
    }

    /**
     * \brief Like parallel_for_br but for loops whose per-item cost varies a lot.
     *
     * Uses fine-grained blocks with a work-stealing/dynamic schedule so that a
     * few very expensive items do not stall a whole worker.
     */
    template <typename F>
    void parallel_for_br_dynamic(const ptrdiff_t start, const ptrdiff_t end, F fun) {
#ifdef SCCD_ENABLE_TBB
        // simple_partitioner + an explicit small grain size keeps the range split
        // down to `grain`, which is what makes stealing effective for skewed work.
        // (The default auto_partitioner stops subdividing far too early here.)
        static const ptrdiff_t GRAIN = 8;
        tbb::parallel_for(
            tbb::blocked_range<ptrdiff_t>(start, end, GRAIN),
            [&](const tbb::blocked_range<ptrdiff_t>& r) { fun(r.begin(), r.end()); },
            tbb::simple_partitioner());
#else
        static const ptrdiff_t TILE_SIZE = 8;
#pragma omp parallel for schedule(dynamic, 8)
        for (ptrdiff_t i = start; i < end; i += TILE_SIZE) {
            ptrdiff_t iend = min(i + TILE_SIZE, end);
            fun(i, iend);
        }
#endif
    }

    namespace detail {
        /**
         * \brief In-place parallel inclusive scan under an arbitrary associative op.
         *
         * Three phases: per-block local scan, serial scan of the block totals,
         * then a parallel fixup of every block but the first. Falls back to the
         * serial scan for short ranges or a single worker.
         */
        template <typename T, typename Op>
        static void parallel_inclusive_scan(T* const SCCD_RESTRICT begin, const ptrdiff_t len, Op op) {
#ifdef _OPENMP
            const ptrdiff_t MIN_PARALLEL_LEN = 4096;
            const int nthreads = omp_get_max_threads();

            if (nthreads > 1 && len >= MIN_PARALLEL_LEN) {
                const ptrdiff_t block = (len + nthreads - 1) / nthreads;
                const int nblocks = static_cast<int>((len + block - 1) / block);
                std::vector<T> block_total(nblocks);

#pragma omp parallel for schedule(static)
                for (int b = 0; b < nblocks; ++b) {
                    const ptrdiff_t s = static_cast<ptrdiff_t>(b) * block;
                    const ptrdiff_t e = sccd::min(s + block, len);
                    T acc = begin[s];
                    for (ptrdiff_t i = s + 1; i < e; ++i) {
                        acc = op(acc, begin[i]);
                        begin[i] = acc;
                    }
                    block_total[b] = acc;
                }

                // Turn the block totals into the offset each block must absorb.
                // Block 0 needs no fixup, so no identity element is required.
                T carry = block_total[0];
                for (int b = 1; b < nblocks; ++b) {
                    const T total = block_total[b];
                    block_total[b] = carry;
                    carry = op(carry, total);
                }

#pragma omp parallel for schedule(static)
                for (int b = 1; b < nblocks; ++b) {
                    const ptrdiff_t s = static_cast<ptrdiff_t>(b) * block;
                    const ptrdiff_t e = sccd::min(s + block, len);
                    const T offset = block_total[b];
                    for (ptrdiff_t i = s; i < e; ++i) {
                        begin[i] = op(offset, begin[i]);
                    }
                }
                return;
            }
#endif  // _OPENMP
            T acc = begin[0];
            for (ptrdiff_t i = 1; i < len; ++i) {
                acc = op(acc, begin[i]);
                begin[i] = acc;
            }
        }

#if !defined(SCCD_ENABLE_TBB) && defined(_OPENMP)
        /**
         * \brief Below this length the merge sort is not worth a parallel region.
         *
         * Shared with parallel_sort so that a short range never opens one: the
         * fork and join alone cost more than sorting the range serially, and the
         * broad phase sorts three lists per step, some of them short.
         */
        static const ptrdiff_t SORT_CUTOFF = 1 << 15;

        template <typename T, typename F>
        static void parallel_sort_impl(T* const begin, T* const end, F fun) {
            const ptrdiff_t len = end - begin;

            if (len < SORT_CUTOFF) {
                std::sort(begin, end, fun);
                return;
            }

            T* const mid = begin + len / 2;
#pragma omp task shared(fun)
            parallel_sort_impl(begin, mid, fun);
            parallel_sort_impl(mid, end, fun);
#pragma omp taskwait
            std::inplace_merge(begin, mid, end, fun);
        }
#endif
    }  // namespace detail

    template <typename T, typename F>
    void parallel_sort(T* const begin, T* const end, F fun) {
#if defined(SCCD_ENABLE_TBB)
        tbb::parallel_sort(begin, end, fun);
#elif defined(_OPENMP)
        if (end - begin < detail::SORT_CUTOFF) {
            std::sort(begin, end, fun);
            return;
        }
#pragma omp parallel
        {
#pragma omp single nowait
            detail::parallel_sort_impl(begin, end, fun);
        }
#else
        std::sort(begin, end, fun);
#endif
    }

    /**
     * \brief Deterministic tiled reduction over [start, end).
     *
     * The range is cut into fixed-size tiles, each tile is reduced by \p tile_fun
     * into its own accumulator, and the accumulators are joined in tile order.
     * Two properties follow, and both are needed here. The tiling depends on the
     * range and not on the worker count, so a floating-point sum gives the same
     * answer at one thread and at seventy-two; and the join order is the index
     * order, so it also matches across runs.
     *
     * \param tile_fun Called as tile_fun(lo, hi) and returns that tile's value.
     * \param join_fun Called as join_fun(acc, value), left-associatively.
     */
    template <typename Acc, typename TileFun, typename JoinFun>
    Acc parallel_tiled_reduce(const ptrdiff_t start,
                              const ptrdiff_t end,
                              TileFun tile_fun,
                              JoinFun join_fun) {
        static const ptrdiff_t REDUCE_TILE = 1024;
        // Sixteen tiles is where the partial array and the fork start paying for
        // themselves; below it the whole range is one tile and no region opens.
        static const ptrdiff_t MIN_TILED_LEN = 16 * REDUCE_TILE;

        const ptrdiff_t len = end - start;
        if (len < MIN_TILED_LEN) {
            // One tile: no partials, no parallel region. The decision depends on
            // the range length alone, so the arithmetic a given input receives
            // is still the same at every worker count.
            return tile_fun(start, end);
        }

        const ptrdiff_t ntiles = (len + REDUCE_TILE - 1) / REDUCE_TILE;
        std::vector<Acc> partial((size_t)ntiles);

        parallel_for_chunks(0, ntiles, [&](const ptrdiff_t t) {
            const ptrdiff_t lo = start + t * REDUCE_TILE;
            const ptrdiff_t hi = min(lo + REDUCE_TILE, end);
            partial[(size_t)t] = tile_fun(lo, hi);
        });

        Acc acc = partial[0];
        for (ptrdiff_t t = 1; t < ntiles; ++t) {
            acc = join_fun(acc, partial[(size_t)t]);
        }
        return acc;
    }

    template <typename T>
    void parallel_cum_sum_br(T* const begin, T* const end) {
        const ptrdiff_t len = end - begin;
        if (len <= 0) {
            return;
        }

#ifdef SCCD_ENABLE_TBB
        tbb::parallel_scan(
            tbb::blocked_range<ptrdiff_t>(0, len),
            T{},
            [=](const tbb::blocked_range<ptrdiff_t>& r, T sum, bool is_final_scan) -> T {
                if (!is_final_scan) {
                    T temp = sum;
                    for (ptrdiff_t i = r.begin(); i < r.end(); ++i) {
                        temp = temp + begin[i];
                    }
                    return temp;
                } else {
                    begin[r.begin()] += sum;
                    for (ptrdiff_t i = r.begin() + 1; i < r.end(); ++i) {
                        begin[i] += begin[i - 1];
                    }

                    return begin[r.end() - 1];
                }
            },
            [](T left, T right) { return left + right; });
#else
        detail::parallel_inclusive_scan<T>(begin, len, [](const T a, const T b) { return a + b; });
#endif  // SCCD_ENABLE_TBB
    }

    /**
     * \brief Monotonically lower a shared minimum.
     *
     * Relaxed ordering is sufficient: the value is a pruning bound, not a
     * synchronization edge, and it is only ever decreased. Nothing is published
     * alongside it, so no acquire/release pairing is needed. Under the default
     * seq_cst this sat on the DFS inner loop and serialized every worker on one
     * cache line.
     *
     * \return The value the atomic held when the exchange settled, i.e. the
     *         prior minimum.
     */
    template <typename T>
    T atomic_min(std::atomic<T>& atm, T val) {
        T expected = atm.load(std::memory_order_relaxed);
        while (expected > val) {
            // compare_exchange_weak refreshes 'expected' on failure, so the loop
            // re-tests and exits as soon as another thread has published a
            // smaller value.
            if (atm.compare_exchange_weak(expected, val, std::memory_order_relaxed, std::memory_order_relaxed)) {
                break;
            }
        }
        return expected;
    }

}  // namespace sccd

#endif  // SCCD_PARALLEL_HPP
