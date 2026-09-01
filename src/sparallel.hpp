#ifndef SCCD_PARALLEL_HPP
#define SCCD_PARALLEL_HPP

// sccd_base.hpp pulls in the generated sccd_config.hpp, which is what defines
// SCCD_ENABLE_TBB/SCCD_ENABLE_OPENMP. It must come before the checks below,
// otherwise the TBB backend silently compiles out.
#include "sccd_base.hpp"

#ifdef SCCD_ENABLE_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_scan.h>
#include <tbb/parallel_sort.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cstddef>
#include <vector>

#include "smath.hpp"

namespace sccd {
    template <typename F>
    void parallel_for_br(const ptrdiff_t start, const ptrdiff_t end, F fun) {
#ifdef SCCD_ENABLE_TBB
        tbb::parallel_for(tbb::blocked_range<ptrdiff_t>(start, end),
                          [&](const tbb::blocked_range<ptrdiff_t>& r) { fun(r.begin(), r.end()); });
#else
        static const ptrdiff_t TILE_SIZE = 128;
#pragma omp parallel for
        for (ptrdiff_t i = start; i < end; i += TILE_SIZE) {
            ptrdiff_t iend = min(i + TILE_SIZE, end);
            fun(i, iend);
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

    namespace sccd_detail {
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
        template <typename T, typename F>
        static void parallel_sort_impl(T* const begin, T* const end, F fun) {
            const ptrdiff_t len = end - begin;
            static const ptrdiff_t SORT_CUTOFF = 1 << 15;

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
    }  // namespace sccd_detail

    template <typename T, typename F>
    void parallel_sort(T* const begin, T* const end, F fun) {
#if defined(SCCD_ENABLE_TBB)
        tbb::parallel_sort(begin, end, fun);
#elif defined(_OPENMP)
#pragma omp parallel
        {
#pragma omp single nowait
            sccd_detail::parallel_sort_impl(begin, end, fun);
        }
#else
        std::sort(begin, end, fun);
#endif
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
        sccd_detail::parallel_inclusive_scan<T>(begin, len, [](const T a, const T b) { return a + b; });
#endif  // SCCD_ENABLE_TBB
    }

    template <typename T>
    void parallel_cum_max_br(T* const begin, T* const end) {
        const ptrdiff_t len = end - begin;
        if (len <= 0) {
            return;
        }

#ifdef SCCD_ENABLE_TBB
        tbb::parallel_scan(
            tbb::blocked_range<ptrdiff_t>(0, len),
            begin[0],
            [=](const tbb::blocked_range<ptrdiff_t>& r, T acc, bool is_final_scan) -> T {
                if (!is_final_scan) {
                    T temp = acc;
                    for (ptrdiff_t i = r.begin(); i < r.end(); ++i) {
                        temp = sccd::max(temp, begin[i]);
                    }
                    return temp;
                } else {
                    begin[r.begin()] = sccd::max(begin[r.begin()], acc);
                    for (ptrdiff_t i = r.begin() + 1; i < r.end(); ++i) {
                        begin[i] = sccd::max(begin[i], begin[i - 1]);
                    }

                    return begin[r.end() - 1];
                }
            },
            [](T left, T right) { return sccd::max(left, right); });
#else
        sccd_detail::parallel_inclusive_scan<T>(begin, len, [](const T a, const T b) { return sccd::max(a, b); });
#endif  // SCCD_ENABLE_TBB
    }
}  // namespace sccd

#endif  // SCCD_PARALLEL_HPP
