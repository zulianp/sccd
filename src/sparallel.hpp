#ifndef SCCD_PARALLEL_HPP
#define SCCD_PARALLEL_HPP

#define SCCD_ENABLE_TBB

#ifdef SCCD_ENABLE_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_scan.h>
#include <tbb/parallel_sort.h>
#endif

#include <cstddef>

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
            ptrdiff_t iend = min(TILE_SIZE, end - i);
            fun(i, iend);
        }
#endif
    }

    template <typename T, typename F>
    void parallel_sort(T* const begin, T* const end, F fun) {
#ifdef SCCD_ENABLE_TBB
        tbb::parallel_sort(begin, end, fun);
#else
        // FIXME: Implement parallel sort without TBB
        std::sort(begin, end, fun);
#endif
    }

    template <typename T>
    void parallel_cum_sum_br(T* const begin, T* const end) {
        ptrdiff_t len = end - begin;

#ifdef SCCD_ENABLE_TBB
        tbb::parallel_scan(
            tbb::blocked_range<ptrdiff_t>(0, len),
            0,
            [=](const tbb::blocked_range<ptrdiff_t>& r, T sum, bool is_final_scan) -> T {
                if (!is_final_scan) {
                    T temp = sum;
                    for (int i = r.begin(); i < r.end(); ++i) {
                        temp = temp + begin[i];
                    }
                    return temp;
                } else {
                    begin[r.begin()] += sum;
                    for (int i = r.begin() + 1; i < r.end(); ++i) {
                        begin[i] += begin[i - 1];
                    }

                    return begin[r.end() - 1];
                }
            },
            [](T left, T right) { return left + right; });
#else
        T acc = 0;
#pragma omp parallel for reduction(inscan, + : acc)
        for (ptrdiff_t i = 0; i < len; i++) {
            acc += begin[i];

#pragma omp scan inclusive(acc)
            begin[i] = acc;
        }
#endif  // SCCD_ENABLE_TBB
    }

    template <typename T>
    void parallel_cum_max_br(T* const begin, T* const end) {
        ptrdiff_t len = end - begin;

#ifdef SCCD_ENABLE_TBB
        tbb::parallel_scan(
            tbb::blocked_range<ptrdiff_t>(0, len),
            begin[0],
            [=](const tbb::blocked_range<ptrdiff_t>& r, T acc, bool is_final_scan) -> T {
                if (!is_final_scan) {
                    T temp = acc;
                    for (int i = r.begin(); i < r.end(); ++i) {
                        temp = sccd::max(temp, begin[i]);
                    }
                    return temp;
                } else {
                    begin[r.begin()] = sccd::max(begin[r.begin()], acc);
                    for (int i = r.begin() + 1; i < r.end(); ++i) {
                        begin[i] = sccd::max(begin[i], begin[i - 1]);
                    }

                    return begin[r.end() - 1];
                }
            },
            [](T left, T right) { return sccd::max(left, right); });
#else
        T acc = begin[0];
#pragma omp parallel for reduction(inscan, max : acc)
        for (ptrdiff_t i = 0; i < len; i++) {
            acc = sccd::max(acc, begin[i]);

#pragma omp scan inclusive(acc)
            begin[i] = acc;
        }
#endif  // SCCD_ENABLE_TBB
    }
}  // namespace sccd

#endif  // SCCD_PARALLEL_HPP
