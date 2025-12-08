#ifndef SCCD_PARALLEL_HPP
#define SCCD_PARALLEL_HPP

#define SCCD_ENABLE_TBB

#ifdef SCCD_ENABLE_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
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
}  // namespace sccd

#endif  // SCCD_PARALLEL_HPP
