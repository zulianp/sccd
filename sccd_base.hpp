#ifndef SCCD_BASE_HPP
#define SCCD_BASE_HPP

#ifndef SFEM_READ_ENV
#define SFEM_READ_ENV(name, conversion)                                        \
    do {                                                                       \
        char* var = getenv(#name);                                             \
        if (var) {                                                             \
            name = conversion(var);                                            \
        }                                                                      \
    } while (0)
#endif

#ifndef SCCD_RESTRICT
#ifndef _WIN32
#define SCCD_RESTRICT __restrict__
#else
#define SCCD_RESTRICT __restrict
#endif
#endif

/// Geometry scalar type used for coordinates and AABB values.
using geom_t = float;
/// Integer type used for element/vertex indices.
using idx_t = int;
/// Integer type used for counts and prefix sums.
using count_t = int;


#endif
