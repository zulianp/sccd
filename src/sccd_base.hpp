#ifndef SCCD_BASE_HPP
#define SCCD_BASE_HPP

#include "sccd_config.hpp"

#ifndef SCCD_READ_ENV
#define SCCD_READ_ENV(name, conversion)                                        \
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
#endif
