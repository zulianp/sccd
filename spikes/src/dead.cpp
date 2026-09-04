// Keeps spikes/src/dead.hpp compiling, so it cannot quietly rot into something
// that no longer parses. Nothing calls these -- that is what makes them dead --
// so the only thing asserted here is that they still build.
//
// The templates are `static inline`, which cannot be explicitly instantiated, so
// this odr-uses one from each family instead. Add a line here when you add a
// function there.

#include "dead.hpp"

namespace sccd {
    namespace dead {

        void typecheck_dead_code() {
            const double p[3] = {0, 0, 0};
            double g[3] = {0, 0, 0};
            double h[9] = {0};
            double f = 0;

            vf_gradient<double>(p, p, p, p, p, p, p, p, p, p, 0.5, 0.25, 0.25, g);
            vf_hessian<double>(p, p, p, p, p, p, p, p, p, p, 0.5, 0.25, 0.25, h);
            ee_objective<double>(p, p, p, p, p, p, p, p, p, p, 0.5, 0.25, 0.25, &f);
            ee_gradient<double>(p, p, p, p, p, p, p, p, p, p, 0.5, 0.25, 0.25, g);

            const double a[4] = {1, 2, 3, 4};
            (void)array_min<double>(4, a);
            (void)array_max<double>(4, a);
            (void)f;
        }

    }  // namespace dead
}  // namespace sccd
