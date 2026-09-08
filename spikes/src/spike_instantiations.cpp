// Forces the template-heavy spikes to be compiled, not merely parsed.
//
// The generated spike_headers.cpp includes every header-only spike, which is
// enough to catch a syntax error or a missing include. It is not enough for a
// template: a call inside an uninstantiated function template is parsed but not
// resolved against the declaration it names, so a spike can go on calling a
// shipped function whose signature has changed and still compile.
//
// That is not hypothetical. `uniform_split.hpp` called `accept_grid_root_vf`
// with thirteen arguments and a `refine` flag for as long as the shipped
// signature had them; when the shipped one became five arguments, the spike kept
// compiling because nothing instantiated the template holding the call.
//
// Each entry point below is instantiated once, for double, which is all it takes
// to have the compiler check the whole chain beneath it.

#include "uniform_split.hpp"

#include <vector>

namespace {

    void instantiate_uniform_split() {
        const double p[3] = {0, 0, 0};
        double t = 1, u = 0, v = 0;
        std::vector<sccd::Box<double>> stack;

        (void)sccd::find_root_grid_uniform_split_vf<double>(
            8, 1e-8, p, p, p, p, p, p, p, p, t, u, v, stack);

        (void)sccd::find_root_grid_uniform_split_ee<double>(
            8, 1e-8, p, p, p, p, p, p, p, p, t, u, v, stack);
    }

}  // namespace

// Referenced so the function is not merely defined and discarded.
void* sccd_spike_instantiation_anchor() {
    return (void*)&instantiate_uniform_split;
}
