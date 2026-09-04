// Simple C exports for SCCD narrowphase vertex-face root finder
// #include "sccd_narrowphase.hpp"

#include "sccd_rootfinder.hpp"

#ifdef SCCD_NP_COUNT_BOXES
namespace sccd {
    unsigned long long g_np_host_boxes = 0;
    unsigned long long g_np_host_hist[24] = {0};
    unsigned long long g_np_host_level[80] = {0};
    unsigned long long g_np_host_depth_accept = 0;
}
#endif

extern "C" {

int sccd_find_root_bisection_vf_f(int max_iter,
                                  float tol,
                                  const float sv[3],
                                  const float s1[3],
                                  const float s2[3],
                                  const float s3[3],
                                  const float ev[3],
                                  const float e1[3],
                                  const float e2[3],
                                  const float e3[3],
                                  float* t,
                                  float* u,
                                  float* v) {
    return sccd::find_root_bisection<float>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, *t, *u, *v);
}

int sccd_find_root_bisection_vf_d(int max_iter,
                                  double tol,
                                  const double sv[3],
                                  const double s1[3],
                                  const double s2[3],
                                  const double s3[3],
                                  const double ev[3],
                                  const double e1[3],
                                  const double e2[3],
                                  const double e3[3],
                                  double* t,
                                  double* u,
                                  double* v) {
    return sccd::find_root_bisection<double>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, *t, *u, *v);
}

int sccd_find_root_vf_f(int max_iter,
                        float tol,
                        const float sv[3],
                        const float s1[3],
                        const float s2[3],
                        const float s3[3],
                        const float ev[3],
                        const float e1[3],
                        const float e2[3],
                        const float e3[3],
                        float* t,
                        float* u,
                        float* v) {
    std::vector<sccd::Box<float>> stack;
    stack.reserve(1024);


    float tols[3];
    compute_face_vertex_tolerance<float>(tol, sv, s1, s2, s3, ev, e1, e2, e3, tols);
    float numerical_error[3];
    sccd::numerical_error_bound<true, float>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);
    const float codomain_widths[3] = {1, 1, 1};

    bool found = false;
    stack.push_back(sccd::unit_domain_box<float>());
    while (!stack.empty()) {
        const sccd::Box<float> box = stack.back();
        stack.pop_back();

        found |= sccd::find_root_grid_adaptive_split_vf<float>(
            max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, box, *t, *u, *v, stack);
    }
    return found;
}

int sccd_find_root_vf_d(int max_iter,
                        double tol,
                        const double sv[3],
                        const double s1[3],
                        const double s2[3],
                        const double s3[3],
                        const double ev[3],
                        const double e1[3],
                        const double e2[3],
                        const double e3[3],
                        double* t,
                        double* u,
                        double* v) {

    std::vector<sccd::Box<double>> stack;
    stack.reserve(1024);

    double tols[3];
    compute_face_vertex_tolerance<double>(tol, sv, s1, s2, s3, ev, e1, e2, e3, tols);
    double numerical_error[3];
    sccd::numerical_error_bound<true, double>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);
    const double codomain_widths[3] = {1, 1, 1};

    bool found = false;
    stack.push_back(sccd::unit_domain_box<double>());
    while (!stack.empty()) {
        const sccd::Box<double> box = stack.back();
        stack.pop_back();

        found |= sccd::find_root_grid_adaptive_split_vf<double>(
            max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, box, *t, *u, *v, stack);
    }
    return found;
}

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
int sccd_find_root_tight_inclusion_vf_d(int max_iter,
                                        double tol,
                                        const double sv[3],
                                        const double s1[3],
                                        const double s2[3],
                                        const double s3[3],
                                        const double ev[3],
                                        const double e1[3],
                                        const double e2[3],
                                        const double e3[3],
                                        double* t,
                                        double* u,
                                        double* v) {
    return sccd::find_root_tight_inclusion_vf<double>(max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, *t, *u, *v);
}

int sccd_find_root_tight_inclusion_ee_d(int max_iter,
                                        double tol,
                                        const double s0[3],
                                        const double s1[3],
                                        const double s2[3],
                                        const double s3[3],
                                        const double e0[3],
                                        const double e1[3],
                                        const double e2[3],
                                        const double e3[3],
                                        double* t,
                                        double* u,
                                        double* v) {
    return sccd::find_root_tight_inclusion_ee<double>(max_iter, tol, s0, s1, s2, s3, e0, e1, e2, e3, *t, *u, *v);
}
#endif

int sccd_find_root_ee_d(int max_iter,
                        double tol,
                        const double s0[3],
                        const double s1[3],
                        const double s2[3],
                        const double s3[3],
                        const double e0[3],
                        const double e1[3],
                        const double e2[3],
                        const double e3[3],
                        double* t,
                        double* u,
                        double* v) {
    std::vector<sccd::Box<double>> stack;
    stack.reserve(1024);


    double tols[3];
    compute_edge_edge_tolerance<double>(tol, s0, s1, s2, s3, e0, e1, e2, e3, tols);
    double numerical_error[3];
    sccd::numerical_error_bound<false, double>(s0, s1, s2, s3, e0, e1, e2, e3, numerical_error);
    const double codomain_widths[3] = {1, 1, 1};

    bool found = false;
    stack.push_back(sccd::unit_domain_box<double>());
    while (!stack.empty()) {
        const sccd::Box<double> box = stack.back();
        stack.pop_back();

        found |= sccd::find_root_grid_adaptive_split_ee<double>(
            max_iter, tol, tols, numerical_error, codomain_widths, s0, s1, s2, s3, e0, e1, e2, e3, box, *t, *u, *v, stack);
    }
    return found;
}

}  // extern "C"
