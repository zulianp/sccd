#include "broadphase.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_test.hpp"

#include <algorithm>
#include <utility>
#include <vector>

template <typename T, ptrdiff_t N>
static void bind_aabb(T (&storage)[6][N], T* (&ptrs)[6]) {
    for (int d = 0; d < 6; ++d) {
        ptrs[d] = storage[d];
    }
}

template <typename I, int NXE, ptrdiff_t N>
static void bind_elements(I (&storage)[NXE][N], I* (&ptrs)[NXE]) {
    for (int d = 0; d < NXE; ++d) {
        ptrs[d] = storage[d];
    }
}

template <typename I>
static std::vector<std::pair<I, I>> sorted_pairs(const I* first, const I* second, const ptrdiff_t n) {
    std::vector<std::pair<I, I>> pairs;
    pairs.reserve(n);
    for (ptrdiff_t i = 0; i < n; ++i) {
        pairs.emplace_back(first[i], second[i]);
    }
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

template <int FaceNXE>
static int run_face_vertex_case(const std::vector<std::pair<int, int>>& expected_pairs) {
    using T = double;
    using I = int;

    constexpr ptrdiff_t n_faces = 4;
    constexpr ptrdiff_t n_vertices = 5;

    T face_storage[6][n_faces] = {
        {4.0, 0.0, 2.0, 7.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {5.0, 1.0, 3.0, 8.0},
        {1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0},
    };
    T vertex_storage[6][n_vertices] = {
        {0.5, 2.5, 4.5, 7.5, 9.0},
        {0.5, 0.5, 0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5, 0.5, 0.5},
        {0.5, 2.5, 4.5, 7.5, 9.0},
        {0.5, 0.5, 0.5, 0.5, 0.5},
        {0.5, 0.5, 0.5, 0.5, 0.5},
    };

    I face_element_storage[FaceNXE][64] = {};
    for (int v = 0; v < FaceNXE; ++v) {
        for (ptrdiff_t f = 0; f < n_faces; ++f) {
            face_element_storage[v][f] = 100 + static_cast<I>(FaceNXE * f + v);
        }
    }
    face_element_storage[1][0] = 2;
    face_element_storage[0][1] = 0;

    T* face_aabb[6];
    T* vertex_aabb[6];
    I* face_elements[FaceNXE];
    bind_aabb(face_storage, face_aabb);
    bind_aabb(vertex_storage, vertex_aabb);
    bind_elements(face_element_storage, face_elements);

    I face_idx[n_faces];
    I vertex_idx[n_vertices];
    constexpr ptrdiff_t scratch_size = n_faces > n_vertices ? n_faces : n_vertices;
    T scratch[scratch_size];
    const int sort_axis = sccd::choose_axis(n_vertices, vertex_aabb);
    SMESH_TEST_EQ(sort_axis, 0);

    sccd::sort_along_axis(n_faces, sort_axis, face_aabb, face_idx, scratch);
    sccd::sort_along_axis(n_vertices, sort_axis, vertex_aabb, vertex_idx, scratch);

    T cumulative_max[n_vertices];
    ptrdiff_t ccdptr[n_faces + 1];
    sccd::cummax(n_vertices, vertex_aabb[3 + sort_axis], cumulative_max);

    const bool has_overlaps = sccd::count_overlaps<FaceNXE, 1, T, I>(sort_axis,
                                                                    n_faces,
                                                                    face_aabb,
                                                                    face_idx,
                                                                    1,
                                                                    face_elements,
                                                                    n_vertices,
                                                                    vertex_aabb,
                                                                    vertex_idx,
                                                                    0,
                                                                    nullptr,
                                                                    ccdptr,
                                                                    cumulative_max);
    SMESH_TEST_ASSERT(has_overlaps);
    SMESH_TEST_EQ(ccdptr[n_faces], static_cast<ptrdiff_t>(expected_pairs.size()));

    I face_overlap[8] = {};
    I vertex_overlap[8] = {};
    sccd::collect_overlaps<FaceNXE, 1, T, I>(sort_axis,
                                             n_faces,
                                             face_aabb,
                                             face_idx,
                                             1,
                                             face_elements,
                                             n_vertices,
                                             vertex_aabb,
                                             vertex_idx,
                                             0,
                                             nullptr,
                                             ccdptr,
                                             cumulative_max,
                                             face_overlap,
                                             vertex_overlap);

    const auto actual_pairs = sorted_pairs(face_overlap, vertex_overlap, ccdptr[n_faces]);
    SMESH_TEST_ASSERT(actual_pairs == expected_pairs);
    return SMESH_TEST_SUCCESS;
}

static int test_broadphase_vertex_triangle() {
    return run_face_vertex_case<3>({{2, 1}, {3, 3}});
}

static int test_broadphase_vertex_quad() {
    if (run_face_vertex_case<4>({{2, 1}, {3, 3}}) != SMESH_TEST_SUCCESS) {
        return SMESH_TEST_FAILURE;
    }

    using T = double;
    using I = int;
    constexpr ptrdiff_t n_edges = 5;

    T edge_storage[6][n_edges] = {
        {2.5, 0.0, 2.0, 7.0, 9.0},
        {0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0},
        {4.5, 2.4, 3.0, 8.0, 10.0},
        {1.0, 1.0, 1.0, 1.0, 1.0},
        {1.0, 1.0, 1.0, 1.0, 1.0},
    };
    I edge_element_storage[2][64] = {
        {40, 10, 20, 30, 50},
        {41, 11, 21, 31, 51},
    };
    edge_element_storage[0][2] = 41;

    T* edge_aabb[6];
    I* edge_elements[2];
    bind_aabb(edge_storage, edge_aabb);
    bind_elements(edge_element_storage, edge_elements);

    I edge_idx[n_edges];
    T scratch[n_edges];
    const int axes_expected[3] = {0, 1, 2};
    int axes[3];
    sccd::largest_variance_axes_sort(n_edges, edge_aabb, axes);
    SMESH_ASSERT_ARRAY_EQ(3, axes_expected, axes);

    sccd::sort_along_axis(n_edges, 0, edge_aabb, edge_idx, scratch);

    ptrdiff_t ccdptr[n_edges + 1];
    const bool has_overlaps =
        sccd::count_self_overlaps<2>(0, n_edges, edge_aabb, edge_idx, 1, edge_elements, ccdptr);
    SMESH_TEST_ASSERT(has_overlaps);
    SMESH_TEST_EQ(ccdptr[n_edges], static_cast<ptrdiff_t>(1));

    I edge0_overlap[2] = {};
    I edge1_overlap[2] = {};
    sccd::collect_self_overlaps<2>(
        0, n_edges, edge_aabb, edge_idx, 1, edge_elements, ccdptr, edge0_overlap, edge1_overlap);

    SMESH_TEST_EQ(edge0_overlap[0], 1);
    SMESH_TEST_EQ(edge1_overlap[0], 2);
    return SMESH_TEST_SUCCESS;
}

int main(int argc, char** argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);

    {
        SMESH_RUN_TEST(test_broadphase_vertex_triangle);
        SMESH_RUN_TEST(test_broadphase_vertex_quad);
    }

    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
