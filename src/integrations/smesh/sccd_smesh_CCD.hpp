#ifndef SCCD_SMESH_CCD_HPP
#define SCCD_SMESH_CCD_HPP

#include "broadphase.hpp"
#include "narrowphase.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_tracer.hpp"

#if defined(SCCD_ENABLE_CUDA)
#include "sccd_broadphase.cuh"
#include "sccd_narrowphase.cuh"
#include "sccd_vaabb.cuh"

#include <cuda_runtime_api.h>
#endif

namespace sccd {
    template <typename scalar_t>
    class CCD {
    public:
        CCD(const std::shared_ptr<smesh::Mesh>& mesh,
            const smesh::ExecutionSpace execution_space = smesh::EXECUTION_SPACE_HOST)
            : mesh_(mesh), execution_space_(execution_space) {}

        static std::shared_ptr<CCD> create(const std::shared_ptr<smesh::Mesh>& mesh,
                                           const smesh::ExecutionSpace execution_space = smesh::EXECUTION_SPACE_HOST) {
#if !defined(SCCD_ENABLE_CUDA)
            SMESH_ASSERT(execution_space == smesh::EXECUTION_SPACE_HOST);
            auto ccd = std::make_shared<CCD>(mesh, smesh::EXECUTION_SPACE_HOST);
#else
            auto ccd = std::make_shared<CCD>(mesh, execution_space);
#endif
            return ccd;
        }

        int find_earliest_impact_time(const smesh::SharedBuffer<scalar_t*>& points_t0,
                                      const smesh::SharedBuffer<scalar_t*>& points_t1,
                                      scalar_t& toi) {
            points_t0_ = points_t0;
            points_t1_ = points_t1;

            int err = find_impact_times_impl(0);

            if (err == SCCD_SUCCESS) {
                auto temp = smesh::to_host(ee_tois_);
                toi = temp->data()[0];
            }

            return err;
        }

        int find_impact_times(const smesh::SharedBuffer<scalar_t*>& points_t0,
                              const smesh::SharedBuffer<scalar_t*>& points_t1,
                              smesh::SharedBuffer<smesh::idx_t>& v_overlap,
                              smesh::SharedBuffer<smesh::idx_t>& f_overlap,
                              smesh::SharedBuffer<scalar_t>& vf_tois,
                              smesh::SharedBuffer<smesh::idx_t>& e0_overlap,
                              smesh::SharedBuffer<smesh::idx_t>& e1_overlap,
                              smesh::SharedBuffer<scalar_t>& ee_tois) {
            points_t0_ = points_t0;
            points_t1_ = points_t1;

            int err = find_impact_times_impl(1);

            v_overlap = v_overlap_;
            f_overlap = f_overlap_;
            vf_tois = vf_tois_;
            e0_overlap = e0_overlap_;
            e1_overlap = e1_overlap_;
            ee_tois = ee_tois_;

            return err;
        }

    private:
        std::shared_ptr<smesh::Mesh> mesh_;
        smesh::ExecutionSpace execution_space_{smesh::EXECUTION_SPACE_HOST};

        smesh::SharedBuffer<smesh::idx_t*> faces_;
        smesh::SharedBuffer<smesh::idx_t> e0_;
        smesh::SharedBuffer<smesh::idx_t> e1_;
        smesh::SharedBuffer<smesh::idx_t*> edges_;

        smesh::SharedBuffer<scalar_t*> points_t0_;
        smesh::SharedBuffer<scalar_t*> points_t1_;

        smesh::SharedBuffer<scalar_t*> vaabb_;
        smesh::SharedBuffer<scalar_t*> faabb_;
        smesh::SharedBuffer<scalar_t*> eaabb_;

        smesh::SharedBuffer<smesh::idx_t> vidx_;
        smesh::SharedBuffer<smesh::idx_t> fidx_;
        smesh::SharedBuffer<smesh::idx_t> eidx_;

        smesh::SharedBuffer<scalar_t> scratch_;
        smesh::SharedBuffer<ptrdiff_t> ccdptr_;

        smesh::SharedBuffer<smesh::idx_t> f_overlap_;
        smesh::SharedBuffer<smesh::idx_t> v_overlap_;
        smesh::SharedBuffer<smesh::idx_t> e0_overlap_;
        smesh::SharedBuffer<smesh::idx_t> e1_overlap_;

        smesh::SharedBuffer<scalar_t> vf_tois_;
        smesh::SharedBuffer<scalar_t> ee_tois_;

        smesh::SharedBuffer<scalar_t> cumulative_max_;

        int find_impact_times_impl_broadphase_host(const ptrdiff_t toi_stride) {
            SMESH_UNUSED(toi_stride);

            if (!vaabb_) {
                init();
            }

            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("Broadphase");

                {
                    SMESH_TRACE_SCOPE("Broadphase: AABB");
                    sccd::compute_aabbs(
                        dim, n_nodes, points_t0_->data(), points_t1_->data(), vaabb_->data(), vaabb_->data() + dim);

                    sccd::compute_aabbs(mesh_->block(0)->n_nodes_per_element(),
                                        n_faces,
                                        faces_->data(),
                                        dim,
                                        points_t0_->data(),
                                        points_t1_->data(),
                                        faabb_->data(),
                                        faabb_->data() + dim);

                    sccd::compute_aabbs(2,
                                        n_edges,
                                        edges_->data(),
                                        dim,
                                        points_t0_->data(),
                                        points_t1_->data(),
                                        eaabb_->data(),
                                        eaabb_->data() + dim);
                }

                int sort_axis = 0;
                {
                    SMESH_TRACE_SCOPE("Sorting AABBs (host)");

                    sort_axis = sccd::choose_axis(n_nodes, vaabb_->data());

                    sccd::sort_along_axis(n_nodes, sort_axis, vaabb_->data(), vidx_->data(), scratch_->data());
                    sccd::sort_along_axis(n_faces, sort_axis, faabb_->data(), fidx_->data(), scratch_->data());
                    sccd::sort_along_axis(n_edges, sort_axis, eaabb_->data(), eidx_->data(), scratch_->data());
                }

                {
                    SMESH_TRACE_SCOPE("Broadphase: F2V");

                    {
                        SMESH_TRACE_SCOPE("cummax");
                        sccd::cummax(n_nodes, vaabb_->data()[dim + sort_axis], cumulative_max_->data());
                    }

                    {
                        SMESH_TRACE_SCOPE("count_overlaps");
                        sccd::count_overlaps<3, 1, scalar_t, smesh::idx_t>(sort_axis,
                                                                           n_faces,
                                                                           faabb_->data(),
                                                                           fidx_->data(),
                                                                           1,
                                                                           faces_->data(),
                                                                           n_nodes,
                                                                           vaabb_->data(),
                                                                           vidx_->data(),
                                                                           0,
                                                                           nullptr,
                                                                           ccdptr_->data(),
                                                                           cumulative_max_->data());
                    }

                    const ptrdiff_t n_vertex_face_overlaps = ccdptr_->data()[n_faces];

                    {
                        SMESH_TRACE_SCOPE("f2v allocations");
                        f_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
                        v_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
                    }

                    {
                        SMESH_TRACE_SCOPE("collect_overlaps");
                        sccd::collect_overlaps<3, 1, scalar_t, smesh::idx_t>(sort_axis,
                                                                             n_faces,
                                                                             faabb_->data(),
                                                                             fidx_->data(),
                                                                             1,
                                                                             faces_->data(),
                                                                             n_nodes,
                                                                             vaabb_->data(),
                                                                             vidx_->data(),
                                                                             0,
                                                                             nullptr,
                                                                             ccdptr_->data(),
                                                                             cumulative_max_->data(),
                                                                             f_overlap_->data(),
                                                                             v_overlap_->data());
                    }
                }

                {
                    SMESH_TRACE_SCOPE("Broadphase: E2E");

                    {
                        SMESH_TRACE_SCOPE("count_self_overlaps");
                        sccd::count_self_overlaps<2>(
                            sort_axis, n_edges, eaabb_->data(), eidx_->data(), 1, edges_->data(), ccdptr_->data());
                    }

                    const ptrdiff_t n_edge_overlaps = ccdptr_->data()[n_edges];

                    {
                        SMESH_TRACE_SCOPE("e2e allocations");
                        e0_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
                        e1_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
                    }

                    {
                        SMESH_TRACE_SCOPE("collect_self_overlaps");
                        sccd::collect_self_overlaps<2>(sort_axis,
                                                       n_edges,
                                                       eaabb_->data(),
                                                       eidx_->data(),
                                                       1,
                                                       edges_->data(),
                                                       ccdptr_->data(),
                                                       e0_overlap_->data(),
                                                       e1_overlap_->data());
                    }
                }
            }

            return SCCD_SUCCESS;
        }

        int find_impact_times_impl_narrowphase_host(const ptrdiff_t toi_stride) {
            const auto toi_storage_size = [toi_stride](const size_t noverlaps) {
                return toi_stride == 0 ? size_t(1) : noverlaps;
            };

            {
                SMESH_TRACE_SCOPE("Narrow phase");

                scalar_t toi = 1;
                {
                    SMESH_TRACE_SCOPE("Narrow phase: F2V");

                    vf_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(v_overlap_->size()), execution_space_);
                    sccd::narrow_phase_vf<3, scalar_t, smesh::idx_t>(v_overlap_->size(),
                                                                     v_overlap_->data(),
                                                                     f_overlap_->data(),
                                                                     points_t0_->data(),
                                                                     points_t1_->data(),
                                                                     1,
                                                                     faces_->data(),
                                                                     toi,
                                                                     vf_tois_->data(),
                                                                     toi_stride);

                    if (toi_stride == 0 && v_overlap_->size() != 0) {
                        toi = vf_tois_->data()[0];
                    }
                }

                {
                    SMESH_TRACE_SCOPE("Narrow phase: E2E");

                    ee_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(e0_overlap_->size()), execution_space_);
                    sccd::narrow_phase_ee<scalar_t, smesh::idx_t>(e0_overlap_->size(),
                                                                  e0_overlap_->data(),
                                                                  e1_overlap_->data(),
                                                                  points_t0_->data(),
                                                                  points_t1_->data(),
                                                                  1,
                                                                  edges_->data(),
                                                                  toi,
                                                                  ee_tois_->data(),
                                                                  toi_stride);
                }
            }

            return SCCD_SUCCESS;
        }

        int find_impact_times_impl_host(const ptrdiff_t toi_stride) {
            int err = find_impact_times_impl_broadphase_host(toi_stride);
            err |= find_impact_times_impl_narrowphase_host(toi_stride);
            return err;
        }

        int find_impact_times_impl_broadphase_device(const ptrdiff_t toi_stride) {
#if defined(SCCD_ENABLE_CUDA)
            if (!vaabb_) {
                init();
            }

            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("Broadphase");

                {
                    SMESH_TRACE_SCOPE("Broadphase: AABB");
                    sccd::device::compute_aabbs(
                        dim, mesh_->n_nodes(), points_t0_->data(), points_t1_->data(), vaabb_->data());

                    sccd::device::compute_aabbs(mesh_->block(0)->n_nodes_per_element(),
                                                n_faces,
                                                faces_->data(),
                                                dim,
                                                points_t0_->data(),
                                                points_t1_->data(),
                                                faabb_->data());

                    sccd::device::compute_aabbs(
                        2, n_edges, edges_->data(), dim, points_t0_->data(), points_t1_->data(), eaabb_->data());
                }

                int sort_axis = 0;
                {
                    SMESH_TRACE_SCOPE("Sorting AABBs (device)");

                    sort_axis = sccd::device::choose_axis(dim, n_nodes, vaabb_->data());

                    sccd::device::sort_along_axis(
                        dim, n_nodes, sort_axis, vaabb_->data(), vidx_->data(), scratch_->data());

                    sccd::device::sort_along_axis(
                        dim, n_faces, sort_axis, faabb_->data(), fidx_->data(), scratch_->data());

                    sccd::device::sort_along_axis(
                        dim, n_edges, sort_axis, eaabb_->data(), eidx_->data(), scratch_->data());
                }

                {
                    SMESH_TRACE_SCOPE("Broadphase: F2V");

                    scalar_t* vaabb_max_axis = nullptr;
                    cudaMemcpy(&vaabb_max_axis,
                               vaabb_->data() + dim + sort_axis,
                               sizeof(vaabb_max_axis),
                               cudaMemcpyDeviceToHost);

                    {
                        SMESH_TRACE_SCOPE("cummax");

                        sccd::device::cummax(n_nodes, vaabb_max_axis, cumulative_max_->data());
                    }

                    {
                        SMESH_TRACE_SCOPE("count_overlaps");

                        sccd::device::count_overlaps<3, 1>(sort_axis,
                                                           n_faces,
                                                           faabb_->data(),
                                                           fidx_->data(),
                                                           1,
                                                           faces_->data(),
                                                           n_nodes,
                                                           vaabb_->data(),
                                                           vidx_->data(),
                                                           0,
                                                           (smesh::idx_t**)nullptr,
                                                           ccdptr_->data(),
                                                           cumulative_max_->data());
                    }

                    ptrdiff_t n_vertex_face_overlaps = 0;
                    cudaMemcpy(&n_vertex_face_overlaps,
                               ccdptr_->data() + n_faces,
                               sizeof(n_vertex_face_overlaps),
                               cudaMemcpyDeviceToHost);

                    {
                        SMESH_TRACE_SCOPE("f2v allocations");

                        f_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
                        v_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
                    }

                    {
                        SMESH_TRACE_SCOPE("collect_overlaps");

                        sccd::device::collect_overlaps<3, 1>(sort_axis,
                                                             n_faces,
                                                             faabb_->data(),
                                                             fidx_->data(),
                                                             1,
                                                             faces_->data(),
                                                             n_nodes,
                                                             vaabb_->data(),
                                                             vidx_->data(),
                                                             0,
                                                             (smesh::idx_t**)nullptr,
                                                             ccdptr_->data(),
                                                             cumulative_max_->data(),
                                                             f_overlap_->data(),
                                                             v_overlap_->data());
                    }
                }

                {
                    SMESH_TRACE_SCOPE("Broadphase: E2E");

                    {
                        SMESH_TRACE_SCOPE("count_self_overlaps");
                        sccd::device::count_self_overlaps<2>(
                            sort_axis, n_edges, eaabb_->data(), eidx_->data(), 1, edges_->data(), ccdptr_->data());
                    }

                    ptrdiff_t n_edge_overlaps = 0;
                    cudaMemcpy(
                        &n_edge_overlaps, ccdptr_->data() + n_edges, sizeof(n_edge_overlaps), cudaMemcpyDeviceToHost);

                    {
                        SMESH_TRACE_SCOPE("e2e allocations");

                        e0_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
                        e1_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
                    }

                    {
                        SMESH_TRACE_SCOPE("collect_self_overlaps");

                        sccd::device::collect_self_overlaps<2>(sort_axis,
                                                               n_edges,
                                                               eaabb_->data(),
                                                               eidx_->data(),
                                                               1,
                                                               edges_->data(),
                                                               ccdptr_->data(),
                                                               e0_overlap_->data(),
                                                               e1_overlap_->data());
                    }
                }
            }

            return SCCD_SUCCESS;
#else
            SMESH_UNUSED(toi_stride);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int find_impact_times_impl_narrowphase_device(const ptrdiff_t toi_stride) {
#if defined(SCCD_ENABLE_CUDA)
            const auto toi_storage_size = [toi_stride](const size_t noverlaps) {
                return toi_stride == 0 ? size_t(1) : noverlaps;
            };

            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("Narrow phase");

                scalar_t toi = 1;
                {
                    SMESH_TRACE_SCOPE("Narrow phase: F2V");

                    vf_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(v_overlap_->size()), execution_space_);
                    sccd::device::narrow_phase_vf<3>(v_overlap_->size(),
                                                     v_overlap_->data(),
                                                     f_overlap_->data(),
                                                     points_t0_->data(),
                                                     points_t1_->data(),
                                                     1,
                                                     faces_->data(),
                                                     toi,
                                                     vf_tois_->data(),
                                                     toi_stride);

                    if (toi_stride == 0 && v_overlap_->size() != 0) {
                        auto temp = smesh::to_host(vf_tois_);
                        toi = temp->data()[0];
                    }
                }

                {
                    SMESH_TRACE_SCOPE("Narrow phase: E2E");

                    ee_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(e0_overlap_->size()), execution_space_);
                    sccd::device::narrow_phase_ee(e0_overlap_->size(),
                                                  e0_overlap_->data(),
                                                  e1_overlap_->data(),
                                                  points_t0_->data(),
                                                  points_t1_->data(),
                                                  1,
                                                  edges_->data(),
                                                  toi,
                                                  ee_tois_->data(),
                                                  toi_stride);
                }
            }
#else
            SMESH_ERROR("Not implemented");
#endif
            return SCCD_SUCCESS;
        }

        int find_impact_times_impl_device(const ptrdiff_t toi_stride) {
            int err = find_impact_times_impl_broadphase_device(toi_stride);
            err |= find_impact_times_impl_narrowphase_device(toi_stride);
            return err;
        }

        int find_impact_times_impl(const ptrdiff_t toi_stride) {
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                return find_impact_times_impl_host(toi_stride);
            } else {
                return find_impact_times_impl_device(toi_stride);
            }
        }

        void init() {
            SMESH_TRACE_SCOPE("CCD::init");

            auto n2n_crs = mesh_->edge_graph();
            auto row_idx_temp = smesh::create_host_buffer<smesh::idx_t>(n2n_crs->nnz());
            smesh::crs_to_coo(mesh_->n_nodes(), n2n_crs->rowptr()->data(), row_idx_temp->data());

            const int dim = mesh_->spatial_dimension();
            SMESH_ASSERT(dim == 3);

            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const ptrdiff_t n_edges = n2n_crs->nnz();

            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                e0_ = row_idx_temp;
                e1_ = n2n_crs->colidx();
                faces_ = mesh_->block(0)->elements();
            } else {
                e0_ = smesh::to_device(row_idx_temp);
                e1_ = smesh::to_device(n2n_crs->colidx());
                faces_ = mesh_->block(0)->device_elements_SoA();
            }

            edges_ = smesh::create_2d(std::vector{e0_, e1_});

            vaabb_ = smesh::create_buffer<scalar_t>(2 * dim, n_nodes, execution_space_);
            faabb_ = smesh::create_buffer<scalar_t>(2 * dim, n_faces, execution_space_);
            eaabb_ = smesh::create_buffer<scalar_t>(2 * dim, n_edges, execution_space_);
            cumulative_max_ = smesh::create_buffer<scalar_t>(n_nodes, execution_space_);

            vidx_ = smesh::create_buffer<smesh::idx_t>(n_nodes, execution_space_);
            fidx_ = smesh::create_buffer<smesh::idx_t>(n_faces, execution_space_);
            eidx_ = smesh::create_buffer<smesh::idx_t>(n_edges, execution_space_);
            scratch_ = smesh::create_buffer<scalar_t>(std::max(n_nodes, std::max(n_faces, n_edges)), execution_space_);
            ccdptr_ = smesh::create_buffer<ptrdiff_t>(std::max(n_faces, n_edges) + 1, execution_space_);
        }
    };
}  // namespace sccd

#endif  // SCCD_SMESH_CCD_HPP
