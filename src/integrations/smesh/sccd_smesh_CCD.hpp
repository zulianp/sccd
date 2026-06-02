#ifndef SCCD_SMESH_CCD_HPP
#define SCCD_SMESH_CCD_HPP

#include "sccd_config.hpp"

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
                                      scalar_t& toi,
                                      const int max_depth,
                                      const scalar_t tol) {
            points_t0_ = points_t0;
            points_t1_ = points_t1;

            int err = find_impact_times_impl(0, max_depth, tol);

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
                              smesh::SharedBuffer<scalar_t>& ee_tois,
                              const int max_depth,
                              const scalar_t tol) {
            points_t0_ = points_t0;
            points_t1_ = points_t1;

            int err = find_impact_times_impl(1, max_depth, tol);

            v_overlap = v_overlap_;
            f_overlap = f_overlap_;
            vf_tois = vf_tois_;
            e0_overlap = e0_overlap_;
            e1_overlap = e1_overlap_;
            ee_tois = ee_tois_;

            return err;
        }

        // FV + EE broad_phase (prep, then face–vertex overlaps, then edge–edge).
        int broad_phase(const smesh::SharedBuffer<scalar_t*>& points_t0,
                        const smesh::SharedBuffer<scalar_t*>& points_t1,
                        smesh::SharedBuffer<smesh::idx_t>& v_overlap,
                        smesh::SharedBuffer<smesh::idx_t>& f_overlap,
                        smesh::SharedBuffer<smesh::idx_t>& e0_overlap,
                        smesh::SharedBuffer<smesh::idx_t>& e1_overlap) {
            points_t0_ = points_t0;
            points_t1_ = points_t1;

            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                err |= broad_phase_host(0);
            } else {
                err |= broad_phase_device(0);
            }

            v_overlap = v_overlap_;
            f_overlap = f_overlap_;
            e0_overlap = e0_overlap_;
            e1_overlap = e1_overlap_;
            return err;
        }

        // FV + EE narrow_phase; `max_toi` is updated after FV when toi_stride == 0.
        int narrow_phase(scalar_t& max_toi,
                         smesh::SharedBuffer<scalar_t>& vf_tois,
                         smesh::SharedBuffer<scalar_t>& ee_tois,
                         const int max_depth,
                         const scalar_t tol,
                         const ptrdiff_t toi_stride = 0) {
            int err = SCCD_SUCCESS;
            err |= narrow_phase_fv(max_toi, vf_tois, max_depth, tol, toi_stride);
            err |= narrow_phase_ee(max_toi, ee_tois, max_depth, tol, toi_stride);
            return err;
        }

        const smesh::SharedBuffer<smesh::idx_t*>& edges() const { return edges_; }

        void set_safe_inflate(const bool safe_inflate) { safe_inflate_ = safe_inflate; }

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

        // Separating axis chosen from the vertex AABBs and reused across the
        // F2V and E2E sorts/scans within the same broad_phase invocation.
        int sort_axis_{0};
        bool safe_inflate_{false};

    public:
        int broad_phase_fv(const smesh::SharedBuffer<scalar_t*>& points_t0,
                           const smesh::SharedBuffer<scalar_t*>& points_t1,
                           smesh::SharedBuffer<smesh::idx_t>& v_overlap,
                           smesh::SharedBuffer<smesh::idx_t>& f_overlap) {
            int err = SCCD_SUCCESS;
            err |= broad_phase_prep(points_t0, points_t1);
            err |= broad_phase_fv_step(v_overlap, f_overlap);
            return err;
        }

        int broad_phase_ee(const smesh::SharedBuffer<scalar_t*>& points_t0,
                           const smesh::SharedBuffer<scalar_t*>& points_t1,
                           smesh::SharedBuffer<smesh::idx_t>& e0_overlap,
                           smesh::SharedBuffer<smesh::idx_t>& e1_overlap) {
            int err = SCCD_SUCCESS;
            err |= broad_phase_prep(points_t0, points_t1);
            err |= broad_phase_ee_step(e0_overlap, e1_overlap);
            return err;
        }

        int broad_phase_prep(const smesh::SharedBuffer<scalar_t*>& points_t0,
                             const smesh::SharedBuffer<scalar_t*>& points_t1) {
            points_t0_ = points_t0;
            points_t1_ = points_t1;
            if (!vaabb_) {
                init();
            }

            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                SMESH_TRACE_SCOPE("Broad_phase prep");
                err |= broad_phase_prep_host_();
            } else {
                SMESH_TRACE_SCOPE("Broad_phase prep");
                err |= broad_phase_prep_device_();
            }
            return err;
        }

        int broad_phase_fv_step(smesh::SharedBuffer<smesh::idx_t>& v_overlap,
                                smesh::SharedBuffer<smesh::idx_t>& f_overlap) {
            if (!vaabb_) {
                init();
            }

            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                SMESH_TRACE_SCOPE("Broad_phase (FV)");
                err |= broad_phase_fv_step_host_();
            } else {
                SMESH_TRACE_SCOPE("Broad_phase (FV)");
                err |= broad_phase_fv_step_device_();
            }
            v_overlap = v_overlap_;
            f_overlap = f_overlap_;
            return err;
        }

        int broad_phase_ee_step(smesh::SharedBuffer<smesh::idx_t>& e0_overlap,
                                smesh::SharedBuffer<smesh::idx_t>& e1_overlap) {
            if (!vaabb_) {
                init();
            }

            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                SMESH_TRACE_SCOPE("Broad_phase (EE)");
                err |= broad_phase_ee_step_host_();
            } else {
                SMESH_TRACE_SCOPE("Broad_phase (EE)");
                err |= broad_phase_ee_step_device_();
            }

            e0_overlap = e0_overlap_;
            e1_overlap = e1_overlap_;
            return err;
        }

        // Run only the FV narrow_phase using the latest FV broad_phase result.
        // `max_toi` upper-bounds the search; on `toi_stride == 0` (shared scalar
        // output) it is updated in place to the earliest TOI observed.
        int narrow_phase_fv(scalar_t& max_toi,
                            smesh::SharedBuffer<scalar_t>& vf_tois,
                            const int max_depth,
                            const scalar_t tol,
                            const ptrdiff_t toi_stride = 0) {
            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                SMESH_TRACE_SCOPE("Narrow phase (FV)");
                err |= narrow_phase_fv_step_host_(max_toi, max_depth, tol, toi_stride);
            } else {
                SMESH_TRACE_SCOPE("Narrow phase (FV)");
                err |= narrow_phase_fv_step_device_(max_toi, max_depth, tol, toi_stride);
            }
            vf_tois = vf_tois_;
            return err;
        }

        // Run only the EE narrow_phase using the latest EE broad_phase result.
        int narrow_phase_ee(const scalar_t max_toi,
                            smesh::SharedBuffer<scalar_t>& ee_tois,
                            const int max_depth,
                            const scalar_t tol,
                            const ptrdiff_t toi_stride = 0) {
            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                SMESH_TRACE_SCOPE("Narrow phase (EE)");
                err |= narrow_phase_ee_step_host_(max_toi, max_depth, tol, toi_stride);
            } else {
                SMESH_TRACE_SCOPE("Narrow phase (EE)");
                err |= narrow_phase_ee_step_device_(max_toi, max_depth, tol, toi_stride);
            }
            ee_tois = ee_tois_;
            return err;
        }

    private:
        // ---- Host stage helpers ----------------------------------------------

        int broad_phase_prep_host_() {
            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("Broad_phase: AABB");
                sccd::compute_aabbs(dim,
                                    n_nodes,
                                    points_t0_->data(),
                                    points_t1_->data(),
                                    vaabb_->data(),
                                    vaabb_->data() + dim,
                                    safe_inflate_);

                sccd::compute_aabbs(mesh_->block(0)->n_nodes_per_element(),
                                    n_faces,
                                    faces_->data(),
                                    dim,
                                    points_t0_->data(),
                                    points_t1_->data(),
                                    faabb_->data(),
                                    faabb_->data() + dim,
                                    safe_inflate_);

                sccd::compute_aabbs(2,
                                    n_edges,
                                    edges_->data(),
                                    dim,
                                    points_t0_->data(),
                                    points_t1_->data(),
                                    eaabb_->data(),
                                    eaabb_->data() + dim,
                                    safe_inflate_);
            }

            {
                SMESH_TRACE_SCOPE("Sorting AABBs (host)");

                sort_axis_ = sccd::choose_axis(n_nodes, vaabb_->data());

                sccd::sort_along_axis(n_nodes, sort_axis_, vaabb_->data(), vidx_->data(), scratch_->data());
                sccd::sort_along_axis(n_faces, sort_axis_, faabb_->data(), fidx_->data(), scratch_->data());
                sccd::sort_along_axis(n_edges, sort_axis_, eaabb_->data(), eidx_->data(), scratch_->data());
            }

            return SCCD_SUCCESS;
        }

        int broad_phase_fv_step_host_() {
            SMESH_TRACE_SCOPE("Broad_phase: F2V");

            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();

            {
                SMESH_TRACE_SCOPE("cummax");
                sccd::cummax(n_nodes, vaabb_->data()[dim + sort_axis_], cumulative_max_->data());
            }

            {
                SMESH_TRACE_SCOPE("count_overlaps");
                sccd::count_overlaps<3, 1, scalar_t, smesh::idx_t>(sort_axis_,
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
                sccd::collect_overlaps<3, 1, scalar_t, smesh::idx_t>(sort_axis_,
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

            return SCCD_SUCCESS;
        }

        int broad_phase_ee_step_host_() {
            SMESH_TRACE_SCOPE("Broad_phase: E2E");

            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("count_self_overlaps");
                sccd::count_self_overlaps<2>(
                    sort_axis_, n_edges, eaabb_->data(), eidx_->data(), 1, edges_->data(), ccdptr_->data());
            }

            const ptrdiff_t n_edge_overlaps = ccdptr_->data()[n_edges];

            {
                SMESH_TRACE_SCOPE("e2e allocations");
                e0_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
                e1_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
            }

            {
                SMESH_TRACE_SCOPE("collect_self_overlaps");
                sccd::collect_self_overlaps<2>(sort_axis_,
                                               n_edges,
                                               eaabb_->data(),
                                               eidx_->data(),
                                               1,
                                               edges_->data(),
                                               ccdptr_->data(),
                                               e0_overlap_->data(),
                                               e1_overlap_->data());
            }

            return SCCD_SUCCESS;
        }

        int narrow_phase_fv_step_host_(scalar_t& max_toi,
                                       const int max_depth,
                                       const scalar_t tol,
                                       const ptrdiff_t toi_stride) {
            SMESH_TRACE_SCOPE("Narrow phase: F2V");

            const auto toi_storage_size = [toi_stride](const size_t noverlaps) {
                return toi_stride == 0 ? size_t(1) : noverlaps;
            };

            vf_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(v_overlap_->size()), execution_space_);
            sccd::narrow_phase_vf<3, scalar_t, smesh::idx_t>(v_overlap_->size(),
                                                             v_overlap_->data(),
                                                             f_overlap_->data(),
                                                             points_t0_->data(),
                                                             points_t1_->data(),
                                                             1,
                                                             faces_->data(),
                                                             max_toi,
                                                             vf_tois_->data(),
                                                             max_depth,
                                                             tol,
                                                             toi_stride);

            if (toi_stride == 0 && v_overlap_->size() != 0) {
                max_toi = vf_tois_->data()[0];
            }

            return SCCD_SUCCESS;
        }

        int narrow_phase_ee_step_host_(const scalar_t max_toi,
                                       const int max_depth,
                                       const scalar_t tol,
                                       const ptrdiff_t toi_stride) {
            SMESH_TRACE_SCOPE("Narrow phase: E2E");

            const auto toi_storage_size = [toi_stride](const size_t noverlaps) {
                return toi_stride == 0 ? size_t(1) : noverlaps;
            };

            ee_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(e0_overlap_->size()), execution_space_);
            sccd::narrow_phase_ee<scalar_t, smesh::idx_t>(e0_overlap_->size(),
                                                          e0_overlap_->data(),
                                                          e1_overlap_->data(),
                                                          points_t0_->data(),
                                                          points_t1_->data(),
                                                          1,
                                                          edges_->data(),
                                                          max_toi,
                                                          ee_tois_->data(),
                                                          max_depth,
                                                          tol,
                                                          toi_stride);

            return SCCD_SUCCESS;
        }

        // ---- Device stage helpers --------------------------------------------

        int broad_phase_prep_device_() {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("AABBs");
            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("Broad_phase: AABB");
                sccd::device::compute_aabbs(
                    dim, n_nodes, points_t0_->data(), points_t1_->data(), vaabb_->data(), safe_inflate_);

                sccd::device::compute_aabbs(mesh_->block(0)->n_nodes_per_element(),
                                            n_faces,
                                            faces_->data(),
                                            dim,
                                            points_t0_->data(),
                                            points_t1_->data(),
                                            faabb_->data(),
                                            safe_inflate_);

                sccd::device::compute_aabbs(2,
                                            n_edges,
                                            edges_->data(),
                                            dim,
                                            points_t0_->data(),
                                            points_t1_->data(),
                                            eaabb_->data(),
                                            safe_inflate_);
            }

            {
                SMESH_TRACE_SCOPE("Sorting AABBs (device)");

                sort_axis_ = sccd::device::choose_axis(dim, n_nodes, vaabb_->data());

                sccd::device::sort_along_axis(
                    dim, n_nodes, sort_axis_, vaabb_->data(), vidx_->data(), scratch_->data());

                sccd::device::sort_along_axis(
                    dim, n_faces, sort_axis_, faabb_->data(), fidx_->data(), scratch_->data());

                sccd::device::sort_along_axis(
                    dim, n_edges, sort_axis_, eaabb_->data(), eidx_->data(), scratch_->data());
            }
            return SCCD_SUCCESS;
#else
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int broad_phase_fv_step_device_() {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Broad_phase: F2V");

            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();

            scalar_t* vaabb_max_axis = nullptr;
            cudaMemcpy(
                &vaabb_max_axis, vaabb_->data() + dim + sort_axis_, sizeof(vaabb_max_axis), cudaMemcpyDeviceToHost);

            {
                SMESH_TRACE_SCOPE("cummax");
                sccd::device::cummax(n_nodes, vaabb_max_axis, cumulative_max_->data());
            }

            {
                SMESH_TRACE_SCOPE("count_overlaps");
                sccd::device::count_overlaps<3, 1>(sort_axis_,
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
                sccd::device::collect_overlaps<3, 1>(sort_axis_,
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

            return SCCD_SUCCESS;
#else
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int broad_phase_ee_step_device_() {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Broad_phase: E2E");

            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("count_self_overlaps");
                sccd::device::count_self_overlaps<2>(
                    sort_axis_, n_edges, eaabb_->data(), eidx_->data(), 1, edges_->data(), ccdptr_->data());
            }

            ptrdiff_t n_edge_overlaps = 0;
            cudaMemcpy(&n_edge_overlaps, ccdptr_->data() + n_edges, sizeof(n_edge_overlaps), cudaMemcpyDeviceToHost);

            {
                SMESH_TRACE_SCOPE("e2e allocations");
                e0_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
                e1_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
            }

            {
                SMESH_TRACE_SCOPE("collect_self_overlaps");
                sccd::device::collect_self_overlaps<2>(sort_axis_,
                                                       n_edges,
                                                       eaabb_->data(),
                                                       eidx_->data(),
                                                       1,
                                                       edges_->data(),
                                                       ccdptr_->data(),
                                                       e0_overlap_->data(),
                                                       e1_overlap_->data());
            }

            return SCCD_SUCCESS;
#else
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int narrow_phase_fv_step_device_(scalar_t& max_toi,
                                         const int max_depth,
                                         const scalar_t tol,
                                         const ptrdiff_t toi_stride) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Narrow phase: F2V");

            const auto toi_storage_size = [toi_stride](const size_t noverlaps) {
                return toi_stride == 0 ? size_t(1) : noverlaps;
            };

            vf_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(v_overlap_->size()), execution_space_);
            sccd::device::narrow_phase_vf<3>(v_overlap_->size(),
                                             v_overlap_->data(),
                                             f_overlap_->data(),
                                             points_t0_->data(),
                                             points_t1_->data(),
                                             1,
                                             faces_->data(),
                                             max_toi,
                                             vf_tois_->data(),
                                             max_depth,
                                             tol,
                                             toi_stride);

            if (toi_stride == 0 && v_overlap_->size() != 0) {
                auto temp = smesh::to_host(vf_tois_);
                max_toi = temp->data()[0];
            }

            return SCCD_SUCCESS;
#else
            SMESH_UNUSED(max_toi);
            SMESH_UNUSED(max_depth);
            SMESH_UNUSED(tol);
            SMESH_UNUSED(toi_stride);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int narrow_phase_ee_step_device_(const scalar_t max_toi,
                                         const int max_depth,
                                         const scalar_t tol,
                                         const ptrdiff_t toi_stride) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Narrow phase: E2E");

            const auto toi_storage_size = [toi_stride](const size_t noverlaps) {
                return toi_stride == 0 ? size_t(1) : noverlaps;
            };

            ee_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(e0_overlap_->size()), execution_space_);
            sccd::device::narrow_phase_ee(e0_overlap_->size(),
                                          e0_overlap_->data(),
                                          e1_overlap_->data(),
                                          points_t0_->data(),
                                          points_t1_->data(),
                                          1,
                                          edges_->data(),
                                          max_toi,
                                          ee_tois_->data(),
                                          max_depth,
                                          tol,
                                          toi_stride);

            return SCCD_SUCCESS;
#else
            SMESH_UNUSED(max_toi);
            SMESH_UNUSED(max_depth);
            SMESH_UNUSED(tol);
            SMESH_UNUSED(toi_stride);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        // ---- Combined-pipeline glue (matches previous behaviour) -------------

        int broad_phase_host(const ptrdiff_t toi_stride) {
            SMESH_UNUSED(toi_stride);

            if (!vaabb_) {
                init();
            }

            SMESH_TRACE_SCOPE("Broad_phase");
            int err = SCCD_SUCCESS;
            err |= broad_phase_prep_host_();
            err |= broad_phase_fv_step_host_();
            err |= broad_phase_ee_step_host_();
            return err;
        }

        int narrow_phase_host(const ptrdiff_t toi_stride, const int max_depth, const scalar_t tol) {
            SMESH_TRACE_SCOPE("Narrow phase");
            scalar_t toi = 1;
            int err = SCCD_SUCCESS;
            err |= narrow_phase_fv_step_host_(toi, max_depth, tol, toi_stride);
            err |= narrow_phase_ee_step_host_(toi, max_depth, tol, toi_stride);
            return err;
        }

        int find_impact_times_impl_host(const ptrdiff_t toi_stride, const int max_depth, const scalar_t tol) {
            SMESH_TRACE_SCOPE("CCD CPU");
            int err = broad_phase_host(toi_stride);
            err |= narrow_phase_host(toi_stride, max_depth, tol);
            return err;
        }

        int broad_phase_device(const ptrdiff_t toi_stride) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_UNUSED(toi_stride);

            if (!vaabb_) {
                init();
            }

            SMESH_TRACE_SCOPE("Broad_phase");
            int err = SCCD_SUCCESS;
            err |= broad_phase_prep_device_();
            err |= broad_phase_fv_step_device_();
            err |= broad_phase_ee_step_device_();
            return err;
#else
            SMESH_UNUSED(toi_stride);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int narrow_phase_device(const ptrdiff_t toi_stride, const int max_depth, const scalar_t tol) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Narrow phase");
            scalar_t toi = 1;
            int err = SCCD_SUCCESS;
            err |= narrow_phase_fv_step_device_(toi, max_depth, tol, toi_stride);
            err |= narrow_phase_ee_step_device_(toi, max_depth, tol, toi_stride);
            return err;
#else
            SMESH_UNUSED(toi_stride);
            SMESH_UNUSED(max_depth);
            SMESH_UNUSED(tol);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int find_impact_times_impl_device(const ptrdiff_t toi_stride, const int max_depth, const scalar_t tol) {
            SMESH_TRACE_SCOPE("CCD GPU");
            int err = broad_phase_device(toi_stride);
            err |= narrow_phase_device(toi_stride, max_depth, tol);
            return err;
        }

        int find_impact_times_impl(const ptrdiff_t toi_stride, const int max_depth, const scalar_t tol) {
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                return find_impact_times_impl_host(toi_stride, max_depth, tol);
            } else {
                return find_impact_times_impl_device(toi_stride, max_depth, tol);
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
