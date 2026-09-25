#ifndef SCCD_SMESH_CCD_HPP
#define SCCD_SMESH_CCD_HPP

#include "sccd_config.hpp"

#include "sccd_broadphase_sweep.hpp"
#include "sccd_broadphase_strategy.hpp"
#include "sccd_broadphase_cell2d.hpp"
#include "sccd_narrowphase.hpp"
#include "sccd_narrowphase_quad.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_elem_type.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_tracer.hpp"

#if defined(SCCD_ENABLE_CUDA)
#include "sccd_broadphase.cuh"
#include "sccd_cell2d_broadphase.cuh"
#include "sccd_narrowphase.cuh"
#include "sccd_narrowphase_vq.cuh"
#include "sccd_vaabb.cuh"

#include <cuda_runtime_api.h>
#endif

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace sccd {
    template <typename scalar_t>
    class CCD {
    public:
        CCD(const std::shared_ptr<smesh::Mesh>& mesh,
            const smesh::ExecutionSpace execution_space = smesh::EXECUTION_SPACE_HOST)
            : mesh_(mesh),
              execution_space_(execution_space),
              face_element_type_(canonical_face_element_type(mesh)) {}

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

            int err = find_impact_times_impl(sccd::ToiOutput::Earliest, max_depth, tol);

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

            int err = find_impact_times_impl(sccd::ToiOutput::PerPair, max_depth, tol);

            v_overlap = v_overlap_;
            f_overlap = f_overlap_;
            vf_tois = vf_tois_;
            e0_overlap = e0_overlap_;
            e1_overlap = e1_overlap_;
            ee_tois = ee_tois_;

            return err;
        }

        // FV + EE broad_phase (prep, then face–vertex overlaps, then edge–edge).

        // FV + EE narrow_phase; `max_toi` is updated after FV when toi_output == sccd::ToiOutput::Earliest.
        int narrow_phase(scalar_t& max_toi,
                         smesh::SharedBuffer<scalar_t>& vf_tois,
                         smesh::SharedBuffer<scalar_t>& ee_tois,
                         const int max_depth,
                         const scalar_t tol,
                         const sccd::ToiOutput toi_output = sccd::ToiOutput::Earliest) {
            int err = SCCD_SUCCESS;
            err |= narrow_phase_vf(max_toi, vf_tois, max_depth, tol, toi_output);
            err |= narrow_phase_ee(max_toi, ee_tois, max_depth, tol, toi_output);
            return err;
        }

        const smesh::SharedBuffer<smesh::idx_t*>& edges() const { return edges_; }

        void set_box_rounding(const sccd::BoxRounding rounding) { rounding_ = rounding; }

    private:
        static smesh::ElemType canonical_face_element_type(const std::shared_ptr<smesh::Mesh>& mesh) {
            if (!mesh || mesh->n_blocks() != 1) {
                SMESH_ERROR("SCCD requires a non-null, single-block surface mesh\n");
            }

            const smesh::ElemType element_type = mesh->block(0)->element_type();
            if (element_type == smesh::TRI3 || element_type == smesh::TRISHELL3) {
                return smesh::TRISHELL3;
            }
            if (element_type == smesh::QUAD4 || element_type == smesh::QUADSHELL4) {
                return smesh::QUADSHELL4;
            }

            SMESH_ERROR("Unsupported SCCD face element type: %s\n", smesh::type_to_string(element_type));
            return smesh::INVALID;
        }

        std::shared_ptr<smesh::Mesh> mesh_;
        smesh::ExecutionSpace execution_space_{smesh::EXECUTION_SPACE_HOST};
        const smesh::ElemType face_element_type_;

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

        // Sized by the longest of the three lists, so it is reserved once with
        // the rest and never reallocated by a step.
        sccd::SortScratch<scalar_t, smesh::idx_t> sort_scratch_;

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

        // Which broad phase this call is using. The first two are first class:
        // the choice is made per call from the geometry, because neither wins
        // everywhere and the same mesh can change character between frames.
        // SCCD_BROADPHASE forces it to sweep, cell2d or cell2dmin for
        // measurement; cell2dmin is under evaluation and is only ever asked for.
        //
        // cell2dmin shares everything with cell2d but the edge-edge query: the
        // vertex grid and the face-vertex step are the same code on both paths,
        // and only the edge binning and its query differ.
        bool use_cell2d_{false};
        bool use_cell2d_min_{false};
        bool use_cell2d_min_sorted_{false};

        // Broad-phase strategy race. The tuner decides; these carry the pending
        // measurement from the steps back to the next prep, which is where it is
        // reported. See stamp_broad_phase_elapsed_ and broadphase_strategy.hpp.
        sccd::BroadPhaseAutoTuner tuner_;
        sccd::BroadPhaseStrategy timed_strategy_{sccd::BroadPhaseStrategy::Auto};
        double broad_phase_t0_{0.0};
        double broad_phase_pending_ms_{0.0};
        bool broad_phase_pending_{false};

        // Elapsed from the start of prep to the end of whichever step just ran.
        // Called after every step, so the value that survives to the next prep
        // covers the whole broad phase however the caller staged it.
        void stamp_broad_phase_elapsed_() {
            if (broad_phase_t0_ <= 0.0) return;
            broad_phase_pending_ms_ = sccd::broadphase_now_ms() - broad_phase_t0_;
            broad_phase_pending_ = true;
        }
        sccd::Cell2DGrid<scalar_t> v_grid_;
        sccd::Cell2DPartition v_part_;
        std::vector<ptrdiff_t> v_cellptr_;
        std::vector<ptrdiff_t> v_cursor_;
        std::vector<smesh::idx_t> v_cellidx_;

        sccd::Cell2DGrid<scalar_t> e_grid_;
        sccd::Cell2DPartition e_part_;
        std::vector<ptrdiff_t> e_cellptr_;
        std::vector<ptrdiff_t> e_cursor_;
        std::vector<smesh::idx_t> e_cellidx_;

        // The bounds the minimum-corner query prunes with, one pair of scalars
        // per cell. They cost 2 * ncells where the binning they replace saves
        // about 3n cell entries, so the two roughly cancel; the query is what the
        // change is for. Held here rather than allocated per step, like
        // everything else on this path.
        std::vector<scalar_t> e_row_prefix_;
        std::vector<scalar_t> e_cell_hi1_;

        // The ordering on the third axis, when it is asked for: one key per cell
        // entry and one bound per cell.
        std::vector<scalar_t> e_cell_key_;
        std::vector<scalar_t> e_cell_hi2_;

        /** \brief Size the grid, bin one list by minimum corner, and take its bounds. */
        static void bin_min_host_(const ptrdiff_t n,
                                  scalar_t** const SCCD_RESTRICT aabb,
                                  sccd::Cell2DGrid<scalar_t>& grid,
                                  sccd::Cell2DPartition& part,
                                  std::vector<ptrdiff_t>& cellptr,
                                  std::vector<smesh::idx_t>& cellidx,
                                  std::vector<ptrdiff_t>& cursor,
                                  std::vector<scalar_t>& row_prefix,
                                  std::vector<scalar_t>& cell_hi1) {
            sccd::cell2d_setup<scalar_t>(n, aabb, grid);
            sccd::cell2dmin_partition<scalar_t>(n, aabb, grid, part);

            cellptr.resize((size_t)grid.ncells() + 1);
            sccd::cell2dmin_count<scalar_t>(n, aabb, grid, part, cellptr.data());

            cellidx.resize((size_t)cellptr[grid.ncells()]);
            cursor.resize((size_t)grid.ncells());
            sccd::cell2dmin_fill<scalar_t, smesh::idx_t>(
                n, aabb, grid, part, cellptr.data(), cellidx.data(), cursor.data());

            row_prefix.resize((size_t)grid.ncells());
            cell_hi1.resize((size_t)grid.ncells());
            sccd::cell2dmin_bounds<scalar_t>(n, aabb, grid, part, row_prefix.data(), cell_hi1.data());
        }

        /** \brief Order the cells on the third axis and take their bound there. */
        static void sort_cells_host_(scalar_t** const SCCD_RESTRICT aabb,
                                     const sccd::Cell2DGrid<scalar_t>& grid,
                                     const std::vector<ptrdiff_t>& cellptr,
                                     std::vector<smesh::idx_t>& cellidx,
                                     std::vector<scalar_t>& cell_key,
                                     std::vector<scalar_t>& cell_hi2) {
            cell_key.resize(cellidx.size());
            cell_hi2.resize((size_t)grid.ncells());
            sccd::cell2dmin_sort_cells<scalar_t, smesh::idx_t>(
                aabb, grid, cellptr.data(), cellidx.data(), cell_key.data(), cell_hi2.data());
        }

#if defined(SCCD_ENABLE_CUDA)
        /**
         * \brief A plain cudaMalloc array owned by this object.
         *
         * The cell list's row pointer, cursor and range arrays are read and
         * atomically updated by the device, and they have to be real device
         * allocations for that to be reliable: taken from smesh's buffer
         * factory they landed outside the device allocator's range, and the
         * scatter pass then read stale counts and misplaced a few thousand of
         * three hundred thousand writes, differently on every run. Holding them
         * here keeps the guarantee local and explicit.
         */
        template <typename T>
        struct DeviceArray {
            T* ptr{nullptr};
            ptrdiff_t count{0};

            DeviceArray() = default;
            DeviceArray(const DeviceArray&) = delete;
            DeviceArray& operator=(const DeviceArray&) = delete;
            ~DeviceArray() {
                if (ptr) cudaFree(ptr);
            }

            /** \brief At least \p n elements, never shrinking. */
            T* reserve(const ptrdiff_t n) {
                if (n > count) {
                    if (ptr) cudaFree(ptr);
                    ptr = nullptr;
                    if (cudaMalloc((void**)&ptr, sizeof(T) * (size_t)(n > 0 ? n : 1)) != cudaSuccess) {
                        SMESH_ERROR("sccd: out of device memory for %ld cell-list entries\n", (long)n);
                    }
                    count = n > 0 ? n : 1;
                }
                return ptr;
            }

            T* get() const { return ptr; }
            ptrdiff_t size() const { return count; }
        };
#endif

#if defined(SCCD_ENABLE_CUDA)
        // The device cell list, mirroring the five host members above. cellptr
        // and cursor are sized once from the element count, because the setup
        // caps the grid at 4n cells, so a step allocates only cellidx, whose
        // length is the span count the setup reports back.
        sccd::device::Cell2DGridD<scalar_t> v_grid_d_;
        DeviceArray<ptrdiff_t> v_cellptr_d_;
        DeviceArray<ptrdiff_t> v_cursor_d_;
        DeviceArray<int> v_ranges_d_;
        DeviceArray<smesh::idx_t> v_cellidx_d_;

        sccd::device::Cell2DGridD<scalar_t> e_grid_d_;
        DeviceArray<ptrdiff_t> e_cellptr_d_;
        DeviceArray<ptrdiff_t> e_cursor_d_;
        DeviceArray<int> e_ranges_d_;
        DeviceArray<smesh::idx_t> e_cellidx_d_;
#endif


        sccd::BoxRounding rounding_{sccd::BoxRounding::Exact};

    public:
        /**
         * \brief Build the edge graph and allocate the working buffers.
         *
         * Everything here is a function of the mesh: its topology and its
         * element counts. It happens once for a mesh and never again, however
         * many steps are simulated, so it is not part of what a step costs.
         *
         * A caller never has to call this, because every entry point below does
         * it on first use. It is public so that a caller who measures the
         * per-step cost can pay it up front instead of inside the first step it
         * times, and so that a caller who cares when allocation happens can
         * choose the moment.
         */
        void initialize() {
            if (!vaabb_) {
                init();
            }
        }

        // fv here, vf in the narrow phase, and both are right: the broad phase
        // pairs the face list against the vertex list, in that order, while the
        // narrow phase answers a vertex-against-face query. The core spells the
        // latter narrow_phase_vf, and this class now matches it.
        int broad_phase_fv(const smesh::SharedBuffer<scalar_t*>& points_t0,
                           const smesh::SharedBuffer<scalar_t*>& points_t1,
                           smesh::SharedBuffer<smesh::idx_t>& v_overlap,
                           smesh::SharedBuffer<smesh::idx_t>& f_overlap) {
            int err = SCCD_SUCCESS;
            err |= broad_phase_prep(points_t0, points_t1);
            err |= broad_phase_fv_step(v_overlap, f_overlap);
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
            stamp_broad_phase_elapsed_();
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

            stamp_broad_phase_elapsed_();
            e0_overlap = e0_overlap_;
            e1_overlap = e1_overlap_;
            return err;
        }

        // Run only the FV narrow_phase using the latest FV broad_phase result.
        // `max_toi` upper-bounds the search; on `toi_output == sccd::ToiOutput::Earliest` (shared scalar
        // output) it is updated in place to the earliest TOI observed.
        int narrow_phase_vf(scalar_t& max_toi,
                            smesh::SharedBuffer<scalar_t>& vf_tois,
                            const int max_depth,
                            const scalar_t tol,
                            const sccd::ToiOutput toi_output = sccd::ToiOutput::Earliest) {
            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                SMESH_TRACE_SCOPE("Narrow phase (FV)");
                err |= narrow_phase_vf_step_host_(max_toi, max_depth, tol, toi_output);
            } else {
                SMESH_TRACE_SCOPE("Narrow phase (FV)");
                err |= narrow_phase_vf_step_device_(max_toi, max_depth, tol, toi_output);
            }
            vf_tois = vf_tois_;
            return err;
        }

        /**
         * \brief Run only the EE narrow phase, using the latest EE broad-phase result.
         *
         * `max_toi` is in/out exactly as in `narrow_phase_vf`: it bounds the
         * search on the way in and, for `toi_output == sccd::ToiOutput::Earliest`,
         * comes back holding the earliest time of impact. Take it by reference in
         * any wrapper -- passing it by value leaves the caller's variable
         * untouched while `ee_tois` still fills, so the miss is silent.
         */
        int narrow_phase_ee(scalar_t& max_toi,
                            smesh::SharedBuffer<scalar_t>& ee_tois,
                            const int max_depth,
                            const scalar_t tol,
                            const sccd::ToiOutput toi_output = sccd::ToiOutput::Earliest) {
            int err = SCCD_SUCCESS;
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                SMESH_TRACE_SCOPE("Narrow phase (EE)");
                err |= narrow_phase_ee_step_host_(max_toi, max_depth, tol, toi_output);
            } else {
                SMESH_TRACE_SCOPE("Narrow phase (EE)");
                err |= narrow_phase_ee_step_device_(max_toi, max_depth, tol, toi_output);
            }
            ee_tois = ee_tois_;
            return err;
        }

    private:
        // ---- Host stage helpers ----------------------------------------------

        /** \brief idx[i] = i, in parallel: it is a full pass over the list. */
        static void fill_identity_(smesh::idx_t* const SCCD_RESTRICT idx, const ptrdiff_t n) {
            sccd::parallel_for_br(0, n, [&](const ptrdiff_t begin, const ptrdiff_t end) {
                for (ptrdiff_t i = begin; i < end; ++i) idx[i] = (smesh::idx_t)i;
            });
        }

        /**
         * \brief Size a grid for one list and bin it, count then fill.
         *
         * The buffers are members and survive the call, so a steady state
         * reallocates nothing. cell2d_count writes every entry of cellptr before
         * reading any, so the resize deliberately does not ask for zeros.
         */
        static void bin_host_(const ptrdiff_t n,
                              scalar_t** const SCCD_RESTRICT aabb,
                              sccd::Cell2DGrid<scalar_t>& grid,
                              sccd::Cell2DPartition& part,
                              std::vector<ptrdiff_t>& cellptr,
                              std::vector<smesh::idx_t>& cellidx,
                              std::vector<ptrdiff_t>& cursor) {
            sccd::cell2d_setup<scalar_t>(n, aabb, grid);
            sccd::cell2d_partition<scalar_t>(n, aabb, grid, part);

            cellptr.resize((size_t)grid.ncells() + 1);
            sccd::cell2d_count<scalar_t>(n, aabb, grid, part, cellptr.data());

            cellidx.resize((size_t)cellptr[grid.ncells()]);
            cursor.resize((size_t)grid.ncells());
            sccd::cell2d_fill<scalar_t, smesh::idx_t>(
                n, aabb, grid, part, cellptr.data(), cellidx.data(), cursor.data());
        }

#if defined(SCCD_ENABLE_CUDA)
        /**
         * \brief Size a grid for one list on the device and bin it.
         *
         * The device counterpart of bin_host_, and the same count-then-fill
         * shape. cellptr and cursor are already sized for the grid's 4n cell
         * cap, so only cellidx depends on the step: the setup returns the span
         * count and the buffer grows to it, never shrinking, so a steady state
         * allocates nothing.
         */
        void bin_device_(const ptrdiff_t n,
                         scalar_t** const SCCD_RESTRICT aabb,
                         sccd::device::Cell2DGridD<scalar_t>& grid,
                         DeviceArray<ptrdiff_t>& cellptr,
                         DeviceArray<smesh::idx_t>& cellidx,
                         DeviceArray<ptrdiff_t>& cursor,
                         DeviceArray<int>& ranges) {
            ptrdiff_t spans = 0;
            sccd::device::cell2d_setup_and_count<scalar_t, smesh::idx_t>(
                n, aabb, grid, cellptr.get(), ranges.get(), &spans);

            ptrdiff_t rejected = 0;
            sccd::device::cell2d_fill<scalar_t, smesh::idx_t>(
                n, aabb, grid, cellptr.get(), ranges.get(),
                cellidx.reserve(spans > 0 ? spans : 1), cellidx.size(),
                cursor.get(), &rejected);
            if (rejected) {
                // The counting pass reserves one slot per span, so this cannot
                // happen while the two passes agree. Saying so is better than
                // running the query over a cell list that is quietly missing
                // entries, which would drop candidate pairs.
                SMESH_ERROR("sccd: device cell list rejected %ld of %ld spans over %ld boxes "
                            "(grid %dx%d)\n",
                            (long)rejected, (long)spans, (long)n, grid.n0, grid.n1);
            }
        }
#endif

        /**
         * \brief Pick the strategy for this call and start the broad-phase clock.
         *
         * Both processors run this, because both implement both strategies. The
         * tuner races the two rather than predicting a winner; see
         * broadphase_strategy.hpp for the four heuristics that failed and the
         * constant that failed after them. The clock starts here and stops when
         * the last of the two steps has run, so what is timed is a whole broad
         * phase and not just the part that differs.
         *
         * A measurement is reported here rather than at the end of a step,
         * because the steps are separately callable: a caller may run
         * face-vertex only, edge-edge only, or both in either order, and there
         * is no step that is reliably "last". Reporting the previous broad phase
         * when the next one starts covers every case, and the only measurement
         * lost is the final one of a run.
         */
        void choose_strategy_(const ptrdiff_t n_nodes) {
            if (broad_phase_pending_) {
                tuner_.record(timed_strategy_, broad_phase_pending_ms_);
                broad_phase_pending_ = false;
            }

            const sccd::BroadPhaseStrategy chosen = tuner_.next();
            use_cell2d_min_sorted_ = (chosen == sccd::BroadPhaseStrategy::Cell2DMinSort);
            use_cell2d_min_ =
                (chosen == sccd::BroadPhaseStrategy::Cell2DMin) || use_cell2d_min_sorted_;
            use_cell2d_ = (chosen == sccd::BroadPhaseStrategy::Cell2D) || use_cell2d_min_;
            timed_strategy_ = chosen;

            if (getenv("SCCD_BROADPHASE_VERBOSE")) {
                // The statistics read the vertex boxes, which live on the device
                // when that is where the run is, so the shape of the geometry is
                // only reported for a host run. The choice itself is reported
                // either way, which is what the diagnostic is mainly for.
                if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                    sccd::BroadPhaseStats<scalar_t> stats =
                        sccd::broadphase_stats<scalar_t>(n_nodes, vaabb_->data());
                    fprintf(stderr,
                            "sccd broad phase: %s (lambda %.1f/%.1f/%.1f, anisotropy %.2f, n %ld; "
                            "sweep %.3f ms, cell2d %.3f ms)\n",
                            sccd::broadphase_strategy_name(chosen),
                            (double)stats.lambda[0],
                            (double)stats.lambda[1],
                            (double)stats.lambda[2],
                            stats.anisotropy(),
                            (long)n_nodes,
                            tuner_.sweep_ms(),
                            tuner_.cell2d_ms());
                } else {
                    fprintf(stderr,
                            "sccd broad phase: %s (device, n %ld; sweep %.3f ms, cell2d %.3f ms)\n",
                            sccd::broadphase_strategy_name(chosen),
                            (long)n_nodes,
                            tuner_.sweep_ms(),
                            tuner_.cell2d_ms());
                }
            }

            // The clock starts after the logging, not before it. Computing the
            // statistics and writing an unbuffered line to stderr took about a
            // millisecond, which on a small mesh was most of the "measurement"
            // -- a diagnostic that changes the number it is diagnosing.
            broad_phase_t0_ = sccd::broadphase_now_ms();
        }

        int broad_phase_prep_host_() {
            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const ptrdiff_t n_edges = e0_->size();

            {
                SMESH_TRACE_SCOPE("Broad_phase: AABB");
                sccd::compute_aabbs(n_nodes,
                                    points_t0_->data(),
                                    points_t1_->data(),
                                    vaabb_->data(),
                                    rounding_);

                sccd::compute_aabbs(mesh_->block(0)->n_nodes_per_element(),
                                    n_faces,
                                    faces_->data(),
                                                                        points_t0_->data(),
                                    points_t1_->data(),
                                    faabb_->data(),
                                    rounding_);

                sccd::compute_aabbs(2,
                                    n_edges,
                                    edges_->data(),
                                                                        points_t0_->data(),
                                    points_t1_->data(),
                                    eaabb_->data(),
                                    rounding_);
            }

            choose_strategy_(n_nodes);

            if (use_cell2d_) {
                // No sorting at all on this path: vertices, faces and edges are
                // all binned. Removing the sort is most of the reason to prefer a
                // cell list, so a variant that kept one would have given up the
                // point of it.
                SMESH_TRACE_SCOPE("Cell list (host)");

                fill_identity_(vidx_->data(), n_nodes);
                fill_identity_(fidx_->data(), n_faces);
                fill_identity_(eidx_->data(), n_edges);

                bin_host_(n_nodes, vaabb_->data(), v_grid_, v_part_, v_cellptr_, v_cellidx_, v_cursor_);
                if (use_cell2d_min_) {
                    bin_min_host_(n_edges, eaabb_->data(), e_grid_, e_part_, e_cellptr_, e_cellidx_,
                                  e_cursor_, e_row_prefix_, e_cell_hi1_);
                    if (use_cell2d_min_sorted_) {
                        sort_cells_host_(eaabb_->data(), e_grid_, e_cellptr_, e_cellidx_,
                                         e_cell_key_, e_cell_hi2_);
                    }
                } else {
                    bin_host_(n_edges, eaabb_->data(), e_grid_, e_part_, e_cellptr_, e_cellidx_, e_cursor_);
                }
            } else {
                SMESH_TRACE_SCOPE("Sorting AABBs (host)");

                // Decided on the vertex boxes and used for all three lists: they
                // come from the same geometry, so their anisotropy agrees, and one
                // decision keeps the passes consistent with each other. The cell
                // list needs no sort axis, so it does not pay for this pass.
                sort_axis_ = sccd::choose_axis(n_nodes, vaabb_->data());

                sccd::sort_along_axis(
                    n_nodes, sort_axis_, vaabb_->data(), vidx_->data(), scratch_->data(), &sort_scratch_);
                sccd::sort_along_axis(
                    n_faces, sort_axis_, faabb_->data(), fidx_->data(), scratch_->data(), &sort_scratch_);
                sccd::sort_along_axis(
                    n_edges, sort_axis_, eaabb_->data(), eidx_->data(), scratch_->data(), &sort_scratch_);
            }

            return SCCD_SUCCESS;
        }

        /**
         * \brief The cell-list face-vertex step, for a face with \p nxe vertices.
         *
         * Split out only so the element type picks an instantiation, the way the
         * sweep does below. The cell list itself is over vertices, so nothing here
         * depends on the face topology beyond the vertex count the pair collector
         * needs in order to skip a face's own vertices.
         */
        template <int nxe>
        int cell2d_fv_step_host_(const ptrdiff_t n_faces) {
            sccd::cell2d_count_overlaps<nxe, 1, scalar_t, smesh::idx_t>(n_faces,
                                                                        faabb_->data(),
                                                                        fidx_->data(),
                                                                        1,
                                                                        faces_->data(),
                                                                        vaabb_->data(),
                                                                        vidx_->data(),
                                                                        0,
                                                                        nullptr,
                                                                        v_grid_,
                                                                        v_cellptr_.data(),
                                                                        v_cellidx_.data(),
                                                                        ccdptr_->data());

            const ptrdiff_t n_pairs = ccdptr_->data()[n_faces];
            f_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);
            v_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);

            sccd::cell2d_fill_overlaps<nxe, 1, scalar_t, smesh::idx_t>(n_faces,
                                                                       faabb_->data(),
                                                                       fidx_->data(),
                                                                       1,
                                                                       faces_->data(),
                                                                       vaabb_->data(),
                                                                       vidx_->data(),
                                                                       0,
                                                                       nullptr,
                                                                       v_grid_,
                                                                       v_cellptr_.data(),
                                                                       v_cellidx_.data(),
                                                                       ccdptr_->data(),
                                                                       f_overlap_->data(),
                                                                       v_overlap_->data());
            return SCCD_SUCCESS;
        }

        /**
         * \brief The edge-edge step over the minimum-corner binning.
         *
         * \tparam sorted Whether the cells were ordered on the third axis. The
         * one branch that picks it is in the caller, so neither instantiation
         * carries a test the other needs.
         */
        template <bool sorted>
        int cell2dmin_ee_step_host_(const ptrdiff_t n_edges) {
            const scalar_t* const key = sorted ? e_cell_key_.data() : nullptr;
            const scalar_t* const hi2 = sorted ? e_cell_hi2_.data() : nullptr;

            sccd::cell2dmin_count_self_overlaps<2, sorted, scalar_t, smesh::idx_t>(n_edges,
                                                                                   eaabb_->data(),
                                                                                   eidx_->data(),
                                                                                   1,
                                                                                   edges_->data(),
                                                                                   e_grid_,
                                                                                   e_cellptr_.data(),
                                                                                   e_cellidx_.data(),
                                                                                   e_row_prefix_.data(),
                                                                                   e_cell_hi1_.data(),
                                                                                   key,
                                                                                   hi2,
                                                                                   ccdptr_->data());

            const ptrdiff_t n_pairs = ccdptr_->data()[n_edges];
            e0_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);
            e1_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);

            sccd::cell2dmin_fill_self_overlaps<2, sorted, scalar_t, smesh::idx_t>(n_edges,
                                                                                  eaabb_->data(),
                                                                                  eidx_->data(),
                                                                                  1,
                                                                                  edges_->data(),
                                                                                  e_grid_,
                                                                                  e_cellptr_.data(),
                                                                                  e_cellidx_.data(),
                                                                                  e_row_prefix_.data(),
                                                                                  e_cell_hi1_.data(),
                                                                                  key,
                                                                                  hi2,
                                                                                  ccdptr_->data(),
                                                                                  e0_overlap_->data(),
                                                                                  e1_overlap_->data());
            return SCCD_SUCCESS;
        }

        int broad_phase_fv_step_host_() {
            SMESH_TRACE_SCOPE("Broad_phase: F2V");

            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const auto element_type = face_element_type_;

            if (use_cell2d_) {
                SMESH_TRACE_SCOPE("cell2d f2v");
                if (element_type == smesh::TRISHELL3) {
                    return cell2d_fv_step_host_<3>(n_faces);
                } else if (element_type == smesh::QUADSHELL4) {
                    return cell2d_fv_step_host_<4>(n_faces);
                } else {
                    SMESH_ERROR("Unsupported CCD face element type: %s\n", smesh::type_to_string(element_type));
                    return SCCD_FAILURE;
                }
            }

            {
                SMESH_TRACE_SCOPE("cummax");
                sccd::cummax(n_nodes, vaabb_->data()[dim + sort_axis_], cumulative_max_->data());
            }

            {
                SMESH_TRACE_SCOPE("count_overlaps");
                if (element_type == smesh::TRISHELL3) {
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
                } else if (element_type == smesh::QUADSHELL4) {
                    sccd::count_overlaps<4, 1, scalar_t, smesh::idx_t>(sort_axis_,
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
                } else {
                    SMESH_ERROR("Unsupported CCD face element type: %s\n", smesh::type_to_string(element_type));
                }
            }

            const ptrdiff_t n_vertex_face_overlaps = ccdptr_->data()[n_faces];

            {
                SMESH_TRACE_SCOPE("f2v allocations");
                f_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
                v_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
            }

            {
                SMESH_TRACE_SCOPE("collect_overlaps");
                if (element_type == smesh::TRISHELL3) {
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
                } else if (element_type == smesh::QUADSHELL4) {
                    sccd::collect_overlaps<4, 1, scalar_t, smesh::idx_t>(sort_axis_,
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

            return SCCD_SUCCESS;
        }

        int broad_phase_ee_step_host_() {
            SMESH_TRACE_SCOPE("Broad_phase: E2E");

            const ptrdiff_t n_edges = e0_->size();

            if (use_cell2d_min_) {
                SMESH_TRACE_SCOPE("cell2dmin e2e");
                return use_cell2d_min_sorted_ ? cell2dmin_ee_step_host_<true>(n_edges)
                                              : cell2dmin_ee_step_host_<false>(n_edges);
            }

            if (use_cell2d_) {
                SMESH_TRACE_SCOPE("cell2d e2e");

                sccd::cell2d_count_self_overlaps<2, scalar_t, smesh::idx_t>(n_edges,
                                                                            eaabb_->data(),
                                                                            eidx_->data(),
                                                                            1,
                                                                            edges_->data(),
                                                                            e_grid_,
                                                                            e_cellptr_.data(),
                                                                            e_cellidx_.data(),
                                                                            ccdptr_->data());

                const ptrdiff_t n_pairs = ccdptr_->data()[n_edges];
                e0_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);
                e1_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);

                sccd::cell2d_fill_self_overlaps<2, scalar_t, smesh::idx_t>(n_edges,
                                                                           eaabb_->data(),
                                                                           eidx_->data(),
                                                                           1,
                                                                           edges_->data(),
                                                                           e_grid_,
                                                                           e_cellptr_.data(),
                                                                           e_cellidx_.data(),
                                                                           ccdptr_->data(),
                                                                           e0_overlap_->data(),
                                                                           e1_overlap_->data());
                return SCCD_SUCCESS;
            }

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

        int narrow_phase_vf_step_host_(scalar_t& max_toi,
                                       const int max_depth,
                                       const scalar_t tol,
                                       const sccd::ToiOutput toi_output) {
            SMESH_TRACE_SCOPE("Narrow phase: F2V");

            const auto toi_storage_size = [toi_output](const size_t noverlaps) {
                return toi_output == sccd::ToiOutput::Earliest ? size_t(1) : noverlaps;
            };

            vf_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(v_overlap_->size()), execution_space_);
            const auto element_type = face_element_type_;
            if (element_type == smesh::TRISHELL3) {
                sccd::narrow_phase_vf<scalar_t, smesh::idx_t>(v_overlap_->size(),
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
                                                                 toi_output);
            } else if (element_type == smesh::QUADSHELL4) {
                sccd::narrow_phase_vq<scalar_t, smesh::idx_t>(v_overlap_->size(),
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
                                                                 toi_output);
            } else {
                SMESH_ERROR("Unsupported CCD face element type: %s\n", smesh::type_to_string(element_type));
            }

            if (toi_output == sccd::ToiOutput::Earliest && v_overlap_->size() != 0) {
                max_toi = vf_tois_->data()[0];
            }

            return SCCD_SUCCESS;
        }

        int narrow_phase_ee_step_host_(scalar_t& max_toi,
                                       const int max_depth,
                                       const scalar_t tol,
                                       const sccd::ToiOutput toi_output) {
            SMESH_TRACE_SCOPE("Narrow phase: E2E");

            const auto toi_storage_size = [toi_output](const size_t noverlaps) {
                return toi_output == sccd::ToiOutput::Earliest ? size_t(1) : noverlaps;
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
                                                          toi_output);

            if (toi_output == sccd::ToiOutput::Earliest) {
                max_toi = ee_tois_->data()[0];
            }
            return SCCD_SUCCESS;
        }

        // ---- Device stage helpers --------------------------------------------

#if defined(SCCD_ENABLE_CUDA)
        // The device broad phase launches asynchronously and synchronises only at
        // its read-backs, so a read-back is where a fault from any kernel since
        // the last one actually lands. Its status has to be checked: dropping it
        // leaves the destination at whatever it held before -- zero, for an
        // overlap count -- which then sizes an empty output buffer that the
        // collect kernel is launched over anyway, scattering through null at the
        // offsets the count kernel wrote. The original fault is then invisible and
        // the sticky error is reported against whichever call comes next.
        static void check_read_back(const cudaError_t err, const char* what) {
            if (err != cudaSuccess) {
                SMESH_ERROR("sccd: device read-back of %s failed: %s\n", what, cudaGetErrorString(err));
            }
        }
#endif

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
                    n_nodes, points_t0_->data(), points_t1_->data(), vaabb_->data(), rounding_);

                sccd::device::compute_aabbs(mesh_->block(0)->n_nodes_per_element(),
                                            n_faces,
                                            faces_->data(),
                                                                                points_t0_->data(),
                                            points_t1_->data(),
                                            faabb_->data(),
                                            rounding_);

                sccd::device::compute_aabbs(2,
                                            n_edges,
                                            edges_->data(),
                                                                                points_t0_->data(),
                                            points_t1_->data(),
                                            eaabb_->data(),
                                            rounding_);
            }

            choose_strategy_(n_nodes);

            if (use_cell2d_min_) {
                // The minimum-corner edge-edge query is a host kernel so far.
                // Saying so beats running the extent-binned one under its name,
                // which would put a timing against a label that did not produce
                // it.
                SMESH_ERROR("sccd: SCCD_BROADPHASE=%s has no device implementation yet\n",
                            sccd::broadphase_strategy_name(timed_strategy_));
                return SCCD_FAILURE;
            }

            if (use_cell2d_) {
                // No sorting at all on this path, as on the host: vertices and
                // edges are binned and faces are queried against the vertex
                // grid. The index arrays carry the sweep's permutation if the
                // previous step took the other strategy, so they are reset.
                SMESH_TRACE_SCOPE("Cell list (device)");

                sccd::device::fill_identity<smesh::idx_t>(n_nodes, vidx_->data());
                sccd::device::fill_identity<smesh::idx_t>(n_faces, fidx_->data());
                sccd::device::fill_identity<smesh::idx_t>(n_edges, eidx_->data());

                bin_device_(n_nodes, vaabb_->data(), v_grid_d_, v_cellptr_d_, v_cellidx_d_,
                            v_cursor_d_, v_ranges_d_);
                bin_device_(n_edges, eaabb_->data(), e_grid_d_, e_cellptr_d_, e_cellidx_d_,
                            e_cursor_d_, e_ranges_d_);
            } else {
                SMESH_TRACE_SCOPE("Sorting AABBs (device)");

                sort_axis_ = sccd::device::choose_axis(n_nodes, vaabb_->data());

                sccd::device::sort_along_axis(
                    n_nodes, sort_axis_, vaabb_->data(), vidx_->data(), scratch_->data());

                sccd::device::sort_along_axis(
                    n_faces, sort_axis_, faabb_->data(), fidx_->data(), scratch_->data());

                sccd::device::sort_along_axis(
                    n_edges, sort_axis_, eaabb_->data(), eidx_->data(), scratch_->data());
            }
            return SCCD_SUCCESS;
#else
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

#if defined(SCCD_ENABLE_CUDA)
        /** \brief The device cell list's face-vertex step, the counterpart of
         *         cell2d_fv_step_host_. Faces are queried against the vertex grid. */
        template <int nxe>
        int cell2d_fv_step_device_(const ptrdiff_t n_faces) {
            {
                SMESH_TRACE_SCOPE("cell2d count_overlaps");
                sccd::device::cell2d_count_overlaps<nxe, 1, scalar_t, smesh::idx_t>(n_faces,
                                                                                    faabb_->data(),
                                                                                    fidx_->data(),
                                                                                    1,
                                                                                    faces_->data(),
                                                                                    vaabb_->data(),
                                                                                    vidx_->data(),
                                                                                    0,
                                                                                    (smesh::idx_t**)nullptr,
                                                                                    v_grid_d_,
                                                                                    v_cellptr_d_.get(),
                                                                                    v_cellidx_d_.get(),
                                                                                    ccdptr_->data());
            }

            ptrdiff_t n_pairs = 0;
            check_read_back(
                cudaMemcpy(&n_pairs, ccdptr_->data() + n_faces, sizeof(n_pairs), cudaMemcpyDeviceToHost),
                "the face-vertex overlap count");

            {
                SMESH_TRACE_SCOPE("f2v allocations");
                f_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);
                v_overlap_ = smesh::create_buffer<smesh::idx_t>(n_pairs, execution_space_);
            }

            // An empty output buffer has a null data pointer, so the collect is
            // only safe to launch when there is something to write.
            if (n_pairs > 0) {
                SMESH_TRACE_SCOPE("cell2d collect_overlaps");
                sccd::device::cell2d_collect_overlaps<nxe, 1, scalar_t, smesh::idx_t>(n_faces,
                                                                                      faabb_->data(),
                                                                                      fidx_->data(),
                                                                                      1,
                                                                                      faces_->data(),
                                                                                      vaabb_->data(),
                                                                                      vidx_->data(),
                                                                                      0,
                                                                                      (smesh::idx_t**)nullptr,
                                                                                      v_grid_d_,
                                                                                      v_cellptr_d_.get(),
                                                                                      v_cellidx_d_.get(),
                                                                                      ccdptr_->data(),
                                                                                      f_overlap_->data(),
                                                                                      v_overlap_->data());
            }
            return SCCD_SUCCESS;
        }
#endif

        int broad_phase_fv_step_device_() {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Broad_phase: F2V");

            const int dim = mesh_->spatial_dimension();
            const ptrdiff_t n_nodes = mesh_->n_nodes();
            const ptrdiff_t n_faces = mesh_->block(0)->n_elements();
            const auto element_type = face_element_type_;

            if (use_cell2d_) {
                SMESH_TRACE_SCOPE("cell2d f2v");
                if (element_type == smesh::TRISHELL3) {
                    return cell2d_fv_step_device_<3>(n_faces);
                } else if (element_type == smesh::QUADSHELL4) {
                    return cell2d_fv_step_device_<4>(n_faces);
                } else {
                    SMESH_ERROR("Unsupported CCD face element type: %s\n", smesh::type_to_string(element_type));
                    return SCCD_FAILURE;
                }
            }

            scalar_t* vaabb_max_axis = nullptr;
            check_read_back(cudaMemcpy(&vaabb_max_axis,
                                       vaabb_->data() + dim + sort_axis_,
                                       sizeof(vaabb_max_axis),
                                       cudaMemcpyDeviceToHost),
                            "the vertex AABB row pointer");

            {
                SMESH_TRACE_SCOPE("cummax");
                sccd::device::cummax(n_nodes, vaabb_max_axis, cumulative_max_->data());
            }

            {
                SMESH_TRACE_SCOPE("count_overlaps");
                if (element_type == smesh::TRISHELL3) {
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
                } else if (element_type == smesh::QUADSHELL4) {
                    sccd::device::count_overlaps<4, 1>(sort_axis_,
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
                } else {
                    SMESH_ERROR("Unsupported CCD face element type: %s\n", smesh::type_to_string(element_type));
                }
            }

            ptrdiff_t n_vertex_face_overlaps = 0;
            check_read_back(cudaMemcpy(&n_vertex_face_overlaps,
                                       ccdptr_->data() + n_faces,
                                       sizeof(n_vertex_face_overlaps),
                                       cudaMemcpyDeviceToHost),
                            "the face-vertex overlap count");

            {
                SMESH_TRACE_SCOPE("f2v allocations");
                f_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
                v_overlap_ = smesh::create_buffer<smesh::idx_t>(n_vertex_face_overlaps, execution_space_);
            }

            // An empty output buffer has a null data pointer, so launching the
            // collect over it is only safe when there is nothing to write.
            if (n_vertex_face_overlaps > 0) {
                SMESH_TRACE_SCOPE("collect_overlaps");
                if (element_type == smesh::TRISHELL3) {
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
                } else if (element_type == smesh::QUADSHELL4) {
                    sccd::device::collect_overlaps<4, 1>(sort_axis_,
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

            if (use_cell2d_) {
                SMESH_TRACE_SCOPE("cell2d count_self_overlaps");
                sccd::device::cell2d_count_self_overlaps<2, scalar_t, smesh::idx_t>(n_edges,
                                                                                    eaabb_->data(),
                                                                                    eidx_->data(),
                                                                                    1,
                                                                                    edges_->data(),
                                                                                    e_grid_d_,
                                                                                    e_cellptr_d_.get(),
                                                                                    e_cellidx_d_.get(),
                                                                                    ccdptr_->data());
            } else {
                SMESH_TRACE_SCOPE("count_self_overlaps");
                sccd::device::count_self_overlaps<2>(
                    sort_axis_, n_edges, eaabb_->data(), eidx_->data(), 1, edges_->data(), ccdptr_->data());
            }

            ptrdiff_t n_edge_overlaps = 0;
            check_read_back(
                cudaMemcpy(
                    &n_edge_overlaps, ccdptr_->data() + n_edges, sizeof(n_edge_overlaps), cudaMemcpyDeviceToHost),
                "the edge-edge overlap count");

            {
                SMESH_TRACE_SCOPE("e2e allocations");
                e0_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
                e1_overlap_ = smesh::create_buffer<smesh::idx_t>(n_edge_overlaps, execution_space_);
            }

            if (n_edge_overlaps > 0) {
                if (use_cell2d_) {
                    SMESH_TRACE_SCOPE("cell2d collect_self_overlaps");
                    sccd::device::cell2d_collect_self_overlaps<2, scalar_t, smesh::idx_t>(n_edges,
                                                                                          eaabb_->data(),
                                                                                          eidx_->data(),
                                                                                          1,
                                                                                          edges_->data(),
                                                                                          e_grid_d_,
                                                                                          e_cellptr_d_.get(),
                                                                                          e_cellidx_d_.get(),
                                                                                          ccdptr_->data(),
                                                                                          e0_overlap_->data(),
                                                                                          e1_overlap_->data());
                } else {
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
            }

            return SCCD_SUCCESS;
#else
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int narrow_phase_vf_step_device_(scalar_t& max_toi,
                                         const int max_depth,
                                         const scalar_t tol,
                                         const sccd::ToiOutput toi_output) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Narrow phase: F2V");

            const auto toi_storage_size = [toi_output](const size_t noverlaps) {
                return toi_output == sccd::ToiOutput::Earliest ? size_t(1) : noverlaps;
            };

            vf_tois_ = smesh::create_buffer<scalar_t>(toi_storage_size(v_overlap_->size()), execution_space_);
            const auto element_type = face_element_type_;
            if (element_type == smesh::QUADSHELL4) {
                // Quads have their own device kernel: the triangle one is built
                // around a four-point query on a barycentric domain, and a
                // vertex-quad query is five points on a domain where u and v range
                // independently. See src/cuda/sccd_narrowphase_vq.cuh.
                sccd::device::narrow_phase_vq<scalar_t, smesh::idx_t>(v_overlap_->size(),
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
                                                                      toi_output);
                if (toi_output == sccd::ToiOutput::Earliest) {
                    auto host_toi = smesh::to_host(vf_tois_);
                    max_toi = host_toi->data()[0];
                }
                return SCCD_SUCCESS;
            }
            if (element_type != smesh::TRISHELL3) {
                SMESH_ERROR("Unsupported CCD face element type: %s\n", smesh::type_to_string(element_type));
            }
            sccd::device::narrow_phase_vf(v_overlap_->size(),
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
                                             toi_output);

            if (toi_output == sccd::ToiOutput::Earliest && v_overlap_->size() != 0) {
                auto temp = smesh::to_host(vf_tois_);
                max_toi = temp->data()[0];
            }

            return SCCD_SUCCESS;
#else
            SMESH_UNUSED(max_toi);
            SMESH_UNUSED(max_depth);
            SMESH_UNUSED(tol);
            SMESH_UNUSED(toi_output);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int narrow_phase_ee_step_device_(scalar_t& max_toi,
                                         const int max_depth,
                                         const scalar_t tol,
                                         const sccd::ToiOutput toi_output) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Narrow phase: E2E");

            const auto toi_storage_size = [toi_output](const size_t noverlaps) {
                return toi_output == sccd::ToiOutput::Earliest ? size_t(1) : noverlaps;
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
                                          toi_output);

            if (toi_output == sccd::ToiOutput::Earliest) {
                auto host_toi = smesh::to_host(ee_tois_);
                max_toi = host_toi->data()[0];
            }
            return SCCD_SUCCESS;
#else
            SMESH_UNUSED(max_toi);
            SMESH_UNUSED(max_depth);
            SMESH_UNUSED(tol);
            SMESH_UNUSED(toi_output);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        // ---- Combined-pipeline glue (matches previous behaviour) -------------

        int broad_phase_host(const sccd::ToiOutput toi_output) {
            SMESH_UNUSED(toi_output);

            if (!vaabb_) {
                init();
            }

            SMESH_TRACE_SCOPE("Broad_phase");
            int err = SCCD_SUCCESS;
            err |= broad_phase_prep_host_();
            err |= broad_phase_fv_step_host_();
            err |= broad_phase_ee_step_host_();
            // The one-shot path calls the private steps directly, so it needs its
            // own stamp; the public staged wrappers stamp themselves. Putting it in
            // the private steps instead would mean covering several early returns.
            stamp_broad_phase_elapsed_();
            return err;
        }

        int narrow_phase_host(const sccd::ToiOutput toi_output, const int max_depth, const scalar_t tol) {
            SMESH_TRACE_SCOPE("Narrow phase");
            scalar_t toi = 1;
            int err = SCCD_SUCCESS;
            err |= narrow_phase_vf_step_host_(toi, max_depth, tol, toi_output);
            err |= narrow_phase_ee_step_host_(toi, max_depth, tol, toi_output);
            return err;
        }

        int find_impact_times_impl_host(const sccd::ToiOutput toi_output, const int max_depth, const scalar_t tol) {
            SMESH_TRACE_SCOPE("CCD CPU");
            int err = broad_phase_host(toi_output);
            err |= narrow_phase_host(toi_output, max_depth, tol);
            return err;
        }

        int broad_phase_device(const sccd::ToiOutput toi_output) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_UNUSED(toi_output);

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
            SMESH_UNUSED(toi_output);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int narrow_phase_device(const sccd::ToiOutput toi_output, const int max_depth, const scalar_t tol) {
#if defined(SCCD_ENABLE_CUDA)
            SMESH_TRACE_SCOPE("Narrow phase");
            scalar_t toi = 1;
            int err = SCCD_SUCCESS;
            err |= narrow_phase_vf_step_device_(toi, max_depth, tol, toi_output);
            err |= narrow_phase_ee_step_device_(toi, max_depth, tol, toi_output);
            return err;
#else
            SMESH_UNUSED(toi_output);
            SMESH_UNUSED(max_depth);
            SMESH_UNUSED(tol);
            SMESH_ERROR("Not implemented");
            return SCCD_SUCCESS;
#endif
        }

        int find_impact_times_impl_device(const sccd::ToiOutput toi_output, const int max_depth, const scalar_t tol) {
            SMESH_TRACE_SCOPE("CCD GPU");
            int err = broad_phase_device(toi_output);
            err |= narrow_phase_device(toi_output, max_depth, tol);
            return err;
        }

        int find_impact_times_impl(const sccd::ToiOutput toi_output, const int max_depth, const scalar_t tol) {
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                return find_impact_times_impl_host(toi_output, max_depth, tol);
            } else {
                return find_impact_times_impl_device(toi_output, max_depth, tol);
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
            if (execution_space_ == smesh::EXECUTION_SPACE_HOST) {
                sort_scratch_.reserve(std::max(n_nodes, std::max(n_faces, n_edges)));
            }
#if defined(SCCD_ENABLE_CUDA)
            else {
                // The cell list's row pointer and cursor, sized for the 4n cell
                // cap the grid setup enforces so that no step reallocates them.
                // One-time per mesh, like everything else here, and so outside
                // what a step is timed for.
                // Sized for the grid's 4n cell cap, so a step reallocates
                // nothing. One-time per mesh, and so outside what a step is
                // timed for.
                v_cellptr_d_.reserve(4 * n_nodes + 2);
                v_cursor_d_.reserve(4 * n_nodes + 2);
                v_ranges_d_.reserve(4 * n_nodes + 4);
                e_cellptr_d_.reserve(4 * n_edges + 2);
                e_cursor_d_.reserve(4 * n_edges + 2);
                e_ranges_d_.reserve(4 * n_edges + 4);
            }
#endif
        }
    };
}  // namespace sccd

#endif  // SCCD_SMESH_CCD_HPP
