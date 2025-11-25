#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP

#include "sbisect.hpp"
#include "stuv.hpp"
#include "snumerr.hpp"
#include "snumtol.hpp"

// TODO reimplement ccd_kernel: 
// - Use the above headers for the helper functions
// - No CUDA
// - No Queue no structs, no classes, no dynamic memory allocation
// - implement a purely serial code of the ccd_kernel by contracting the tuvs in place

// namespace sccd {

// struct Interval {
//   __device__ Interval() = default;
//   __device__ Interval(const Scalar &l, const Scalar &u) : lower(l), upper(u) {}
//   Scalar lower;
//   Scalar upper;
// };

// struct SplitInterval {
//   __device__ SplitInterval(const Interval &interval) {
//     const Scalar mid = (interval.lower + interval.upper) / 2;
//     first = Interval(interval.lower, mid);
//     second = Interval(mid, interval.upper);
//   }

//   Interval first;
//   Interval second;
// };

// class CCDDomain {
// public:
//   __device__ void init(int i) {
//     tuv[0] = Interval(0, 1);
//     tuv[1] = Interval(0, 1);
//     tuv[2] = Interval(0, 1);
//     query_id = i;
//   }

//   /// @brief The intervals for the t, u, and v parameters
//   Interval tuv[3];
//   /// @brief The query id
//   int query_id;
// };

// // this is to calculate the vertices of the inclusion function
// struct DomainCorner {
//   /// @brief Update the t, u, and v parameters based on the corner.
//   /// @param domain Domain intervals
//   /// @param corner The corner to use (the first bit is for the t parameter,
//   /// second for the u parameter, and third for the v parameter)
//   __device__ void update_tuv(const CCDDomain &domain, const uint8_t corner) {
//     t = (corner & 1) ? domain.tuv[0].upper : domain.tuv[0].lower;
//     u = (corner & 2) ? domain.tuv[1].upper : domain.tuv[1].lower;
//     v = (corner & 4) ? domain.tuv[2].upper : domain.tuv[2].lower;
//   }

//   Scalar t; ///< @brief The t parameter at the corner
//   Scalar u; ///< @brief The u parameter at the corner
//   Scalar v; ///< @brief The v parameter at the corner
// };

// struct CCDBuffer {
//   /// @brief Check if the buffer is full
//   __device__ bool is_full() const {
//     return (m_tail + 1) % m_capacity == m_head;
//   }

//   /// @brief Check if the buffer is empty
//   __device__ bool is_empty() const { return m_tail == m_head; }

//   // __device__ CCDDomain pop()
//   // {
//   //     if (!is_empty()) {
//   //         return m_data[atomicInc(&m_head, m_capacity - 1)];
//   //     }
//   //     assert(false);
//   // }

//   __device__ CCDDomain &push(const CCDDomain &val) {
//     if (is_full()) {
//       atomicCAS(&m_overflow_flag, 0, 1);
//       return m_data[m_tail]; // Return a dummy value
//     }
//     const int i = atomicInc(&m_tail, m_capacity - 1);
//     m_data[i] = val;
//     return m_data[i];
//   }

//   __device__ CCDDomain &operator[](const int i) {
//     return m_data[(i + m_head) % m_capacity];
//   }

//   __device__ void shift_queue_start() {
//     // Update the head to the new starting position (assuming all starting
//     // elements were consumed)
//     m_head = (m_head + m_starting_size) % m_capacity;
//     // Calculate the new starting size
//     if (m_head <= m_tail) {
//       m_starting_size = m_tail - m_head;
//     } else {
//       m_starting_size = m_capacity - m_head + m_tail;
//     }
//   }

//   __device__ unsigned starting_size() const { return m_starting_size; }
//   __device__ unsigned capacity() const { return m_capacity; }
//   __device__ unsigned head() const { return m_head; }
//   __device__ unsigned tail() const { return m_tail; }
//   __device__ int overflow_flag() const { return m_overflow_flag; }

//   // These should be private, but we need to access them from the kernel.
//   // private:
//   CCDDomain *m_data;
//   unsigned m_starting_size = 0;
//   unsigned m_capacity = 0;
//   unsigned m_head = 0;
//   unsigned m_tail = 0;
//   int m_overflow_flag = 0;
// };

// __global__ void initialize_buffer(CCDBuffer *buffer) {
//   int tx = threadIdx.x + blockIdx.x * blockDim.x;
//   if (tx >= buffer->starting_size())
//     return;
//   assert(buffer->head() == 0 && buffer->tail() >= tx);
//   buffer->m_data[tx].init(tx);
// }

// __global__ void shift_queue_start(CCDBuffer *buffer) {
//   assert(threadIdx.x == 0 && blockDim.x == 1 && gridDim.x == 1);
//   buffer->shift_queue_start();
// }

// struct CCDData {
//   Vector3 v0s;
//   Vector3 v1s;
//   Vector3 v2s;
//   Vector3 v3s;
//   Vector3 v0e;
//   Vector3 v1e;
//   Vector3 v2e;
//   Vector3 v3e;
//   Array3 err; // error bound of each query, calculated from each scene
//   Array3 tol; // domain tolerance to help decide which dimension to split
//   Scalar ms;  // minimum separation
// #ifdef SCALABLE_CCD_TOI_PER_QUERY
//   Scalar toi;
//   int aid;
//   int bid;
// #endif
//   int nbr_checks = 0;
// };

// __device__ inline bool sum_less_than_one(const Scalar &num1,
//                                          const Scalar &num2) {
// #ifdef SCALABLE_CCD_USE_DOUBLE
//   return num1 + num2 <= 1 / (1 - DBL_EPSILON);
// #else
//   return num1 + num2 <= 1 / (1 - FLT_EPSILON);
// #endif
// }

// __device__ inline Scalar max_Linf_4(const Vector3 &p1, const Vector3 &p2,
//                                     const Vector3 &p3, const Vector3 &p4,
//                                     const Vector3 &p1e, const Vector3 &p2e,
//                                     const Vector3 &p3e, const Vector3 &p4e) {
//   return max(max((p1e - p1).lpNorm<Eigen::Infinity>(),
//                  (p2e - p2).lpNorm<Eigen::Infinity>()),
//              max((p3e - p3).lpNorm<Eigen::Infinity>(),
//                  (p4e - p4).lpNorm<Eigen::Infinity>()));
// }

// __device__ void compute_face_vertex_tolerance(CCDData &data_in) {
//   const Vector3 p000 = data_in.v0s - data_in.v1s;
//   const Vector3 p001 = data_in.v0s - data_in.v3s;
//   const Vector3 p011 = data_in.v0s - (data_in.v2s + data_in.v3s - data_in.v1s);
//   const Vector3 p010 = data_in.v0s - data_in.v2s;
//   const Vector3 p100 = data_in.v0e - data_in.v1e;
//   const Vector3 p101 = data_in.v0e - data_in.v3e;
//   const Vector3 p111 = data_in.v0e - (data_in.v2e + data_in.v3e - data_in.v1e);
//   const Vector3 p110 = data_in.v0e - data_in.v2e;

//   data_in.tol[0] =
//       CONFIG.co_domain_tolerance /
//       (3 * max_Linf_4(p000, p001, p011, p010, p100, p101, p111, p110));
//   data_in.tol[1] =
//       CONFIG.co_domain_tolerance /
//       (3 * max_Linf_4(p000, p100, p101, p001, p010, p110, p111, p011));
//   data_in.tol[2] =
//       CONFIG.co_domain_tolerance /
//       (3 * max_Linf_4(p000, p100, p110, p010, p001, p101, p111, p011));
// }

// __device__ void compute_edge_edge_tolerance(CCDData &data_in) {
//   // WARNING: This differs from the original implementation in
//   // Tight-Inclusion! This may be a bug, but results in better accuracy.
//   const Vector3 p000 = data_in.v0s - data_in.v2s;
//   const Vector3 p001 = data_in.v0s - data_in.v3s;
//   const Vector3 p010 = data_in.v1s - data_in.v2s;
//   const Vector3 p011 = data_in.v1s - data_in.v3s;
//   const Vector3 p100 = data_in.v0e - data_in.v2e;
//   const Vector3 p101 = data_in.v0e - data_in.v3e;
//   const Vector3 p110 = data_in.v1e - data_in.v2e;
//   const Vector3 p111 = data_in.v1e - data_in.v3e;

//   data_in.tol[0] =
//       CONFIG.co_domain_tolerance /
//       (3 * max_Linf_4(p000, p001, p011, p010, p100, p101, p111, p110));
//   data_in.tol[1] =
//       CONFIG.co_domain_tolerance /
//       (3 * max_Linf_4(p000, p001, p011, p010, p100, p101, p111, p110));
//   data_in.tol[2] =
//       CONFIG.co_domain_tolerance /
//       (3 * max_Linf_4(p000, p100, p101, p001, p010, p110, p111, p011));
// }

// template <bool is_vf>
// __device__ __host__ void get_numerical_error(CCDData &data_in, bool use_ms) {
//   Scalar filter;
//   if (!use_ms) {
// #ifdef SCALABLE_CCD_USE_DOUBLE
//     if constexpr (is_vf) {
//       filter = 6.661338147750939e-15;
//     } else {
//       filter = 6.217248937900877e-15;
//     }
// #else
//     if constexpr (is_vf) {
//       filter = 3.576279e-06;
//     } else {
//       filter = 3.337861e-06;
//     }
// #endif
//   } else {
// #ifdef SCALABLE_CCD_USE_DOUBLE
//     if constexpr (is_vf) {
//       filter = 7.549516567451064e-15;
//     } else {
//       filter = 7.105427357601002e-15;
//     }
// #else
//     if constexpr (is_vf) {
//       filter = 4.053116e-06;
//     } else {
//       filter = 3.814698e-06;
//     }
// #endif
//   }

//   const Vector3 max = data_in.v0s.cwiseAbs()
//                           .cwiseMax(data_in.v1s.cwiseAbs())
//                           .cwiseMax(data_in.v2s.cwiseAbs())
//                           .cwiseMax(data_in.v3s.cwiseAbs())
//                           .cwiseMax(data_in.v0e.cwiseAbs())
//                           .cwiseMax(data_in.v1e.cwiseAbs())
//                           .cwiseMax(data_in.v2e.cwiseAbs())
//                           .cwiseMax(data_in.v3e.cwiseAbs())
//                           .cwiseMax(Vector3::Ones());

//   data_in.err = max.array() * max.array() * max.array() * filter;
// }

// __device__ Vector3 calculate_vf(const CCDData &data_in,
//                                 const DomainCorner &tuv) {
//   const Vector3 v = (data_in.v0e - data_in.v0s) * tuv.t + data_in.v0s;
//   const Vector3 t0 = (data_in.v1e - data_in.v1s) * tuv.t + data_in.v1s;
//   const Vector3 t1 = (data_in.v2e - data_in.v2s) * tuv.t + data_in.v2s;
//   const Vector3 t2 = (data_in.v3e - data_in.v3s) * tuv.t + data_in.v3s;
//   return v - (t1 - t0) * tuv.u - (t2 - t0) * tuv.v - t0;
// }

// __device__ Vector3 calculate_ee(const CCDData &data_in,
//                                 const DomainCorner &tuv) {
//   const Vector3 ea0 = (data_in.v0e - data_in.v0s) * tuv.t + data_in.v0s;
//   const Vector3 ea1 = (data_in.v1e - data_in.v1s) * tuv.t + data_in.v1s;
//   const Vector3 eb0 = (data_in.v2e - data_in.v2s) * tuv.t + data_in.v2s;
//   const Vector3 eb1 = (data_in.v3e - data_in.v3s) * tuv.t + data_in.v3s;
//   return ((ea1 - ea0) * tuv.u + ea0) - ((eb1 - eb0) * tuv.v + eb0);
// }

// template <bool is_vf>
// __device__ bool origin_in_inclusion_function(const CCDData &data_in,
//                                              const CCDDomain &domain,
//                                              Scalar &true_tol, bool &box_in) {
//   Array3 codomain_min = Array3::Constant(SCALAR_MAX);
//   Array3 codomain_max = Array3::Constant(-SCALAR_MAX);

//   DomainCorner domain_corner;
//   for (uint8_t corner = 0; corner < 8; corner++) {
//     domain_corner.update_tuv(domain, corner);

//     Vector3 codomain_corner;
//     if constexpr (is_vf) {
//       codomain_corner = calculate_vf(data_in, domain_corner);
//     } else {
//       codomain_corner = calculate_ee(data_in, domain_corner);
//     }

//     codomain_min = codomain_min.min(codomain_corner.array());
//     codomain_max = codomain_max.max(codomain_corner.array());
//   }

//   // get the min and max in one dimension
//   true_tol = max(0.0, (codomain_max - codomain_min).maxCoeff());

//   box_in = true;

//   if ((codomain_min - data_in.ms > data_in.err).any() ||
//       (codomain_max + data_in.ms < -data_in.err).any()) {
//     return false;
//   }

//   if ((codomain_min + data_in.ms < -data_in.err).any() ||
//       (codomain_max - data_in.ms > data_in.err).any()) {
//     box_in = false;
//   }

//   return true;
// }

// __device__ int split_dimension(const CCDData &data, Array3 width) {
//   const Array3 res = width / data.tol;
//   if (res[0] >= res[1] && res[0] >= res[2]) {
//     return 0;
//   } else if (res[1] >= res[0] && res[1] >= res[2]) {
//     return 1;
//   } else {
//     assert(res[2] >= res[0] && res[2] >= res[1]);
//     return 2;
//   }
// }

// template <bool is_vf>
// __device__ inline bool bisect(const CCDDomain &domain, const int split,
//                               const Scalar *const toi,
//                               CCDBuffer *const buffer) {
//   const SplitInterval halves(domain.tuv[split]);

//   if (halves.first.lower >= halves.first.upper ||
//       halves.second.lower >= halves.second.upper) {
//     return true;
//   }

//   buffer->push(domain).tuv[split] = halves.first;

//   if (split == 0) {
//     if (halves.second.lower <= *toi) {
//       buffer->push(domain).tuv[0] = halves.second;
//     }
//   } else {
//     if constexpr (is_vf) {
//       if (split == 1) {
//         // check if u+v<=1
//         if (sum_less_than_one(halves.second.lower, domain.tuv[2].lower)) {
//           buffer->push(domain).tuv[1] = halves.second;
//         }
//       } else if (split == 2) {
//         // check if u+v<=1
//         if (sum_less_than_one(halves.second.lower, domain.tuv[1].lower)) {
//           buffer->push(domain).tuv[2] = halves.second;
//         }
//       }
//     } else {
//       buffer->push(domain).tuv[split] = halves.second;
//     }
//   }

//   return false;
// }

// } // namespace sccd

// // === the memory pool method =================================================

// template <bool is_vf>
// __global__ void compute_tolerance(CCDData *data, const int query_size) {
//   int tx = threadIdx.x + blockIdx.x * blockDim.x;
//   if (tx >= query_size)
//     return;

//   if constexpr (is_vf) {
//     compute_face_vertex_tolerance(data[tx]);
//   } else {
//     compute_edge_edge_tolerance(data[tx]);
//   }

//   data[tx].nbr_checks = 0;
//   get_numerical_error<is_vf>(data[tx], CONFIG.use_ms);
// }

// template <bool is_vf>
// __global__ void ccd_kernel(CCDBuffer *const buffer, CCDData *const data,
//                            Scalar *const toi) {
//   const int tx = threadIdx.x + blockIdx.x * blockDim.x;
//   if (tx >= buffer->starting_size())
//     return;

//   // Get the tx element from the buffer without advancing the head.
//   // This allows for better coalescing of memory access compared to pop().
//   const CCDDomain domain_in = (*buffer)[tx];
//   const CCDData data_in = data[domain_in.query_id];
//   atomicAdd(&data[domain_in.query_id].nbr_checks, 1);

//   const Scalar min_t = domain_in.tuv[0].lower; // the time of this domain

// #ifndef SCALABLE_CCD_TOI_PER_QUERY
//   // if the time is larger than toi, return
//   if (min_t >= *toi) {
// #else
//   if (min_t >= data_in.toi) {
// #endif
//     return;
//   }

//   // Check if exceeds max checks
//   if (CONFIG.max_iter >= 0 && data_in.nbr_checks > CONFIG.max_iter) {
//     return;
//   }
//   // else if (buffer.mp_remaining > buffer->capacity() / 2) { // overflow
//   //     atomicCAS(&buffer->overflow_flag, 0, 1);
//   //     return;
//   // }

//   Scalar true_tol = 0; // set by origin_in_inclusion_function
//   bool box_in;         // set by origin_in_inclusion_function
//   if (origin_in_inclusion_function<is_vf>(data_in, domain_in, true_tol,
//                                           box_in)) {

//     const Array3 widths(domain_in.tuv[0].upper - domain_in.tuv[0].lower,
//                         domain_in.tuv[1].upper - domain_in.tuv[1].lower,
//                         domain_in.tuv[2].upper - domain_in.tuv[2].lower);

//     // Condition 1: the domain is smaller than the tolerance.
//     if ((widths <= data_in.tol).all()) {
//       atomicMin(toi, min_t);
// #ifdef SCALABLE_CCD_TOI_PER_QUERY
//       atomicMin(&data[domain_in.query_id].toi, min_t);
// #endif
//       return;
//     }

//     // Condition 2: the box is inside the epsilon box
//     if (box_in && (CONFIG.allow_zero_toi || min_t > 0)) {
//       atomicMin(toi, min_t);
// #ifdef SCALABLE_CCD_TOI_PER_QUERY
//       atomicMin(&data[domain_in.query_id].toi, min_t);
// #endif
//       return;
//     }

//     // Condition 3: real tolerance is smaller than the int tolerance
//     if (true_tol <= CONFIG.co_domain_tolerance &&
//         (CONFIG.allow_zero_toi || min_t > 0)) {
//       atomicMin(toi, min_t);
// #ifdef SCALABLE_CCD_TOI_PER_QUERY
//       atomicMin(&data[domain_in.query_id].toi, min_t);
// #endif
//       return;
//     }

//     // Get the next dimension to split
//     const int split = split_dimension(data_in, widths);

//     const bool sure_in = bisect<is_vf>(domain_in, split,
// #ifdef SCALABLE_CCD_TOI_PER_QUERY
//                                        &data_in.toi,
// #else
//                                        toi,
// #endif
//                                        buffer);

//     // Condition 4 (rare): the interval is too small that overflow happens
//     if (sure_in) {
//       atomicMin(toi, min_t);
// #ifdef SCALABLE_CCD_TOI_PER_QUERY
//       atomicMin(&data[domain_in.query_id].toi, min_t);
// #endif
//       return;
//     }
//   }
// }

// template <bool is_vf>
// bool ccd(thrust::device_vector<CCDData> &d_data,
//          const std::shared_ptr<MemoryHandler> memory_handler,
//          const int parallel_nbr, const int max_iter, const Scalar tol,
//          const bool use_ms, const bool allow_zero_toi, Scalar &toi) {
//   const int nbr = d_data.size();

//   CCDBuffer *d_buffer;
//   {
//     const size_t unit_size = memory_handler->MAX_UNIT_SIZE;

//     logger().trace("CCD Buffer of size {:d} ({:g} GB)", unit_size,
//                    sizeof(CCDDomain) * unit_size / 1e9);

//     CCDBuffer buffer;
//     gpuErrchk(cudaMalloc(&buffer.m_data, sizeof(CCDDomain) * unit_size));
//     buffer.m_starting_size = nbr;
//     buffer.m_capacity = unit_size;
//     buffer.m_head = 0;
//     buffer.m_tail = nbr;
//     buffer.m_overflow_flag = 0;

//     gpuErrchk(cudaMalloc(&d_buffer, sizeof(CCDBuffer)));
//     gpuErrchk(cudaMemcpy(d_buffer, &buffer, sizeof(CCDBuffer),
//                          cudaMemcpyHostToDevice));

//     initialize_buffer<<<nbr / parallel_nbr + 1, parallel_nbr>>>(d_buffer);
//     gpuErrchk(cudaDeviceSynchronize());
//   }

//   // Initialize the global configuration variable
//   {
//     CCDConfig config;
//     config.co_domain_tolerance = tol;
//     config.use_ms = use_ms;
//     config.allow_zero_toi = allow_zero_toi;
//     config.max_iter = max_iter;
//     gpuErrchk(cudaMemcpyToSymbol(CONFIG, &config, sizeof(CCDConfig)));
//   }

//   DeviceVariable d_toi(toi);

//   // ---

//   compute_tolerance<is_vf><<<nbr / parallel_nbr + 1, parallel_nbr>>>(
//       thrust::raw_pointer_cast(d_data.data()), nbr);
//   gpuErrchk(cudaDeviceSynchronize());

//   logger().trace("Max queries: {:d}", memory_handler->MAX_QUERIES);

//   int nbr_per_loop = nbr;
//   logger().trace("Queue size t0: {:d}", nbr_per_loop);
//   while (nbr_per_loop > 0) {
//     ccd_kernel<is_vf><<<nbr_per_loop / parallel_nbr + 1, parallel_nbr>>>(
//         d_buffer, thrust::raw_pointer_cast(d_data.data()), &d_toi);
//     gpuErrchk(cudaDeviceSynchronize());

//     // Update the starting index to reflect that the starting_size was
//     // procesed.
//     shift_queue_start<<<1, 1>>>(d_buffer);
//     gpuErrchk(cudaDeviceSynchronize());

//     gpuErrchk(cudaMemcpy(&nbr_per_loop, &(d_buffer->m_starting_size),
//                          sizeof(unsigned), cudaMemcpyDeviceToHost));

//     logger().trace("Queue size: {:d}", nbr_per_loop);
//     logger().trace("toi={:g}", *d_toi);
//   }
//   gpuErrchk(cudaDeviceSynchronize());

//   toi = d_toi;

//   int overflow;
//   gpuErrchk(cudaMemcpy(&overflow, &(d_buffer->m_overflow_flag), sizeof(int),
//                        cudaMemcpyDeviceToHost));
//   return overflow;
}

#endif // NARROWPHASE_HPP
