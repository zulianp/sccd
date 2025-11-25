// VF
/* Auto-generated with SymPy (CSE + pow expansion) */

template <typename T>
static inline void vf_objective(const T sv0[3], const T s1[3], const T s2[3],
                                const T s3[3], const T s4[3], const T ev0[3],
                                const T e1[3], const T e2[3], const T e3[3],
                                const T e4[3], const T t, const T u, const T v,
                                T *out_f) {
  const T ssa0 = t - 1;
  const T ssa1 = u + v - 1;
  *out_f =
      sccd::pow2<T>(
          -ev[0] * t + ssa0 * sv[0] - ssa1 * (e1[0] * t - s1[0] * ssa0) +
          u * (e2[0] * t - s2[0] * ssa0) + v * (e3[0] * t - s3[0] * ssa0)) +
      sccd::pow2<T>(
          -ev[1] * t + ssa0 * sv[1] - ssa1 * (e1[1] * t - s1[1] * ssa0) +
          u * (e2[1] * t - s2[1] * ssa0) + v * (e3[1] * t - s3[1] * ssa0)) +
      sccd::pow2<T>(
          -ev[2] * t + ssa0 * sv[2] - ssa1 * (e1[2] * t - s1[2] * ssa0) +
          u * (e2[2] * t - s2[2] * ssa0) + v * (e3[2] * t - s3[2] * ssa0));
}

template <typename T>
static inline void vf_gradient(const T sv0[3], const T s1[3], const T s2[3],
                               const T s3[3], const T s4[3], const T ev0[3],
                               const T e1[3], const T e2[3], const T e3[3],
                               const T e4[3], const T t, const T u, const T v,
                               T out_g[3]) {
  const T ssa0 = u + v - 1;
  const T ssa1 = t - 1;
  const T ssa2 = e2[0] * t;
  const T ssa3 = s2[0] * ssa1;
  const T ssa4 = e3[0] * t;
  const T ssa5 = s3[0] * ssa1;
  const T ssa6 = e1[0] * t - s1[0] * ssa1;
  const T ssa7 = -2 * ev[0] * t - 2 * ssa0 * ssa6 + 2 * ssa1 * sv[0] +
                 2 * u * (ssa2 - ssa3) + 2 * v * (ssa4 - ssa5);
  const T ssa8 = e2[1] * t;
  const T ssa9 = s2[1] * ssa1;
  const T ssa10 = e3[1] * t;
  const T ssa11 = s3[1] * ssa1;
  const T ssa12 = e1[1] * t - s1[1] * ssa1;
  const T ssa13 = -2 * ev[1] * t - 2 * ssa0 * ssa12 + 2 * ssa1 * sv[1] +
                  2 * u * (ssa8 - ssa9) + 2 * v * (ssa10 - ssa11);
  const T ssa14 = e2[2] * t;
  const T ssa15 = s2[2] * ssa1;
  const T ssa16 = e3[2] * t;
  const T ssa17 = s3[2] * ssa1;
  const T ssa18 = e1[2] * t - s1[2] * ssa1;
  const T ssa19 = -2 * ev[2] * t - 2 * ssa0 * ssa18 + 2 * ssa1 * sv[2] +
                  2 * u * (ssa14 - ssa15) + 2 * v * (ssa16 - ssa17);
  out_g[0] = ssa13 * (-ev[1] - ssa0 * (e1[1] - s1[1]) + sv[1] +
                      u * (e2[1] - s2[1]) + v * (e3[1] - s3[1])) +
             ssa19 * (-ev[2] - ssa0 * (e1[2] - s1[2]) + sv[2] +
                      u * (e2[2] - s2[2]) + v * (e3[2] - s3[2])) +
             ssa7 * (-ev[0] - ssa0 * (e1[0] - s1[0]) + sv[0] +
                     u * (e2[0] - s2[0]) + v * (e3[0] - s3[0]));
  out_g[1] = -ssa13 * (ssa12 - ssa8 + ssa9) - ssa19 * (-ssa14 + ssa15 + ssa18) -
             ssa7 * (-ssa2 + ssa3 + ssa6);
  out_g[2] = -ssa13 * (-ssa10 + ssa11 + ssa12) -
             ssa19 * (-ssa16 + ssa17 + ssa18) - ssa7 * (-ssa4 + ssa5 + ssa6);
}

template <typename T>
static inline void vf_hessian(const T sv0[3], const T s1[3], const T s2[3],
                              const T s3[3], const T s4[3], const T ev0[3],
                              const T e1[3], const T e2[3], const T e3[3],
                              const T e4[3], const T t, const T u, const T v,
                              T out_H[9]) {
  const T ssa0 = u + v - 1;
  const T ssa1 = e1[0] - s1[0];
  const T ssa2 =
      -ev[0] - ssa0 * ssa1 + sv[0] + u * (e2[0] - s2[0]) + v * (e3[0] - s3[0]);
  const T ssa3 = e1[1] - s1[1];
  const T ssa4 =
      -ev[1] - ssa0 * ssa3 + sv[1] + u * (e2[1] - s2[1]) + v * (e3[1] - s3[1]);
  const T ssa5 = e1[2] - s1[2];
  const T ssa6 =
      -ev[2] - ssa0 * ssa5 + sv[2] + u * (e2[2] - s2[2]) + v * (e3[2] - s3[2]);
  const T ssa7 = t - 1;
  const T ssa8 = s2[0] * ssa7;
  const T ssa9 = e2[0] * t;
  const T ssa10 = e1[0] * t - s1[0] * ssa7;
  const T ssa11 = ssa10 + ssa8 - ssa9;
  const T ssa12 = 2 * ssa2;
  const T ssa13 = s2[1] * ssa7;
  const T ssa14 = e2[1] * t;
  const T ssa15 = e1[1] * t - s1[1] * ssa7;
  const T ssa16 = ssa13 - ssa14 + ssa15;
  const T ssa17 = 2 * ssa4;
  const T ssa18 = s2[2] * ssa7;
  const T ssa19 = e2[2] * t;
  const T ssa20 = e1[2] * t - s1[2] * ssa7;
  const T ssa21 = ssa18 - ssa19 + ssa20;
  const T ssa22 = 2 * ssa6;
  const T ssa23 = e3[0] * t;
  const T ssa24 = s3[0] * ssa7;
  const T ssa25 = -2 * ev[0] * t - 2 * ssa0 * ssa10 + 2 * ssa7 * sv[0] +
                  2 * u * (-ssa8 + ssa9) + 2 * v * (ssa23 - ssa24);
  const T ssa26 = e3[1] * t;
  const T ssa27 = s3[1] * ssa7;
  const T ssa28 = -2 * ev[1] * t - 2 * ssa0 * ssa15 + 2 * ssa7 * sv[1] +
                  2 * u * (-ssa13 + ssa14) + 2 * v * (ssa26 - ssa27);
  const T ssa29 = e3[2] * t;
  const T ssa30 = s3[2] * ssa7;
  const T ssa31 = -2 * ev[2] * t - 2 * ssa0 * ssa20 + 2 * ssa7 * sv[2] +
                  2 * u * (-ssa18 + ssa19) + 2 * v * (ssa29 - ssa30);
  const T ssa32 = -ssa11 * ssa12 - ssa16 * ssa17 - ssa21 * ssa22 -
                  ssa25 * (-e2[0] + s2[0] + ssa1) -
                  ssa28 * (-e2[1] + s2[1] + ssa3) -
                  ssa31 * (-e2[2] + s2[2] + ssa5);
  const T ssa33 = ssa10 - ssa23 + ssa24;
  const T ssa34 = ssa15 - ssa26 + ssa27;
  const T ssa35 = ssa20 - ssa29 + ssa30;
  const T ssa36 = -ssa12 * ssa33 - ssa17 * ssa34 - ssa22 * ssa35 -
                  ssa25 * (-e3[0] + s3[0] + ssa1) -
                  ssa28 * (-e3[1] + s3[1] + ssa3) -
                  ssa31 * (-e3[2] + s3[2] + ssa5);
  const T ssa37 = 2 * ssa11 * ssa33 + 2 * ssa16 * ssa34 + 2 * ssa21 * ssa35;
  out_H[0] = 2 * sccd::pow2<T>(ssa2) + 2 * sccd::pow2<T>(ssa4) +
             2 * sccd::pow2<T>(ssa6);
  out_H[1] = ssa32;
  out_H[2] = ssa36;
  out_H[3] = ssa32;
  out_H[4] = 2 * sccd::pow2<T>(ssa11) + 2 * sccd::pow2<T>(ssa16) +
             2 * sccd::pow2<T>(ssa21);
  out_H[5] = ssa37;
  out_H[6] = ssa36;
  out_H[7] = ssa37;
  out_H[8] = 2 * sccd::pow2<T>(ssa33) + 2 * sccd::pow2<T>(ssa34) +
             2 * sccd::pow2<T>(ssa35);
}

template <typename T>
static inline void
vf_objective_dir(const T sv0[3], const T s1[3], const T s2[3], const T s3[3],
                 const T s4[3], const T ev0[3], const T e1[3], const T e2[3],
                 const T e3[3], const T e4[3], const T t, const T u, const T v,
                 T *out_f, T out_p[3]) {
  const T ssa0 = t - 1;
  const T ssa1 = e2[0] * t;
  const T ssa2 = s2[0] * ssa0;
  const T ssa3 = e3[0] * t;
  const T ssa4 = s3[0] * ssa0;
  const T ssa5 = u + v - 1;
  const T ssa6 = e1[0] * t - s1[0] * ssa0;
  const T ssa7 = -ev[0] * t + ssa0 * sv[0] - ssa5 * ssa6 + u * (ssa1 - ssa2) +
                 v * (ssa3 - ssa4);
  const T ssa8 = e2[1] * t;
  const T ssa9 = s2[1] * ssa0;
  const T ssa10 = e3[1] * t;
  const T ssa11 = s3[1] * ssa0;
  const T ssa12 = e1[1] * t - s1[1] * ssa0;
  const T ssa13 = -ev[1] * t + ssa0 * sv[1] - ssa12 * ssa5 + u * (ssa8 - ssa9) +
                  v * (ssa10 - ssa11);
  const T ssa14 = e2[2] * t;
  const T ssa15 = s2[2] * ssa0;
  const T ssa16 = e3[2] * t;
  const T ssa17 = s3[2] * ssa0;
  const T ssa18 = e1[2] * t - s1[2] * ssa0;
  const T ssa19 = -ev[2] * t + ssa0 * sv[2] - ssa18 * ssa5 +
                  u * (ssa14 - ssa15) + v * (ssa16 - ssa17);
  const T ssa20 = e1[0] - s1[0];
  const T ssa21 =
      -ev[0] - ssa20 * ssa5 + sv[0] + u * (e2[0] - s2[0]) + v * (e3[0] - s3[0]);
  const T ssa22 = e1[1] - s1[1];
  const T ssa23 =
      -ev[1] - ssa22 * ssa5 + sv[1] + u * (e2[1] - s2[1]) + v * (e3[1] - s3[1]);
  const T ssa24 = e1[2] - s1[2];
  const T ssa25 =
      -ev[2] - ssa24 * ssa5 + sv[2] + u * (e2[2] - s2[2]) + v * (e3[2] - s3[2]);
  const T ssa26 = 1.0 / (2 * sccd::pow2<T>(ssa21) + 2 * sccd::pow2<T>(ssa23) +
                         2 * sccd::pow2<T>(ssa25));
  const T ssa27 = -ssa1 + ssa2 + ssa6;
  const T ssa28 = 2 * ssa21;
  const T ssa29 = ssa12 - ssa8 + ssa9;
  const T ssa30 = 2 * ssa23;
  const T ssa31 = -ssa14 + ssa15 + ssa18;
  const T ssa32 = 2 * ssa25;
  const T ssa33 = 2 * ssa7;
  const T ssa34 = 2 * ssa13;
  const T ssa35 = 2 * ssa19;
  const T ssa36 = -ssa27 * ssa28 - ssa29 * ssa30 - ssa31 * ssa32 -
                  ssa33 * (-e2[0] + s2[0] + ssa20) -
                  ssa34 * (-e2[1] + s2[1] + ssa22) -
                  ssa35 * (-e2[2] + s2[2] + ssa24);
  const T ssa37 =
      1.0 / (-ssa26 * sccd::pow2<T>(ssa36) + 2 * sccd::pow2<T>(ssa27) +
             2 * sccd::pow2<T>(ssa29) + 2 * sccd::pow2<T>(ssa31));
  const T ssa38 = -ssa3 + ssa4 + ssa6;
  const T ssa39 = -ssa10 + ssa11 + ssa12;
  const T ssa40 = -ssa16 + ssa17 + ssa18;
  const T ssa41 = -ssa28 * ssa38 - ssa30 * ssa39 - ssa32 * ssa40 -
                  ssa33 * (-e3[0] + s3[0] + ssa20) -
                  ssa34 * (-e3[1] + s3[1] + ssa22) -
                  ssa35 * (-e3[2] + s3[2] + ssa24);
  const T ssa42 = ssa26 * ssa36;
  const T ssa43 =
      2 * ssa27 * ssa38 + 2 * ssa29 * ssa39 + 2 * ssa31 * ssa40 - ssa41 * ssa42;
  const T ssa44 = ssa21 * ssa33 + ssa23 * ssa34 + ssa25 * ssa35;
  const T ssa45 = ssa27 * ssa33 + ssa29 * ssa34 + ssa31 * ssa35 + ssa42 * ssa44;
  const T ssa46 = (-ssa26 * ssa41 * ssa44 - ssa33 * ssa38 - ssa34 * ssa39 -
                   ssa35 * ssa40 + ssa37 * ssa43 * ssa45) /
                  (-ssa26 * sccd::pow2<T>(ssa41) -
                   ssa37 * sccd::pow2<T>(ssa43) + 2 * sccd::pow2<T>(ssa38) +
                   2 * sccd::pow2<T>(ssa39) + 2 * sccd::pow2<T>(ssa40));
  const T ssa47 = ssa37 * (-ssa43 * ssa46 - ssa45);
  *out_f = sccd::pow2<T>(ssa13) + sccd::pow2<T>(ssa19) + sccd::pow2<T>(ssa7);
  out_p[0] = ssa26 * (-ssa36 * ssa47 - ssa41 * ssa46 + ssa44);
  out_p[1] = ssa47;
  out_p[2] = ssa46;
}

template <typename T>
static inline void vf_all(const T sv0[3], const T s1[3], const T s2[3],
                          const T s3[3], const T s4[3], const T ev0[3],
                          const T e1[3], const T e2[3], const T e3[3],
                          const T e4[3], const T t, const T u, const T v,
                          T *out_f, T out_g[3], T out_H[9]) {
  const T ssa0 = t - 1;
  const T ssa1 = e2[0] * t;
  const T ssa2 = s2[0] * ssa0;
  const T ssa3 = e3[0] * t;
  const T ssa4 = s3[0] * ssa0;
  const T ssa5 = u + v - 1;
  const T ssa6 = e1[0] * t - s1[0] * ssa0;
  const T ssa7 = -ev[0] * t + ssa0 * sv[0] - ssa5 * ssa6 + u * (ssa1 - ssa2) +
                 v * (ssa3 - ssa4);
  const T ssa8 = e2[1] * t;
  const T ssa9 = s2[1] * ssa0;
  const T ssa10 = e3[1] * t;
  const T ssa11 = s3[1] * ssa0;
  const T ssa12 = e1[1] * t - s1[1] * ssa0;
  const T ssa13 = -ev[1] * t + ssa0 * sv[1] - ssa12 * ssa5 + u * (ssa8 - ssa9) +
                  v * (ssa10 - ssa11);
  const T ssa14 = e2[2] * t;
  const T ssa15 = s2[2] * ssa0;
  const T ssa16 = e3[2] * t;
  const T ssa17 = s3[2] * ssa0;
  const T ssa18 = e1[2] * t - s1[2] * ssa0;
  const T ssa19 = -ev[2] * t + ssa0 * sv[2] - ssa18 * ssa5 +
                  u * (ssa14 - ssa15) + v * (ssa16 - ssa17);
  const T ssa20 = e1[0] - s1[0];
  const T ssa21 =
      -ev[0] - ssa20 * ssa5 + sv[0] + u * (e2[0] - s2[0]) + v * (e3[0] - s3[0]);
  const T ssa22 = 2 * ssa7;
  const T ssa23 = e1[1] - s1[1];
  const T ssa24 =
      -ev[1] - ssa23 * ssa5 + sv[1] + u * (e2[1] - s2[1]) + v * (e3[1] - s3[1]);
  const T ssa25 = 2 * ssa13;
  const T ssa26 = e1[2] - s1[2];
  const T ssa27 =
      -ev[2] - ssa26 * ssa5 + sv[2] + u * (e2[2] - s2[2]) + v * (e3[2] - s3[2]);
  const T ssa28 = 2 * ssa19;
  const T ssa29 = -ssa1 + ssa2 + ssa6;
  const T ssa30 = ssa12 - ssa8 + ssa9;
  const T ssa31 = -ssa14 + ssa15 + ssa18;
  const T ssa32 = -ssa3 + ssa4 + ssa6;
  const T ssa33 = -ssa10 + ssa11 + ssa12;
  const T ssa34 = -ssa16 + ssa17 + ssa18;
  const T ssa35 = 2 * ssa21;
  const T ssa36 = 2 * ssa24;
  const T ssa37 = 2 * ssa27;
  const T ssa38 = -ssa22 * (-e2[0] + s2[0] + ssa20) -
                  ssa25 * (-e2[1] + s2[1] + ssa23) -
                  ssa28 * (-e2[2] + s2[2] + ssa26) - ssa29 * ssa35 -
                  ssa30 * ssa36 - ssa31 * ssa37;
  const T ssa39 = -ssa22 * (-e3[0] + s3[0] + ssa20) -
                  ssa25 * (-e3[1] + s3[1] + ssa23) -
                  ssa28 * (-e3[2] + s3[2] + ssa26) - ssa32 * ssa35 -
                  ssa33 * ssa36 - ssa34 * ssa37;
  const T ssa40 = 2 * ssa29 * ssa32 + 2 * ssa30 * ssa33 + 2 * ssa31 * ssa34;
  *out_f = sccd::pow2<T>(ssa13) + sccd::pow2<T>(ssa19) + sccd::pow2<T>(ssa7);
  out_g[0] = ssa21 * ssa22 + ssa24 * ssa25 + ssa27 * ssa28;
  out_g[1] = -ssa22 * ssa29 - ssa25 * ssa30 - ssa28 * ssa31;
  out_g[2] = -ssa22 * ssa32 - ssa25 * ssa33 - ssa28 * ssa34;
  out_H[0] = 2 * sccd::pow2<T>(ssa21) + 2 * sccd::pow2<T>(ssa24) +
             2 * sccd::pow2<T>(ssa27);
  out_H[1] = ssa38;
  out_H[2] = ssa39;
  out_H[3] = ssa38;
  out_H[4] = 2 * sccd::pow2<T>(ssa29) + 2 * sccd::pow2<T>(ssa30) +
             2 * sccd::pow2<T>(ssa31);
  out_H[5] = ssa40;
  out_H[6] = ssa39;
  out_H[7] = ssa40;
  out_H[8] = 2 * sccd::pow2<T>(ssa32) + 2 * sccd::pow2<T>(ssa33) +
             2 * sccd::pow2<T>(ssa34);
}
// EE
/* Auto-generated with SymPy (CSE + pow expansion) */

template <typename T>
static inline void ee_objective(const T sv0[3], const T s1[3], const T s2[3],
                                const T s3[3], const T s4[3], const T ev0[3],
                                const T e1[3], const T e2[3], const T e3[3],
                                const T e4[3], const T t, const T u, const T v,
                                T *out_f) {
  const T ssa0 = t - 1;
  const T ssa1 = -s3[0] * ssa0;
  const T ssa2 = v - 1;
  const T ssa3 = u - 1;
  const T ssa4 = -s3[1] * ssa0;
  const T ssa5 = -s3[2] * ssa0;
  *out_f =
      sccd::pow2<T>(ssa2 * (e3[0] * t + ssa1) -
                    ssa3 * (e1[0] * t - s1[0] * ssa0) +
                    u * (e2[0] * t - s2[0] * ssa0) - v * (e4[0] * t + ssa1)) +
      sccd::pow2<T>(ssa2 * (e3[1] * t + ssa4) -
                    ssa3 * (e1[1] * t - s1[1] * ssa0) +
                    u * (e2[1] * t - s2[1] * ssa0) - v * (e4[1] * t + ssa4)) +
      sccd::pow2<T>(ssa2 * (e3[2] * t + ssa5) -
                    ssa3 * (e1[2] * t - s1[2] * ssa0) +
                    u * (e2[2] * t - s2[2] * ssa0) - v * (e4[2] * t + ssa5));
}

template <typename T>
static inline void ee_gradient(const T sv0[3], const T s1[3], const T s2[3],
                               const T s3[3], const T s4[3], const T ev0[3],
                               const T e1[3], const T e2[3], const T e3[3],
                               const T e4[3], const T t, const T u, const T v,
                               T out_g[3]) {
  const T ssa0 = -s3[0];
  const T ssa1 = v - 1;
  const T ssa2 = u - 1;
  const T ssa3 = e2[0] * t;
  const T ssa4 = t - 1;
  const T ssa5 = s2[0] * ssa4;
  const T ssa6 = -s3[0] * ssa4;
  const T ssa7 = e1[0] * t - s1[0] * ssa4;
  const T ssa8 = ssa1 * (e3[0] * t + ssa6) - ssa2 * ssa7 + u * (ssa3 - ssa5) -
                 v * (e4[0] * t + ssa6);
  const T ssa9 = 2 * ssa8;
  const T ssa10 = -s3[1];
  const T ssa11 = e2[1] * t;
  const T ssa12 = s2[1] * ssa4;
  const T ssa13 = -s3[1] * ssa4;
  const T ssa14 = e1[1] * t - s1[1] * ssa4;
  const T ssa15 = ssa1 * (e3[1] * t + ssa13) - ssa14 * ssa2 +
                  u * (ssa11 - ssa12) - v * (e4[1] * t + ssa13);
  const T ssa16 = 2 * ssa15;
  const T ssa17 = -s3[2];
  const T ssa18 = e2[2] * t;
  const T ssa19 = s2[2] * ssa4;
  const T ssa20 = -s3[2] * ssa4;
  const T ssa21 = e1[2] * t - s1[2] * ssa4;
  const T ssa22 = ssa1 * (e3[2] * t + ssa20) - ssa2 * ssa21 +
                  u * (ssa18 - ssa19) - v * (e4[2] * t + ssa20);
  const T ssa23 = 2 * ssa22;
  out_g[0] = ssa16 * (ssa1 * (e3[1] + ssa10) - ssa2 * (e1[1] - s1[1]) +
                      u * (e2[1] - s2[1]) - v * (e4[1] + ssa10)) +
             ssa23 * (ssa1 * (e3[2] + ssa17) - ssa2 * (e1[2] - s1[2]) +
                      u * (e2[2] - s2[2]) - v * (e4[2] + ssa17)) +
             ssa9 * (ssa1 * (e3[0] + ssa0) - ssa2 * (e1[0] - s1[0]) +
                     u * (e2[0] - s2[0]) - v * (e4[0] + ssa0));
  out_g[1] = -ssa16 * (-ssa11 + ssa12 + ssa14) -
             ssa23 * (-ssa18 + ssa19 + ssa21) - ssa9 * (-ssa3 + ssa5 + ssa7);
  out_g[2] = 2 * t *
             (ssa15 * (e3[1] - e4[1]) + ssa22 * (e3[2] - e4[2]) +
              ssa8 * (e3[0] - e4[0]));
}

template <typename T>
static inline void ee_hessian(const T sv0[3], const T s1[3], const T s2[3],
                              const T s3[3], const T s4[3], const T ev0[3],
                              const T e1[3], const T e2[3], const T e3[3],
                              const T e4[3], const T t, const T u, const T v,
                              T out_H[9]) {
  const T ssa0 = -s3[0];
  const T ssa1 = v - 1;
  const T ssa2 = u - 1;
  const T ssa3 = e1[0] - s1[0];
  const T ssa4 = ssa1 * (e3[0] + ssa0) - ssa2 * ssa3 + u * (e2[0] - s2[0]) -
                 v * (e4[0] + ssa0);
  const T ssa5 = -s3[1];
  const T ssa6 = e1[1] - s1[1];
  const T ssa7 = ssa1 * (e3[1] + ssa5) - ssa2 * ssa6 + u * (e2[1] - s2[1]) -
                 v * (e4[1] + ssa5);
  const T ssa8 = -s3[2];
  const T ssa9 = e1[2] - s1[2];
  const T ssa10 = ssa1 * (e3[2] + ssa8) - ssa2 * ssa9 + u * (e2[2] - s2[2]) -
                  v * (e4[2] + ssa8);
  const T ssa11 = t - 1;
  const T ssa12 = s2[0] * ssa11;
  const T ssa13 = e2[0] * t;
  const T ssa14 = e1[0] * t - s1[0] * ssa11;
  const T ssa15 = ssa12 - ssa13 + ssa14;
  const T ssa16 = s2[1] * ssa11;
  const T ssa17 = e2[1] * t;
  const T ssa18 = e1[1] * t - s1[1] * ssa11;
  const T ssa19 = ssa16 - ssa17 + ssa18;
  const T ssa20 = s2[2] * ssa11;
  const T ssa21 = e2[2] * t;
  const T ssa22 = e1[2] * t - s1[2] * ssa11;
  const T ssa23 = ssa20 - ssa21 + ssa22;
  const T ssa24 = -s3[0] * ssa11;
  const T ssa25 = 2 * ssa1 * (e3[0] * t + ssa24) - 2 * ssa14 * ssa2 +
                  2 * u * (-ssa12 + ssa13) - 2 * v * (e4[0] * t + ssa24);
  const T ssa26 = -s3[1] * ssa11;
  const T ssa27 = 2 * ssa1 * (e3[1] * t + ssa26) - 2 * ssa18 * ssa2 +
                  2 * u * (-ssa16 + ssa17) - 2 * v * (e4[1] * t + ssa26);
  const T ssa28 = -s3[2] * ssa11;
  const T ssa29 = 2 * ssa1 * (e3[2] * t + ssa28) - 2 * ssa2 * ssa22 +
                  2 * u * (-ssa20 + ssa21) - 2 * v * (e4[2] * t + ssa28);
  const T ssa30 = -2 * ssa10 * ssa23 - 2 * ssa15 * ssa4 - 2 * ssa19 * ssa7 -
                  ssa25 * (-e2[0] + s2[0] + ssa3) -
                  ssa27 * (-e2[1] + s2[1] + ssa6) -
                  ssa29 * (-e2[2] + s2[2] + ssa9);
  const T ssa31 = e3[0] - e4[0];
  const T ssa32 = sccd::pow2<T>(t);
  const T ssa33 = e3[1] - e4[1];
  const T ssa34 = e3[2] - e4[2];
  const T ssa35 = ssa10 * ssa32 * ssa34 + ssa25 * ssa31 + ssa27 * ssa33 +
                  ssa29 * ssa34 + ssa31 * ssa32 * ssa4 + ssa32 * ssa33 * ssa7;
  const T ssa36 = ssa32 * (-ssa15 * ssa31 - ssa19 * ssa33 - ssa23 * ssa34);
  out_H[0] = 2 * sccd::pow2<T>(ssa10) + 2 * sccd::pow2<T>(ssa4) +
             2 * sccd::pow2<T>(ssa7);
  out_H[1] = ssa30;
  out_H[2] = ssa35;
  out_H[3] = ssa30;
  out_H[4] = 2 * sccd::pow2<T>(ssa15) + 2 * sccd::pow2<T>(ssa19) +
             2 * sccd::pow2<T>(ssa23);
  out_H[5] = ssa36;
  out_H[6] = ssa35;
  out_H[7] = ssa36;
  out_H[8] = (2 * sccd::pow2<T>(ssa31) + 2 * sccd::pow2<T>(ssa33) +
              2 * sccd::pow2<T>(ssa34)) *
             (t * t);
}

template <typename T>
static inline void
ee_objective_dir(const T sv0[3], const T s1[3], const T s2[3], const T s3[3],
                 const T s4[3], const T ev0[3], const T e1[3], const T e2[3],
                 const T e3[3], const T e4[3], const T t, const T u, const T v,
                 T *out_f, T out_p[3]) {
  const T ssa0 = e2[0] * t;
  const T ssa1 = t - 1;
  const T ssa2 = s2[0] * ssa1;
  const T ssa3 = -s3[0] * ssa1;
  const T ssa4 = v - 1;
  const T ssa5 = u - 1;
  const T ssa6 = e1[0] * t - s1[0] * ssa1;
  const T ssa7 = ssa4 * (e3[0] * t + ssa3) - ssa5 * ssa6 + u * (ssa0 - ssa2) -
                 v * (e4[0] * t + ssa3);
  const T ssa8 = e2[1] * t;
  const T ssa9 = s2[1] * ssa1;
  const T ssa10 = -s3[1] * ssa1;
  const T ssa11 = e1[1] * t - s1[1] * ssa1;
  const T ssa12 = -ssa11 * ssa5 + ssa4 * (e3[1] * t + ssa10) +
                  u * (ssa8 - ssa9) - v * (e4[1] * t + ssa10);
  const T ssa13 = e2[2] * t;
  const T ssa14 = s2[2] * ssa1;
  const T ssa15 = -s3[2] * ssa1;
  const T ssa16 = e1[2] * t - s1[2] * ssa1;
  const T ssa17 = -ssa16 * ssa5 + ssa4 * (e3[2] * t + ssa15) +
                  u * (ssa13 - ssa14) - v * (e4[2] * t + ssa15);
  const T ssa18 = -s3[0];
  const T ssa19 = e1[0] - s1[0];
  const T ssa20 = -ssa19 * ssa5 + ssa4 * (e3[0] + ssa18) + u * (e2[0] - s2[0]) -
                  v * (e4[0] + ssa18);
  const T ssa21 = -s3[1];
  const T ssa22 = e1[1] - s1[1];
  const T ssa23 = -ssa22 * ssa5 + ssa4 * (e3[1] + ssa21) + u * (e2[1] - s2[1]) -
                  v * (e4[1] + ssa21);
  const T ssa24 = -s3[2];
  const T ssa25 = e1[2] - s1[2];
  const T ssa26 = -ssa25 * ssa5 + ssa4 * (e3[2] + ssa24) + u * (e2[2] - s2[2]) -
                  v * (e4[2] + ssa24);
  const T ssa27 = 1.0 / (2 * sccd::pow2<T>(ssa20) + 2 * sccd::pow2<T>(ssa23) +
                         2 * sccd::pow2<T>(ssa26));
  const T ssa28 = -ssa0 + ssa2 + ssa6;
  const T ssa29 = ssa11 - ssa8 + ssa9;
  const T ssa30 = -ssa13 + ssa14 + ssa16;
  const T ssa31 = 2 * ssa7;
  const T ssa32 = 2 * ssa12;
  const T ssa33 = 2 * ssa17;
  const T ssa34 = -2 * ssa20 * ssa28 - 2 * ssa23 * ssa29 - 2 * ssa26 * ssa30 -
                  ssa31 * (-e2[0] + s2[0] + ssa19) -
                  ssa32 * (-e2[1] + s2[1] + ssa22) -
                  ssa33 * (-e2[2] + s2[2] + ssa25);
  const T ssa35 =
      1.0 / (-ssa27 * sccd::pow2<T>(ssa34) + 2 * sccd::pow2<T>(ssa28) +
             2 * sccd::pow2<T>(ssa29) + 2 * sccd::pow2<T>(ssa30));
  const T ssa36 = e3[0] - e4[0];
  const T ssa37 = e3[1] - e4[1];
  const T ssa38 = e3[2] - e4[2];
  const T ssa39 = sccd::pow2<T>(t);
  const T ssa40 = ssa36 * ssa7;
  const T ssa41 = ssa12 * ssa37;
  const T ssa42 = ssa17 * ssa38;
  const T ssa43 = ssa20 * ssa36 * ssa39 + ssa23 * ssa37 * ssa39 +
                  ssa26 * ssa38 * ssa39 + 2 * ssa40 + 2 * ssa41 + 2 * ssa42;
  const T ssa44 = ssa27 * ssa34;
  const T ssa45 =
      ssa39 * (-ssa28 * ssa36 - ssa29 * ssa37 - ssa30 * ssa38) - ssa43 * ssa44;
  const T ssa46 = ssa20 * ssa31 + ssa23 * ssa32 + ssa26 * ssa33;
  const T ssa47 = ssa28 * ssa31 + ssa29 * ssa32 + ssa30 * ssa33 + ssa44 * ssa46;
  const T ssa48 = -ssa27 * ssa43 * ssa46 + ssa35 * ssa45 * ssa47 +
                  2 * t * (ssa40 + ssa41 + ssa42);
  const T ssa49 =
      -1.0 / (ssa27 * sccd::pow2<T>(ssa43) + ssa35 * sccd::pow2<T>(ssa45) -
              2 *
                  (sccd::pow2<T>(ssa36) + sccd::pow2<T>(ssa37) +
                   sccd::pow2<T>(ssa38)) *
                  t * t);
  const T ssa50 = -ssa45 * ssa48 * ssa49 - ssa47;
  *out_f = sccd::pow2<T>(ssa12) + sccd::pow2<T>(ssa17) + sccd::pow2<T>(ssa7);
  out_p[0] = ssa27 * (-ssa34 * ssa35 * ssa50 - ssa43 * ssa48 * ssa49 + ssa46);
  out_p[1] = ssa35 * ssa50;
  out_p[2] = ssa48 * ssa49;
}

template <typename T>
static inline void ee_all(const T sv0[3], const T s1[3], const T s2[3],
                          const T s3[3], const T s4[3], const T ev0[3],
                          const T e1[3], const T e2[3], const T e3[3],
                          const T e4[3], const T t, const T u, const T v,
                          T *out_f, T out_g[3], T out_H[9]) {
  const T ssa0 = e2[0] * t;
  const T ssa1 = t - 1;
  const T ssa2 = s2[0] * ssa1;
  const T ssa3 = -s3[0] * ssa1;
  const T ssa4 = v - 1;
  const T ssa5 = u - 1;
  const T ssa6 = e1[0] * t - s1[0] * ssa1;
  const T ssa7 = ssa4 * (e3[0] * t + ssa3) - ssa5 * ssa6 + u * (ssa0 - ssa2) -
                 v * (e4[0] * t + ssa3);
  const T ssa8 = e2[1] * t;
  const T ssa9 = s2[1] * ssa1;
  const T ssa10 = -s3[1] * ssa1;
  const T ssa11 = e1[1] * t - s1[1] * ssa1;
  const T ssa12 = -ssa11 * ssa5 + ssa4 * (e3[1] * t + ssa10) +
                  u * (ssa8 - ssa9) - v * (e4[1] * t + ssa10);
  const T ssa13 = e2[2] * t;
  const T ssa14 = s2[2] * ssa1;
  const T ssa15 = -s3[2] * ssa1;
  const T ssa16 = e1[2] * t - s1[2] * ssa1;
  const T ssa17 = -ssa16 * ssa5 + ssa4 * (e3[2] * t + ssa15) +
                  u * (ssa13 - ssa14) - v * (e4[2] * t + ssa15);
  const T ssa18 = -s3[0];
  const T ssa19 = e1[0] - s1[0];
  const T ssa20 = -ssa19 * ssa5 + ssa4 * (e3[0] + ssa18) + u * (e2[0] - s2[0]) -
                  v * (e4[0] + ssa18);
  const T ssa21 = 2 * ssa7;
  const T ssa22 = -s3[1];
  const T ssa23 = e1[1] - s1[1];
  const T ssa24 = -ssa23 * ssa5 + ssa4 * (e3[1] + ssa22) + u * (e2[1] - s2[1]) -
                  v * (e4[1] + ssa22);
  const T ssa25 = 2 * ssa12;
  const T ssa26 = -s3[2];
  const T ssa27 = e1[2] - s1[2];
  const T ssa28 = -ssa27 * ssa5 + ssa4 * (e3[2] + ssa26) + u * (e2[2] - s2[2]) -
                  v * (e4[2] + ssa26);
  const T ssa29 = 2 * ssa17;
  const T ssa30 = -ssa0 + ssa2 + ssa6;
  const T ssa31 = ssa11 - ssa8 + ssa9;
  const T ssa32 = -ssa13 + ssa14 + ssa16;
  const T ssa33 = e3[0] - e4[0];
  const T ssa34 = ssa33 * ssa7;
  const T ssa35 = e3[1] - e4[1];
  const T ssa36 = ssa12 * ssa35;
  const T ssa37 = e3[2] - e4[2];
  const T ssa38 = ssa17 * ssa37;
  const T ssa39 = sccd::pow2<T>(t);
  const T ssa40 = -2 * ssa20 * ssa30 - ssa21 * (-e2[0] + s2[0] + ssa19) -
                  2 * ssa24 * ssa31 - ssa25 * (-e2[1] + s2[1] + ssa23) -
                  2 * ssa28 * ssa32 - ssa29 * (-e2[2] + s2[2] + ssa27);
  const T ssa41 = ssa20 * ssa33 * ssa39 + ssa24 * ssa35 * ssa39 +
                  ssa28 * ssa37 * ssa39 + 2 * ssa34 + 2 * ssa36 + 2 * ssa38;
  const T ssa42 = ssa39 * (-ssa30 * ssa33 - ssa31 * ssa35 - ssa32 * ssa37);
  *out_f = sccd::pow2<T>(ssa12) + sccd::pow2<T>(ssa17) + sccd::pow2<T>(ssa7);
  out_g[0] = ssa20 * ssa21 + ssa24 * ssa25 + ssa28 * ssa29;
  out_g[1] = -ssa21 * ssa30 - ssa25 * ssa31 - ssa29 * ssa32;
  out_g[2] = ssa39 * (ssa34 + ssa36 + ssa38);
  out_H[0] = 2 * sccd::pow2<T>(ssa20) + 2 * sccd::pow2<T>(ssa24) +
             2 * sccd::pow2<T>(ssa28);
  out_H[1] = ssa40;
  out_H[2] = ssa41;
  out_H[3] = ssa40;
  out_H[4] = 2 * sccd::pow2<T>(ssa30) + 2 * sccd::pow2<T>(ssa31) +
             2 * sccd::pow2<T>(ssa32);
  out_H[5] = ssa42;
  out_H[6] = ssa41;
  out_H[7] = ssa42;
  out_H[8] = (2 * sccd::pow2<T>(ssa33) + 2 * sccd::pow2<T>(ssa35) +
              2 * sccd::pow2<T>(ssa37)) *
             (t * t);
}
