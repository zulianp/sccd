template <typename T>
static inline void sccd_calculate_vf_soa(
    const T t, const T u, const T v, const T v0sx, const T v0sy, const T v0sz,
    const T v1sx, const T v1sy, const T v1sz, const T v2sx, const T v2sy,
    const T v2sz, const T v3sx, const T v3sy, const T v3sz, const T v0ex,
    const T v0ey, const T v0ez, const T v1ex, const T v1ey, const T v1ez,
    const T v2ex, const T v2ey, const T v2ez, const T v3ex, const T v3ey,
    const T v3ez, T *const SFEM_RESTRICT rx, T *const SFEM_RESTRICT ry,
    T *const SFEM_RESTRICT rz) {
  const T ssa0 = -v1sx;
  const T ssa1 = t * (ssa0 + v1ex);
  const T ssa2 = -v2sx;
  const T ssa3 = ssa1 + v1sx;
  const T ssa4 = -v3sx;
  const T ssa5 = -v1sy;
  const T ssa6 = t * (ssa5 + v1ey);
  const T ssa7 = -v2sy;
  const T ssa8 = ssa6 + v1sy;
  const T ssa9 = -v3sy;
  const T ssa10 = -v1sz;
  const T ssa11 = t * (ssa10 + v1ez);
  const T ssa12 = -v2sz;
  const T ssa13 = ssa11 + v1sz;
  const T ssa14 = -v3sz;
  *rx = ssa0 - ssa1 + t * (v0ex - v0sx) +
        u * (ssa2 + ssa3 - t * (ssa2 + v2ex)) +
        v * (ssa3 + ssa4 - t * (ssa4 + v3ex)) + v0sx;
  *ry = ssa5 - ssa6 + t * (v0ey - v0sy) +
        u * (ssa7 + ssa8 - t * (ssa7 + v2ey)) +
        v * (ssa8 + ssa9 - t * (ssa9 + v3ey)) + v0sy;
  *rz = ssa10 - ssa11 + t * (v0ez - v0sz) +
        u * (ssa12 + ssa13 - t * (ssa12 + v2ez)) +
        v * (ssa13 + ssa14 - t * (ssa14 + v3ez)) + v0sz;
}

template <typename T>
static inline void sccd_calculate_ee_soa(
    const T t, const T u, const T v, const T v0sx, const T v0sy, const T v0sz,
    const T v1sx, const T v1sy, const T v1sz, const T v2sx, const T v2sy,
    const T v2sz, const T v3sx, const T v3sy, const T v3sz, const T v0ex,
    const T v0ey, const T v0ez, const T v1ex, const T v1ey, const T v1ez,
    const T v2ex, const T v2ey, const T v2ez, const T v3ex, const T v3ey,
    const T v3ez, T *const SFEM_RESTRICT rx, T *const SFEM_RESTRICT ry,
    T *const SFEM_RESTRICT rz) {
  const T ssa0 = -v2sx;
  const T ssa1 = -v3sx;
  const T ssa2 = t * (ssa0 + v2ex);
  const T ssa3 = -v1sx;
  const T ssa4 = t * (v0ex - v0sx) + v0sx;
  const T ssa5 = -v2sy;
  const T ssa6 = -v3sy;
  const T ssa7 = t * (ssa5 + v2ey);
  const T ssa8 = -v1sy;
  const T ssa9 = t * (v0ey - v0sy) + v0sy;
  const T ssa10 = -v2sz;
  const T ssa11 = -v3sz;
  const T ssa12 = t * (ssa10 + v2ez);
  const T ssa13 = -v1sz;
  const T ssa14 = t * (v0ez - v0sz) + v0sz;
  *rx = ssa0 - ssa2 + ssa4 - u * (ssa3 + ssa4 - t * (ssa3 + v1ex)) +
        v * (ssa1 + ssa2 - t * (ssa1 + v3ex) + v2sx);
  *ry = ssa5 - ssa7 + ssa9 - u * (ssa8 + ssa9 - t * (ssa8 + v1ey)) +
        v * (ssa6 + ssa7 - t * (ssa6 + v3ey) + v2sy);
  *rz = ssa10 - ssa12 + ssa14 - u * (ssa13 + ssa14 - t * (ssa13 + v1ez)) +
        v * (ssa11 + ssa12 - t * (ssa11 + v3ez) + v2sz);
}

template <typename T>
static inline int sccd_origin_in_inclusion_vf_soa(
    const T t_l, const T t_u, const T u_l, const T u_u, const T v_l,
    const T v_u, const T v0sx, const T v0sy, const T v0sz, const T v1sx,
    const T v1sy, const T v1sz, const T v2sx, const T v2sy, const T v2sz,
    const T v3sx, const T v3sy, const T v3sz, const T v0ex, const T v0ey,
    const T v0ez, const T v1ex, const T v1ey, const T v1ez, const T v2ex,
    const T v2ey, const T v2ez, const T v3ex, const T v3ey, const T v3ez,
    const T msx, const T msy, const T msz, const T errx, const T erry,
    const T errz, T *const SFEM_RESTRICT true_tol,
    int *const SFEM_RESTRICT box_in) {
  const T ssa0 = -v3sx;
  const T ssa1 = ssa0 + v3ex;
  const T ssa2 = -v1sx;
  const T ssa3 = ssa2 + v1ex;
  const T ssa4 = ssa3 * t_l;
  const T ssa5 = ssa4 + v1sx;
  const T ssa6 = ssa0 - ssa1 * t_l + ssa5;
  const T ssa7 = ssa6 * v_l;
  const T ssa8 = -v2sx;
  const T ssa9 = ssa8 + v2ex;
  const T ssa10 = ssa5 + ssa8 - ssa9 * t_l;
  const T ssa11 = -v0sx;
  const T ssa12 = ssa11 + v0ex;
  const T ssa13 = ssa2 + v0sx;
  const T ssa14 = ssa12 * t_l + ssa13 - ssa4;
  const T ssa15 = ssa10 * u_l + ssa14;
  const T ssa16 = ssa15 + ssa7;
  const T ssa17 = ssa6 * v_u;
  const T ssa18 = ssa15 + ssa17;
  const T ssa19 = ssa10 * u_u + ssa14;
  const T ssa20 = ssa19 + ssa7;
  const T ssa21 = ssa17 + ssa19;
  const T ssa22 = ssa3 * t_u;
  const T ssa23 = ssa22 + v1sx;
  const T ssa24 = ssa0 - ssa1 * t_u + ssa23;
  const T ssa25 = ssa24 * v_l;
  const T ssa26 = ssa23 + ssa8 - ssa9 * t_u;
  const T ssa27 = ssa12 * t_u + ssa13 - ssa22;
  const T ssa28 = ssa26 * u_l + ssa27;
  const T ssa29 = ssa25 + ssa28;
  const T ssa30 = ssa24 * v_u;
  const T ssa31 = ssa28 + ssa30;
  const T ssa32 = ssa26 * u_u + ssa27;
  const T ssa33 = ssa25 + ssa32;
  const T ssa34 = ssa30 + ssa32;
  const T ssa35 = sccd::max<T>(
      ssa16,
      sccd::max<T>(
          ssa18,
          sccd::max<T>(
              ssa20,
              sccd::max<T>(
                  ssa21, sccd::max<T>(
                             ssa29, sccd::max<T>(
                                        ssa31, sccd::max<T>(ssa33, ssa34)))))));
  const T ssa36 = sccd::min<T>(
      ssa16,
      sccd::min<T>(
          ssa18,
          sccd::min<T>(
              ssa20,
              sccd::min<T>(
                  ssa21, sccd::min<T>(
                             ssa29, sccd::min<T>(
                                        ssa31, sccd::min<T>(ssa33, ssa34)))))));
  const T ssa37 = -v3sy;
  const T ssa38 = ssa37 + v3ey;
  const T ssa39 = -v1sy;
  const T ssa40 = ssa39 + v1ey;
  const T ssa41 = ssa40 * t_l;
  const T ssa42 = ssa41 + v1sy;
  const T ssa43 = ssa37 - ssa38 * t_l + ssa42;
  const T ssa44 = ssa43 * v_l;
  const T ssa45 = -v2sy;
  const T ssa46 = ssa45 + v2ey;
  const T ssa47 = ssa42 + ssa45 - ssa46 * t_l;
  const T ssa48 = -v0sy;
  const T ssa49 = ssa48 + v0ey;
  const T ssa50 = ssa39 + v0sy;
  const T ssa51 = -ssa41 + ssa49 * t_l + ssa50;
  const T ssa52 = ssa47 * u_l + ssa51;
  const T ssa53 = ssa44 + ssa52;
  const T ssa54 = ssa43 * v_u;
  const T ssa55 = ssa52 + ssa54;
  const T ssa56 = ssa47 * u_u + ssa51;
  const T ssa57 = ssa44 + ssa56;
  const T ssa58 = ssa54 + ssa56;
  const T ssa59 = ssa40 * t_u;
  const T ssa60 = ssa59 + v1sy;
  const T ssa61 = ssa37 - ssa38 * t_u + ssa60;
  const T ssa62 = ssa61 * v_l;
  const T ssa63 = ssa45 - ssa46 * t_u + ssa60;
  const T ssa64 = ssa49 * t_u + ssa50 - ssa59;
  const T ssa65 = ssa63 * u_l + ssa64;
  const T ssa66 = ssa62 + ssa65;
  const T ssa67 = ssa61 * v_u;
  const T ssa68 = ssa65 + ssa67;
  const T ssa69 = ssa63 * u_u + ssa64;
  const T ssa70 = ssa62 + ssa69;
  const T ssa71 = ssa67 + ssa69;
  const T ssa72 = sccd::max<T>(
      ssa53,
      sccd::max<T>(
          ssa55,
          sccd::max<T>(
              ssa57,
              sccd::max<T>(
                  ssa58, sccd::max<T>(
                             ssa66, sccd::max<T>(
                                        ssa68, sccd::max<T>(ssa70, ssa71)))))));
  const T ssa73 = sccd::min<T>(
      ssa53,
      sccd::min<T>(
          ssa55,
          sccd::min<T>(
              ssa57,
              sccd::min<T>(
                  ssa58, sccd::min<T>(
                             ssa66, sccd::min<T>(
                                        ssa68, sccd::min<T>(ssa70, ssa71)))))));
  const T ssa74 = -v3sz;
  const T ssa75 = ssa74 + v3ez;
  const T ssa76 = -v1sz;
  const T ssa77 = ssa76 + v1ez;
  const T ssa78 = ssa77 * t_l;
  const T ssa79 = ssa78 + v1sz;
  const T ssa80 = ssa74 - ssa75 * t_l + ssa79;
  const T ssa81 = ssa80 * v_l;
  const T ssa82 = -v2sz;
  const T ssa83 = ssa82 + v2ez;
  const T ssa84 = ssa79 + ssa82 - ssa83 * t_l;
  const T ssa85 = -v0sz;
  const T ssa86 = ssa85 + v0ez;
  const T ssa87 = ssa76 + v0sz;
  const T ssa88 = -ssa78 + ssa86 * t_l + ssa87;
  const T ssa89 = ssa84 * u_l + ssa88;
  const T ssa90 = ssa81 + ssa89;
  const T ssa91 = ssa80 * v_u;
  const T ssa92 = ssa89 + ssa91;
  const T ssa93 = ssa84 * u_u + ssa88;
  const T ssa94 = ssa81 + ssa93;
  const T ssa95 = ssa91 + ssa93;
  const T ssa96 = ssa77 * t_u;
  const T ssa97 = ssa96 + v1sz;
  const T ssa98 = ssa74 - ssa75 * t_u + ssa97;
  const T ssa99 = ssa98 * v_l;
  const T ssa100 = ssa82 - ssa83 * t_u + ssa97;
  const T ssa101 = ssa86 * t_u + ssa87 - ssa96;
  const T ssa102 = ssa100 * u_l + ssa101;
  const T ssa103 = ssa102 + ssa99;
  const T ssa104 = ssa98 * v_u;
  const T ssa105 = ssa102 + ssa104;
  const T ssa106 = ssa100 * u_u + ssa101;
  const T ssa107 = ssa106 + ssa99;
  const T ssa108 = ssa104 + ssa106;
  const T ssa109 = sccd::max<T>(
      ssa103,
      sccd::max<T>(
          ssa105,
          sccd::max<T>(
              ssa107,
              sccd::max<T>(
                  ssa108,
                  sccd::max<T>(
                      ssa90,
                      sccd::max<T>(ssa92, sccd::max<T>(ssa94, ssa95)))))));
  const T ssa110 = sccd::min<T>(
      ssa103,
      sccd::min<T>(
          ssa105,
          sccd::min<T>(
              ssa107,
              sccd::min<T>(
                  ssa108,
                  sccd::min<T>(
                      ssa90,
                      sccd::min<T>(ssa92, sccd::min<T>(ssa94, ssa95)))))));
  const T ssa111 = -errx;
  const T ssa112 = -erry;
  const T ssa113 = -errz;
  const T ssa114 = -ssa6;
  const T ssa115 = ssa114 * v_l;
  const T ssa116 = -ssa10;
  const T ssa117 = ssa11 - ssa12 * t_l + ssa5;
  const T ssa118 = ssa116 * u_l + ssa117;
  const T ssa119 = -ssa115 - ssa118;
  const T ssa120 = ssa114 * v_u;
  const T ssa121 = -ssa118 - ssa120;
  const T ssa122 = ssa116 * u_u + ssa117;
  const T ssa123 = -ssa115 - ssa122;
  const T ssa124 = -ssa120 - ssa122;
  const T ssa125 = -ssa24;
  const T ssa126 = ssa125 * v_l;
  const T ssa127 = -ssa26;
  const T ssa128 = ssa11 - ssa12 * t_u + ssa23;
  const T ssa129 = ssa127 * u_l + ssa128;
  const T ssa130 = -ssa126 - ssa129;
  const T ssa131 = ssa125 * v_u;
  const T ssa132 = -ssa129 - ssa131;
  const T ssa133 = ssa127 * u_u + ssa128;
  const T ssa134 = -ssa126 - ssa133;
  const T ssa135 = -ssa131 - ssa133;
  const T ssa136 = -ssa43;
  const T ssa137 = ssa136 * v_l;
  const T ssa138 = -ssa47;
  const T ssa139 = ssa42 + ssa48 - ssa49 * t_l;
  const T ssa140 = ssa138 * u_l + ssa139;
  const T ssa141 = -ssa137 - ssa140;
  const T ssa142 = ssa136 * v_u;
  const T ssa143 = -ssa140 - ssa142;
  const T ssa144 = ssa138 * u_u + ssa139;
  const T ssa145 = -ssa137 - ssa144;
  const T ssa146 = -ssa142 - ssa144;
  const T ssa147 = -ssa61;
  const T ssa148 = ssa147 * v_l;
  const T ssa149 = -ssa63;
  const T ssa150 = ssa48 - ssa49 * t_u + ssa60;
  const T ssa151 = ssa149 * u_l + ssa150;
  const T ssa152 = -ssa148 - ssa151;
  const T ssa153 = ssa147 * v_u;
  const T ssa154 = -ssa151 - ssa153;
  const T ssa155 = ssa149 * u_u + ssa150;
  const T ssa156 = -ssa148 - ssa155;
  const T ssa157 = -ssa153 - ssa155;
  const T ssa158 = -ssa80;
  const T ssa159 = ssa158 * v_l;
  const T ssa160 = -ssa84;
  const T ssa161 = ssa79 + ssa85 - ssa86 * t_l;
  const T ssa162 = ssa160 * u_l + ssa161;
  const T ssa163 = -ssa159 - ssa162;
  const T ssa164 = ssa158 * v_u;
  const T ssa165 = -ssa162 - ssa164;
  const T ssa166 = ssa160 * u_u + ssa161;
  const T ssa167 = -ssa159 - ssa166;
  const T ssa168 = -ssa164 - ssa166;
  const T ssa169 = -ssa98;
  const T ssa170 = ssa169 * v_l;
  const T ssa171 = -ssa100;
  const T ssa172 = ssa85 - ssa86 * t_u + ssa97;
  const T ssa173 = ssa171 * u_l + ssa172;
  const T ssa174 = -ssa170 - ssa173;
  const T ssa175 = ssa169 * v_u;
  const T ssa176 = -ssa173 - ssa175;
  const T ssa177 = ssa171 * u_u + ssa172;
  const T ssa178 = -ssa170 - ssa177;
  const T ssa179 = -ssa175 - ssa177;
  *true_tol =
      sccd::max<T>(0, sccd::max<T>(ssa109 - ssa110,
                                   sccd::max<T>(ssa35 - ssa36, ssa72 - ssa73)));
  if (-msx + sccd::min<T>(
                 ssa119,
                 sccd::min<T>(
                     ssa121,
                     sccd::min<T>(
                         ssa123,
                         sccd::min<T>(
                             ssa124,
                             sccd::min<T>(
                                 ssa130, sccd::min<T>(
                                             ssa132, sccd::min<T>(
                                                         ssa134, ssa135))))))) >
          errx ||
      -msy + sccd::min<T>(
                 ssa141,
                 sccd::min<T>(
                     ssa143,
                     sccd::min<T>(
                         ssa145,
                         sccd::min<T>(
                             ssa146,
                             sccd::min<T>(
                                 ssa152, sccd::min<T>(
                                             ssa154, sccd::min<T>(
                                                         ssa156, ssa157))))))) >
          erry ||
      -msz + sccd::min<T>(
                 ssa163,
                 sccd::min<T>(
                     ssa165,
                     sccd::min<T>(
                         ssa167,
                         sccd::min<T>(
                             ssa168,
                             sccd::min<T>(
                                 ssa174, sccd::min<T>(
                                             ssa176, sccd::min<T>(
                                                         ssa178, ssa179))))))) >
          errz ||
      msx + ssa35 < ssa111 || msy + ssa72 < ssa112 || msz + ssa109 < ssa113) {
    *box_in = 0;
    return 0;
  }
  *box_in =
      ((ssa111 > msx + ssa36 || ssa112 > msy + ssa73 || ssa113 > msz + ssa110 ||
        errx <
            -msx + sccd::max<T>(
                       ssa119,
                       sccd::max<T>(
                           ssa121,
                           sccd::max<T>(
                               ssa123,
                               sccd::max<T>(
                                   ssa124,
                                   sccd::max<T>(
                                       ssa130,
                                       sccd::max<T>(
                                           ssa132,
                                           sccd::max<T>(ssa134, ssa135))))))) ||
        erry <
            -msy + sccd::max<T>(
                       ssa141,
                       sccd::max<T>(
                           ssa143,
                           sccd::max<T>(
                               ssa145,
                               sccd::max<T>(
                                   ssa146,
                                   sccd::max<T>(
                                       ssa152,
                                       sccd::max<T>(
                                           ssa154,
                                           sccd::max<T>(ssa156, ssa157))))))) ||
        errz < -msz + sccd::max<T>(
                          ssa163,
                          sccd::max<T>(
                              ssa165,
                              sccd::max<T>(
                                  ssa167,
                                  sccd::max<T>(
                                      ssa168,
                                      sccd::max<T>(
                                          ssa174,
                                          sccd::max<T>(
                                              ssa176, sccd::max<T>(
                                                          ssa178, ssa179))))))))
           ? (0)
           : (1));
  return 1;
}

template <typename T>
static inline int sccd_origin_in_inclusion_ee_soa(
    const T t_l, const T t_u, const T u_l, const T u_u, const T v_l,
    const T v_u, const T v0sx, const T v0sy, const T v0sz, const T v1sx,
    const T v1sy, const T v1sz, const T v2sx, const T v2sy, const T v2sz,
    const T v3sx, const T v3sy, const T v3sz, const T v0ex, const T v0ey,
    const T v0ez, const T v1ex, const T v1ey, const T v1ez, const T v2ex,
    const T v2ey, const T v2ez, const T v3ex, const T v3ey, const T v3ez,
    const T msx, const T msy, const T msz, const T errx, const T erry,
    const T errz, T *const SFEM_RESTRICT true_tol,
    int *const SFEM_RESTRICT box_in) {
  const T ssa0 = -v2sx;
  const T ssa1 = ssa0 + v2ex;
  const T ssa2 = ssa1 * t_l;
  const T ssa3 = -v3sx;
  const T ssa4 = ssa3 + v3ex;
  const T ssa5 = ssa3 + v2sx;
  const T ssa6 = ssa2 - ssa4 * t_l + ssa5;
  const T ssa7 = ssa6 * v_l;
  const T ssa8 = v0ex - v0sx;
  const T ssa9 = ssa8 * t_l;
  const T ssa10 = -v1sx;
  const T ssa11 = ssa10 + v1ex;
  const T ssa12 = ssa10 + v0sx;
  const T ssa13 = -ssa11 * t_l + ssa12 + ssa9;
  const T ssa14 = ssa0 + v0sx;
  const T ssa15 = ssa14 - ssa2 + ssa9;
  const T ssa16 = -ssa13 * u_l + ssa15;
  const T ssa17 = ssa16 + ssa7;
  const T ssa18 = ssa6 * v_u;
  const T ssa19 = ssa16 + ssa18;
  const T ssa20 = -ssa13 * u_u + ssa15;
  const T ssa21 = ssa20 + ssa7;
  const T ssa22 = ssa18 + ssa20;
  const T ssa23 = ssa1 * t_u;
  const T ssa24 = ssa23 - ssa4 * t_u + ssa5;
  const T ssa25 = ssa24 * v_l;
  const T ssa26 = ssa8 * t_u;
  const T ssa27 = -ssa11 * t_u + ssa12 + ssa26;
  const T ssa28 = ssa14 - ssa23 + ssa26;
  const T ssa29 = -ssa27 * u_l + ssa28;
  const T ssa30 = ssa25 + ssa29;
  const T ssa31 = ssa24 * v_u;
  const T ssa32 = ssa29 + ssa31;
  const T ssa33 = -ssa27 * u_u + ssa28;
  const T ssa34 = ssa25 + ssa33;
  const T ssa35 = ssa31 + ssa33;
  const T ssa36 = sccd::max<T>(
      ssa17,
      sccd::max<T>(
          ssa19,
          sccd::max<T>(
              ssa21,
              sccd::max<T>(
                  ssa22, sccd::max<T>(
                             ssa30, sccd::max<T>(
                                        ssa32, sccd::max<T>(ssa34, ssa35)))))));
  const T ssa37 = sccd::min<T>(
      ssa17,
      sccd::min<T>(
          ssa19,
          sccd::min<T>(
              ssa21,
              sccd::min<T>(
                  ssa22, sccd::min<T>(
                             ssa30, sccd::min<T>(
                                        ssa32, sccd::min<T>(ssa34, ssa35)))))));
  const T ssa38 = -ssa37;
  const T ssa39 = -v2sy;
  const T ssa40 = ssa39 + v2ey;
  const T ssa41 = ssa40 * t_l;
  const T ssa42 = -v3sy;
  const T ssa43 = ssa42 + v3ey;
  const T ssa44 = ssa42 + v2sy;
  const T ssa45 = ssa41 - ssa43 * t_l + ssa44;
  const T ssa46 = ssa45 * v_l;
  const T ssa47 = v0ey - v0sy;
  const T ssa48 = ssa47 * t_l;
  const T ssa49 = -v1sy;
  const T ssa50 = ssa49 + v1ey;
  const T ssa51 = ssa49 + v0sy;
  const T ssa52 = ssa48 - ssa50 * t_l + ssa51;
  const T ssa53 = ssa39 + v0sy;
  const T ssa54 = -ssa41 + ssa48 + ssa53;
  const T ssa55 = -ssa52 * u_l + ssa54;
  const T ssa56 = ssa46 + ssa55;
  const T ssa57 = ssa45 * v_u;
  const T ssa58 = ssa55 + ssa57;
  const T ssa59 = -ssa52 * u_u + ssa54;
  const T ssa60 = ssa46 + ssa59;
  const T ssa61 = ssa57 + ssa59;
  const T ssa62 = ssa40 * t_u;
  const T ssa63 = -ssa43 * t_u + ssa44 + ssa62;
  const T ssa64 = ssa63 * v_l;
  const T ssa65 = ssa47 * t_u;
  const T ssa66 = -ssa50 * t_u + ssa51 + ssa65;
  const T ssa67 = ssa53 - ssa62 + ssa65;
  const T ssa68 = -ssa66 * u_l + ssa67;
  const T ssa69 = ssa64 + ssa68;
  const T ssa70 = ssa63 * v_u;
  const T ssa71 = ssa68 + ssa70;
  const T ssa72 = -ssa66 * u_u + ssa67;
  const T ssa73 = ssa64 + ssa72;
  const T ssa74 = ssa70 + ssa72;
  const T ssa75 = sccd::max<T>(
      ssa56,
      sccd::max<T>(
          ssa58,
          sccd::max<T>(
              ssa60,
              sccd::max<T>(
                  ssa61, sccd::max<T>(
                             ssa69, sccd::max<T>(
                                        ssa71, sccd::max<T>(ssa73, ssa74)))))));
  const T ssa76 = sccd::min<T>(
      ssa56,
      sccd::min<T>(
          ssa58,
          sccd::min<T>(
              ssa60,
              sccd::min<T>(
                  ssa61, sccd::min<T>(
                             ssa69, sccd::min<T>(
                                        ssa71, sccd::min<T>(ssa73, ssa74)))))));
  const T ssa77 = -ssa76;
  const T ssa78 = -v2sz;
  const T ssa79 = ssa78 + v2ez;
  const T ssa80 = ssa79 * t_l;
  const T ssa81 = -v3sz;
  const T ssa82 = ssa81 + v3ez;
  const T ssa83 = ssa81 + v2sz;
  const T ssa84 = ssa80 - ssa82 * t_l + ssa83;
  const T ssa85 = ssa84 * v_l;
  const T ssa86 = v0ez - v0sz;
  const T ssa87 = ssa86 * t_l;
  const T ssa88 = -v1sz;
  const T ssa89 = ssa88 + v1ez;
  const T ssa90 = ssa88 + v0sz;
  const T ssa91 = ssa87 - ssa89 * t_l + ssa90;
  const T ssa92 = ssa78 + v0sz;
  const T ssa93 = -ssa80 + ssa87 + ssa92;
  const T ssa94 = -ssa91 * u_l + ssa93;
  const T ssa95 = ssa85 + ssa94;
  const T ssa96 = ssa84 * v_u;
  const T ssa97 = ssa94 + ssa96;
  const T ssa98 = -ssa91 * u_u + ssa93;
  const T ssa99 = ssa85 + ssa98;
  const T ssa100 = ssa96 + ssa98;
  const T ssa101 = ssa79 * t_u;
  const T ssa102 = ssa101 - ssa82 * t_u + ssa83;
  const T ssa103 = ssa102 * v_l;
  const T ssa104 = ssa86 * t_u;
  const T ssa105 = ssa104 - ssa89 * t_u + ssa90;
  const T ssa106 = -ssa101 + ssa104 + ssa92;
  const T ssa107 = -ssa105 * u_l + ssa106;
  const T ssa108 = ssa103 + ssa107;
  const T ssa109 = ssa102 * v_u;
  const T ssa110 = ssa107 + ssa109;
  const T ssa111 = -ssa105 * u_u + ssa106;
  const T ssa112 = ssa103 + ssa111;
  const T ssa113 = ssa109 + ssa111;
  const T ssa114 = sccd::max<T>(
      ssa100,
      sccd::max<T>(
          ssa108,
          sccd::max<T>(
              ssa110,
              sccd::max<T>(
                  ssa112,
                  sccd::max<T>(
                      ssa113,
                      sccd::max<T>(ssa95, sccd::max<T>(ssa97, ssa99)))))));
  const T ssa115 = sccd::min<T>(
      ssa100,
      sccd::min<T>(
          ssa108,
          sccd::min<T>(
              ssa110,
              sccd::min<T>(
                  ssa112,
                  sccd::min<T>(
                      ssa113,
                      sccd::min<T>(ssa95, sccd::min<T>(ssa97, ssa99)))))));
  const T ssa116 = -ssa115;
  const T ssa117 = -errx;
  const T ssa118 = -erry;
  const T ssa119 = -errz;
  *true_tol =
      sccd::max<T>(0, sccd::max<T>(ssa114 + ssa116,
                                   sccd::max<T>(ssa36 + ssa38, ssa75 + ssa77)));
  if (-msx - ssa38 > errx || -msy - ssa77 > erry || -msz - ssa116 > errz ||
      msx + ssa36 < ssa117 || msy + ssa75 < ssa118 || msz + ssa114 < ssa119) {
    *box_in = 0;
    return 0;
  }
  *box_in =
      ((ssa117 > msx + ssa37 || ssa118 > msy + ssa76 || ssa119 > msz + ssa115 ||
        errx < -msx + ssa36 || erry < -msy + ssa75 || errz < -msz + ssa114)
           ? (0)
           : (1));
  return 1;
}
