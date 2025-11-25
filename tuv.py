#!/usr/bin/env python3

import sympy as sp


# Reference code
# struct Interval {
#   __device__ Interval() = default;
#   __device__ Interval(const Scalar &l, const Scalar &u) : lower(l), upper(u) {}
#   Scalar lower;
#   Scalar upper;
# };

# struct SplitInterval {
#   __device__ SplitInterval(const Interval &interval) {
#     const Scalar mid = (interval.lower + interval.upper) / 2;
#     first = Interval(interval.lower, mid);
#     second = Interval(mid, interval.upper);
#   }

#   Interval first;
#   Interval second;
# };

# class CCDDomain {
# public:
#   __device__ void init(int i) {
#     tuv[0] = Interval(0, 1);
#     tuv[1] = Interval(0, 1);
#     tuv[2] = Interval(0, 1);
#     query_id = i;
#   }

#   /// @brief The intervals for the t, u, and v parameters
#   Interval tuv[3];
#   /// @brief The query id
#   int query_id;
# };

# // this is to calculate the vertices of the inclusion function
# struct DomainCorner {
#   /// @brief Update the t, u, and v parameters based on the corner.
#   /// @param domain Domain intervals
#   /// @param corner The corner to use (the first bit is for the t parameter,
#   /// second for the u parameter, and third for the v parameter)
#   __device__ void update_tuv(const CCDDomain &domain, const uint8_t corner) {
#     t = (corner & 1) ? domain.tuv[0].upper : domain.tuv[0].lower;
#     u = (corner & 2) ? domain.tuv[1].upper : domain.tuv[1].lower;
#     v = (corner & 4) ? domain.tuv[2].upper : domain.tuv[2].lower;
#   }

#   Scalar t; ///< @brief The t parameter at the corner
#   Scalar u; ///< @brief The u parameter at the corner
#   Scalar v; ///< @brief The v parameter at the corner
# };

# __device__ Vector3 calculate_vf(const CCDData &data_in,
#                                 const DomainCorner &tuv) {
#   const Vector3 v = (data_in.v0e - data_in.v0s) * tuv.t + data_in.v0s;
#   const Vector3 t0 = (data_in.v1e - data_in.v1s) * tuv.t + data_in.v1s;
#   const Vector3 t1 = (data_in.v2e - data_in.v2s) * tuv.t + data_in.v2s;
#   const Vector3 t2 = (data_in.v3e - data_in.v3s) * tuv.t + data_in.v3s;
#   return v - (t1 - t0) * tuv.u - (t2 - t0) * tuv.v - t0;
# }

# __device__ Vector3 calculate_ee(const CCDData &data_in,
#                                 const DomainCorner &tuv) {
#   const Vector3 ea0 = (data_in.v0e - data_in.v0s) * tuv.t + data_in.v0s;
#   const Vector3 ea1 = (data_in.v1e - data_in.v1s) * tuv.t + data_in.v1s;
#   const Vector3 eb0 = (data_in.v2e - data_in.v2s) * tuv.t + data_in.v2s;
#   const Vector3 eb1 = (data_in.v3e - data_in.v3s) * tuv.t + data_in.v3s;
#   return ((ea1 - ea0) * tuv.u + ea0) - ((eb1 - eb0) * tuv.v + eb0);


# template <bool is_vf>
# __device__ bool origin_in_inclusion_function(const CCDData &data_in,
#                                              const CCDDomain &domain,
#                                              Scalar &true_tol, bool &box_in) {
#   Array3 codomain_min = Array3::Constant(SCALAR_MAX);
#   Array3 codomain_max = Array3::Constant(-SCALAR_MAX);

#   DomainCorner domain_corner;
#   for (uint8_t corner = 0; corner < 8; corner++) {
#     domain_corner.update_tuv(domain, corner);

#     Vector3 codomain_corner;
#     if constexpr (is_vf) {
#       codomain_corner = calculate_vf(data_in, domain_corner);
#     } else {
#       codomain_corner = calculate_ee(data_in, domain_corner);
#     }

#     codomain_min = codomain_min.min(codomain_corner.array());
#     codomain_max = codomain_max.max(codomain_corner.array());
#   }

#   // get the min and max in one dimension
#   true_tol = max(0.0, (codomain_max - codomain_min).maxCoeff());

#   box_in = true;

#   if ((codomain_min - data_in.ms > data_in.err).any() ||
#       (codomain_max + data_in.ms < -data_in.err).any()) {
#     return false;
#   }

#   if ((codomain_min + data_in.ms < -data_in.err).any() ||
#       (codomain_max - data_in.ms > data_in.err).any()) {
#     box_in = false;
#   }

#   return true;
# }

# --- SymPy C code generator for TUV kernels ---
import sympy as sp

# Common codegen helper: replace C math with sccd::* intrinsics
def ccode_custom(expr):
    s = sp.ccode(expr)
    s = s.replace("fmaxf(", "sccd::max<T>(").replace("fmax(", "sccd::max<T>(")
    s = s.replace("fminf(", "sccd::min<T>(").replace("fmin(", "sccd::min<T>(")
    s = s.replace("fabsf(", "sccd::abs<T>(").replace("fabs(", "sccd::abs<T>(")
    return s

# Geometry symbols (SoA)
def geom_symbols():
    v0sx, v0sy, v0sz = sp.symbols("v0sx v0sy v0sz", real=True)
    v1sx, v1sy, v1sz = sp.symbols("v1sx v1sy v1sz", real=True)
    v2sx, v2sy, v2sz = sp.symbols("v2sx v2sy v2sz", real=True)
    v3sx, v3sy, v3sz = sp.symbols("v3sx v3sy v3sz", real=True)
    v0ex, v0ey, v0ez = sp.symbols("v0ex v0ey v0ez", real=True)
    v1ex, v1ey, v1ez = sp.symbols("v1ex v1ey v1ez", real=True)
    v2ex, v2ey, v2ez = sp.symbols("v2ex v2ey v2ez", real=True)
    v3ex, v3ey, v3ez = sp.symbols("v3ex v3ey v3ez", real=True)
    v0s = (v0sx, v0sy, v0sz); v1s = (v1sx, v1sy, v1sz)
    v2s = (v2sx, v2sy, v2sz); v3s = (v3sx, v3sy, v3sz)
    v0e = (v0ex, v0ey, v0ez); v1e = (v1ex, v1ey, v1ez)
    v2e = (v2ex, v2ey, v2ez); v3e = (v3ex, v3ey, v3ez)
    arg_names = [
        "v0sx", "v0sy", "v0sz",
        "v1sx", "v1sy", "v1sz",
        "v2sx", "v2sy", "v2sz",
        "v3sx", "v3sy", "v3sz",
        "v0ex", "v0ey", "v0ez",
        "v1ex", "v1ey", "v1ez",
        "v2ex", "v2ey", "v2ez",
        "v3ex", "v3ey", "v3ez",
    ]
    return (v0s, v1s, v2s, v3s, v0e, v1e, v2e, v3e, arg_names)

def vadd(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def vsub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def vscale(a, s): return (a[0]*s, a[1]*s, a[2]*s)

# Linear interpolation over time t in [0,1]
def lerp3(ps, pe, t):
    return vadd(ps, vscale(vsub(pe, ps), t))

# Symbolic calculate_vf
def vf_expr(tt, uu, vv):
    v0s,v1s,v2s,v3s,v0e,v1e,v2e,v3e,_ = geom_symbols()
    v = lerp3(v0s, v0e, tt)
    t0 = lerp3(v1s, v1e, tt)
    t1 = lerp3(v2s, v2e, tt)
    t2 = lerp3(v3s, v3e, tt)
    term1 = vsub(t1, t0)  # t1 - t0
    term2 = vsub(t2, t0)  # t2 - t0
    r = vsub(v, vadd(vadd(vscale(term1, uu), vscale(term2, vv)), t0))
    return r

# Symbolic calculate_ee
def ee_expr(tt, uu, vv):
    v0s,v1s,v2s,v3s,v0e,v1e,v2e,v3e,_ = geom_symbols()
    ea0 = lerp3(v0s, v0e, tt)
    ea1 = lerp3(v1s, v1e, tt)
    eb0 = lerp3(v2s, v2e, tt)
    eb1 = lerp3(v3s, v3e, tt)
    a = vadd(vscale(vsub(ea1, ea0), uu), ea0)
    b = vadd(vscale(vsub(eb1, eb0), vv), eb0)
    r = vsub(a, b)
    return r

def gen_vec_function(name, comp_exprs, extra_scalar_params):
    # Build signature (T for all scalars; SoA geometry; rx,ry,rz)
    _,_,_,_,_,_,_,_, geom_names = geom_symbols()
    args = (
        [f"T {p}" for p in extra_scalar_params]
        + [f"T {n}" for n in geom_names]
        + [
            "T * const SFEM_RESTRICT rx",
            "T * const SFEM_RESTRICT ry",
            "T * const SFEM_RESTRICT rz",
        ]
    )
    body = []
    body.append(f"template<typename T> void {name}({', '.join(args)}) {{")
    reps, reduced = sp.cse(list(comp_exprs), symbols=sp.numbered_symbols("ssa"), optimizations="basic")
    for sym, expr in reps:
        body.append(f"  const T {str(sym)} = {ccode_custom(expr)};")
    body.append(f"  *rx = {ccode_custom(reduced[0])};")
    body.append(f"  *ry = {ccode_custom(reduced[1])};")
    body.append(f"  *rz = {ccode_custom(reduced[2])};")
    body.append("}")
    return "\n".join(body)

def gen_origin_function(name, is_vf):
    # Domain interval scalars
    tl, tu, ul, uu, vl, vu = sp.symbols("t_l t_u u_l u_u v_l v_u", real=True)
    # ms and err (per component)
    msx, msy, msz = sp.symbols("msx msy msz", real=True)
    ex, ey, ez = sp.symbols("errx erry errz", real=True)
    # Geometry
    v0s,v1s,v2s,v3s,v0e,v1e,v2e,v3e, geom_names = geom_symbols()
    # Evaluate inclusion function at 8 corners
    corners = []
    for t_ in (tl, tu):
        for u_ in (ul, uu):
            for v_ in (vl, vu):
                r = vf_expr(t_, u_, v_) if is_vf else ee_expr(t_, u_, v_)
                corners.append(r)
    xs = [c[0] for c in corners]
    ys = [c[1] for c in corners]
    zs = [c[2] for c in corners]
    cminx = sp.Min(*xs); cmaxx = sp.Max(*xs)
    cminy = sp.Min(*ys); cmaxy = sp.Max(*ys)
    cminz = sp.Min(*zs); cmaxz = sp.Max(*zs)
    # true_tol = max(0, max((cmax - cmin) per axis))
    true_tol = sp.Max(
        0,
        sp.Max(cmaxx - cminx, cmaxy - cminy, cmaxz - cminz)
    )
    # early out condition
    cond_false = sp.Or(
        sp.Or(cminx - msx > ex, cminy - msy > ey, cminz - msz > ez),
        sp.Or(cmaxx + msx < -ex, cmaxy + msy < -ey, cmaxz + msz < -ez),
    )
    # box_in false condition
    cond_not_box = sp.Or(
        sp.Or(cminx + msx < -ex, cminy + msy < -ey, cminz + msz < -ez),
        sp.Or(cmaxx - msx > ex, cmaxy - msy > ey, cmaxz - msz > ez),
    )
    # Build signature
    args = (
        [
            "T t_l", "T t_u",
            "T u_l", "T u_u",
            "T v_l", "T v_u",
        ]
        + [f"T {n}" for n in geom_names]
        + [
            "T msx", "T msy", "T msz",
            "T errx", "T erry", "T errz",
            "T * const SFEM_RESTRICT true_tol",
            "int * const SFEM_RESTRICT box_in",
        ]
    )
    fn = []
    fn.append(f"template<typename T> int {name}({', '.join(args)}) {{")
    # CSE across all needed expressions
    reps, reduced = sp.cse(
        [true_tol, cond_false, cond_not_box],
        symbols=sp.numbered_symbols("ssa"),
        optimizations="basic",
    )
    for sym, expr in reps:
        fn.append(f"  const T {str(sym)} = {ccode_custom(expr)};")
    # Unpack reduced
    rt, rfalse, rnotbox = reduced
    fn.append(f"  *true_tol = {ccode_custom(rt)};")
    fn.append(f"  if ({ccode_custom(rfalse)}) {{")
    fn.append(f"    *box_in = 0;")
    fn.append(f"    return 0;")
    fn.append(f"  }}")
    fn.append(f"  *box_in = {ccode_custom(sp.Piecewise((0, rnotbox), (1, True)))};")
    fn.append(f"  return 1;")
    fn.append("}")
    return "\n".join(fn)

if __name__ == "__main__":
    # calculate_vf
    tt, uu, vv = sp.symbols("t u v", real=True)
    rx, ry, rz = vf_expr(tt, uu, vv)
    print(gen_vec_function("sccd_calculate_vf_soa", (rx, ry, rz), ["t", "u", "v"]))
    print()
    # calculate_ee
    rx2, ry2, rz2 = ee_expr(tt, uu, vv)
    print(gen_vec_function("sccd_calculate_ee_soa", (rx2, ry2, rz2), ["t", "u", "v"]))
    print()
    # origin_in_inclusion (vf)
    print(gen_origin_function("sccd_origin_in_inclusion_vf_soa", is_vf=True))
    print()
    # origin_in_inclusion (ee)
    print(gen_origin_function("sccd_origin_in_inclusion_ee_soa", is_vf=False))