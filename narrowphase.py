#!/usr/bin/env python3

import sympy as sp

# --- SymPy C code generator for tolerance functions ---

# This script generates C code (no classes/structs) for:
#   - compute_face_vertex_tolerance_soa
#   - compute_edge_edge_tolerance_soa
# using minimal SoA-style scalar inputs (xyz per point at start/end) and
# pointer outputs (tol0, tol1, tol2).


def symbols_inputs():
    codomain_tol = sp.symbols("codomain_tol", real=True)
    v0sx, v0sy, v0sz = sp.symbols("v0sx v0sy v0sz", real=True)
    v1sx, v1sy, v1sz = sp.symbols("v1sx v1sy v1sz", real=True)
    v2sx, v2sy, v2sz = sp.symbols("v2sx v2sy v2sz", real=True)
    v3sx, v3sy, v3sz = sp.symbols("v3sx v3sy v3sz", real=True)
    v0ex, v0ey, v0ez = sp.symbols("v0ex v0ey v0ez", real=True)
    v1ex, v1ey, v1ez = sp.symbols("v1ex v1ey v1ez", real=True)
    v2ex, v2ey, v2ez = sp.symbols("v2ex v2ey v2ez", real=True)
    v3ex, v3ey, v3ez = sp.symbols("v3ex v3ey v3ez", real=True)

    v0s = (v0sx, v0sy, v0sz)
    v1s = (v1sx, v1sy, v1sz)
    v2s = (v2sx, v2sy, v2sz)
    v3s = (v3sx, v3sy, v3sz)
    v0e = (v0ex, v0ey, v0ez)
    v1e = (v1ex, v1ey, v1ez)
    v2e = (v2ex, v2ey, v2ez)
    v3e = (v3ex, v3ey, v3ez)

    names = [
        "v0sx", "v0sy", "v0sz",
        "v1sx", "v1sy", "v1sz",
        "v2sx", "v2sy", "v2sz",
        "v3sx", "v3sy", "v3sz",
        "v0ex", "v0ey", "v0ez",
        "v1ex", "v1ey", "v1ez",
        "v2ex", "v2ey", "v2ez",
        "v3ex", "v3ey", "v3ez",
    ]
    return codomain_tol, v0s, v1s, v2s, v3s, v0e, v1e, v2e, v3e, names


def vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def linf_norm(v):
    return sp.Max(sp.Abs(v[0]), sp.Abs(v[1]), sp.Abs(v[2]))


def max_linf_4(p1s, p2s, p3s, p4s, p1e, p2e, p3e, p4e):
    return sp.Max(
        linf_norm(vsub(p1e, p1s)),
        linf_norm(vsub(p2e, p2s)),
        linf_norm(vsub(p3e, p3s)),
        linf_norm(vsub(p4e, p4s)),
    )


def vf_tolerances():
    codomain_tol, v0s, v1s, v2s, v3s, v0e, v1e, v2e, v3e, arg_names = symbols_inputs()
    # p*** as in CUDA reference
    p000 = vsub(v0s, v1s)
    p001 = vsub(v0s, v3s)
    p011 = vsub(v0s, (v2s[0] + v3s[0] - v1s[0], v2s[1] + v3s[1] - v1s[1], v2s[2] + v3s[2] - v1s[2]))
    p010 = vsub(v0s, v2s)
    p100 = vsub(v0e, v1e)
    p101 = vsub(v0e, v3e)
    p111 = vsub(v0e, (v2e[0] + v3e[0] - v1e[0], v2e[1] + v3e[1] - v1e[1], v2e[2] + v3e[2] - v1e[2]))
    p110 = vsub(v0e, v2e)

    den0 = max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110)
    den1 = max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011)
    den2 = max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011)

    t0 = codomain_tol / (3 * den0)
    t1 = codomain_tol / (3 * den1)
    t2 = codomain_tol / (3 * den2)
    return (t0, t1, t2), arg_names


def ee_tolerances():
    codomain_tol, v0s, v1s, v2s, v3s, v0e, v1e, v2e, v3e, arg_names = symbols_inputs()
    # p*** as in CUDA reference
    p000 = vsub(v0s, v2s)
    p001 = vsub(v0s, v3s)
    p010 = vsub(v1s, v2s)
    p011 = vsub(v1s, v3s)
    p100 = vsub(v0e, v2e)
    p101 = vsub(v0e, v3e)
    p110 = vsub(v1e, v2e)
    p111 = vsub(v1e, v3e)

    den0 = max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110)
    den1 = max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110)  # same as den0
    den2 = max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011)

    t0 = codomain_tol / (3 * den0)
    t1 = codomain_tol / (3 * den1)
    t2 = codomain_tol / (3 * den2)
    return (t0, t1, t2), arg_names


def ccode_custom(expr):
    # Replace standard fabs/fmax with custom sccd::abs/sccd::max
    s = sp.ccode(expr)
    s = s.replace("fmaxf(", "sccd::max<T>(").replace("fmax(", "sccd::max<T>(")
    s = s.replace("fabsf(", "sccd::abs<T>(").replace("fabs(", "sccd::abs<T>(")
    return s


def gen_c_function(name, tol_exprs, arg_names):
    args = ", ".join(
        ["const T codomain_tol"]
        + [f"const T {n}" for n in arg_names]
        + [
            "T * const SFEM_RESTRICT tol0",
            "T * const SFEM_RESTRICT tol1",
            "T * const SFEM_RESTRICT tol2",
        ]
    )
    body = []
    body.append(f"template<typename T> static inline void {name}({args}) {{")
    reps, reduced = sp.cse(list(tol_exprs), symbols=sp.numbered_symbols("ssa"), optimizations="basic")
    for sym, expr in reps:
        body.append(f"  const T {str(sym)} = {ccode_custom(expr)};")
    body.append(f"  *tol0 = {ccode_custom(reduced[0])};")
    body.append(f"  *tol1 = {ccode_custom(reduced[1])};")
    body.append(f"  *tol2 = {ccode_custom(reduced[2])};")
    body.append("}")
    return "\n".join(body)


if __name__ == "__main__":
    (vf_t0, vf_t1, vf_t2), vf_arg_names = vf_tolerances()
    vf_code = gen_c_function(
        "compute_face_vertex_tolerance_soa",
        (vf_t0, vf_t1, vf_t2),
        vf_arg_names,
    )

    (ee_t0, ee_t1, ee_t2), ee_arg_names = ee_tolerances()
    ee_code = gen_c_function(
        "compute_edge_edge_tolerance_soa",
        (ee_t0, ee_t1, ee_t2),
        ee_arg_names,
    )

    print("/* Generated by SymPy; requires <math.h> for fabs */")
    print(vf_code)
    print()
    print(ee_code)

