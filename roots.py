#!/usr/bin/env python3

import sympy as sp
from sympy.codegen.rewriting import optimize, create_expand_pow_optimization
from sympy.printing.c import C99CodePrinter

class SCCDPrinter(C99CodePrinter):
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        if exp.is_Integer:
            n = int(exp)
            if 2 <= n <= 9:
                return f"sccd::pow{n}<T>({self._print(base)})"
        return super()._print_Pow(expr)

    def _print_Abs(self, expr):
        return f"sccd::abs<T>({self._print(expr.args[0])})"

    def _print_Max(self, expr):
        args = [self._print(a) for a in expr.args]
        acc = args[0]
        for a in args[1:]:
            acc = f"sccd::max<T>({acc}, {a})"
        return acc


_PRINTER = SCCDPrinter()

def simplify(expr):
    return expr

sv = sp.Matrix(3, 1, list(sp.symbols("sv[0:3]")))
s1 = sp.Matrix(3, 1, list(sp.symbols("s1[0:3]")))
s2 = sp.Matrix(3, 1, list(sp.symbols("s2[0:3]")))
s3 = sp.Matrix(3, 1, list(sp.symbols("s3[0:3]")))
s4 = sp.Matrix(3, 1, list(sp.symbols("s3[0:3]")))

ev = sp.Matrix(3, 1, list(sp.symbols("ev[0:3]")))
e1 = sp.Matrix(3, 1, list(sp.symbols("e1[0:3]")))
e2 = sp.Matrix(3, 1, list(sp.symbols("e2[0:3]")))
e3 = sp.Matrix(3, 1, list(sp.symbols("e3[0:3]")))
e4 = sp.Matrix(3, 1, list(sp.symbols("e4[0:3]")))

t, u, v = sp.symbols('t u v')
theta = sp.Matrix(3, 1, [t, u, v])

p0 = (1 - t) * sv + t * ev

p1 = (1 - t) * s1 + t * e1
p2 = (1 - t) * s2 + t * e2
p3 = (1 - t) * s3 + t * e3

p4 = (1 - t) * s4 + t * e4


# ----------------------------------------------
# VF
# ----------------------------------------------

# Point
pt = p0

# Face
v1t = p1
v2t = p2
v3t = p3

# Root for vertex face collision
vf_F = pt - ((1 - u - v) * v1t + u * v2t + v * v3t)
vf_dFdt = vf_F.jacobian(theta)

# Least-square solution
vf_objective = simplify(vf_F.T * vf_F)
vf_gradient  = simplify(vf_objective.jacobian(theta))
vf_Hessian   = simplify(vf_gradient.jacobian(theta))

# 1) Find minimum
# 2) check if it is a zero
# 3) check if t,u,v are valid
# 4) if 2 and 3 are true t is the local impact time



# ----------------------------------------------
# EE
# ----------------------------------------------

# Edge 1
p1t = p1
p2t = p2

# Edge 2
p3t = p3
p4t = p4

# Root for the edge-edge collision
ee_F = ((1 - u) * p1t + u * p2t) - ((1 - v) * p3t + v * p4t)
ee_dFdt = ee_F.jacobian(theta)

# Least-square solution
ee_objective = simplify(ee_F.T * ee_F)
ee_gradient  = simplify(ee_objective.jacobian(theta))
ee_Hessian   = simplify(ee_gradient.jacobian(theta))

def generate_c_code(name,objective, gradient, Hessian):
    """
    Generate C code for objective (scalar), gradient (3), and Hessian (3x3)
    using SymPy CSE and simple power-expansion optimizations.
    Returns a single C translation unit as a string.
    """
    def _as_scalar(expr_or_mat):
        if isinstance(expr_or_mat, sp.MatrixBase):
            return sp.simplify(expr_or_mat[0, 0])
        return sp.simplify(expr_or_mat)

    def _as_vector(expr_or_mat):
        if isinstance(expr_or_mat, sp.MatrixBase):
            flat = list(expr_or_mat)
            return [sp.simplify(e) for e in flat]
        return [sp.simplify(expr_or_mat)]

    def _as_matrix(expr_or_mat, rows=3, cols=3):
        if isinstance(expr_or_mat, sp.MatrixBase):
            M = expr_or_mat
        else:
            M = sp.Matrix(expr_or_mat)
        return [[sp.simplify(M[i, j]) for j in range(cols)] for i in range(rows)]

    def _pow_opt(exprs):
        opt = [create_expand_pow_optimization(4)]
        return [optimize(e, opt) for e in exprs]

    def _emit_function(signature, exprs, out_assigns):
        # exprs: list[Expr] corresponding to each assignment in out_assigns sequence
        # out_assigns: iterable of strings like "out_f[0] = {expr};"
        c_lines = []
        c_lines.append(signature)
        c_lines.append("{")
        # CSE on all exprs
        exprs_opt = _pow_opt(exprs)
        tmp_syms = sp.numbered_symbols(prefix="ssa")
        cse_repls, cse_exprs = sp.cse(exprs_opt, symbols=tmp_syms, order="canonical")
        # Decls for temporaries
        for symb, rhs in cse_repls:
            c_lines.append(f"  const T {_PRINTER.doprint(symb)} = {_PRINTER.doprint(rhs)};")
        # Final assignments
        for dst, e in zip(out_assigns, cse_exprs):
            c_lines.append(f"  {dst} {_PRINTER.doprint(e)};")
        c_lines.append("}")
        return "\n".join(c_lines)

    # Normalize inputs
    f_expr = _as_scalar(objective)
    g_exprs = _as_vector(gradient)
    # Flatten Hessian row-major
    H_mat = _as_matrix(Hessian, rows=3, cols=3)
    H_exprs = [H_mat[i][j] for i in range(3) for j in range(3)]

    # Common header and typedef
    header = []
    header.append("/* Auto-generated with SymPy (CSE + pow expansion) */")
    header.append("#include \"smath.hpp\"")

    # Argument list shared by all functions
    args = (
        "const T sv[3], const T s1[3], const T s2[3], "
        "const T s3[3], const T s4[3], "
        "const T ev[3], const T e1[3], const T e2[3], "
        "const T e3[3], const T e4[3], "
        "const T t, const T u, const T v"
    )

    # Objective
    obj_signature = f"template<typename T> static inline void {name}_objective({args}, T *out_f)"
    obj_exprs = [f_expr]
    obj_assigns = ["*out_f ="]
    obj_body = _emit_function(obj_signature, obj_exprs, obj_assigns)

    # Gradient (size 3, row-vector in input)
    grad_signature = f"template<typename T> static inline void {name}_gradient({args}, T out_g[3])"
    grad_exprs = g_exprs
    grad_assigns = [f"out_g[{i}] =" for i in range(3)]
    grad_body = _emit_function(grad_signature, grad_exprs, grad_assigns)

    # Hessian (3x3, row-major)
    hess_signature = f"template<typename T> static inline void {name}_hessian({args}, T out_H[9])"
    hess_exprs = H_exprs
    hess_assigns = [f"out_H[{i*3 + j}] =" for i in range(3) for j in range(3)]
    hess_body = _emit_function(hess_signature, hess_exprs, hess_assigns)

    # Objective + Newton direction p = inv(H) * g (solve H p = g)
    p_exprs = list(sp.Matrix(H_mat).LUsolve(sp.Matrix(g_exprs)))
    objdir_signature = f"template<typename T> static inline void {name}_objective_dir({args}, T *out_f, T out_p[3])"
    objdir_exprs = [f_expr] + p_exprs
    objdir_assigns = ["*out_f ="] + [f"out_p[{i}] =" for i in range(3)]
    objdir_body = _emit_function(objdir_signature, objdir_exprs, objdir_assigns)

    # Fused computation (shared CSE across all three)
    fused_signature = (
        f"template<typename T> static inline void {name}_all({args}, T *out_f, T out_g[3], T out_H[9])"
    )
    fused_lines = []
    fused_lines.append(fused_signature)
    fused_lines.append("{")
    all_exprs_opt = _pow_opt([f_expr] + g_exprs + H_exprs)
    tmp_syms = sp.numbered_symbols(prefix="ssa")
    cse_repls, cse_exprs = sp.cse(all_exprs_opt, symbols=tmp_syms, order="canonical")
    # temporaries
    for symb, rhs in cse_repls:
        fused_lines.append(f"  const T {_PRINTER.doprint(symb)} = {_PRINTER.doprint(rhs)};")
    # outputs
    idx = 0
    fused_lines.append(f"  *out_f = {_PRINTER.doprint(cse_exprs[idx])};")
    idx += 1
    for i in range(3):
        fused_lines.append(f"  out_g[{i}] = {_PRINTER.doprint(cse_exprs[idx])};")
        idx += 1
    for i in range(3):
        for j in range(3):
            fused_lines.append(f"  out_H[{i*3 + j}] = {_PRINTER.doprint(cse_exprs[idx])};")
            idx += 1
    fused_lines.append("}")
    fused_body = "\n".join(fused_lines)

    return "\n\n".join(header + [obj_body, "", grad_body, "", hess_body, "", objdir_body, "", fused_body])

if __name__ == "__main__":
    print("// VF")
    print(generate_c_code("vf", vf_objective, vf_gradient, vf_Hessian))
    
    print("// EE")
    print(generate_c_code("ee", ee_objective, ee_gradient, ee_Hessian))



# # ------------------------------------------------------------
# # 2D CCD (point vs triangle) utilities and visualization
# # ------------------------------------------------------------
# from math import floor

# # Keep chunk size similar to C++ heuristic
# ROOT_FINDING_CHUNK_SIZE = 4096


# def _vf_F_2d(sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d, tt, uu, vv):
#     """
#     2D version of the vertex-face residual:
#       F(t,u,v) = V(t) - BaryFace(t,u,v)
#     All inputs are 2D sequences (len==2). Returns tuple (Fx, Fy).
#     """
#     t0 = 1.0 - tt
#     t1 = tt
#     one_minus_uv = 1.0 - uu - vv

#     vx = t0 * sv_2d[0] + t1 * ev_2d[0]
#     vy = t0 * sv_2d[1] + t1 * ev_2d[1]

#     fx = t0 * (one_minus_uv * s1_2d[0] + uu * s2_2d[0] + vv * s3_2d[0]) \
#        + t1 * (one_minus_uv * e1_2d[0] + uu * e2_2d[0] + vv * e3_2d[0])
#     fy = t0 * (one_minus_uv * s1_2d[1] + uu * s2_2d[1] + vv * s3_2d[1]) \
#        + t1 * (one_minus_uv * e1_2d[1] + uu * e2_2d[1] + vv * e3_2d[1])

#     return (vx - fx, vy - fy)


# def _sample_Fvf_component_2d(
#     n_t, n_u, n_v,
#     t_stride, u_stride,
#     t_start, u_start, v_start,
#     t_step, u_step, v_step,
#     comp,  # 0 for x, 1 for y
#     sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d
# ):
#     """
#     Sample a single component (x or y) of F_vf on a regular 3D grid over (t,u,v).
#     Returns a flat list/array of length (n_t+1)*(n_u+1)*(n_v+1).
#     """
#     total = (n_t + 1) * (n_u + 1) * (n_v + 1)
#     F = [0.0] * total
#     for ti in range(n_t + 1):
#         t_val = t_start + ti * t_step
#         base_t = ti * t_stride
#         for ui in range(n_u + 1):
#             u_val = u_start + ui * u_step
#             base_u = base_t + ui * u_stride
#             for vi in range(n_v + 1):
#                 v_val = v_start + vi * v_step
#                 idx = base_u + vi
#                 Fx, Fy = _vf_F_2d(sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d, t_val, u_val, v_val)
#                 F[idx] = Fx if comp == 0 else Fy
#     return F


# def _detect_zero_cells_2d(n_t, n_u, n_v, t_stride, u_stride, Fx, Fy, tol, contains_zero):
#     """
#     For each 3D cell (8 corners) check whether 0 is contained in the interval
#     of both Fx and Fy components: [min(corners), max(corners)] overlapping 0.
#     contains_zero is a list sized n_t*n_u*n_v (one per cell) that gets updated in-place.
#     """
#     for ti in range(n_t):
#         for ui in range(n_u):
#             for vi in range(n_v):
#                 cell_idx = ti * n_u * n_v + ui * n_v + vi
#                 i0 = ti * t_stride + ui * u_stride + vi
#                 i1 = i0 + 1
#                 i2 = i0 + u_stride
#                 i3 = i2 + 1
#                 i4 = i0 + t_stride
#                 i5 = i1 + t_stride
#                 i6 = i2 + t_stride
#                 i7 = i3 + t_stride

#                 # Fx mins/maxs over 8 corners
#                 fxs = (Fx[i0], Fx[i1], Fx[i2], Fx[i3], Fx[i4], Fx[i5], Fx[i6], Fx[i7])
#                 fys = (Fy[i0], Fy[i1], Fy[i2], Fy[i3], Fy[i4], Fy[i5], Fy[i6], Fy[i7])

#                 fx_min, fx_max = min(fxs), max(fxs)
#                 fy_min, fy_max = min(fys), max(fys)

#                 # Check interval containment of zero with tolerance
#                 has_zero_fx = (fx_min <= tol) and (fx_max >= -tol)
#                 has_zero_fy = (fy_min <= tol) and (fy_max >= -tol)
#                 contains_zero[cell_idx] &= int(has_zero_fx and has_zero_fy)


# def find_root_2d(
#     _max_iter,
#     tol,
#     sv_2d, s1_2d, s2_2d, s3_2d,  # 2D arrays-like
#     ev_2d, e1_2d, e2_2d, e3_2d   # 2D arrays-like
# ):
#     """
#     Broad-phase root search for 2D vertex-face CCD using interval tests on Fx,Fy.
#     Returns (found, t, u, v). If found=False, the remaining values are None.
#     """
#     # Coarse 2x2x2 sampling to estimate ranges
#     Fx_coarse = _sample_Fvf_component_2d(
#         1, 1, 1,  # n_t, n_u, n_v -> 2 samples in each dim
#         t_stride=4, u_stride=2,
#         t_start=0.0, u_start=0.0, v_start=0.0,
#         t_step=1.0, u_step=1.0, v_step=1.0,
#         comp=0, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
#     )
#     Fy_coarse = _sample_Fvf_component_2d(
#         1, 1, 1,
#         t_stride=4, u_stride=2,
#         t_start=0.0, u_start=0.0, v_start=0.0,
#         t_step=1.0, u_step=1.0, v_step=1.0,
#         comp=1, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
#     )

#     fx_min, fx_max = min(Fx_coarse), max(Fx_coarse)
#     fy_min, fy_max = min(Fy_coarse), max(Fy_coarse)

#     # Quick AABB reject
#     fx_intersect = (fx_min <= tol) and (fx_max >= -tol)
#     fy_intersect = (fy_min <= tol) and (fy_max >= -tol)
#     if not (fx_intersect and fy_intersect):
#         return (False, None, None, None)

#     # Allocate more resolution to dimensions that have smaller value ranges,
#     # similar to C++ heuristic.
#     def safe_inv_range(vmin, vmax):
#         rng = max(1e-5, (vmax - vmin))
#         return max(2.0, 1.0 / rng)

#     # For 2D we only have two components; combine them conservatively
#     inv_x = safe_inv_range(fx_min, fx_max)
#     inv_y = safe_inv_range(fy_min, fy_max)
#     # Use the product to allocate per-dimension counts
#     total_inv = inv_x * inv_y

#     # Start balanced, then redistribute proportionally
#     # We have three parameters (t,u,v); split effort roughly evenly,
#     # then bias towards larger total_inv.
#     base = ROOT_FINDING_CHUNK_SIZE
#     # A simple proportional allocation that respects the total budget
#     # and keeps at least 2 samples per dimension.
#     # We bias t slightly to find the earliest impact.
#     bias_t, bias_u, bias_v = 1.2, 1.0, 1.0
#     weight_sum = bias_t + bias_u + bias_v
#     Nt = max(2, int(floor((base ** (1/3)) * (bias_t / weight_sum) * total_inv ** (1/3))))
#     Nu = max(2, int(floor((base ** (1/3)) * (bias_u / weight_sum) * total_inv ** (1/3))))
#     Nv = max(2, int(floor((base ** (1/3)) * (bias_v / weight_sum) * total_inv ** (1/3))))
#     # Ensure product under budget; adjust down if necessary
#     while (Nt * Nu * Nv) > ROOT_FINDING_CHUNK_SIZE:
#         # reduce the largest
#         if Nt >= Nu and Nt >= Nv and Nt > 2:
#             Nt -= 1
#         elif Nu >= Nt and Nu >= Nv and Nu > 2:
#             Nu -= 1
#         elif Nv > 2:
#             Nv -= 1
#         else:
#             break

#     # Convert to cell counts
#     t_n = Nt - 1
#     u_n = Nu - 1
#     v_n = Nv - 1
#     t_min = 0.0
#     u_min = 0.0
#     v_min = 0.0
#     t_max = 1.0
#     u_max = 1.0
#     v_max = 1.0
#     t_h = (t_max - t_min) / t_n
#     u_h = (u_max - u_min) / u_n
#     v_h = (v_max - v_min) / v_n

#     t_stride = (u_n + 1) * (v_n + 1)
#     u_stride = (v_n + 1)

#     Fx = _sample_Fvf_component_2d(
#         t_n, u_n, v_n,
#         t_stride, u_stride,
#         t_min, u_min, v_min,
#         t_h, u_h, v_h,
#         comp=0, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
#     )
#     Fy = _sample_Fvf_component_2d(
#         t_n, u_n, v_n,
#         t_stride, u_stride,
#         t_min, u_min, v_min,
#         t_h, u_h, v_h,
#         comp=1, sv_2d=sv_2d, s1_2d=s1_2d, s2_2d=s2_2d, s3_2d=s3_2d, ev_2d=ev_2d, e1_2d=e1_2d, e2_2d=e2_2d, e3_2d=e3_2d
#     )

#     contains_zero = [1] * (t_n * u_n * v_n)
#     _detect_zero_cells_2d(t_n, u_n, v_n, t_stride, u_stride, Fx, Fy, tol, contains_zero)

#     found = False
#     best_t = None
#     best_u = None
#     best_v = None
#     # Scan in t-major order to get earliest t in this discretization
#     for ti in range(t_n):
#         if found:
#             break
#         for ui in range(u_n):
#             for vi in range(v_n):
#                 cell_idx = ti * u_n * v_n + ui * v_n + vi
#                 if contains_zero[cell_idx]:
#                     t_current = t_min + ti * t_h
#                     u_current = u_min + ui * u_h
#                     v_current = v_min + vi * v_h
#                     if (u_h + v_h) > 1.0 + tol:
#                         # coarse cell entirely outside feasible region; skip
#                         continue
#                     # Accept this coarse root; optional refinement could be added
#                     best_t, best_u, best_v = t_current, u_current, v_current
#                     found = True
#                     break
#     if not found:
#         return (False, None, None, None)
#     return (True, best_t, best_u, best_v)


# def visualize_Fvf_xy_2d(
#     sv_2d, s1_2d, s2_2d, s3_2d, ev_2d, e1_2d, e2_2d, e3_2d,
#     t_value=0.5,
#     resolution=256,
#     out_path="vf2d_Fxy.png"
# ):
#     """
#     Save two heatmaps (Fx and Fy) of F_vf(t,u,v) evaluated at a fixed t over (u,v) in [0,1], u+v<=1.
#     """
#     import numpy as np
#     import matplotlib
#     matplotlib.use("Agg")
#     import matplotlib.pyplot as plt

#     sv_arr = np.asarray(sv_2d, dtype=float)
#     s1_arr = np.asarray(s1_2d, dtype=float)
#     s2_arr = np.asarray(s2_2d, dtype=float)
#     s3_arr = np.asarray(s3_2d, dtype=float)
#     ev_arr = np.asarray(ev_2d, dtype=float)
#     e1_arr = np.asarray(e1_2d, dtype=float)
#     e2_arr = np.asarray(e2_2d, dtype=float)
#     e3_arr = np.asarray(e3_2d, dtype=float)

#     U, V = np.meshgrid(
#         np.linspace(0.0, 1.0, num=resolution),
#         np.linspace(0.0, 1.0, num=resolution),
#         indexing="xy"
#     )
#     mask = (U + V) <= 1.0
#     O = 1.0 - U - V
#     t0 = 1.0 - t_value
#     t1 = t_value

#     vx = t0 * sv_arr[0] + t1 * ev_arr[0]
#     vy = t0 * sv_arr[1] + t1 * ev_arr[1]

#     fx = t0 * (O * s1_arr[0] + U * s2_arr[0] + V * s3_arr[0]) + t1 * (O * e1_arr[0] + U * e2_arr[0] + V * e3_arr[0])
#     fy = t0 * (O * s1_arr[1] + U * s2_arr[1] + V * s3_arr[1]) + t1 * (O * e1_arr[1] + U * e2_arr[1] + V * e3_arr[1])

#     Fx = vx - fx
#     Fy = vy - fy

#     # Mask outside simplex
#     Fx_plot = np.full_like(Fx, np.nan)
#     Fy_plot = np.full_like(Fy, np.nan)
#     Fx_plot[mask] = Fx[mask]
#     Fy_plot[mask] = Fy[mask]

#     fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
#     im0 = axs[0].imshow(
#         Fx_plot.T, origin="lower", extent=(0, 1, 0, 1), cmap="coolwarm"
#     )
#     axs[0].plot([0, 1], [1, 0], 'k--', lw=0.8)  # u+v=1 boundary
#     axs[0].set_title(f"F_vf.x at t={t_value:.3f}")
#     axs[0].set_xlabel("u"); axs[0].set_ylabel("v")
#     fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

#     im1 = axs[1].imshow(
#         Fy_plot.T, origin="lower", extent=(0, 1, 0, 1), cmap="coolwarm"
#     )
#     axs[1].plot([0, 1], [1, 0], 'k--', lw=0.8)
#     axs[1].set_title(f"F_vf.y at t={t_value:.3f}")
#     axs[1].set_xlabel("u"); axs[1].set_ylabel("v")
#     fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

#     fig.suptitle("2D Vertex-Face Residual Components")
#     fig.savefig(out_path, dpi=150)
#     plt.close(fig)
