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

print("// VF")
print(generate_c_code("vf", vf_objective, vf_gradient, vf_Hessian))

print("// EE")
print(generate_c_code("ee", ee_objective, ee_gradient, ee_Hessian))


