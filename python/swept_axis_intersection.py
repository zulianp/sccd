#!/usr/bin/env python3
"""
Align a swept triangle to the x-axis and extract intersection roots.
Steps:
 1) Translate so sv is at the origin.
 2) Rotate so ev - sv aligns with +x.
 3) Scale so ev lands at x = 1.
 4) Solve the swept triangle (s1->e1, s2->e2, s3->e3) against the x-axis
    by enforcing y = z = 0.
"""

import sympy as sp
from sympy.polys.polyerrors import PolynomialError

Float = sp.Float


def _skew(v):
    """Skew-symmetric matrix for cross products."""
    return sp.Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0],
    ])


def rotation_to_x(d_hat, eps_ang=Float("1e-12")):
    """
    Rotation matrix that maps d_hat to +x, with stable handling of
    near-parallel and anti-parallel cases.
    """
    e1 = sp.Matrix([1, 0, 0])
    v = d_hat.cross(e1)
    s2 = sp.simplify(v.dot(v))
    c = sp.simplify(d_hat.dot(e1))
    I = sp.eye(3)

    vx = _skew(v)
    generic = sp.simplify(I + vx + vx * vx * ((1 - c) / sp.Max(s2, eps_ang)))

    axis_seed = sp.Matrix([0, 1, 0]) if sp.Abs(d_hat[0]) > sp.Abs(d_hat[1]) else sp.Matrix([0, 0, 1])
    axis = d_hat.cross(axis_seed)
    axis = sp.simplify(axis / sp.sqrt(sp.Max(axis.dot(axis), eps_ang)))
    vx_pi = _skew(axis)
    half_turn = sp.simplify(I + 2 * vx_pi * vx_pi)

    small = sp.Lt(sp.Abs(s2), eps_ang)
    return sp.Piecewise(
        (I, sp.And(small, c >= 0)),
        (half_turn, sp.And(small, c < 0)),
        (generic, True),
    )


def build_transform(sv, ev, eps_len=Float("1e-12"), eps_ang=Float("1e-12")):
    """Return transform(p) plus metadata for the translate/rotate/scale."""
    sv = sp.Matrix(sv)
    ev = sp.Matrix(ev)
    d = ev - sv
    dlen2 = sp.simplify(d.dot(d))
    dlen = sp.sqrt(sp.Max(dlen2, 0))
    too_short = sp.Lt(dlen, eps_len)

    safe_len = sp.Max(dlen, eps_len)
    d_hat = sp.Piecewise((sp.Matrix([1, 0, 0]), too_short), (sp.simplify(d / safe_len), True))
    scale = sp.Piecewise((1, too_short), (1 / safe_len, True))
    R = rotation_to_x(sp.simplify(d_hat), eps_ang)

    def transform(p):
        return sp.simplify(scale * (R * (sp.Matrix(p) - sv)))

    meta = {"scale": scale, "rotation": R, "dlen": dlen, "too_short": too_short}
    return transform, meta


def swept_triangle_axis_system(sv, s1, s2, s3, ev, e1, e2, e3,
                               eps_len=Float("1e-12"), eps_ang=Float("1e-12")):
    """
    Build symbolic expressions for the swept triangle vs x-axis intersection.
    Returns expressions for t, u(t), v(t), x(t), determinant, and barycentric numerators.
    """
    transform, meta = build_transform(sv, ev, eps_len, eps_ang)
    S1, S2, S3 = map(transform, (s1, s2, s3))
    E1, E2, E3 = map(transform, (e1, e2, e3))

    t, u, v = sp.symbols("t u v", real=True)
    S = S1 + u * (S2 - S1) + v * (S3 - S1)
    E = E1 + u * (E2 - E1) + v * (E3 - E1)
    P = (1 - t) * S + t * E

    A = sp.Matrix([
        [(1 - t) * (S2[1] - S1[1]) + t * (E2[1] - E1[1]),
         (1 - t) * (S3[1] - S1[1]) + t * (E3[1] - E1[1])],
        [(1 - t) * (S2[2] - S1[2]) + t * (E2[2] - E1[2]),
         (1 - t) * (S3[2] - S1[2]) + t * (E3[2] - E1[2])],
    ])
    b = -sp.Matrix([
        (1 - t) * S1[1] + t * E1[1],
        (1 - t) * S1[2] + t * E1[2],
    ])

    detA = sp.simplify(A.det())
    u_expr = sp.simplify(sp.Matrix([[b[0], A[0, 1]], [b[1], A[1, 1]]]).det() / detA)
    v_expr = sp.simplify(sp.Matrix([[A[0, 0], b[0]], [A[1, 0], b[1]]]).det() / detA)

    x_expr = sp.simplify(P[0].subs({u: u_expr, v: v_expr}))

    num_u, den_u = sp.fraction(sp.together(u_expr))
    num_v, den_v = sp.fraction(sp.together(v_expr))
    num_sum, den_sum = sp.fraction(sp.together(u_expr + v_expr - 1))

    return {
        "t": t,
        "u_expr": u_expr,
        "v_expr": v_expr,
        "x_expr": x_expr,
        "det_expr": detA,
        "num_u": sp.simplify(num_u),
        "num_v": sp.simplify(num_v),
        "num_sum": sp.simplify(num_sum),
        "den_u": sp.simplify(den_u),
        "den_v": sp.simplify(den_v),
        "den_sum": sp.simplify(den_sum),
        "meta": meta,
    }


def _as_float_matrix(vec):
    return sp.Matrix([Float(v) for v in vec])


def solve_axis_hits_numeric(sv, s1, s2, s3, ev, e1, e2, e3,
                            eps_len=1e-12, eps_ang=1e-12, eps_det=1e-12, tol=1e-9):
    """
    Numeric helper: returns a list of intersections dicts {t,u,v,x} with 0<=t<=1.
    All inputs are iterables of length 3 (numbers).
    """
    sv = _as_float_matrix(sv)
    s1, s2, s3 = map(_as_float_matrix, (s1, s2, s3))
    ev = _as_float_matrix(ev)
    e1, e2, e3 = map(_as_float_matrix, (e1, e2, e3))

    system = swept_triangle_axis_system(
        sv, s1, s2, s3, ev, e1, e2, e3,
        eps_len=Float(str(eps_len)), eps_ang=Float(str(eps_ang)),
    )
    t = system["t"]
    u_expr = system["u_expr"]
    v_expr = system["v_expr"]
    x_expr = system["x_expr"]
    det_expr = system["det_expr"]

    # Detect continuous contact (u,v constant for all t and valid barycentric coords)
    if t not in u_expr.free_symbols and t not in v_expr.free_symbols:
        u_const = float(sp.N(u_expr))
        v_const = float(sp.N(v_expr))
        if u_const >= -tol and v_const >= -tol and (u_const + v_const) <= 1 + tol:
            try:
                det_poly = sp.Poly(sp.simplify(det_expr), t)
                det_roots = det_poly.nroots()
            except PolynomialError:
                det_roots = []
            det_root_in_range = any(
                abs(sp.im(r)) < 1e-10 and (-tol <= float(sp.re(r)) <= 1 + tol) for r in det_roots
            )
            det_samples = [det_expr.subs(t, s) for s in (0, 0.5, 1)]
            det_ok = all(abs(float(sp.N(ds))) >= eps_det for ds in det_samples)
            if det_ok and not det_root_in_range:
                x0 = float(sp.N(x_expr.subs(t, 0)))
                x1 = float(sp.N(x_expr.subs(t, 1)))
                # Represent the whole interval contact; keep endpoints for convenience.
                return [
                    {"t": 0.0, "u": u_const, "v": v_const, "x": x0, "segment": (0.0, 1.0)},
                    {"t": 1.0, "u": u_const, "v": v_const, "x": x1, "segment": (0.0, 1.0)},
                ]

    candidates = [0.0, 1.0]
    for expr in (system["num_u"], system["num_v"], system["num_sum"], det_expr):
        expr = sp.simplify(sp.N(expr))
        if expr == 0:
            continue
        try:
            poly = sp.Poly(expr, t)
            roots = poly.nroots()
        except PolynomialError:
            roots = []
        for r in roots:
            if abs(sp.im(r)) < 1e-10:
                candidates.append(float(sp.re(r)))

    hits = []
    for tv in candidates:
        if tv < -tol or tv > 1 + tol:
            continue
        det_val = complex(sp.N(det_expr.subs(t, tv)))
        if abs(det_val) < eps_det:
            continue
        u_val = sp.N(u_expr.subs(t, tv))
        v_val = sp.N(v_expr.subs(t, tv))
        if not u_val.is_real or not v_val.is_real:
            continue
        u_f = float(u_val)
        v_f = float(v_val)
        if u_f < -tol or v_f < -tol or (u_f + v_f) > 1 + tol:
            continue
        x_f = float(sp.N(x_expr.subs(t, tv)))
        hits.append({"t": tv, "u": u_f, "v": v_f, "x": x_f})

    hits.sort(key=lambda h: h["t"])
    return hits


def solve_vf_ccd_numeric(sv, ev, s1, s2, s3, e1, e2, e3,
                         tol=1e-9, eps_det=1e-12, eps_normal=1e-12):
    """
    Vertex-face CCD: find all (t,u,v) with 0<=t<=1 where the moving point
    (sv->ev) lies inside the moving triangle (s1..s3 -> e1..e3).
    Returns a list of dicts {t,u,v}.
    """
    sv = _as_float_matrix(sv)
    ev = _as_float_matrix(ev)
    s1, s2, s3 = map(_as_float_matrix, (s1, s2, s3))
    e1, e2, e3 = map(_as_float_matrix, (e1, e2, e3))

    t = sp.symbols("t", real=True)
    p = (1 - t) * sv + t * ev
    v1t = (1 - t) * s1 + t * e1
    v2t = (1 - t) * s2 + t * e2
    v3t = (1 - t) * s3 + t * e3

    e1t = v2t - v1t
    e2t = v3t - v1t
    n = sp.simplify(sp.Matrix(e1t).cross(sp.Matrix(e2t)))
    coplanar = sp.simplify((p - v1t).dot(n))

    hits = []
    try:
        poly = sp.Poly(coplanar, t)
        roots = poly.nroots()
    except PolynomialError:
        roots = []

    for r in roots:
        if abs(sp.im(r)) > 1e-10:
            continue
        tr = float(sp.re(r))
        if tr < -tol or tr > 1 + tol:
            continue

        v1v = v1t.subs(t, tr)
        e1v = e1t.subs(t, tr)
        e2v = e2t.subs(t, tr)
        n_val = sp.Matrix(e1v).cross(sp.Matrix(e2v))
        if float(sp.N(n_val.norm())) < eps_normal:
            continue

        A = sp.Matrix.hstack(e1v, e2v)
        b = p.subs(t, tr) - v1v
        detA = float(sp.N(A.det()))
        if abs(detA) < eps_det:
            continue
        uv = A.LUsolve(b)
        uval = float(sp.N(uv[0]))
        vval = float(sp.N(uv[1]))
        if uval < -tol or vval < -tol or (uval + vval) > 1 + tol:
            continue
        hits.append({"t": tr, "u": uval, "v": vval})

    hits.sort(key=lambda h: h["t"])
    return hits


if __name__ == "__main__":
    # Minimal smoke test with dummy points (replace with real data)
    sv = [0, 0.2, 0.2]
    ev = [1, 0.2, 0.2]
    s1, s2, s3 = [0, 0, 0], [0, 1, 0], [0, 0, 1]
    e1, e2, e3 = [1, 0, 0], [1, 1, 0], [1, 0, 1]

    print("axis hits:", solve_axis_hits_numeric(sv, s1, s2, s3, ev, e1, e2, e3))
    print("vf hits:", solve_vf_ccd_numeric(sv, ev, s1, s2, s3, e1, e2, e3))
