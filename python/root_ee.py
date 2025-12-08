from __future__ import annotations

import json
from typing import Any, Dict, List, Sequence, Union

from sympy import (
    And,
    Eq,
    I,
    Interval,
    S,
    Symbol,
    Union as SymUnion,
    simplify,
    sympify,
)
from sympy.core.expr import Expr
from sympy.sets.sets import FiniteSet, Set
from sympy.solvers.inequalities import solve_univariate_inequality
from sympy.solvers.solveset import linsolve, nonlinsolve


def _to_expr(x: Any) -> Expr:
    if isinstance(x, Expr):
        return x
    return sympify(x, evaluate=True)


def _row_ratios(row: Sequence[Any]) -> List[Expr]:
    # Expect row to have 6 entries: [n0, d0, n1, d1, n2, d2]
    # Returns 3 simplified ratios.
    ratios: List[Expr] = []
    for i in range(3):
        num = _to_expr(row[2 * i])
        den = _to_expr(row[2 * i + 1])
        ratios.append(simplify(num / den))
    return ratios


def _set_measure(possible_t: Set) -> Union[int, float]:
    # Lebesgue measure for sets of reals; FiniteSet has measure 0.
    try:
        return float(possible_t.measure)
    except (AttributeError, TypeError, NotImplementedError):
        # If measure isn't available, conservatively fallback to non-emptiness
        # but treat FiniteSet as zero-measure.
        if isinstance(possible_t, FiniteSet):
            return 0.0
        return 0.0


def _pick_representative_t(possible_t: Set) -> Union[float, None]:
    # Choose a representative t in the feasible set, prefer interval interior.
    if possible_t.is_EmptySet:
        return None
    if isinstance(possible_t, Interval):
        if possible_t.measure == 0:
            # Single point (closed-degenerate interval)
            return float(possible_t.start)
        mid = (possible_t.start + possible_t.end) / 2
        return float(mid)
    if isinstance(possible_t, SymUnion):
        for subset in possible_t.args:
            val = _pick_representative_t(subset)
            if val is not None:
                return val
        return None
    if isinstance(possible_t, FiniteSet):
        # Pick the first point
        return float(list(possible_t)[0])
    # Fallback: not a known set shape
    return None


def roots(data: Sequence[Sequence[Any]], outfile: str = "") -> bool:
    """
    Convert the Mathematica routine to SymPy:
    - data: 8 rows, each with 6 entries: [n0, d0, n1, d1, n2, d2]
      rows 0..3 -> a0s, a1s, b0s, b1s (start)
      rows 4..7 -> a0e, a1e, b0e, b1e (end)
    - Solve for t, la, lb such that point A(t, la) == B(t, lb)
    - Accept solutions where t, la, lb are real and in [0, 1]
      If parameterized by t, accept if feasible t-set in [0,1] has positive measure.
    - If outfile is non-empty and solutions exist, write JSON array of dicts.
    """
    if len(data) < 8:
        raise ValueError("Expected at least 8 rows in data.")
    for r in range(8):
        if len(data[r]) < 6:
            raise ValueError(f"Row {r} must have at least 6 entries.")

    a0s = _row_ratios(data[0])
    a1s = _row_ratios(data[1])
    b0s = _row_ratios(data[2])
    b1s = _row_ratios(data[3])

    a0e = _row_ratios(data[4])
    a1e = _row_ratios(data[5])
    b0e = _row_ratios(data[6])
    b1e = _row_ratios(data[7])

    t: Symbol = Symbol("t", real=True)
    la: Symbol = Symbol("la", real=True)
    lb: Symbol = Symbol("lb", real=True)

    las = [simplify((1 - la) * a0s[i] + la * a1s[i]) for i in range(3)]
    lae = [simplify((1 - la) * a0e[i] + la * a1e[i]) for i in range(3)]

    lbs = [simplify((1 - lb) * b0s[i] + lb * b1s[i]) for i in range(3)]
    lbe = [simplify((1 - lb) * b0e[i] + lb * b1e[i]) for i in range(3)]

    lla = [simplify((1 - t) * las[i] + t * lae[i]) for i in range(3)]
    llb = [simplify((1 - t) * lbs[i] + t * lbe[i]) for i in range(3)]

    eqs = [Eq(lla[i], llb[i]) for i in range(3)]

    # Try linear solve first; fall back to nonlinear solve if needed.
    sols: List[Dict[Symbol, Expr]] = []
    try:
        lsol = list(linsolve(eqs, (t, la, lb)))
        for s in lsol:
            sols.append({t: simplify(s[0]), la: simplify(s[1]), lb: simplify(s[2])})
    except (TypeError, ValueError, NotImplementedError):
        pass

    if not sols:
        try:
            nsol = list(nonlinsolve(eqs, (t, la, lb)))
            for s in nsol:
                sols.append({t: simplify(s[0]), la: simplify(s[1]), lb: simplify(s[2])})
        except (TypeError, ValueError, NotImplementedError):
            pass

    valid_roots: List[Dict[str, Union[float, str]]] = []

    for sol in sols:
        t_expr = simplify(sol.get(t, t))
        la_expr = simplify(sol.get(la, la))
        lb_expr = simplify(sol.get(lb, lb))

        # Reject complex-valued solutions
        if any(e.has(I) for e in (t_expr, la_expr, lb_expr)):
            continue

        # Case 1: All three are numeric constants
        if all(e.free_symbols == set() for e in (t_expr, la_expr, lb_expr)):
            try:
                t_val = float(t_expr.evalf())
                la_val = float(la_expr.evalf())
                lb_val = float(lb_expr.evalf())
            except (TypeError, ValueError):
                continue
            if 0.0 <= t_val <= 1.0 and 0.0 <= la_val <= 1.0 and 0.0 <= lb_val <= 1.0:
                valid_roots.append({"t": t_val, "a": la_val, "b": lb_val})
            continue

        # Case 2: Parameterized by t. Determine feasible t-set in [0,1].
        cond = And(
            t >= 0,
            t <= 1,
            la_expr >= 0,
            la_expr <= 1,
            lb_expr >= 0,
            lb_expr <= 1,
        )
        try:
            feasible_t: Set = solve_univariate_inequality(cond, t, relational=False)
            # Intersect with [0, 1] explicitly (some solvers already include it)
            feasible_t = feasible_t.intersect(Interval(0, 1))
        except (TypeError, ValueError, NotImplementedError):
            feasible_t = S.EmptySet

        if feasible_t.is_EmptySet:
            continue

        measure = _set_measure(feasible_t)
        # Accept if set has positive measure (interval of non-zero length)
        if measure > 0.0:
            t_pick = _pick_representative_t(feasible_t)
            if t_pick is None:
                continue
            la_val = float(la_expr.subs(t, t_pick).evalf())
            lb_val = float(lb_expr.subs(t, t_pick).evalf())
            if 0.0 <= la_val <= 1.0 and 0.0 <= lb_val <= 1.0:
                valid_roots.append({"t": float(t_pick), "a": la_val, "b": lb_val})
        else:
            # Zero-measure feasible set (e.g., isolated points) -> match Mathematica's
            # behavior where "simpleSol = True" was required for direct acceptance.
            # We won't accept zero-measure sets here.
            pass

    if valid_roots and outfile:
        with open(outfile, "w", encoding="utf-8") as f:
            json.dump(valid_roots, f)

    return len(valid_roots) > 0
