from __future__ import annotations

import re
import tarfile
from typing import Any, Dict, Iterable, List, Tuple

from sympy import Interval, RootOf, S, Symbol, sympify
from sympy.solvers.inequalities import solve_univariate_inequality


_ROOT_PATTERN = re.compile(r"Root\[(.*?&),\s*(\d+),\s*0\]", re.DOTALL)
_QNUM_PATTERN = re.compile(r"_q(\d+)_")


def _root_to_sympy_str(expr_str: str) -> str:
    def repl(match: re.Match[str]) -> str:
        poly = match.group(1).strip()
        if poly.endswith("&"):
            poly = poly[:-1].strip()
        poly = poly.replace("#1", "x").replace("^", "**")
        idx = int(match.group(2)) - 1
        return f"RootOf({poly}, {idx})"

    expr_str = expr_str.replace("^", "**").replace("#1", "x")
    return _ROOT_PATTERN.sub(repl, expr_str)


def _to_float(val: Any) -> float | None:
    try:
        from wolframclient.language.expression import WLFunction
    except Exception:  # pragma: no cover - optional dependency already required
        WLFunction = ()  # type: ignore

    if isinstance(val, (int, float)):
        return float(val)

    if isinstance(val, WLFunction):
        try:
            expr = _wl_to_sympy(val)
            return float(expr.evalf())
        except Exception:
            return None

    try:
        return float(val)
    except Exception:
        pass
    try:
        expr = sympify(_root_to_sympy_str(str(val)))
        return float(expr.evalf())
    except Exception:
        return None


def _iter_rules(expr: Any) -> Iterable[Tuple[str, Any]]:
    if isinstance(expr, dict):
        for k, v in expr.items():
            yield str(k), v
        return

    for item in expr:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            yield str(item[0]), item[1]
            continue
        head = getattr(item, "head", None)
        if getattr(head, "name", None) in {"Rule", "RuleDelayed"}:
            args = getattr(item, "args", ())
            if len(args) == 2:
                yield str(args[0]), args[1]


def _wl_to_sympy(obj: Any) -> Any:
    try:
        from wolframclient.language.expression import WLFunction
    except Exception:  # pragma: no cover - optional dependency already required
        WLFunction = ()  # type: ignore

    if isinstance(obj, WLFunction):
        head = str(obj.head)
        args = obj.args
        if head == "Plus":
            return sum(_wl_to_sympy(a) for a in args)
        if head == "Times":
            result = 1
            for a in args:
                result *= _wl_to_sympy(a)
            return result
        if head == "Power":
            base, exp = (_wl_to_sympy(args[0]), _wl_to_sympy(args[1]))
            return base ** exp
        if head == "Slot":
            return Symbol(f"x{args[0]}")
        if head == "Rational":
            num, den = (_wl_to_sympy(args[0]), _wl_to_sympy(args[1]))
            return num / den
        if head == "Function":
            return _wl_to_sympy(args[0])
        if head == "Root":
            poly_expr = _wl_to_sympy(args[0])
            idx = int(args[1]) - 1
            return RootOf(poly_expr, idx)

    if isinstance(obj, (list, tuple)):
        return [_wl_to_sympy(o) for o in obj]

    try:
        from wolframclient.language.expression import WLSymbol
    except Exception:  # pragma: no cover - optional dependency already required
        WLSymbol = ()  # type: ignore

    if isinstance(obj, WLSymbol):
        # `Global`t$3013` -> `t_3013`: sympify parses neither the context
        # backtick nor the `$` Mathematica appends when it scopes a variable.
        name = str(obj).split("`")[-1].replace("$", "_")
        return Symbol(name, real=True)

    return sympify(obj)


def _to_expr(val: Any) -> Any:
    """`val` as a sympy expression, or None if it cannot be interpreted."""
    try:
        return _wl_to_sympy(val)
    except Exception:
        pass
    try:
        return sympify(_root_to_sympy_str(str(val)))
    except Exception:
        return None


def _degenerate_family_toi(raw: Dict[str, Any]) -> float | None:
    """
    Earliest time of impact of a contact that persists over an interval.

    When the vertex stays in the face's plane for the whole step, the contact is
    not a point in time but a continuum, and Mathematica returns a *family* of
    solutions: `t` comes back as the unevaluated `First[False]` while `a` and `b`
    are exact rational-linear functions of the scoped time variable it renamed to
    something like `t$3013`. Nine cloth-funnel cases are of this shape --
    verified against the geometry, where every point of 419vf query 0 sits at
    z = -2 at both ends of the step.

    The root is still fully determined by what is in the file. Contact holds
    where the barycentric coordinates are admissible, so the time of impact is

        inf { t in [0, 1] : a(t) >= 0, b(t) >= 0, a(t) + b(t) <= 1 }

    computed in exact rational arithmetic and only then rounded to double. Return
    None when the set is empty, which would mean the file disagrees with itself;
    fabricating a time of impact there would be worse than reporting nothing.
    """
    a_expr = _to_expr(raw.get("a"))
    b_expr = _to_expr(raw.get("b"))
    if a_expr is None or b_expr is None:
        return None

    free = (a_expr.free_symbols | b_expr.free_symbols)
    if len(free) != 1:
        return None
    var = next(iter(free))

    feasible = Interval(0, 1)
    for constraint in (a_expr >= 0, b_expr >= 0, a_expr + b_expr <= 1):
        try:
            solution = solve_univariate_inequality(
                constraint, var, relational=False, domain=S.Reals)
        except Exception:
            return None
        feasible = feasible.intersect(solution)

    if feasible.is_empty:
        return None
    try:
        # `inf` is the earliest admissible time; rounding it to double is the
        # only inexact step, and float() of a Rational rounds to nearest, which
        # can land a hair past the exact infimum. Nudge to the representable
        # value at or below it -- an earlier reference time can only make the
        # late-time-of-impact gate stricter, never let a violation through.
        import math
        value = float(feasible.inf)
        if value > 0.0 and feasible.inf < value:
            value = math.nextafter(value, 0.0)
        return value
    except Exception:
        return None


def read_wxf_roots(archive_path: str) -> Dict[int, Dict[str, float]]:
    """
    Given a .tar.gz archive containing .wxf Mathematica files, read every file,
    interpret the rules for t, a, and b, and return them as numeric values.
    Requires the `wolframclient` package.
    """
    try:
        from wolframclient.deserializers import binary_deserialize
    except Exception as exc:  # pragma: no cover - dependency is external
        raise ImportError(
            "Reading .wxf files requires the 'wolframclient' package."
        ) from exc

    roots_by_query: Dict[int, Dict[str, float]] = {}
    with tarfile.open(archive_path, mode="r:gz") as tar:
        for member in tar:
            if not (member.isfile() and member.name.endswith(".wxf")):
                continue
            # Extract query number from member name: ..._q<number>_...
            m = _QNUM_PATTERN.search(member.name)
            if not m:
                continue
            qnum = int(m.group(1))
            fileobj = tar.extractfile(member)
            if fileobj is None:
                continue
            data = fileobj.read()
            deserialized = binary_deserialize(data)
            roots_in_member: List[Dict[str, float]] = []
            for root_expr in deserialized:
                values: Dict[str, float] = {}
                raw_values: Dict[str, Any] = {}
                for key, val in _iter_rules(root_expr):
                    key_lower = key.lower()
                    if key_lower in {"t", "a", "b"}:
                        raw_values[key_lower] = val
                        num = _to_float(val)
                        if num is not None:
                            values[key_lower] = num
                if "t" not in values:
                    # A contact that holds over an interval rather than at an
                    # instant; the earliest time of it is recoverable exactly.
                    family_t = _degenerate_family_toi(raw_values)
                    if family_t is not None:
                        values["t"] = family_t
                        values.pop("a", None)
                        values.pop("b", None)
                # `t` is the requirement; `a` and `b` are not.
                #
                # Requiring all three discarded a usable time of impact whenever
                # the two parameter coordinates failed to evaluate, which on
                # cloth-funnel is most of them: every one of 229ee's fourteen
                # queries has a convertible `t`, and thirteen carry an `a` and
                # `b` that arrive as unevaluated WLFunction. The whole root went
                # in the bin for the sake of two fields that
                # benchmark/roots_to_raw.py never reads -- it writes
                # `toi[query] = root["t"]` and nothing else. Across the scene
                # that lost 3,572 of 7,552 roots, 47%, and every one of them
                # then read as "no collision" downstream.
                if "t" in values:
                    roots_in_member.append(values)
            if roots_in_member:
                # Pick earliest t root for this query
                best = min(roots_in_member, key=lambda r: r["t"])
                roots_by_query[qnum] = best
    return roots_by_query


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python sccd_read_roots.py <path_to_wxf_file>")
        sys.exit(1)

    roots = read_wxf_roots(sys.argv[1])
    for q in sorted(roots.keys()):
        root = roots[q]
        print(f'q{q}) {root["t"]}, {root.get("a", "-")}, {root.get("b", "-")}')
        print()
