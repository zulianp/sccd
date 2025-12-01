from __future__ import annotations

import re
import tarfile
from typing import Any, Dict, Iterable, List, Tuple

from sympy import RootOf, Symbol, sympify


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

    return sympify(obj)


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
                for key, val in _iter_rules(root_expr):
                    key_lower = key.lower()
                    if key_lower in {"t", "a", "b"}:
                        num = _to_float(val)
                        if num is not None:
                            values[key_lower] = num
                if {"t", "a", "b"} <= set(values):
                    roots_in_member.append(values)
            if roots_in_member:
                # Pick earliest t root for this query
                best = min(roots_in_member, key=lambda r: r["t"])
                roots_by_query[qnum] = best
    return roots_by_query


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python read_wxf.py <path_to_wxf_file>")
        sys.exit(1)

    roots = read_wxf_roots(sys.argv[1])
    for q in sorted(roots.keys()):
        root = roots[q]
        print(f'q{q}) {root["t"]}, {root["a"]}, {root["b"]}')
        print()
