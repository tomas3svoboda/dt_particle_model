"""Chebyshev construction and evaluation primitives for the study lane.

STUDY LANE ONLY.  Nothing here is imported by the certified lane; the module
never executes unless :mod:`dtdc_simulator.core2.study_lane.activation` has
already refused-or-approved an explicit opt-in.

Two operations are needed and both are deliberately dependency-free at
evaluation time (plain floats, no numpy), because evaluation happens inside
the tower's inner loops where a numpy round trip would cost more than the
polynomial itself:

* ``fit_1d`` / ``fit_2d`` — Chebyshev interpolation of a reference callable at
  Chebyshev-Gauss nodes, coefficients by the exact discrete cosine sums.  Used
  once, offline, by ``scripts/build_study_surrogates.py``.
* ``eval_1d`` / ``eval_2d`` — Clenshaw / basis-recurrence evaluation of the
  serialized coefficients.  Used at run time.

Interpolation at the Chebyshev-Gauss nodes ``x_k = cos(pi (k + 1/2) / n)`` is
near-minimax for analytic integrands: the error decays geometrically in the
degree, which is why a degree ~10-16 series reaches the reference solver's own
convergence floor on the smooth saturation and density branches.  Where the
reference is NOT analytic (a bisection's staircase, a Newton guard that trips)
the decay stalls, the measured verification error stays large, and the build
script refuses to adopt that surrogate.  The stall is the detector; nothing
here assumes smoothness.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

__all__ = [
    "eval_1d",
    "eval_2d",
    "fit_1d",
    "fit_2d",
    "gauss_nodes",
]


def gauss_nodes(lo: float, hi: float, count: int) -> tuple[float, ...]:
    """The ``count`` Chebyshev-Gauss nodes of ``[lo, hi]``, ascending."""

    if count < 1:
        raise ValueError("node count must be positive")
    if not (hi > lo):
        raise ValueError(f"degenerate interval [{lo}, {hi}]")
    half_span = 0.5 * (hi - lo)
    centre = 0.5 * (hi + lo)
    nodes = [
        centre + half_span * math.cos(math.pi * (k + 0.5) / count) for k in range(count)
    ]
    nodes.reverse()
    return tuple(nodes)


def fit_1d(
    reference: Callable[[float], float], lo: float, hi: float, degree: int
) -> tuple[float, ...]:
    """Chebyshev coefficients of ``reference`` on ``[lo, hi]`` at ``degree``.

    Returns ``degree + 1`` coefficients in the convention consumed by
    :func:`eval_1d` (``c[0]`` already halved).
    """

    count = degree + 1
    nodes = gauss_nodes(lo, hi, count)
    values = [float(reference(x)) for x in nodes]
    # ``gauss_nodes`` reversed the natural cosine order, so undo that here to
    # keep the transform's index alignment exact.
    values.reverse()
    coefficients = []
    for j in range(count):
        total = math.fsum(
            values[k] * math.cos(math.pi * j * (k + 0.5) / count) for k in range(count)
        )
        coefficients.append(2.0 * total / count)
    coefficients[0] *= 0.5
    return tuple(coefficients)


def eval_1d(coefficients: Sequence[float], lo: float, hi: float, x: float) -> float:
    """Clenshaw evaluation of a 1-D series produced by :func:`fit_1d`."""

    t = (2.0 * x - (hi + lo)) / (hi - lo)
    t2 = t + t
    b1 = 0.0
    b2 = 0.0
    for index in range(len(coefficients) - 1, 0, -1):
        b1, b2 = t2 * b1 - b2 + coefficients[index], b1
    return t * b1 - b2 + coefficients[0]


def fit_2d(
    reference: Callable[[float, float], float],
    x_lo: float,
    x_hi: float,
    y_lo: float,
    y_hi: float,
    x_degree: int,
    y_degree: int,
) -> tuple[tuple[float, ...], ...]:
    """Tensor-product Chebyshev coefficients ``C[i][j]`` on the given box."""

    nx = x_degree + 1
    ny = y_degree + 1
    x_nodes = list(gauss_nodes(x_lo, x_hi, nx))
    y_nodes = list(gauss_nodes(y_lo, y_hi, ny))
    x_nodes.reverse()
    y_nodes.reverse()
    grid = [[float(reference(x, y)) for y in y_nodes] for x in x_nodes]

    # Transform along y first, then along x.
    partial = []
    for row in grid:
        row_coefficients = []
        for j in range(ny):
            total = math.fsum(
                row[m] * math.cos(math.pi * j * (m + 0.5) / ny) for m in range(ny)
            )
            row_coefficients.append(2.0 * total / ny)
        row_coefficients[0] *= 0.5
        partial.append(row_coefficients)

    coefficients = []
    for i in range(nx):
        column = []
        for j in range(ny):
            total = math.fsum(
                partial[k][j] * math.cos(math.pi * i * (k + 0.5) / nx) for k in range(nx)
            )
            column.append(2.0 * total / nx)
        coefficients.append(column)
    coefficients[0] = [value * 0.5 for value in coefficients[0]]
    return tuple(tuple(column) for column in coefficients)


def eval_2d(
    coefficients: Sequence[Sequence[float]],
    x_lo: float,
    x_hi: float,
    y_lo: float,
    y_hi: float,
    x: float,
    y: float,
) -> float:
    """Evaluate a tensor-product series produced by :func:`fit_2d`.

    Row-wise Clenshaw in ``y`` followed by one Clenshaw in ``x``: the same
    arithmetic count as the explicit double sum, with the numerically stable
    recurrence in both directions.
    """

    tx = (2.0 * x - (x_hi + x_lo)) / (x_hi - x_lo)
    ty = (2.0 * y - (y_hi + y_lo)) / (y_hi - y_lo)
    ty2 = ty + ty
    reduced = []
    for row in coefficients:
        b1 = 0.0
        b2 = 0.0
        for index in range(len(row) - 1, 0, -1):
            b1, b2 = ty2 * b1 - b2 + row[index], b1
        reduced.append(ty * b1 - b2 + row[0])
    tx2 = tx + tx
    b1 = 0.0
    b2 = 0.0
    for index in range(len(reduced) - 1, 0, -1):
        b1, b2 = tx2 * b1 - b2 + reduced[index], b1
    return tx * b1 - b2 + reduced[0]
