r"""Unit-invariant local Jacobian certification for Core2 particle solves.

Raw Jacobian condition numbers depend on the units chosen for residuals and
physical variables.  This module instead certifies the explicitly
dimensionless matrix

``J_star = diag(1 / residual_scale) @ J_physical @ diag(variable_scale)``.

If residual values are converted by positive factors ``a_i`` and variable
values by positive factors ``b_j``, then ``J_physical`` changes as
``diag(a) @ J_physical @ diag(1 / b)``.  Converting the characteristic scales
by the same factors leaves ``J_star`` unchanged.  The reported singular
values, rank, and 2-norm condition number are consequently independent of
those unit choices.

The utility is diagnostic only.  It neither changes nonlinear residuals nor
relaxes an acceptance tolerance.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


class JacobianCertificateError(ValueError):
    """A physical Jacobian cannot receive a full-rank finite certificate."""


@dataclass(frozen=True)
class UnitInvariantJacobianCertificate:
    """Full-rank SVD certificate for one dimensionless square Jacobian."""

    dimensionless_jacobian: tuple[tuple[float, ...], ...]
    singular_values: tuple[float, ...]
    numerical_rank: int
    expected_rank: int
    rank_relative_tolerance: float
    rank_absolute_threshold: float
    condition_number_2: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def full_rank(self) -> bool:
        return self.numerical_rank == self.expected_rank


def certify_unit_invariant_jacobian(
    physical_jacobian: Sequence[Sequence[float]] | np.ndarray,
    residual_characteristic_scales: Sequence[float] | np.ndarray,
    physical_variable_scales: Sequence[float] | np.ndarray,
    *,
    rank_relative_tolerance: float = 1.0e-10,
) -> UnitInvariantJacobianCertificate:
    r"""Dimensionless and certify one physical square Jacobian.

    ``residual_characteristic_scales[i]`` has the units of residual row ``i``.
    ``physical_variable_scales[j]`` has the units of physical variable column
    ``j``.  Both sets of scales are part of the model's declared numerical
    contract; they are not inferred from the raw matrix.

    Rank is evaluated against
    ``rank_relative_tolerance * largest_singular_value``.  A deficient matrix
    fails closed rather than returning an infinite or misleading condition
    number.
    """

    if (
        not isinstance(rank_relative_tolerance, (int, float))
        or isinstance(rank_relative_tolerance, bool)
        or not math.isfinite(rank_relative_tolerance)
        or not 0.0 < rank_relative_tolerance < 1.0
    ):
        raise JacobianCertificateError(
            "rank relative tolerance must be finite and lie strictly in (0, 1)"
        )

    try:
        matrix = np.asarray(physical_jacobian, dtype=float)
    except (TypeError, ValueError) as exc:
        raise JacobianCertificateError(
            "physical Jacobian must be a rectangular array of real values"
        ) from exc
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[0] != matrix.shape[1]:
        raise JacobianCertificateError(
            "physical Jacobian must be a non-empty square matrix"
        )
    if not np.all(np.isfinite(matrix)):
        raise JacobianCertificateError("physical Jacobian must contain only finite values")

    dimension = matrix.shape[0]
    residual_scales = _positive_scale_vector(
        residual_characteristic_scales,
        dimension,
        "residual characteristic scales",
    )
    variable_scales = _positive_scale_vector(
        physical_variable_scales,
        dimension,
        "physical-variable scales",
    )

    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        dimensionless = (matrix / residual_scales[:, None]) * variable_scales[None, :]
    if not np.all(np.isfinite(dimensionless)):
        raise JacobianCertificateError(
            "declared scales overflowed the dimensionless Jacobian"
        )
    try:
        singular_values_array = np.linalg.svd(dimensionless, compute_uv=False)
    except np.linalg.LinAlgError as exc:
        raise JacobianCertificateError(
            "dimensionless Jacobian SVD did not converge"
        ) from exc
    if (
        singular_values_array.shape != (dimension,)
        or not np.all(np.isfinite(singular_values_array))
        or singular_values_array[0] <= 0.0
    ):
        raise JacobianCertificateError(
            "dimensionless Jacobian has no finite positive singular spectrum"
        )

    threshold = float(rank_relative_tolerance * singular_values_array[0])
    numerical_rank = int(np.count_nonzero(singular_values_array > threshold))
    if numerical_rank != dimension:
        raise JacobianCertificateError(
            "dimensionless Jacobian is rank deficient under the declared tolerance: "
            f"rank {numerical_rank} != {dimension}"
        )
    smallest = float(singular_values_array[-1])
    condition = float(singular_values_array[0] / smallest)
    if not math.isfinite(condition) or condition < 1.0:
        raise JacobianCertificateError(
            "dimensionless Jacobian condition number is not finite and valid"
        )

    return UnitInvariantJacobianCertificate(
        dimensionless_jacobian=tuple(
            tuple(float(value) for value in row) for row in dimensionless
        ),
        singular_values=tuple(float(value) for value in singular_values_array),
        numerical_rank=numerical_rank,
        expected_rank=dimension,
        rank_relative_tolerance=float(rank_relative_tolerance),
        rank_absolute_threshold=threshold,
        condition_number_2=condition,
    )


def _positive_scale_vector(
    values: Sequence[float] | np.ndarray,
    expected_length: int,
    name: str,
) -> np.ndarray:
    try:
        vector = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise JacobianCertificateError(f"{name} must be a real vector") from exc
    if vector.ndim != 1 or vector.shape != (expected_length,):
        raise JacobianCertificateError(
            f"{name} must have length {expected_length}"
        )
    if not np.all(np.isfinite(vector)) or np.any(vector <= 0.0):
        raise JacobianCertificateError(
            f"{name} must contain only positive finite values"
        )
    return vector


__all__ = [
    "JacobianCertificateError",
    "UnitInvariantJacobianCertificate",
    "certify_unit_invariant_jacobian",
]
