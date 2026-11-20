"""Nonlinear-interface damping and Anderson mixing utilities."""

from __future__ import annotations

import math

import numpy as np


def _adapt_relaxation(
    value: float,
    residual: float,
    previous_residual: float,
    *,
    minimum: float,
    maximum: float,
) -> float:
    """Safeguarded residual-monotonic damping update for the outer coupling."""
    if not math.isfinite(residual) or residual > 1.20 * previous_residual:
        return max(minimum, 0.7 * value)
    if residual < 0.98 * previous_residual:
        return min(maximum, 1.25 * value)
    if residual < 0.995 * previous_residual:
        return min(maximum, 1.10 * value)
    return value


class _AndersonInterfaceMixer:
    """Safeguarded Type-II Anderson mixing for the small DT interface map.

    The resolved outer problem has only four iterated boundary variables:
    solid temperature into DCZ and water flow, hexane flow, and temperature
    into FTRZ. Independent scalar damping is retained as the fallback. Once
    two map evaluations are available, this mixer forms a multisecant inverse
    of the coupled residual and proposes a common accelerated update.

    The proposal is accepted only when its linear residual prediction improves
    the current residual, has adequate numerical rank, and stays inside a
    componentwise trust region around the proven damped step. Thus acceleration
    changes the nonlinear path only; the physical map and convergence test are
    untouched.
    """

    def __init__(self, *, depth: int = 4, trust_factor: float = 100.0) -> None:
        if depth < 1 or trust_factor <= 1.0:
            raise ValueError("invalid Anderson mixer controls")
        self.depth = depth
        self.trust_factor = trust_factor
        self._states: list[np.ndarray] = []
        self._residuals: list[np.ndarray] = []
        self._previous_norm = math.inf
        self._cooldown = 0

    def clear(self, *, cooldown: int = 0) -> None:
        self._states.clear()
        self._residuals.clear()
        self._previous_norm = math.inf
        self._cooldown = max(cooldown, 0)

    def propose(
        self,
        state: np.ndarray,
        residual: np.ndarray,
        damped_candidate: np.ndarray,
        *,
        tolerance: float,
    ) -> tuple[np.ndarray, bool]:
        state = np.asarray(state, dtype=float)
        residual = np.asarray(residual, dtype=float)
        damped_candidate = np.asarray(damped_candidate, dtype=float)
        if (
            state.shape != residual.shape
            or state.shape != damped_candidate.shape
            or state.ndim != 1
            or not np.all(np.isfinite(state))
            or not np.all(np.isfinite(residual))
            or not np.all(np.isfinite(damped_candidate))
            or tolerance <= 0.0
        ):
            raise ValueError("invalid Anderson mixer state")

        residual_norm = float(np.max(np.abs(residual)))
        # A discontinuous active-set change invalidates the local secant
        # history. Keep the current point as the first sample of a new model.
        if residual_norm > 1.5 * self._previous_norm:
            self.clear(cooldown=2)
        self._states.append(state.copy())
        self._residuals.append(residual.copy())
        if len(self._states) > self.depth + 1:
            self._states.pop(0)
            self._residuals.pop(0)
        self._previous_norm = residual_norm

        if self._cooldown > 0:
            self._cooldown -= 1
            return damped_candidate, False
        if len(self._states) < 3 or residual_norm <= tolerance:
            return damped_candidate, False

        delta_state = np.column_stack(
            [
                self._states[index + 1] - self._states[index]
                for index in range(len(self._states) - 1)
            ]
        )
        delta_residual = np.column_stack(
            [
                self._residuals[index + 1] - self._residuals[index]
                for index in range(len(self._residuals) - 1)
            ]
        )
        try:
            coefficients, _, rank, singular_values = np.linalg.lstsq(
                delta_residual,
                residual,
                rcond=1.0e-10,
            )
        except np.linalg.LinAlgError:
            return damped_candidate, False
        if (
            rank == 0
            or singular_values.size == 0
            or not np.all(np.isfinite(coefficients))
        ):
            return damped_candidate, False

        predicted_residual = residual - delta_residual @ coefficients
        if (
            float(np.max(np.abs(predicted_residual)))
            >= 0.75 * residual_norm
        ):
            return damped_candidate, False

        # Type-II Anderson update with beta=1:
        # x+ = x + f - (DeltaX + DeltaF) gamma.
        accelerated_step = (
            residual
            - (delta_state + delta_residual) @ coefficients
        )
        if not np.all(np.isfinite(accelerated_step)):
            return damped_candidate, False

        damped_step = damped_candidate - state
        trust_limit = self.trust_factor * np.maximum(
            np.abs(damped_step),
            tolerance,
        )
        nonzero = np.abs(accelerated_step) > 0.0
        scale = 1.0
        if np.any(nonzero):
            scale = min(
                1.0,
                float(
                    np.min(
                        trust_limit[nonzero]
                        / np.abs(accelerated_step[nonzero])
                    )
                ),
            )
        if scale < 0.02:
            return damped_candidate, False
        return state + scale * accelerated_step, True
