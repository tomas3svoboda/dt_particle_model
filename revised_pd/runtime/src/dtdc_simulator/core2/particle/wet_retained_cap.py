r"""Exact wet-side PHY-031 retained-cap normal-cone graph.

The liquid-filled wet particle has no pore gas or internal free-water phase.
Its retained-water loading therefore uses the hard-cap inclusion selected by
GT-PS-2-P1E-06,

``0 <= W_cap - W  perpendicular  0 <= lambda``.

``lambda`` is the dimensionless normal-cone contribution to
``mu_ret/(R*T)``.  It has no mass, volume, enthalpy, or internal energy and is
not interpreted as gas activity or pressure.  Storage and calorics use ``W``
only.  The exact graph is parameterized by one semismooth coordinate, so a
capacity-active wet piece does not add an equation, source, or phase.

This module deliberately does not alter the dry-pore cap.  There the explicit
shared pore gas carries excess water while the retained potential remains
pinned at the evidence cap.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Sequence

from dtdc_simulator.core2.props import sorption as sp


class WetRetainedCapGraphError(ValueError):
    """A state or coordinate is not on the exact wet retained-cap graph."""


class WetRetainedWaterBranch(str, Enum):
    """Exact branches of the liquid-filled wet retained-water inclusion."""

    LUIKOV_SMOOTH = "luikov_smooth"
    RETAINED_CAP_NORMAL_CONE = "retained_cap_normal_cone"


class WetRetainedCapTransition(str, Enum):
    """Exact accepted-step transition of one wet retained-water piece."""

    SMOOTH_CONTINUATION = "smooth_continuation"
    CAP_ENTRY = "cap_entry"
    CAP_CONTINUATION = "cap_continuation"
    CAP_EXIT = "cap_exit"


@dataclass(frozen=True)
class WetRetainedWaterGraphPoint:
    """One exact primal/dual point on the wet retained-water graph."""

    retained_water_loading: float
    capacity_dual_over_rt: float
    branch: WetRetainedWaterBranch
    gap_to_capacity: float
    complementarity_product: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.retained_water_loading,
            self.capacity_dual_over_rt,
            self.gap_to_capacity,
            self.complementarity_product,
        )
        if not all(math.isfinite(value) for value in values):
            raise WetRetainedCapGraphError("wet retained-cap point must be finite")
        if not isinstance(self.branch, WetRetainedWaterBranch):
            raise TypeError("wet retained-cap point needs a named graph branch")
        if self.capacity_dual_over_rt < 0.0 or self.gap_to_capacity < 0.0:
            raise WetRetainedCapGraphError("wet retained-cap primal and dual gaps must be nonnegative")
        if self.complementarity_product != 0.0:
            raise WetRetainedCapGraphError(
                "wet retained-cap complementarity must hold exactly, without a tolerance"
            )
        if self.branch is WetRetainedWaterBranch.LUIKOV_SMOOTH:
            if self.gap_to_capacity <= 0.0 or self.capacity_dual_over_rt != 0.0:
                raise WetRetainedCapGraphError(
                    "smooth wet retained water requires positive cap gap and exact zero dual"
                )
        elif self.gap_to_capacity != 0.0:
            raise WetRetainedCapGraphError(
                "capacity-active wet retained water requires exact cap contact"
            )

    @property
    def is_capacity_active(self) -> bool:
        return self.branch is WetRetainedWaterBranch.RETAINED_CAP_NORMAL_CONE


@dataclass(frozen=True)
class WetRetainedCapCertificate:
    """Immutable exact-graph certificate for all positive-volume wet pieces."""

    points: tuple[WetRetainedWaterGraphPoint, ...]
    maximum_loading: float
    maximum_capacity_dual_over_rt: float
    active_piece_indices: tuple[int, ...]
    storage_uses_primal_loading_only: bool = field(default=True, init=False)
    dual_has_mass_volume_or_energy: bool = field(default=False, init=False)
    introduces_free_water_or_other_phase: bool = field(default=False, init=False)
    exact_graph_without_tolerance: bool = field(default=True, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.points:
            raise WetRetainedCapGraphError("wet retained-cap certificate needs at least one piece")
        if not all(isinstance(point, WetRetainedWaterGraphPoint) for point in self.points):
            raise TypeError("wet retained-cap certificate contains an invalid point")
        expected_maximum_loading = max(
            point.retained_water_loading for point in self.points
        )
        expected_maximum_dual = max(
            point.capacity_dual_over_rt for point in self.points
        )
        expected_active = tuple(
            index for index, point in enumerate(self.points) if point.is_capacity_active
        )
        if self.maximum_loading != expected_maximum_loading:
            raise WetRetainedCapGraphError("wet retained-cap maximum loading is inconsistent")
        if self.maximum_capacity_dual_over_rt != expected_maximum_dual:
            raise WetRetainedCapGraphError("wet retained-cap maximum dual is inconsistent")
        if self.active_piece_indices != expected_active:
            raise WetRetainedCapGraphError("wet retained-cap active-piece index is inconsistent")

    @property
    def active_piece_count(self) -> int:
        return len(self.active_piece_indices)

    @property
    def maximum_complementarity_product(self) -> float:
        return max(point.complementarity_product for point in self.points)


@dataclass(frozen=True)
class WetRetainedCapTransitionCertificate:
    """Exact branch-transition certificate across one accepted same-rank step."""

    piece_transitions: tuple[WetRetainedCapTransition, ...]
    entry_piece_indices: tuple[int, ...]
    continuation_piece_indices: tuple[int, ...]
    exit_piece_indices: tuple[int, ...]
    smooth_piece_indices: tuple[int, ...]
    exact_branch_comparisons_without_tolerance: bool = field(
        default=True,
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.piece_transitions:
            raise WetRetainedCapGraphError(
                "wet retained-cap transition certificate needs at least one piece"
            )
        groups = (
            self.entry_piece_indices,
            self.continuation_piece_indices,
            self.exit_piece_indices,
            self.smooth_piece_indices,
        )
        flattened = tuple(index for group in groups for index in group)
        if tuple(sorted(flattened)) != tuple(range(len(self.piece_transitions))):
            raise WetRetainedCapGraphError(
                "wet retained-cap transition groups must partition every piece"
            )
        if len(set(flattened)) != len(flattened):
            raise WetRetainedCapGraphError(
                "wet retained-cap transition groups must be disjoint"
            )
        expected = {
            WetRetainedCapTransition.CAP_ENTRY: self.entry_piece_indices,
            WetRetainedCapTransition.CAP_CONTINUATION: (
                self.continuation_piece_indices
            ),
            WetRetainedCapTransition.CAP_EXIT: self.exit_piece_indices,
            WetRetainedCapTransition.SMOOTH_CONTINUATION: (
                self.smooth_piece_indices
            ),
        }
        for transition, indices in expected.items():
            if indices != tuple(
                index
                for index, value in enumerate(self.piece_transitions)
                if value is transition
            ):
                raise WetRetainedCapGraphError(
                    "wet retained-cap transition index certificate is inconsistent"
                )

    @property
    def entry_piece_count(self) -> int:
        return len(self.entry_piece_indices)

    @property
    def continuation_piece_count(self) -> int:
        return len(self.continuation_piece_indices)

    @property
    def exit_piece_count(self) -> int:
        return len(self.exit_piece_indices)


@dataclass(frozen=True)
class WetRetainedCapSemismoothChart:
    r"""One-coordinate exact parameterization of the primal/dual graph.

    For coordinate ``q < 0``,

    ``W = W_lo + (W_cap-W_lo)*exp(q)`` and ``lambda=0``.

    For ``q >= 0``, ``W=W_cap`` and ``lambda=kappa*q``.  ``kappa`` matches
    the one-sided derivative of the modified-Luikov log-activity potential at
    the cap, which removes a gratuitous potential-slope jump without smoothing
    the storage kink or its active set.
    """

    lower_loading: float
    luikov: sp.LuikovParams = sp.LuikovParams()
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.lower_loading):
            raise WetRetainedCapGraphError("wet retained-cap chart lower loading must be finite")
        if not self.luikov.W_ref <= self.lower_loading < self.luikov.W_cap:
            raise WetRetainedCapGraphError(
                "wet retained-cap chart lower loading must lie inside the frozen authority"
            )

    @property
    def loading_span(self) -> float:
        return self.luikov.W_cap - self.lower_loading

    @property
    def dual_coordinate_scale(self) -> float:
        derivative_at_cap = self.luikov.A1 / (
            self.luikov.A2 * self.luikov.W_cap * self.luikov.W_cap
        )
        scale = self.loading_span * derivative_at_cap
        if not math.isfinite(scale) or scale <= 0.0:  # pragma: no cover - parameter guard
            raise WetRetainedCapGraphError("wet retained-cap dual scale is not positive")
        return scale

    def encode(self, loading: float, capacity_dual_over_rt: float) -> float:
        """Encode one exact graph point without projection or clipping."""

        point = classify_exact_graph_point(
            loading,
            capacity_dual_over_rt,
            self.luikov,
        )
        if point.branch is WetRetainedWaterBranch.LUIKOV_SMOOTH:
            if not self.lower_loading < loading < self.luikov.W_cap:
                raise WetRetainedCapGraphError(
                    "smooth wet retained loading is outside the active solver chart"
                )
            gap_fraction = (self.luikov.W_cap - loading) / self.loading_span
            if gap_fraction <= 0.5:
                # Preserve the adjacent representable smooth state at the cap;
                # forming ``(W-Wlo)/span`` first can round that fraction to 1.
                coordinate = math.log1p(-gap_fraction)
            else:
                coordinate = math.log(
                    (loading - self.lower_loading) / self.loading_span
                )
            if not math.isfinite(coordinate) or not coordinate < 0.0:
                raise WetRetainedCapGraphError(
                    "smooth wet retained state did not encode on the negative chart side"
                )
            return coordinate
        coordinate = capacity_dual_over_rt / self.dual_coordinate_scale
        if not math.isfinite(coordinate) or coordinate < 0.0:
            raise WetRetainedCapGraphError(
                "capacity-active wet retained state did not encode on the nonnegative side"
            )
        return coordinate

    def decode(self, coordinate: float) -> WetRetainedWaterGraphPoint:
        """Decode exactly onto the smooth or active graph, never between them."""

        if not math.isfinite(coordinate):
            raise WetRetainedCapGraphError("wet retained-cap coordinate must be finite")
        if coordinate < 0.0:
            fraction = math.exp(coordinate)
            loading = self.lower_loading + self.loading_span * fraction
            if not self.lower_loading < loading < self.luikov.W_cap:
                raise WetRetainedCapGraphError(
                    "smooth wet retained coordinate collapsed onto a chart endpoint"
                )
            return classify_exact_graph_point(loading, 0.0, self.luikov)
        return classify_exact_graph_point(
            self.luikov.W_cap,
            self.dual_coordinate_scale * coordinate,
            self.luikov,
        )


def classify_exact_graph_point(
    loading: float,
    capacity_dual_over_rt: float,
    luikov: sp.LuikovParams = sp.LuikovParams(),
) -> WetRetainedWaterGraphPoint:
    """Validate and classify one point using exact branch comparisons."""

    if not all(math.isfinite(value) for value in (loading, capacity_dual_over_rt)):
        raise WetRetainedCapGraphError("wet retained-cap loading and dual must be finite")
    if not luikov.W_ref <= loading <= luikov.W_cap:
        raise WetRetainedCapGraphError(
            "wet retained loading left the frozen positive-moisture/cap authority"
        )
    if capacity_dual_over_rt < 0.0:
        raise WetRetainedCapGraphError("wet retained-cap dual must be nonnegative")
    gap = luikov.W_cap - loading
    if loading < luikov.W_cap:
        if capacity_dual_over_rt != 0.0:
            raise WetRetainedCapGraphError(
                "a positive wet retained-cap dual requires exact cap contact"
            )
        branch = WetRetainedWaterBranch.LUIKOV_SMOOTH
    elif loading == luikov.W_cap:
        branch = WetRetainedWaterBranch.RETAINED_CAP_NORMAL_CONE
        gap = 0.0
    else:  # pragma: no cover - exhaustive after interval validation
        raise WetRetainedCapGraphError("unrecognized wet retained-cap branch")
    product = gap * capacity_dual_over_rt
    return WetRetainedWaterGraphPoint(
        retained_water_loading=loading,
        capacity_dual_over_rt=capacity_dual_over_rt,
        branch=branch,
        gap_to_capacity=gap,
        complementarity_product=product,
    )


def certify_exact_graph(
    loadings: Sequence[float],
    capacity_duals_over_rt: Sequence[float],
    luikov: sp.LuikovParams = sp.LuikovParams(),
) -> WetRetainedCapCertificate:
    """Return an immutable certificate or fail on the first off-graph point."""

    loading_tuple = tuple(loadings)
    dual_tuple = tuple(capacity_duals_over_rt)
    if len(loading_tuple) != len(dual_tuple) or not loading_tuple:
        raise WetRetainedCapGraphError(
            "wet retained-cap loadings and duals need one nonempty aligned field"
        )
    points = tuple(
        classify_exact_graph_point(loading, dual, luikov)
        for loading, dual in zip(loading_tuple, dual_tuple)
    )
    return WetRetainedCapCertificate(
        points=points,
        maximum_loading=max(point.retained_water_loading for point in points),
        maximum_capacity_dual_over_rt=max(
            point.capacity_dual_over_rt for point in points
        ),
        active_piece_indices=tuple(
            index for index, point in enumerate(points) if point.is_capacity_active
        ),
    )


def certify_exact_transition(
    before: WetRetainedCapCertificate,
    after: WetRetainedCapCertificate,
) -> WetRetainedCapTransitionCertificate:
    """Certify exact branch changes for one accepted same-rank wet field."""

    if len(before.points) != len(after.points):
        raise WetRetainedCapGraphError(
            "wet retained-cap transition requires unchanged wet-piece rank"
        )
    transitions = []
    for old_point, new_point in zip(before.points, after.points):
        old_active = old_point.is_capacity_active
        new_active = new_point.is_capacity_active
        if not old_active and not new_active:
            transition = WetRetainedCapTransition.SMOOTH_CONTINUATION
        elif not old_active and new_active:
            transition = WetRetainedCapTransition.CAP_ENTRY
        elif old_active and new_active:
            transition = WetRetainedCapTransition.CAP_CONTINUATION
        else:
            transition = WetRetainedCapTransition.CAP_EXIT
        transitions.append(transition)
    piece_transitions = tuple(transitions)
    return WetRetainedCapTransitionCertificate(
        piece_transitions=piece_transitions,
        entry_piece_indices=tuple(
            index
            for index, transition in enumerate(piece_transitions)
            if transition is WetRetainedCapTransition.CAP_ENTRY
        ),
        continuation_piece_indices=tuple(
            index
            for index, transition in enumerate(piece_transitions)
            if transition is WetRetainedCapTransition.CAP_CONTINUATION
        ),
        exit_piece_indices=tuple(
            index
            for index, transition in enumerate(piece_transitions)
            if transition is WetRetainedCapTransition.CAP_EXIT
        ),
        smooth_piece_indices=tuple(
            index
            for index, transition in enumerate(piece_transitions)
            if transition is WetRetainedCapTransition.SMOOTH_CONTINUATION
        ),
    )


def effective_capacity_duals(
    values: Sequence[float],
    expected: int,
) -> tuple[float, ...]:
    """Expand the backward-compatible empty tuple to exact smooth zeros."""

    result = tuple(values)
    if not result:
        return (0.0,) * expected
    if len(result) != expected:
        raise WetRetainedCapGraphError(
            f"wet retained-cap dual field has length {len(result)}, expected {expected}"
        )
    if not all(math.isfinite(value) and value >= 0.0 for value in result):
        raise WetRetainedCapGraphError("wet retained-cap duals must be finite and nonnegative")
    return result


__all__ = [
    "WetRetainedCapCertificate",
    "WetRetainedCapGraphError",
    "WetRetainedCapSemismoothChart",
    "WetRetainedCapTransition",
    "WetRetainedCapTransitionCertificate",
    "WetRetainedWaterBranch",
    "WetRetainedWaterGraphPoint",
    "certify_exact_graph",
    "certify_exact_transition",
    "classify_exact_graph_point",
    "effective_capacity_duals",
]
