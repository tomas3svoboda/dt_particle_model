r"""Guarded Gate-1g coefficient selection on the frozen G1G-01 equations.

This module does not replace or diagonalize the coupled water/n-hexane model.
It selects its two Onsager coefficients from one local, isothermal
small-signal constraint.  Around a uniform dry-pore equilibrium with
``y = y_h``, the frozen component fluxes are

``N_w = y_w N_t + J + R`` and ``N_h = y_h N_t - J``,

where ``J = -L_b grad(phi_b)`` is the one binary pore-gas flux and
``R = -L_r grad(phi_r)`` is retained-water transport.  Linearizing both
component balances and eliminating the algebraic total Stefan flux gives

``B_theta d(delta y)/dt = div(J + y_h R)``,

with the projected conserved-storage capacity

``B_theta = y_w dC_h/dy - y_h dC_w/dy > 0``.

Consequently, for negative thermodynamic-force slopes ``phi_b'`` and
``phi_r'``, the eliminated apparent diffusivity is

``D_app = [L_b*(-phi_b') + y_h*L_r*(-phi_r')] / B_theta``.

This is deliberately not a sum of two independently assigned diffusivities.
The storage projection and the ``y_h`` factor are consequences of eliminating
``N_t`` from the same two conservative balances used by G1G-01.  Once the
binary mobility is mapped from a supplied pore-effective Cardarelli-scale
``D_b``, the expression is affine in ``L_r`` and has at most one non-negative
solution.

The storage and force slopes are evaluated as the same symmetric,
topology-preserving secants of the frozen equilibrium functions.  Thus the
implemented constraint is a declared *local discrete small-signal* response,
not a claim that a whole-particle uptake experiment has zero surface
resistance.  It is valid only at fixed ``T`` and ``P``, away from an active-set
kink, with no external film in the reduction.  The retained-cap branch has
``phi_r' = 0`` and is rejected as rank deficient; it is never regularized.

The numerical coefficient brackets below preserve the evidence classes in
the Gate-1g audit.  They are model-form/stress bounds, not identified
soybean-DT parameter bands, so every result remains
``physically_qualifying=False``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import transport_coefficients as tc
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import sorption as sp


BINARY_DIFFUSIVITY_MIN_M2_S = 2.0e-10
BINARY_DIFFUSIVITY_MAX_M2_S = 6.0e-10
APPARENT_WATER_DIFFUSIVITY_MIN_M2_S = 1.55e-11
APPARENT_WATER_DIFFUSIVITY_MAX_M2_S = 6.2e-10
PARTICLE_RADIUS_MIN_M = 0.5e-3
PARTICLE_RADIUS_MAX_M = 1.5e-3
RETAINED_WATER_MIN_KG_KG_DRY = 0.023
RETAINED_WATER_MAX_KG_KG_DRY = 0.2356709040
TEMPERATURE_MIN_K = 323.15
TEMPERATURE_DIRECT_EVIDENCE_MAX_K = 343.15
TEMPERATURE_SENSITIVITY_MAX_K = 433.0


class CoefficientManifoldError(ValueError):
    """Base error for an inadmissible local coefficient selection."""


class CoefficientManifoldStateError(CoefficientManifoldError):
    """The requested state or perturbation left the declared fidelity band."""


class CoefficientManifoldRankError(CoefficientManifoldError):
    """The eliminated local response cannot identify retained mobility."""


class NoFeasibleCoefficientPartitionError(CoefficientManifoldError):
    """No non-negative retained mobility can match the supplied target."""


class ZeroRetainedMobilityEndpointError(CoefficientManifoldError):
    """A zero retained-path endpoint cannot enter the positive solver API."""


@dataclass(frozen=True)
class ApparentWaterDiffusivityTarget:
    """One provenance-bearing apparent-water numerical/model-form target."""

    value_m2_s: float
    provenance: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.value_m2_s) or not (
            APPARENT_WATER_DIFFUSIVITY_MIN_M2_S
            <= self.value_m2_s
            <= APPARENT_WATER_DIFFUSIVITY_MAX_M2_S
        ):
            raise ValueError(
                "apparent-water target must lie in the declared numerical/model-form "
                f"bracket [{APPARENT_WATER_DIFFUSIVITY_MIN_M2_S}, "
                f"{APPARENT_WATER_DIFFUSIVITY_MAX_M2_S}] m2/s"
            )
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("apparent-water target needs explicit provenance")


@dataclass(frozen=True)
class LocalCoefficientManifoldRequest:
    """A phase-local coefficient request and its explicit fidelity coordinates.

    ``composition_half_width`` defines the two states ``y_h +/- width`` used
    by the discrete small-signal reduction.  The selector rejects a requested
    width that crosses a phase or retained-water active-set boundary; it does
    not shrink, clip, or otherwise change the experiment.
    """

    temperature_k: float
    pressure_pa: float
    y_hexane: float
    particle_radius_m: float
    binary_diffusivity: tc.EffectiveBinaryDiffusivitySelection
    apparent_water_diffusivity: ApparentWaterDiffusivityTarget
    pore: cp.CoupledPoreParams = cp.CoupledPoreParams()
    composition_half_width: float = 1.0e-6
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(
            self.binary_diffusivity,
            tc.EffectiveBinaryDiffusivitySelection,
        ):
            raise TypeError("binary diffusivity must be an explicit selection")
        if not isinstance(
            self.apparent_water_diffusivity,
            ApparentWaterDiffusivityTarget,
        ):
            raise TypeError("apparent-water diffusivity must be an explicit target")
        if not isinstance(self.pore, cp.CoupledPoreParams):
            raise TypeError("coefficient-manifold pore data require CoupledPoreParams")

        values = (
            self.temperature_k,
            self.pressure_pa,
            self.y_hexane,
            self.particle_radius_m,
            self.composition_half_width,
            self.pore.w_o,
        )
        if not all(math.isfinite(value) for value in values):
            raise CoefficientManifoldStateError("coefficient-manifold coordinates must be finite")
        if not TEMPERATURE_MIN_K <= self.temperature_k <= TEMPERATURE_SENSITIVITY_MAX_K:
            raise CoefficientManifoldStateError(
                "temperature is outside the 323.15--433 K Gate-1g fidelity band"
            )
        if not bg.PROJECT_PRESSURE_MIN_PA <= self.pressure_pa <= bg.PROJECT_PRESSURE_MAX_PA:
            raise CoefficientManifoldStateError(
                "pressure is outside the 101--170 kPa Gate-1g fidelity band"
            )
        if not PARTICLE_RADIUS_MIN_M <= self.particle_radius_m <= PARTICLE_RADIUS_MAX_M:
            raise CoefficientManifoldStateError(
                "particle radius is outside the 0.5--1.5 mm Gate-1g fidelity band"
            )
        if not sp.W_O_MIN <= self.pore.w_o <= sp.W_O_MAX:
            raise CoefficientManifoldStateError(
                "residual-oil label is outside the 0.012--0.030 Gate-1g fidelity band"
            )
        if not 0.0 < self.y_hexane < 1.0:
            raise CoefficientManifoldStateError(
                "local coefficient selection requires 0 < y_hexane < 1"
            )
        if self.composition_half_width <= 0.0:
            raise CoefficientManifoldStateError(
                "composition half-width must be positive and finite"
            )
        binary_value = self.binary_diffusivity.value_m2_s
        if not BINARY_DIFFUSIVITY_MIN_M2_S <= binary_value <= BINARY_DIFFUSIVITY_MAX_M2_S:
            raise CoefficientManifoldStateError(
                "selected pore-effective binary diffusivity is outside the declared "
                "Cardarelli-scale numerical bracket"
            )

    @property
    def uses_retained_water_temperature_extrapolation(self) -> bool:
        """Whether retained-water closure use exceeds its direct 343.15 K evidence."""

        return self.temperature_k > TEMPERATURE_DIRECT_EVIDENCE_MAX_K


@dataclass(frozen=True)
class EliminatedLocalResponse:
    """Auditable secants and projected capacity after eliminating ``N_t``."""

    water_storage_slope_mol_m3_per_y_hexane: float
    hexane_storage_slope_mol_m3_per_y_hexane: float
    projected_capacity_mol_m3_per_y_hexane: float
    gas_force_slope_per_y_hexane: float
    retained_force_slope_per_y_hexane: float
    binary_only_apparent_diffusivity_m2_s: float
    retained_response_gain_m3_mol: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.water_storage_slope_mol_m3_per_y_hexane,
            self.hexane_storage_slope_mol_m3_per_y_hexane,
            self.projected_capacity_mol_m3_per_y_hexane,
            self.gas_force_slope_per_y_hexane,
            self.retained_force_slope_per_y_hexane,
            self.binary_only_apparent_diffusivity_m2_s,
            self.retained_response_gain_m3_mol,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("eliminated local-response quantities must be finite")
        if self.projected_capacity_mol_m3_per_y_hexane <= 0.0:
            raise ValueError("projected conserved-storage capacity must be positive")
        if self.gas_force_slope_per_y_hexane >= 0.0:
            raise ValueError("binary thermodynamic-force slope must be negative")
        if self.retained_force_slope_per_y_hexane >= 0.0:
            raise ValueError("retained-water thermodynamic-force slope must be negative")
        if self.binary_only_apparent_diffusivity_m2_s <= 0.0:
            raise ValueError("binary-only apparent response must be positive")
        if self.retained_response_gain_m3_mol <= 0.0:
            raise ValueError("retained-mobility response gain must be positive")


@dataclass(frozen=True)
class CoefficientPartitionSelection:
    """The unique local non-negative partition matching one apparent target."""

    request: LocalCoefficientManifoldRequest
    base_state: cp.EquilibriumPoreState
    left_state: cp.EquilibriumPoreState
    right_state: cp.EquilibriumPoreState
    binary_mobility: tc.EvaluatedBinaryMobilitySelection
    retained_water_mobility_mol_m_s: float
    response: EliminatedLocalResponse
    recovered_apparent_diffusivity_m2_s: float
    apparent_diffusivity_residual_m2_s: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.retained_water_mobility_mol_m_s,
            self.recovered_apparent_diffusivity_m2_s,
            self.apparent_diffusivity_residual_m2_s,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("coefficient partition must be finite")
        if self.retained_water_mobility_mol_m_s < 0.0:
            raise ValueError("retained-water mobility must be non-negative")
        target = self.request.apparent_water_diffusivity.value_m2_s
        tolerance = 64.0 * math.ulp(max(abs(target), 1.0e-300))
        if abs(self.apparent_diffusivity_residual_m2_s) > tolerance:
            raise ValueError("coefficient partition does not recover its apparent target")

    @property
    def partition_is_unique(self) -> bool:
        """The positive response gain makes the affine inverse one-to-one."""

        return self.response.retained_response_gain_m3_mol > 0.0

    @property
    def uses_retained_water_temperature_extrapolation(self) -> bool:
        return self.request.uses_retained_water_temperature_extrapolation

    @property
    def retained_water_is_zero_endpoint(self) -> bool:
        return self.retained_water_mobility_mol_m_s == 0.0

    def as_binary_mobility_selection(self) -> cp.MobilitySelection:
        """Return the existing G1G-01 binary-mobility solver adapter."""

        return self.binary_mobility.as_coupled_pore_selection()

    def as_retained_water_mobility_selection(self) -> cp.MobilitySelection:
        """Return a positive G1G-01 retained-mobility solver adapter.

        The existing solver type correctly requires a positive interval.  An
        exact zero endpoint remains an auditable limiting case and is rejected
        here instead of being replaced by a small positive floor.
        """

        value = self.retained_water_mobility_mol_m_s
        if value == 0.0:
            raise ZeroRetainedMobilityEndpointError(
                "zero retained mobility is a limiting endpoint; do not floor it for the solver"
            )
        label = (
            "Gate-1g local eliminated-response partition; "
            f"D_app provenance: {self.request.apparent_water_diffusivity.provenance}; "
            "numerical/model-form case, not soybean-DT identification"
        )
        return cp.MobilityInterval(value, value, label).select(0.0)

    def entropy_production_w_m3_k(
        self,
        gas_force_gradient_m_inv: float,
        retained_force_gradient_m_inv: float,
    ) -> float:
        """Return the two uncoupled Onsager contributions, each counted once."""

        gradients = (gas_force_gradient_m_inv, retained_force_gradient_m_inv)
        if not all(math.isfinite(value) for value in gradients):
            raise ValueError("thermodynamic-force gradients must be finite")
        production = bg.R * math.fsum(
            (
                self.binary_mobility.value_mol_m_s * gas_force_gradient_m_inv**2,
                self.retained_water_mobility_mol_m_s
                * retained_force_gradient_m_inv**2,
            )
        )
        if not math.isfinite(production) or production < 0.0:
            raise RuntimeError("selected coefficient partition violates non-negative dissipation")
        return production


def select_local_coefficient_partition(
    request: LocalCoefficientManifoldRequest,
) -> CoefficientPartitionSelection:
    r"""Infer the unique retained-water mobility at one guarded base state.

    The exact affine inverse is

    ``L_r = (D_app*B_theta - L_b*(-phi_b')) / (y_h*(-phi_r'))``.

    A negative numerator means that the supplied binary path alone is faster
    than the apparent target, so the non-negative feasible set is empty.  A
    zero retained-force slope means the map has lost rank.  Both conditions
    fail closed; no coefficient, composition, or target is clamped.
    """

    if not isinstance(request, LocalCoefficientManifoldRequest):
        raise TypeError("local coefficient selection requires a manifold request")

    interval = cp.gas_only_composition_interval(
        request.temperature_k,
        request.pressure_pa,
        request.pore,
    )
    left_y = request.y_hexane - request.composition_half_width
    right_y = request.y_hexane + request.composition_half_width
    if not (
        interval.lower_y_hexane
        < left_y
        < request.y_hexane
        < right_y
        < interval.upper_y_hexane
    ):
        raise CoefficientManifoldStateError(
            "declared small-signal secant crosses the open gas-only composition boundary; "
            "change the base experiment, do not clip the perturbation"
        )

    base = cp.evaluate_equilibrium(
        request.temperature_k,
        request.pressure_pa,
        request.y_hexane,
        request.pore,
    )
    left = cp.evaluate_equilibrium(
        request.temperature_k,
        request.pressure_pa,
        left_y,
        request.pore,
    )
    right = cp.evaluate_equilibrium(
        request.temperature_k,
        request.pressure_pa,
        right_y,
        request.pore,
    )
    for state in (left, base, right):
        _validate_retained_water_loading(state.retained_water_loading)

    active_sets = {
        left.retained_water_active_set,
        base.retained_water_active_set,
        right.retained_water_active_set,
    }
    if active_sets != {cp.RetainedWaterActiveSet.LUIKOV_SMOOTH}:
        raise CoefficientManifoldRankError(
            "retained-cap or active-set-crossing response has no unique retained-mobility "
            "inverse; preserve the named active set and do not regularize it"
        )

    delta_y = right_y - left_y
    force_secants = cp.isothermal_potential_composition_secant_slopes(
        request.temperature_k,
        request.pressure_pa,
        left_y,
        delta_y,
        request.pore,
    )
    if (
        force_secants.crosses_retained_cap
        or force_secants.retained_water_secant_branch
        is not cp.RetainedWaterSecantBranch.LUIKOV_SMOOTH
    ):
        raise CoefficientManifoldRankError(
            "small-signal force secant is not wholly on the smooth retained-water branch"
        )

    water_storage_slope = math.fsum(
        (
            right.total_water_concentration_mol_m3,
            -left.total_water_concentration_mol_m3,
        )
    ) / delta_y
    hexane_storage_slope = math.fsum(
        (
            right.total_hexane_concentration_mol_m3,
            -left.total_hexane_concentration_mol_m3,
        )
    ) / delta_y
    projected_capacity = math.fsum(
        (
            base.y_water * hexane_storage_slope,
            -base.y_hexane * water_storage_slope,
        )
    )
    if not math.isfinite(projected_capacity) or projected_capacity <= 0.0:
        raise CoefficientManifoldRankError(
            "projected component-storage Jacobian is non-positive or rank deficient"
        )

    binary_face = tc.symmetric_binary_face_state(
        left.binary_gas.molar_density_mol_m3,
        left.y_water,
        left.y_hexane,
        right.binary_gas.molar_density_mol_m3,
        right.y_water,
        right.y_hexane,
    )
    binary_mobility = tc.evaluate_binary_mobility_interval(
        request.binary_diffusivity.interval,
        binary_face,
    ).select(request.binary_diffusivity.fraction)
    gas_force_gain = -force_secants.gas_exchange_potential_secant_per_y_hexane
    retained_force_gain = -force_secants.retained_water_potential_secant_per_y_hexane
    retained_response_denominator = base.y_hexane * retained_force_gain
    if (
        not math.isfinite(retained_response_denominator)
        or retained_response_denominator <= 0.0
    ):
        raise CoefficientManifoldRankError(
            "retained-water force map has zero or negative rank at this state"
        )

    binary_conductance = binary_mobility.value_mol_m_s * gas_force_gain
    binary_only_apparent = binary_conductance / projected_capacity
    retained_response_gain = retained_response_denominator / projected_capacity
    target = request.apparent_water_diffusivity.value_m2_s
    retained_numerator = math.fsum(
        (
            target * projected_capacity,
            -binary_conductance,
        )
    )
    if not math.isfinite(retained_numerator):
        raise CoefficientManifoldRankError("coefficient inverse produced a non-finite numerator")
    if retained_numerator < 0.0:
        raise NoFeasibleCoefficientPartitionError(
            "binary pore transport alone exceeds the apparent-water target; "
            "the non-negative retained-mobility feasible set is empty"
        )
    retained_mobility = retained_numerator / retained_response_denominator
    if not math.isfinite(retained_mobility) or retained_mobility < 0.0:
        raise NoFeasibleCoefficientPartitionError(
            "no finite non-negative retained-water mobility matches the apparent target"
        )

    recovered = math.fsum(
        (
            binary_conductance,
            retained_response_denominator * retained_mobility,
        )
    ) / projected_capacity
    residual = math.fsum((recovered, -target))
    response = EliminatedLocalResponse(
        water_storage_slope_mol_m3_per_y_hexane=water_storage_slope,
        hexane_storage_slope_mol_m3_per_y_hexane=hexane_storage_slope,
        projected_capacity_mol_m3_per_y_hexane=projected_capacity,
        gas_force_slope_per_y_hexane=(
            force_secants.gas_exchange_potential_secant_per_y_hexane
        ),
        retained_force_slope_per_y_hexane=(
            force_secants.retained_water_potential_secant_per_y_hexane
        ),
        binary_only_apparent_diffusivity_m2_s=binary_only_apparent,
        retained_response_gain_m3_mol=retained_response_gain,
    )
    return CoefficientPartitionSelection(
        request=request,
        base_state=base,
        left_state=left,
        right_state=right,
        binary_mobility=binary_mobility,
        retained_water_mobility_mol_m_s=retained_mobility,
        response=response,
        recovered_apparent_diffusivity_m2_s=recovered,
        apparent_diffusivity_residual_m2_s=residual,
    )


def _validate_retained_water_loading(loading_kg_kg_dry: float) -> None:
    if not math.isfinite(loading_kg_kg_dry) or not (
        RETAINED_WATER_MIN_KG_KG_DRY
        <= loading_kg_kg_dry
        <= RETAINED_WATER_MAX_KG_KG_DRY
    ):
        raise CoefficientManifoldStateError(
            "retained-water loading is outside the 0.023--0.2356709040 kg/kg-dry "
            "Gate-1g fidelity band; reject rather than clamp"
        )


__all__ = [
    "APPARENT_WATER_DIFFUSIVITY_MAX_M2_S",
    "APPARENT_WATER_DIFFUSIVITY_MIN_M2_S",
    "ApparentWaterDiffusivityTarget",
    "BINARY_DIFFUSIVITY_MAX_M2_S",
    "BINARY_DIFFUSIVITY_MIN_M2_S",
    "CoefficientManifoldError",
    "CoefficientManifoldRankError",
    "CoefficientManifoldStateError",
    "CoefficientPartitionSelection",
    "EliminatedLocalResponse",
    "LocalCoefficientManifoldRequest",
    "NoFeasibleCoefficientPartitionError",
    "PARTICLE_RADIUS_MAX_M",
    "PARTICLE_RADIUS_MIN_M",
    "RETAINED_WATER_MAX_KG_KG_DRY",
    "RETAINED_WATER_MIN_KG_KG_DRY",
    "TEMPERATURE_DIRECT_EVIDENCE_MAX_K",
    "TEMPERATURE_MIN_K",
    "TEMPERATURE_SENSITIVITY_MAX_K",
    "ZeroRetainedMobilityEndpointError",
    "select_local_coefficient_partition",
]
