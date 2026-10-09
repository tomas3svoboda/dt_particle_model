r"""Square cut-cell ALE/RH residual foundation for Gate 1g.

This module is the first simultaneous wet/front/dry numerical assembly on the
fixed material grid.  It deliberately stops at a residual foundation: a later
driver may put a safeguarded nonlinear/event solver around
:func:`assemble_backward_euler`, but this module does not claim a converged
trajectory or physical coefficient qualification.

For a front strictly inside one master cell, positive wet pieces carry
``(T, W)``, positive dry pieces carry ``(T, y_h)``, every dry face carries an
algebraic total Stefan flux ``N_t``, and the front carries ``(z, T_Gamma)``.
The rank is therefore

``2*n_w + 2*n_d + (n_d + 1) + 2 = 2*n_w + 3*n_d + 3``.

It is matched by wet water/energy balances, dry water/hexane/energy balances,
and the three water/hexane/common-datum-energy Rankine--Hugoniot equations.
The same exact swept-volume scalar from :mod:`cut_geometry` is used in all
regional and jump residuals.  No porosity/accessibility fraction is used as a
geometric cut fraction.

Current scope is one within-cell primary-drainage step.  A proposal reaching
or crossing a material face raises :class:`CutTransportStepError` with the
exact input state as rollback and explicitly requests an event restart.  It
never creates an epsilon piece or silently changes equation rank.

The interface n-hexane composition is obtained from the frozen pure-liquid /
binary-gas fugacity equality at ``T_Gamma``.  That saturated ``a_h=1`` trace
is an interface topology and has its own caller-declared root bracket; it is
not forced into the open, gas-only dry-cell primitive band.  Water activity
and the retained-water active set are nevertheless evaluated and validated by
the same coupled-pore authority.

The default mass-force discretization is the explicitly labelled
negligible-Soret branch.  Internal and surface faces use one symmetric
isothermal force temperature.  At the moving RH face, both the saturated
trace and first dry composition are instead evaluated at the physical
``T_Gamma``.  This one-sided isothermal secant is the finite-volume boundary
counterpart of the exact-face tangent: it neither evaluates a saturated trace
at a colder artificial temperature nor constrains the sign of the physical
temperature gradient.  Conduction continues to use the endpoint-temperature
gradient.  An immutable, opt-in complete-potential reference-gauge
counterfactual evaluates endpoint ``mu/(R*T)`` terms from the same frozen
EOS/fugacity authorities.  It changes no balance, mobility, energy-flux,
grid, bound, or solver tolerance and is not an identified Soret coefficient.

The PHY-031 retained-cap branch fixes the retained-matrix potential at the
cap fugacity while leaving the actual pore-gas water activity unchanged above
the evidence endpoint.  It is not an activity clamp.  A cap-to-cap face has
zero retained-force flux but retains its positive binary-gas mobility; a
finite cap/Luikov crossing uses the exact continuous piecewise secant.  The
nonunique derivative at the exact kink is never smoothed.

The n-hexane front is not a water phase boundary.  The retained meal matrix is
continuous, so the wet-side retained-water trace uses the dry-side
local-equilibrium retained fugacity/loading at the shared ``T_Gamma``.  This
is the frozen continuous-matrix/local-equilibrium assumption that closes the
wet water face; it is intentionally visible in :func:`evaluate_interface_state`
rather than hidden as an incidental numerical boundary value.

That trace closes the wet water FACE only.  Since the O9a front-donor ruling
(docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md, ruled
2026-08-21) it is NOT the ALE mass donor: the receding front consumes wet
material at that material's own bulk retained-water loading, and the excess
over the dry-side equilibrium leaves through the unchanged Rankine--Hugoniot
water jump.  Booking the sweep at the trace had implicitly assumed the
arriving retained water pre-drains to the dry-pore equilibrium loading, which
the measured Peclet number 3.4e3 forbids; the correction removes that implicit
assumption rather than adding one.  See :func:`_residual_blocks`.

Nor is it the ALE ENERGY donor.  Since the O10a completion
(docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md, ruled
2026-08-21) the swept material's enthalpy is convected at the SAME (bulk,
old-time) material state as its mass, because the two-phase moving-boundary
Rankine--Hugoniot jump conditions must convect mass and enthalpy of the same
material state - an internal-consistency requirement of the balance laws - and
that state is the bulk, undrained one the O9 traverse measured (Peclet 3.4e3).
The acceptance instrument is gauge covariance: with both donors on the same
state, a shift of the water caloric datum moves ``rh_energy`` by exactly
``cw * rh_water`` again.

The outer Dirichlet trace is only a numerical oracle.  Coefficients and
feasibility bands are caller supplied, no state is clipped, there is no
background reservoir, no Faner ``s(X)``, no fitted coefficient, no contact
fallback, and no added latent source.  Every public result reports
``physically_qualifying=False``.

An optional prescribed-pressure target evaluates old storage with the input
configuration at ``P_n`` and every candidate storage, interface, flux, and RH
term with a pressure-only target configuration at ``P_(n+1)``.  It is a
finite backward-Euler transition, not an instantaneous remap.  The target is
carried by the residual assembly for commitment only after the safeguarded
integrator accepts all nonlinear and conservation ledgers.  This numerical
foundation remains nonphysical until finite-film and bed-hydraulic pressure
authorities are connected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from enum import Enum
from functools import lru_cache
from typing import Sequence

from dtdc_simulator.core2.exact_cache import float_bits_key
from dtdc_simulator.core2.particle import actual_composition_force_path as acfp
from dtdc_simulator.core2.particle import condensed_solvent_reduced_film as csrf
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import front
from dtdc_simulator.core2.particle import nonisothermal_potential as nip
from dtdc_simulator.core2.particle import transport_coefficients as tc
from dtdc_simulator.core2.particle import wet_retained_cap as wrc
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa


class CutTransportTopologyError(ValueError):
    """The candidate does not belong to this fixed-rank partial-front chart."""


class DryThermodynamicForceReconstruction(str, Enum):
    """Auditable temperature rule for one isothermal dry-face force secant."""

    SYMMETRIC_ENDPOINT_TEMPERATURE = "symmetric_endpoint_temperature"
    MOVING_INTERFACE_T_GAMMA = "moving_interface_t_gamma"
    MOVING_INTERFACE_EXACT_TRACE_SECANT = "moving_interface_exact_trace_secant"
    MOVING_INTERFACE_ACTUAL_T_Y_PATH = "moving_interface_actual_t_y_path"
    COMPLETE_ENDPOINT_POTENTIALS = "complete_endpoint_potentials"
    MOVING_INTERFACE_COMPLETE_ENDPOINT_POTENTIALS = "moving_interface_complete_endpoint_potentials"
    MOVING_INTERFACE_COMPLETE_ENDPOINT_POTENTIALS_ACTUAL_STATES = (
        "moving_interface_complete_endpoint_potentials_actual_states"
    )
    REFERENCE_INVARIANT_BINARY_THERMAL_DIFFUSION_FACTOR = (
        "reference_invariant_binary_thermal_diffusion_factor"
    )
    REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS = "reference_invariant_thermal_force_factors"
    CENTER_SYMMETRY = "center_symmetry"


class MovingInterfaceCompositionForceAuthority(str, Enum):
    """Explicit moving-interface state/force reconstruction authority."""

    LEGACY_COMMON_T_GAMMA = "legacy_common_t_gamma"
    ACTUAL_LINEAR_T_Y_PATH = "actual_linear_t_y_path"


class VanishingWetCutThresholdAuthority(str, Enum):
    """Prospectively frozen Amendment-12 cut-cell threshold authorities."""

    ORIGINAL_SAME_CELL_UNSTABILIZED = "original_same_cell_unstabilized"
    P1EF_PRODUCTION_HALF_VOLUME = "p1ef_production_half_volume"
    P1EF_SENSITIVITY_QUARTER_VOLUME_NONPRODUCTION = "p1ef_sensitivity_quarter_volume_nonproduction"

    @property
    def target_wet_volume_fraction(self) -> float | None:
        if self is VanishingWetCutThresholdAuthority.ORIGINAL_SAME_CELL_UNSTABILIZED:
            return None
        if self is VanishingWetCutThresholdAuthority.P1EF_PRODUCTION_HALF_VOLUME:
            return 0.5
        if self is (
            VanishingWetCutThresholdAuthority.P1EF_SENSITIVITY_QUARTER_VOLUME_NONPRODUCTION
        ):
            return 0.25
        raise RuntimeError("unknown vanishing-wet-cut threshold authority")


P1EF_PRODUCTION_VANISHING_WET_CUT_AUTHORITY = (
    VanishingWetCutThresholdAuthority.P1EF_PRODUCTION_HALF_VOLUME
)
ORIGINAL_SAME_CELL_VANISHING_WET_CUT_AUTHORITY = (
    VanishingWetCutThresholdAuthority.ORIGINAL_SAME_CELL_UNSTABILIZED
)
P1EF_SENSITIVITY_VANISHING_WET_CUT_AUTHORITY = (
    VanishingWetCutThresholdAuthority.P1EF_SENSITIVITY_QUARTER_VOLUME_NONPRODUCTION
)


def vanishing_wet_cut_weight(
    wet_volume_fraction: float,
    authority: VanishingWetCutThresholdAuthority = (ORIGINAL_SAME_CELL_VANISHING_WET_CUT_AUTHORITY),
) -> float:
    """Return the frozen C1 Amendment-12 geometric reconstruction weight."""

    if not math.isfinite(wet_volume_fraction) or not 0.0 <= wet_volume_fraction <= 1.0:
        raise ValueError("wet cut-cell volume fraction must lie in [0,1]")
    if not isinstance(authority, VanishingWetCutThresholdAuthority):
        raise TypeError("vanishing-wet-cut treatment needs a frozen threshold authority")
    target = authority.target_wet_volume_fraction
    if target is None:
        return 0.0
    if wet_volume_fraction >= target:
        return 0.0
    coordinate = wet_volume_fraction / target
    # Algebraically identical to 1-3*x**2+2*x**3, but the factored spelling
    # preserves the declared quadratic shutoff instead of cancellation to an
    # O(eps) false weight when x is one representable bit below one.
    return (1.0 - coordinate) ** 2 * (1.0 + 2.0 * coordinate)


class CutTransportStepError(RuntimeError):
    """Rejected residual evaluation carrying the exact immutable rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: "CutTransportState",
        *,
        event_restart_required: bool = False,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.event_restart_required = event_restart_required


@dataclass(frozen=True)
class InterfaceCompositionBracket:
    """Caller-declared bracket for the saturated interface ``y_h`` root."""

    y_hexane_bounds: tuple[float, float]
    label: str
    log_fugacity_tolerance: float = 2.0e-11
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        lower, upper = self.y_hexane_bounds
        if not all(math.isfinite(value) for value in (lower, upper)):
            raise ValueError("interface composition bracket must be finite")
        if not 0.0 < lower < upper < 1.0:
            raise ValueError("interface composition bracket must lie inside (0,1)")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("interface composition bracket needs a provenance label")
        if (
            not math.isfinite(self.log_fugacity_tolerance)
            or not 0.0 < self.log_fugacity_tolerance <= 1.0e-10
        ):
            raise ValueError("interface root tolerance must be in (0, 1e-10]")


@dataclass(frozen=True)
class CutTransportConfig:
    """Shared wet/dry authorities and explicit numerical coefficient choices."""

    dry: ct.FullyDryTransportConfig
    wet: ww.WetWaterModel
    interface_composition: InterfaceCompositionBracket
    mass_force_mode: nip.NonisothermalMassForceMode = (
        nip.NonisothermalMassForceMode.COMMON_FACE_TEMPERATURE_NEGLIGIBLE_SORET
    )
    complete_potential_reference_gauge: nip.CompletePotentialReferenceGauge = (
        nip.DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE
    )
    binary_thermal_diffusion_factor: nip.BinaryThermalDiffusionFactorSelection = (
        nip.ZERO_BINARY_THERMAL_DIFFUSION_FACTOR
    )
    retained_water_thermal_force_factor: nip.RetainedWaterThermalForceFactorSelection = (
        nip.ZERO_RETAINED_WATER_THERMAL_FORCE_FACTOR
    )
    moving_interface_composition_force_authority: MovingInterfaceCompositionForceAuthority = (
        MovingInterfaceCompositionForceAuthority.LEGACY_COMMON_T_GAMMA
    )
    vanishing_wet_cut_threshold_authority: VanishingWetCutThresholdAuthority = (
        ORIGINAL_SAME_CELL_VANISHING_WET_CUT_AUTHORITY
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.dry, ct.FullyDryTransportConfig):
            raise TypeError("cut transport requires the bounded dry configuration")
        if not isinstance(self.wet, ww.WetWaterModel):
            raise TypeError("cut transport requires the wet-water model")
        if not isinstance(self.interface_composition, InterfaceCompositionBracket):
            raise TypeError("cut transport requires an interface composition bracket")
        if not isinstance(self.mass_force_mode, nip.NonisothermalMassForceMode):
            raise TypeError("cut transport requires an immutable mass-force mode")
        if not isinstance(
            self.complete_potential_reference_gauge,
            nip.CompletePotentialReferenceGauge,
        ):
            raise TypeError("cut transport requires an explicit complete-potential gauge")
        if not isinstance(
            self.binary_thermal_diffusion_factor,
            nip.BinaryThermalDiffusionFactorSelection,
        ):
            raise TypeError("cut transport requires an explicit binary alpha_T selection")
        if not isinstance(
            self.retained_water_thermal_force_factor,
            nip.RetainedWaterThermalForceFactorSelection,
        ):
            raise TypeError("cut transport requires an explicit retained-water alpha selection")
        if not isinstance(
            self.moving_interface_composition_force_authority,
            MovingInterfaceCompositionForceAuthority,
        ):
            raise TypeError("cut transport requires an explicit moving-interface force authority")
        if (
            self.moving_interface_composition_force_authority
            is MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
            and not isinstance(self.dry.primitive_band, ct.GasOnlyPrimitiveDomain)
        ):
            raise ValueError(
                "the actual-path moving-interface authority is restricted to the "
                "exact GasOnlyPrimitiveDomain route"
            )
        if not isinstance(
            self.vanishing_wet_cut_threshold_authority,
            VanishingWetCutThresholdAuthority,
        ):
            raise TypeError("cut transport requires a frozen vanishing-wet-cut threshold authority")
        if (
            not nip.reference_invariant_thermal_factor_mode(self.mass_force_mode)
            and self.binary_thermal_diffusion_factor.alpha_t != 0.0
        ):
            raise ValueError(
                "a nonzero binary alpha_T is allowed only in its explicit "
                "reference-invariant counterfactual mode"
            )
        if (
            self.mass_force_mode
            is not nip.NonisothermalMassForceMode.REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS_COUNTERFACTUAL
            and self.retained_water_thermal_force_factor.alpha_ret != 0.0
        ):
            raise ValueError(
                "a nonzero retained-water alpha is allowed only in the explicit "
                "reference-invariant thermal-force-factors counterfactual mode"
            )
        if self.dry.pressure_pa != self.wet.wet.pressure_pa:
            raise ValueError("wet, interface, and dry pressures must be identical")
        if self.dry.pore.rho_dm_p != self.wet.wet.rho_dm_p:
            raise ValueError("wet and dry dry-meal densities must be identical")
        if self.dry.pore.eps_g != self.wet.wet.epsilon_p:
            raise ValueError("wet and dry porosity authorities must be identical")
        if self.dry.retained_water_mobility != self.wet.retained_water_mobility:
            raise ValueError("one retained-water mobility must be shared across the front")
        if self.dry.pore.luikov != self.wet.luikov:
            raise ValueError("wet and dry retained-water state authorities must match")
        if self.dry.pore.retained_water_liquid_volume_fraction != 0.0:
            raise ValueError(
                "cut transport requires the nominal zero retained-water partial-"
                "volume energy convention on both sides of the front; a nonzero "
                "sensitivity requires a matched wet-side implementation"
            )
        shared_state_fields = (
            ("dry-meal heat capacity", self.dry.pore.cp_dry_meal, self.wet.wet.cp_dry_meal),
            ("oil heat capacity", self.dry.pore.cp_oil, self.wet.wet.cp_oil),
            ("solid energy datum", self.dry.pore.T_ref_solid, self.wet.wet.T_ref_solid),
            ("hexane GAB authority", self.dry.pore.gab, self.wet.wet.gab),
            ("oil isotherm authority", self.dry.pore.oil, self.wet.wet.oil),
        )
        for name, dry_value, wet_value in shared_state_fields:
            if dry_value != wet_value:
                raise ValueError(f"wet and dry {name} must be identical")


@dataclass(frozen=True)
class CutTransportLayout:
    """Positive-piece indexing and the explicit square rank."""

    wet_cell_indices: tuple[int, ...]
    dry_cell_indices: tuple[int, ...]
    cut_cell_index: int
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def wet_piece_count(self) -> int:
        return len(self.wet_cell_indices)

    @property
    def dry_piece_count(self) -> int:
        return len(self.dry_cell_indices)

    @property
    def dry_face_count(self) -> int:
        return self.dry_piece_count + 1

    @property
    def unknown_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("wet_temperature", self.wet_piece_count),
            ("wet_retained_water", self.wet_piece_count),
            ("dry_temperature", self.dry_piece_count),
            ("dry_y_hexane", self.dry_piece_count),
            ("dry_total_stefan_flux", self.dry_face_count),
            ("front_z", 1),
            ("interface_temperature", 1),
        )

    @property
    def residual_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("wet_water", self.wet_piece_count),
            ("wet_energy", self.wet_piece_count),
            ("dry_water", self.dry_piece_count),
            ("dry_hexane", self.dry_piece_count),
            ("dry_energy", self.dry_piece_count),
            ("rh_water", 1),
            ("rh_hexane", 1),
            ("rh_energy", 1),
        )

    @property
    def unknown_count(self) -> int:
        return sum(size for _, size in self.unknown_blocks)

    @property
    def residual_count(self) -> int:
        return sum(size for _, size in self.residual_blocks)

    @property
    def is_square(self) -> bool:
        return self.unknown_count == self.residual_count

    @property
    def rank_formula(self) -> int:
        return 2 * self.wet_piece_count + 3 * self.dry_piece_count + 3


@dataclass(frozen=True)
class CutTransportState:
    """Immutable differential state on one strict partial-front cut chart."""

    geometry: cg.CutGeometry
    config: CutTransportConfig
    time_s: float
    wet_temperatures_k: tuple[float, ...]
    wet_retained_water_loadings: tuple[float, ...]
    dry_temperatures_k: tuple[float, ...]
    dry_y_hexane: tuple[float, ...]
    historical_hexane_loadings: tuple[float, ...]
    oil_fraction_labels: tuple[float, ...]
    wet_retained_water_capacity_duals_over_rt: tuple[float, ...] = ()
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        layout = layout_for_geometry(self.geometry)
        if not math.isfinite(self.time_s) or self.time_s < 0.0:
            raise ValueError("cut-transport time must be finite and non-negative")
        _require_length("wet temperatures", self.wet_temperatures_k, layout.wet_piece_count)
        _require_length(
            "wet retained-water loadings",
            self.wet_retained_water_loadings,
            layout.wet_piece_count,
        )
        _require_length("dry temperatures", self.dry_temperatures_k, layout.dry_piece_count)
        _require_length("dry compositions", self.dry_y_hexane, layout.dry_piece_count)
        _require_length(
            "historical n-hexane labels",
            self.historical_hexane_loadings,
            self.geometry.master_grid.n,
        )
        _require_length("oil labels", self.oil_fraction_labels, self.geometry.master_grid.n)
        if not all(
            math.isfinite(value) and value >= 0.0 for value in self.historical_hexane_loadings
        ):
            raise ValueError("historical n-hexane labels must be finite and non-negative")
        if not all(math.isfinite(value) for value in self.oil_fraction_labels):
            raise ValueError("oil labels must be finite")
        wet_duals = wrc.effective_capacity_duals(
            self.wet_retained_water_capacity_duals_over_rt,
            layout.wet_piece_count,
        )
        _evaluate_wet_piece_states(
            layout,
            self.wet_temperatures_k,
            self.wet_retained_water_loadings,
            self.historical_hexane_loadings,
            self.oil_fraction_labels,
            self.config,
            capacity_duals_over_rt=wet_duals,
        )
        wrc.certify_exact_graph(
            self.wet_retained_water_loadings,
            wet_duals,
            self.config.wet.luikov,
        )
        _evaluate_dry_piece_states(
            layout,
            self.dry_temperatures_k,
            self.dry_y_hexane,
            self.oil_fraction_labels,
            self.config,
        )

    @property
    def layout(self) -> CutTransportLayout:
        return layout_for_geometry(self.geometry)

    @property
    def effective_wet_retained_water_capacity_duals_over_rt(self) -> tuple[float, ...]:
        """Return one explicit exact-graph dual per positive-volume wet piece."""

        return wrc.effective_capacity_duals(
            self.wet_retained_water_capacity_duals_over_rt,
            self.layout.wet_piece_count,
        )


@dataclass(frozen=True)
class CutTransportUnknowns:
    """One nonlinear candidate in the documented rank ordering."""

    wet_temperatures_k: tuple[float, ...]
    wet_retained_water_loadings: tuple[float, ...]
    dry_temperatures_k: tuple[float, ...]
    dry_y_hexane: tuple[float, ...]
    dry_total_stefan_fluxes_mol_m2_s: tuple[float, ...]
    front_z: float
    interface_temperature_k: float
    wet_retained_water_capacity_duals_over_rt: tuple[float, ...] = ()
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        wrc.effective_capacity_duals(
            self.wet_retained_water_capacity_duals_over_rt,
            len(self.wet_retained_water_loadings),
        )

    @property
    def effective_wet_retained_water_capacity_duals_over_rt(self) -> tuple[float, ...]:
        return wrc.effective_capacity_duals(
            self.wet_retained_water_capacity_duals_over_rt,
            len(self.wet_retained_water_loadings),
        )

    def vector(self) -> tuple[float, ...]:
        """Return the equation-rank physical ordering described by ``unknown_blocks``.

        Each wet loading and its capacity dual occupy one semismooth graph
        coordinate in the nonlinear integrator.  The dual is therefore branch
        metadata on that one coordinate, not an additional equation-rank
        entry in this physical vector.
        """

        if any(value != 0.0 for value in self.effective_wet_retained_water_capacity_duals_over_rt):
            raise CutTransportTopologyError(
                "the legacy physical vector cannot represent a positive wet retained-cap "
                "dual; use the integrator graph coordinate or lossless_state_vector"
            )
        return self.equation_rank_vector()

    def equation_rank_vector(self) -> tuple[float, ...]:
        """Return the square equation-rank fields, excluding cap branch metadata."""

        return (
            *self.wet_temperatures_k,
            *self.wet_retained_water_loadings,
            *self.dry_temperatures_k,
            *self.dry_y_hexane,
            *self.dry_total_stefan_fluxes_mol_m2_s,
            self.front_z,
            self.interface_temperature_k,
        )

    def lossless_state_vector(self) -> tuple[float, ...]:
        """Return a restart/identity representation that distinguishes cap duals."""

        return (
            *self.equation_rank_vector(),
            *self.effective_wet_retained_water_capacity_duals_over_rt,
        )

    @classmethod
    def from_vector(
        cls, layout: CutTransportLayout, values: Sequence[float]
    ) -> "CutTransportUnknowns":
        vector = tuple(values)
        if len(vector) != layout.unknown_count:
            raise ValueError(
                f"candidate vector has length {len(vector)}, expected {layout.unknown_count}"
            )
        if not all(math.isfinite(value) for value in vector):
            raise ValueError("candidate vector must be finite")
        nw = layout.wet_piece_count
        nd = layout.dry_piece_count
        nf = layout.dry_face_count
        cursor = 0

        def take(count: int) -> tuple[float, ...]:
            nonlocal cursor
            result = vector[cursor : cursor + count]
            cursor += count
            return result

        return cls(
            wet_temperatures_k=take(nw),
            wet_retained_water_loadings=take(nw),
            dry_temperatures_k=take(nd),
            dry_y_hexane=take(nd),
            dry_total_stefan_fluxes_mol_m2_s=take(nf),
            front_z=take(1)[0],
            interface_temperature_k=take(1)[0],
        )


@dataclass(frozen=True)
class InterfaceCompositionRootAudit:
    """Per-root certificate for the representable PHY-039 interface trace."""

    caller_lower_y_hexane: float
    caller_upper_y_hexane: float
    final_lower_y_hexane: float
    final_upper_y_hexane: float
    returned_y_hexane: float
    returned_log_fugacity_residual: float
    returned_product_hexane_activity: float
    actual_oil_fraction_label: float
    log_fugacity_tolerance: float
    bisection_iteration_count: int
    rounding_plateau_encountered: bool
    no_distinct_midpoint_encountered: bool
    final_upper_classification: str
    final_upper_log_fugacity_residual: float
    final_upper_product_hexane_activity: float | None
    full_coupled_pore_certification_passed: bool = field(default=True, init=False)
    numerical_amendment_id: str = field(default="GT-PS-2-P1E-18", init=False)
    postprojection_or_nextafter_used: bool = field(default=False, init=False)
    log_fugacity_equation_changed: bool = field(default=False, init=False)
    log_fugacity_tolerance_changed: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        finite = (
            self.caller_lower_y_hexane,
            self.caller_upper_y_hexane,
            self.final_lower_y_hexane,
            self.final_upper_y_hexane,
            self.returned_y_hexane,
            self.returned_log_fugacity_residual,
            self.returned_product_hexane_activity,
            self.actual_oil_fraction_label,
            self.log_fugacity_tolerance,
            self.final_upper_log_fugacity_residual,
        )
        if not all(math.isfinite(value) for value in finite):
            raise ValueError("interface-composition root audit values must be finite")
        if not (0.0 < self.caller_lower_y_hexane < self.caller_upper_y_hexane < 1.0):
            raise ValueError("interface-composition caller bracket must lie inside (0,1)")
        if not (
            self.caller_lower_y_hexane
            <= self.final_lower_y_hexane
            <= self.returned_y_hexane
            <= self.final_upper_y_hexane
            <= self.caller_upper_y_hexane
        ):
            raise ValueError("interface-composition final bracket lost its caller interval")
        if not 0.0 <= self.returned_log_fugacity_residual <= self.log_fugacity_tolerance:
            raise ValueError("interface-composition returned residual exceeds its tolerance")
        if not 0.0 < self.returned_product_hexane_activity <= 1.0:
            raise ValueError("interface-composition returned activity is supersaturated")
        if (
            isinstance(self.bisection_iteration_count, bool)
            or not isinstance(self.bisection_iteration_count, int)
            or not 0 <= self.bisection_iteration_count <= 120
        ):
            raise ValueError("interface-composition iteration count must lie in [0,120]")
        if not isinstance(self.final_upper_classification, str) or not (
            self.final_upper_classification.strip()
        ):
            raise ValueError("interface-composition upper classification is required")
        allowed_upper_classifications = {
            "negative_log_fugacity_residual_liquid_side",
            "product_activity_supersaturated_on_log_tolerance_plateau",
            "returned_fully_certified_upper_endpoint",
        }
        if self.final_upper_classification not in allowed_upper_classifications:
            raise ValueError("interface-composition upper classification is unknown")
        if self.final_upper_product_hexane_activity is not None and not math.isfinite(
            self.final_upper_product_hexane_activity
        ):
            raise ValueError("interface-composition upper activity must be finite when present")
        if (
            self.final_upper_classification == ("negative_log_fugacity_residual_liquid_side")
            and self.final_upper_log_fugacity_residual >= 0.0
        ):
            raise ValueError("log-liquid upper classification requires a negative residual")
        if self.final_upper_classification == (
            "product_activity_supersaturated_on_log_tolerance_plateau"
        ) and not (
            0.0 <= self.final_upper_log_fugacity_residual <= self.log_fugacity_tolerance
            and self.final_upper_product_hexane_activity is not None
            and self.final_upper_product_hexane_activity > 1.0
        ):
            raise ValueError("product-supersaturated upper classification is inconsistent")
        if self.final_upper_classification == "returned_fully_certified_upper_endpoint" and not (
            self.returned_y_hexane == self.final_upper_y_hexane
            and self.returned_log_fugacity_residual == self.final_upper_log_fugacity_residual
            and self.returned_product_hexane_activity == self.final_upper_product_hexane_activity
        ):
            raise ValueError("returned upper-endpoint classification is inconsistent")
        adjacent_final_endpoints = (
            math.nextafter(self.final_lower_y_hexane, math.inf) == self.final_upper_y_hexane
        )
        if self.rounding_plateau_encountered and not (
            self.no_distinct_midpoint_encountered and adjacent_final_endpoints
        ):
            raise ValueError(
                "an interface-composition rounding plateau must terminate at "
                "adjacent binary64 endpoints with no distinct midpoint"
            )
        if not self.rounding_plateau_encountered and self.no_distinct_midpoint_encountered:
            raise ValueError("a no-midpoint root audit requires a detected rounding plateau")
        object.__setattr__(
            self,
            "log_fugacity_tolerance_changed",
            self.log_fugacity_tolerance != 2.0e-11,
        )
        if not self.full_coupled_pore_certification_passed:
            raise ValueError("a returned interface root must carry full pore certification")


@dataclass(frozen=True)
class CutInterfaceState:
    """One continuous-temperature, fugacity-matched wet/dry trace pair."""

    temperature_k: float
    y_hexane: float
    log_hexane_fugacity_residual: float
    log_retained_water_fugacity_residual: float
    retained_water_trace_capacity_dual_over_rt: float
    dry_actual_water_potential_isothermal: float
    dry: cp.EquilibriumPoreState
    wet: ww.WetWaterCellState
    historical_hexane_loading: float
    oil_fraction_label: float
    interface_composition_root_audit: InterfaceCompositionRootAudit
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class DryFaceFlux:
    """Constitutive and algebraic fluxes at one dry-domain face.

    ``conductive_heat_flux_w_m2`` and ``fourier_heat_flux_w_m2`` are the same
    structural Fourier flux.  In the explicit thermal-factor counterfactual,
    the energy object's historical ``conductive`` slot carries the complete
    non-advective reduced heat ``q_F+q_D``; the separate audit fields prevent
    the reciprocal part from being mislabeled as material conductivity.
    """

    component: cp.ComponentMolarFluxes
    energy: cp.ComponentEnergyFlux
    independent_water_flux_mol_m2_s: float
    retained_water_flux_mol_m2_s: float
    conductive_heat_flux_w_m2: float
    gas_entropy_production_w_m3_k: float
    retained_water_entropy_production_w_m3_k: float
    conduction_entropy_production_w_m3_k: float
    thermodynamic_force_temperature_k: float | None
    thermodynamic_force_reconstruction: DryThermodynamicForceReconstruction
    composition_force_path: acfp.ActualCompositionForcePath | None = None
    mass_force_mode: nip.NonisothermalMassForceMode = (
        nip.NonisothermalMassForceMode.COMMON_FACE_TEMPERATURE_NEGLIGIBLE_SORET
    )
    selected_binary_mobility_mol_m_s: float = 0.0
    gas_chemical_force_gradient_m_inv: float = 0.0
    binary_thermal_force_gradient_m_inv: float = 0.0
    binary_thermal_diffusion_factor: float = 0.0
    retained_water_chemical_force_gradient_m_inv: float = 0.0
    retained_water_thermal_force_gradient_m_inv: float = 0.0
    retained_water_thermal_force_factor: float = 0.0
    log_temperature_gradient_m_inv: float = 0.0
    fourier_heat_flux_w_m2: float = 0.0
    reciprocal_heat_of_transport_flux_w_m2: float = 0.0
    reduced_heat_flux_w_m2: float = 0.0
    discrete_conjugate_temperature_k: float | None = None
    surface_film_audit: ct.SurfaceFilmAudit | None = None
    # GT-PS-2-P1E-20 section 3.4 requires the accepted moving-face audit to
    # serialize enough independent primitives to recompute the actual-endpoint
    # mobility selection.  The scalar above is the selected value only; the
    # symmetric face state it was evaluated from was never carried forward, so
    # a verifier could not recompute it and had to rebuild the state from
    # serialized primitives instead.  That rebuild follows a different
    # floating-point operation order and lands one ULP away, which made the
    # exact-equality contract unsatisfiable.  Retaining the state here is
    # audit-only: it changes no equation, no flux, and no accepted value.
    selected_binary_mobility_face_state: tc.SymmetricBinaryFaceState | None = None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(
            self.thermodynamic_force_reconstruction,
            DryThermodynamicForceReconstruction,
        ):
            raise TypeError("dry-face force reconstruction must be explicit")
        if not isinstance(self.mass_force_mode, nip.NonisothermalMassForceMode):
            raise TypeError("dry-face mass-force mode must be explicit")
        actual_path_reconstructions = (
            DryThermodynamicForceReconstruction.MOVING_INTERFACE_ACTUAL_T_Y_PATH,
            DryThermodynamicForceReconstruction.MOVING_INTERFACE_COMPLETE_ENDPOINT_POTENTIALS_ACTUAL_STATES,
        )
        if self.thermodynamic_force_reconstruction in actual_path_reconstructions:
            if not isinstance(
                self.composition_force_path,
                acfp.ActualCompositionForcePath,
            ):
                raise TypeError("actual-path moving-face reconstruction requires its path audit")
            if self.thermodynamic_force_temperature_k is not None:
                raise ValueError(
                    "an actual-path force has endpoint temperatures, not one face temperature"
                )
        elif self.composition_force_path is not None:
            raise ValueError("a non-path dry-face reconstruction cannot carry an actual-path audit")
        force_values = (
            self.selected_binary_mobility_mol_m_s,
            self.gas_chemical_force_gradient_m_inv,
            self.binary_thermal_force_gradient_m_inv,
            self.binary_thermal_diffusion_factor,
            self.retained_water_chemical_force_gradient_m_inv,
            self.retained_water_thermal_force_gradient_m_inv,
            self.retained_water_thermal_force_factor,
            self.log_temperature_gradient_m_inv,
            self.fourier_heat_flux_w_m2,
            self.reciprocal_heat_of_transport_flux_w_m2,
            self.reduced_heat_flux_w_m2,
        )
        if not all(math.isfinite(value) for value in force_values):
            raise ValueError("dry-face force audit values must be finite")
        if self.selected_binary_mobility_mol_m_s < 0.0:
            raise ValueError("dry-face binary mobility cannot be negative")
        if self.thermodynamic_force_reconstruction is (
            DryThermodynamicForceReconstruction.CENTER_SYMMETRY
        ):
            if self.thermodynamic_force_temperature_k is not None:
                raise ValueError("center symmetry has no thermodynamic-force temperature")
            if self.discrete_conjugate_temperature_k is not None:
                raise ValueError("center symmetry has no conjugate face temperature")
        elif self.thermodynamic_force_reconstruction in actual_path_reconstructions:
            if not (
                self.discrete_conjugate_temperature_k is not None
                and math.isfinite(self.discrete_conjugate_temperature_k)
                and self.discrete_conjugate_temperature_k > 0.0
            ):
                raise ValueError("actual-path dry face needs a positive conjugate temperature")
        elif not (
            self.thermodynamic_force_temperature_k is not None
            and math.isfinite(self.thermodynamic_force_temperature_k)
            and self.thermodynamic_force_temperature_k > 0.0
        ):
            raise ValueError("dry-face force temperature must be positive and finite")
        elif not (
            self.discrete_conjugate_temperature_k is not None
            and math.isfinite(self.discrete_conjugate_temperature_k)
            and self.discrete_conjugate_temperature_k > 0.0
        ):
            raise ValueError("dry-face conjugate temperature must be positive and finite")
        if self.conductive_heat_flux_w_m2 != self.fourier_heat_flux_w_m2:
            raise ValueError("dry-face Fourier heat audit lost its identity")
        expected_reduced_heat = (
            self.fourier_heat_flux_w_m2
            if self.reciprocal_heat_of_transport_flux_w_m2 == 0.0
            else math.fsum(
                (
                    self.fourier_heat_flux_w_m2,
                    self.reciprocal_heat_of_transport_flux_w_m2,
                )
            )
        )
        if self.reduced_heat_flux_w_m2 != expected_reduced_heat:
            raise ValueError("dry-face reciprocal reduced-heat identity was lost")
        if self.energy.conductive_heat_flux_w_m2 != self.reduced_heat_flux_w_m2:
            raise ValueError("dry-face energy does not carry the selected reduced heat")
        if (
            self.thermodynamic_force_reconstruction
            is not DryThermodynamicForceReconstruction.CENTER_SYMMETRY
        ):
            effective_gas_force = math.fsum(
                (
                    self.gas_chemical_force_gradient_m_inv,
                    self.binary_thermal_force_gradient_m_inv,
                )
            )
            expected_hexane_counterflux = (
                self.selected_binary_mobility_mol_m_s * effective_gas_force
            )
            if -self.independent_water_flux_mol_m2_s != expected_hexane_counterflux:
                raise ValueError("dry-face binary mobility/force/counterflux identity was lost")

    @property
    def total_entropy_production_w_m3_k(self) -> float:
        return math.fsum(
            (
                self.gas_entropy_production_w_m3_k,
                self.retained_water_entropy_production_w_m3_k,
                self.conduction_entropy_production_w_m3_k,
            )
        )


@dataclass(frozen=True)
class CutResidualBlocks:
    """The square regional and RH residual blocks, all in extensive rates."""

    wet_water_mol_s: tuple[float, ...]
    wet_energy_w: tuple[float, ...]
    dry_water_mol_s: tuple[float, ...]
    dry_hexane_mol_s: tuple[float, ...]
    dry_energy_w: tuple[float, ...]
    rh_water_mol_s: float
    rh_hexane_mol_s: float
    rh_energy_w: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def vector(self) -> tuple[float, ...]:
        return (
            *self.wet_water_mol_s,
            *self.wet_energy_w,
            *self.dry_water_mol_s,
            *self.dry_hexane_mol_s,
            *self.dry_energy_w,
            self.rh_water_mol_s,
            self.rh_hexane_mol_s,
            self.rh_energy_w,
        )


@dataclass(frozen=True)
class CutTransportLedger:
    """Independent global/telescoping checks for one residual evaluation."""

    pressure_before_pa: float
    pressure_after_pa: float
    finite_pressure_transition: bool
    wet_hexane_sweep_identity_mol_s: float
    water_global_from_blocks_mol_s: float
    water_global_direct_mol_s: float
    water_telescoping_error_mol_s: float
    hexane_global_from_blocks_mol_s: float
    hexane_global_direct_mol_s: float
    hexane_telescoping_error_mol_s: float
    energy_global_from_blocks_w: float
    energy_global_direct_w: float
    energy_telescoping_error_w: float
    geometry_cell_gcl_residual_m3: float
    minimum_entropy_production_w_m3_k: float
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class VanishingWetCutTreatmentAudit:
    """Immutable Amendment-12 local reconstruction and conservation audit."""

    threshold_authority: VanishingWetCutThresholdAuthority
    target_wet_volume_fraction: float | None
    treatment_enabled: bool
    wet_volume_fraction: float
    geometric_weight: float
    applied_weight: float
    full_wet_neighbor_available: bool
    active: bool
    activation_count: int
    cut_master_cell_index: int
    adjacent_wet_piece_position: int | None
    adjacent_master_cell_index: int | None
    adjacent_to_interface_distance_m: float | None
    interface_area_m2: float | None
    original_adjacent_water_residual_mol_s: float | None
    original_cut_water_residual_mol_s: float
    transformed_adjacent_water_residual_mol_s: float | None
    transformed_cut_water_residual_mol_s: float
    original_adjacent_common_energy_residual_w: float | None
    original_cut_common_energy_residual_w: float
    transformed_adjacent_common_energy_residual_w: float | None
    transformed_cut_common_energy_residual_w: float
    original_cut_interface_water_flux_mol_m2_s: float | None
    reconstructed_adjacent_interface_water_flux_mol_m2_s: float | None
    original_cut_interface_common_energy_flux_w_m2: float | None
    reconstructed_adjacent_interface_common_energy_flux_w_m2: float | None
    water_flux_difference_area_weighted_mol_s: float | None
    common_energy_flux_difference_area_weighted_w: float | None
    original_water_pair_sum_mol_s: float | None
    transformed_water_pair_sum_mol_s: float | None
    water_pair_sum_difference_mol_s: float
    original_common_energy_pair_sum_w: float | None
    transformed_common_energy_pair_sum_w: float | None
    common_energy_pair_sum_difference_w: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(
            self.threshold_authority,
            VanishingWetCutThresholdAuthority,
        ):
            raise TypeError("wet-cut audit needs a frozen threshold authority")
        expected_target = self.threshold_authority.target_wet_volume_fraction
        if self.target_wet_volume_fraction != expected_target:
            raise ValueError("wet-cut audit target disagrees with its authority")
        expected_enabled = expected_target is not None
        if self.treatment_enabled is not expected_enabled:
            raise ValueError("wet-cut audit enabled state disagrees with its authority")
        if not 0.0 < self.wet_volume_fraction < 1.0:
            raise ValueError("same-cell wet-cut audit requires a strict partial fraction")
        expected_geometric = vanishing_wet_cut_weight(
            self.wet_volume_fraction,
            self.threshold_authority,
        )
        if self.geometric_weight != expected_geometric:
            raise ValueError("wet-cut audit geometric weight is inconsistent")
        expected_applied = (
            expected_geometric if expected_enabled and self.full_wet_neighbor_available else 0.0
        )
        if self.applied_weight != expected_applied:
            raise ValueError("wet-cut audit applied weight is inconsistent")
        expected_active = expected_applied > 0.0
        if self.active is not expected_active or self.activation_count != int(expected_active):
            raise ValueError("wet-cut audit activation state is inconsistent")
        if (
            isinstance(self.cut_master_cell_index, bool)
            or not isinstance(self.cut_master_cell_index, int)
            or self.cut_master_cell_index < 0
        ):
            raise ValueError("wet-cut audit cut-cell index must be nonnegative")
        if not all(
            math.isfinite(value)
            for value in (
                self.original_cut_water_residual_mol_s,
                self.transformed_cut_water_residual_mol_s,
                self.original_cut_common_energy_residual_w,
                self.transformed_cut_common_energy_residual_w,
                self.water_pair_sum_difference_mol_s,
                self.common_energy_pair_sum_difference_w,
            )
        ):
            raise ValueError("wet-cut audit residual evidence must be finite")

        optional = (
            self.adjacent_wet_piece_position,
            self.adjacent_master_cell_index,
            self.adjacent_to_interface_distance_m,
            self.interface_area_m2,
            self.original_adjacent_water_residual_mol_s,
            self.transformed_adjacent_water_residual_mol_s,
            self.original_adjacent_common_energy_residual_w,
            self.transformed_adjacent_common_energy_residual_w,
            self.original_cut_interface_water_flux_mol_m2_s,
            self.reconstructed_adjacent_interface_water_flux_mol_m2_s,
            self.original_cut_interface_common_energy_flux_w_m2,
            self.reconstructed_adjacent_interface_common_energy_flux_w_m2,
            self.water_flux_difference_area_weighted_mol_s,
            self.common_energy_flux_difference_area_weighted_w,
            self.original_water_pair_sum_mol_s,
            self.transformed_water_pair_sum_mol_s,
            self.original_common_energy_pair_sum_w,
            self.transformed_common_energy_pair_sum_w,
        )
        reconstruction_available = expected_active
        if not reconstruction_available:
            if any(value is not None for value in optional):
                raise ValueError("inactive wet-cut audit cannot carry reconstruction data")
            if (
                self.original_cut_water_residual_mol_s != self.transformed_cut_water_residual_mol_s
                or self.original_cut_common_energy_residual_w
                != self.transformed_cut_common_energy_residual_w
            ):
                raise ValueError("inactive wet-cut audit changed a cut residual row")
            if (
                self.water_pair_sum_difference_mol_s != 0.0
                or self.common_energy_pair_sum_difference_w != 0.0
            ):
                raise ValueError("inactive wet-cut audit cannot change a pair sum")
            return

        if any(value is None for value in optional):
            raise ValueError("eligible wet-cut audit lost reconstruction data")
        assert self.adjacent_wet_piece_position is not None
        assert self.adjacent_master_cell_index is not None
        assert self.adjacent_to_interface_distance_m is not None
        assert self.interface_area_m2 is not None
        if self.adjacent_wet_piece_position < 0 or self.adjacent_master_cell_index < 0:
            raise ValueError("wet-cut adjacent indices must be nonnegative")
        if self.adjacent_to_interface_distance_m <= 0.0 or self.interface_area_m2 <= 0.0:
            raise ValueError("wet-cut reconstruction geometry must be positive")
        numeric = tuple(float(value) for value in optional[2:]) + (
            self.water_pair_sum_difference_mol_s,
            self.common_energy_pair_sum_difference_w,
        )
        if not all(math.isfinite(value) for value in numeric):
            raise ValueError("wet-cut reconstruction audit must be finite")

        def roundoff_limit(*values: float) -> float:
            scale = max((abs(value) for value in values), default=0.0)
            return 128.0 * math.ulp(scale if scale > 0.0 else 1.0)

        assert self.original_water_pair_sum_mol_s is not None
        assert self.transformed_water_pair_sum_mol_s is not None
        assert self.original_common_energy_pair_sum_w is not None
        assert self.transformed_common_energy_pair_sum_w is not None
        assert self.original_adjacent_water_residual_mol_s is not None
        assert self.transformed_adjacent_water_residual_mol_s is not None
        assert self.original_adjacent_common_energy_residual_w is not None
        assert self.transformed_adjacent_common_energy_residual_w is not None
        assert self.original_cut_interface_water_flux_mol_m2_s is not None
        assert self.reconstructed_adjacent_interface_water_flux_mol_m2_s is not None
        assert self.original_cut_interface_common_energy_flux_w_m2 is not None
        assert self.reconstructed_adjacent_interface_common_energy_flux_w_m2 is not None
        assert self.water_flux_difference_area_weighted_mol_s is not None
        assert self.common_energy_flux_difference_area_weighted_w is not None

        expected_water_difference = self.interface_area_m2 * math.fsum(
            (
                self.original_cut_interface_water_flux_mol_m2_s,
                -self.reconstructed_adjacent_interface_water_flux_mol_m2_s,
            )
        )
        expected_energy_difference = self.interface_area_m2 * math.fsum(
            (
                self.original_cut_interface_common_energy_flux_w_m2,
                -self.reconstructed_adjacent_interface_common_energy_flux_w_m2,
            )
        )
        if self.water_flux_difference_area_weighted_mol_s != expected_water_difference:
            raise ValueError("wet-cut water E is not area times its direct flux difference")
        if self.common_energy_flux_difference_area_weighted_w != expected_energy_difference:
            raise ValueError("wet-cut energy E is not area times its direct flux difference")

        expected_adjacent_water = math.fsum(
            (
                self.original_adjacent_water_residual_mol_s,
                self.applied_weight
                * (self.original_cut_water_residual_mol_s - expected_water_difference),
            )
        )
        expected_cut_water = math.fsum(
            (
                (1.0 - self.applied_weight) * self.original_cut_water_residual_mol_s,
                self.applied_weight * expected_water_difference,
            )
        )
        expected_adjacent_energy = math.fsum(
            (
                self.original_adjacent_common_energy_residual_w,
                self.applied_weight
                * (self.original_cut_common_energy_residual_w - expected_energy_difference),
            )
        )
        expected_cut_energy = math.fsum(
            (
                (1.0 - self.applied_weight) * self.original_cut_common_energy_residual_w,
                self.applied_weight * expected_energy_difference,
            )
        )
        if (
            self.transformed_adjacent_water_residual_mol_s != expected_adjacent_water
            or self.transformed_cut_water_residual_mol_s != expected_cut_water
        ):
            raise ValueError("wet-cut transformed water rows do not reproduce")
        if (
            self.transformed_adjacent_common_energy_residual_w != expected_adjacent_energy
            or self.transformed_cut_common_energy_residual_w != expected_cut_energy
        ):
            raise ValueError("wet-cut transformed energy rows do not reproduce")

        expected_original_water_pair = math.fsum(
            (
                self.original_adjacent_water_residual_mol_s,
                self.original_cut_water_residual_mol_s,
            )
        )
        expected_transformed_water_pair = math.fsum(
            (
                self.transformed_adjacent_water_residual_mol_s,
                self.transformed_cut_water_residual_mol_s,
            )
        )
        expected_water_pair_difference = math.fsum(
            (expected_transformed_water_pair, -expected_original_water_pair)
        )
        expected_original_energy_pair = math.fsum(
            (
                self.original_adjacent_common_energy_residual_w,
                self.original_cut_common_energy_residual_w,
            )
        )
        expected_transformed_energy_pair = math.fsum(
            (
                self.transformed_adjacent_common_energy_residual_w,
                self.transformed_cut_common_energy_residual_w,
            )
        )
        expected_energy_pair_difference = math.fsum(
            (expected_transformed_energy_pair, -expected_original_energy_pair)
        )
        if (
            self.original_water_pair_sum_mol_s != expected_original_water_pair
            or self.transformed_water_pair_sum_mol_s != expected_transformed_water_pair
            or self.water_pair_sum_difference_mol_s != expected_water_pair_difference
        ):
            raise ValueError("wet-cut water pair identities do not reproduce")
        if (
            self.original_common_energy_pair_sum_w != expected_original_energy_pair
            or self.transformed_common_energy_pair_sum_w != expected_transformed_energy_pair
            or self.common_energy_pair_sum_difference_w != expected_energy_pair_difference
        ):
            raise ValueError("wet-cut energy pair identities do not reproduce")

        if abs(self.water_pair_sum_difference_mol_s) > roundoff_limit(
            self.original_water_pair_sum_mol_s,
            self.transformed_water_pair_sum_mol_s,
            self.original_adjacent_water_residual_mol_s,
            self.original_cut_water_residual_mol_s,
            self.transformed_adjacent_water_residual_mol_s,
            self.transformed_cut_water_residual_mol_s,
        ):
            raise RuntimeError("wet-cut water pair sum was not preserved to roundoff")
        if abs(self.common_energy_pair_sum_difference_w) > roundoff_limit(
            self.original_common_energy_pair_sum_w,
            self.transformed_common_energy_pair_sum_w,
            self.original_adjacent_common_energy_residual_w,
            self.original_cut_common_energy_residual_w,
            self.transformed_adjacent_common_energy_residual_w,
            self.transformed_cut_common_energy_residual_w,
        ):
            raise RuntimeError("wet-cut energy pair sum was not preserved to roundoff")


@dataclass(frozen=True)
class CutTransportAssembly:
    """Complete nonlinear residual evaluation without state commitment."""

    before: CutTransportState
    target_config: CutTransportConfig
    candidate: CutTransportUnknowns
    candidate_geometry: cg.CutGeometry
    swept_geometry: cg.SweptCutGeometry
    layout: CutTransportLayout
    interface: CutInterfaceState
    wet_cells: tuple[ww.WetWaterCellState, ...]
    dry_cells: tuple[cp.EquilibriumPoreState, ...]
    wet_face_fluxes: tuple[ww.WetWaterFaceFlux, ...]
    dry_face_fluxes: tuple[DryFaceFlux, ...]
    wet_retained_cap_certificate: wrc.WetRetainedCapCertificate
    vanishing_wet_cut_treatment_audit: VanishingWetCutTreatmentAudit
    residuals: CutResidualBlocks
    ledger: CutTransportLedger
    surface_boundary: ct.PoreBoundary
    surface_film_audit: ct.SurfaceFilmAudit | None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.layout.is_square:
            raise RuntimeError("cut transport residual lost its square rank")
        if len(self.candidate.equation_rank_vector()) != self.layout.unknown_count:
            raise RuntimeError("candidate vector does not match declared rank")
        if len(self.residuals.vector) != self.layout.residual_count:
            raise RuntimeError("residual vector does not match declared rank")
        if not isinstance(
            self.wet_retained_cap_certificate,
            wrc.WetRetainedCapCertificate,
        ):
            raise TypeError("cut transport assembly needs an exact wet retained-cap certificate")
        if not isinstance(
            self.vanishing_wet_cut_treatment_audit,
            VanishingWetCutTreatmentAudit,
        ):
            raise TypeError("cut transport assembly needs a wet-cut treatment audit")
        if (
            not math.isfinite(self.ledger.minimum_entropy_production_w_m3_k)
            or self.ledger.minimum_entropy_production_w_m3_k < 0.0
        ):
            raise RuntimeError("cut transport assembly has uncertified local entropy production")


@dataclass(frozen=True)
class _PieceInventories:
    wet_volumes_m3: tuple[float, ...]
    dry_volumes_m3: tuple[float, ...]
    wet_water_concentration_mol_m3: tuple[float, ...]
    wet_hexane_concentration_mol_m3: tuple[float, ...]
    wet_energy_density_j_m3: tuple[float, ...]
    dry_water_concentration_mol_m3: tuple[float, ...]
    dry_hexane_concentration_mol_m3: tuple[float, ...]
    dry_energy_density_j_m3: tuple[float, ...]
    wet_water_mol: tuple[float, ...]
    wet_hexane_mol: tuple[float, ...]
    wet_energy_j: tuple[float, ...]
    dry_water_mol: tuple[float, ...]
    dry_hexane_mol: tuple[float, ...]
    dry_energy_j: tuple[float, ...]

    @property
    def total_water_mol(self) -> float:
        return math.fsum((*self.wet_water_mol, *self.dry_water_mol))

    @property
    def total_hexane_mol(self) -> float:
        return math.fsum((*self.wet_hexane_mol, *self.dry_hexane_mol))

    @property
    def total_energy_j(self) -> float:
        return math.fsum((*self.wet_energy_j, *self.dry_energy_j))


@dataclass(frozen=True)
class _PieceInventoryChanges:
    """Cancellation-safe exact-geometry changes over one fixed cut chart."""

    wet_water_mol: tuple[float, ...]
    wet_hexane_mol: tuple[float, ...]
    wet_energy_j: tuple[float, ...]
    dry_water_mol: tuple[float, ...]
    dry_hexane_mol: tuple[float, ...]
    dry_energy_j: tuple[float, ...]


def layout_for_geometry(geometry: cg.CutGeometry) -> CutTransportLayout:
    """Return the fixed-rank chart for one strict interior cut geometry."""

    if not isinstance(geometry, cg.CutGeometry):
        raise TypeError("cut transport requires a CutGeometry")
    if geometry.front.regime != "partial" or geometry.cut_cell_index is None:
        raise CutTransportTopologyError(
            "this foundation requires a strict positive-volume wet/dry cut cell; "
            "endpoints and master-face events need their separate event chart"
        )
    wet = tuple(cell.index for cell in geometry.cells if cell.wet_volume_m3 > 0.0)
    dry = tuple(cell.index for cell in geometry.cells if cell.dry_volume_m3 > 0.0)
    layout = CutTransportLayout(wet, dry, geometry.cut_cell_index)
    if not layout.is_square or layout.unknown_count != layout.rank_formula:
        raise RuntimeError("cut-cell rank formula is internally inconsistent")
    if wet[-1] != geometry.cut_cell_index or dry[0] != geometry.cut_cell_index:
        raise RuntimeError("cut-cell phase pieces are not contiguous")
    return layout


def initialize_cut_state(
    geometry: cg.CutGeometry,
    config: CutTransportConfig,
    wet_temperatures_k: Sequence[float],
    wet_retained_water_loadings: Sequence[float],
    dry_temperatures_k: Sequence[float],
    dry_y_hexane: Sequence[float],
    historical_hexane_loadings: Sequence[float],
    *,
    oil_fraction_labels: Sequence[float] | None = None,
    wet_retained_water_capacity_duals_over_rt: Sequence[float] = (),
    time_s: float = 0.0,
) -> CutTransportState:
    """Create a validated partial-front state without inventing algebraic fluxes."""

    if oil_fraction_labels is None:
        oil = (config.dry.pore.w_o,) * geometry.master_grid.n
    else:
        oil = tuple(oil_fraction_labels)
    return CutTransportState(
        geometry=geometry,
        config=config,
        time_s=time_s,
        wet_temperatures_k=tuple(wet_temperatures_k),
        wet_retained_water_loadings=tuple(wet_retained_water_loadings),
        dry_temperatures_k=tuple(dry_temperatures_k),
        dry_y_hexane=tuple(dry_y_hexane),
        historical_hexane_loadings=tuple(historical_hexane_loadings),
        oil_fraction_labels=oil,
        wet_retained_water_capacity_duals_over_rt=tuple(wet_retained_water_capacity_duals_over_rt),
    )


def pressure_only_target_config(
    config: CutTransportConfig,
    pressure_pa: float,
) -> CutTransportConfig:
    """Build a matched wet/dry target config without transforming any state."""

    if not isinstance(config, CutTransportConfig):
        raise TypeError("pressure target requires a cut transport config")
    target_dry = ct.pressure_only_target_config(config.dry, pressure_pa)
    target_wet = replace(
        config.wet,
        wet=replace(config.wet.wet, pressure_pa=pressure_pa),
    )
    return replace(config, dry=target_dry, wet=target_wet)


def validate_pressure_only_target_config(
    before: CutTransportConfig,
    target: CutTransportConfig,
) -> None:
    """Reject any cut target change other than the shared prescribed pressure."""

    if not isinstance(before, CutTransportConfig) or not isinstance(
        target,
        CutTransportConfig,
    ):
        raise TypeError("pressure transition requires cut transport configs")
    expected = pressure_only_target_config(before, target.dry.pressure_pa)
    if target != expected:
        raise ValueError(
            "finite-pressure cut target may change pressure only; wet/dry "
            "authorities, coefficients, charts, and solver contracts must be identical"
        )


def candidate_from_state(
    state: CutTransportState,
    *,
    dry_total_stefan_fluxes_mol_m2_s: Sequence[float],
    interface_temperature_k: float,
    front_z: float | None = None,
) -> CutTransportUnknowns:
    """Build an uncommitted candidate using the state's differential primitives."""

    return CutTransportUnknowns(
        wet_temperatures_k=state.wet_temperatures_k,
        wet_retained_water_loadings=state.wet_retained_water_loadings,
        dry_temperatures_k=state.dry_temperatures_k,
        dry_y_hexane=state.dry_y_hexane,
        dry_total_stefan_fluxes_mol_m2_s=tuple(dry_total_stefan_fluxes_mol_m2_s),
        front_z=state.geometry.front.z if front_z is None else front_z,
        interface_temperature_k=interface_temperature_k,
        wet_retained_water_capacity_duals_over_rt=(
            state.effective_wet_retained_water_capacity_duals_over_rt
        ),
    )


def evaluate_interface_state(
    temperature_k: float,
    config: CutTransportConfig,
    historical_hexane_loading: float,
    oil_fraction_label: float,
) -> CutInterfaceState:
    """Evaluate the interface, memoized for bit-exact repeated inputs."""

    return _evaluate_interface_state_cached(
        temperature_k,
        config,
        historical_hexane_loading,
        oil_fraction_label,
        float_bits_key(
            temperature_k,
            historical_hexane_loading,
            oil_fraction_label,
        ),
        hx.caloric_datum_signature(),
        wa.caloric_datum_signature(),
    )


@lru_cache(maxsize=2048, typed=True)
def _evaluate_interface_state_cached(
    temperature_k: float,
    config: CutTransportConfig,
    historical_hexane_loading: float,
    oil_fraction_label: float,
    _scalar_key: tuple[bytes, ...],
    _hexane_datum_key: tuple[bytes, ...],
    _water_datum_key: tuple[bytes, ...],
) -> CutInterfaceState:
    """Evaluate the two fugacity-matched traces without using a dry-cell band.

    The dry local-equilibrium retained loading is reused on the wet trace
    because the retained matrix is continuous and the n-hexane front is not a
    water phase boundary.  On a smooth dry trace its normal-cone multiplier is
    exactly zero.  On a retained-cap dry trace the zero-volume wet trace gets
    the algebraic KKT multiplier ``log(a_w,dry/a_w,cap)`` so its generalized
    retained potential equals the actual dry-gas water potential.  The
    multiplier changes no storage, caloric, phase, or conserved inventory.
    """

    values = (temperature_k, historical_hexane_loading, oil_fraction_label)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("interface temperature and material labels must be finite")
    if historical_hexane_loading < 0.0:
        raise ValueError("interface historical n-hexane loading must be non-negative")
    if not config.wet.wet.T_min <= temperature_k <= config.wet.wet.T_max:
        raise CutTransportTopologyError("interface temperature left the wet caloric bracket")
    pore_lo, pore_hi = config.dry.pore.temperature_bounds_k
    if not pore_lo <= temperature_k <= pore_hi:
        raise CutTransportTopologyError("interface temperature left the pore-state bracket")

    root = _solve_interface_composition_with_audit(
        temperature_k,
        config,
        oil_fraction_label,
    )
    y_hexane = root.y_hexane
    residual = root.log_fugacity_residual
    dry = root.dry
    trace_capacity_dual, dry_actual_water_potential = _retained_water_interface_trace_kkt(
        dry, config
    )
    wet = ww.evaluate_cell(
        temperature_k,
        dry.retained_water_loading,
        historical_hexane_loading,
        oil_fraction_label,
        config.wet,
        retained_water_capacity_dual_over_rt=trace_capacity_dual,
    )
    water_residual = dry_actual_water_potential - wet.retained_water_potential
    if abs(water_residual) > 1.0e-10:
        raise CutTransportTopologyError(
            "wet generalized retained-water potential does not match the actual "
            "dry-gas water potential"
        )
    return CutInterfaceState(
        temperature_k=temperature_k,
        y_hexane=y_hexane,
        log_hexane_fugacity_residual=residual,
        log_retained_water_fugacity_residual=water_residual,
        retained_water_trace_capacity_dual_over_rt=trace_capacity_dual,
        dry_actual_water_potential_isothermal=dry_actual_water_potential,
        dry=dry,
        wet=wet,
        historical_hexane_loading=historical_hexane_loading,
        oil_fraction_label=oil_fraction_label,
        interface_composition_root_audit=root.audit,
    )


def _retained_water_interface_trace_kkt(
    dry: cp.EquilibriumPoreState,
    config: CutTransportConfig,
) -> tuple[float, float]:
    """Return the exact zero-volume wet-trace dual and actual dry potential.

    The helper is shared by the within-cell and exact-face charts so the two
    representations cannot silently use different retained-water interface
    laws.  The returned dual belongs only to the zero-volume wet trace.  It
    carries no storage, volume, energy, source, or additional equation.
    """

    if not isinstance(dry, cp.EquilibriumPoreState):
        raise TypeError("retained-water interface KKT needs an equilibrium pore state")
    if not isinstance(config, CutTransportConfig):
        raise TypeError("retained-water interface KKT needs a cut transport config")
    if not 0.0 < dry.water_activity <= 1.0:
        raise CutTransportTopologyError("interface water activity left (0,1]")
    activity_at_cap = sp.water_activity(config.wet.luikov.W_cap, config.wet.luikov)
    if dry.retained_water_active_set is cp.RetainedWaterActiveSet.RETAINED_CAP:
        trace_capacity_dual = math.log(dry.water_activity / activity_at_cap)
        if trace_capacity_dual < 0.0:
            raise CutTransportTopologyError(
                "cap-active interface produced a negative retained-water KKT multiplier"
            )
    else:
        trace_capacity_dual = 0.0
    dry_actual_water_potential = math.log(
        dry.binary_gas.fugacity_water_pa / cp.REFERENCE_FUGACITY_PA
    )
    return trace_capacity_dual, dry_actual_water_potential


def clear_exact_input_caches() -> None:
    """Clear all bounded exact-input caches used by the cut residual."""

    _evaluate_interface_state_cached.cache_clear()
    cp.clear_equilibrium_cache()
    ww.clear_evaluate_cell_cache()
    bg.clear_state_cache()
    hx.clear_state_tp_cache()
    wa.clear_state_tp_cache()


def exact_input_cache_info() -> dict[str, object]:
    """Return per-process bounded-cache diagnostics for profiling/audits."""

    return {
        "interface_state": _evaluate_interface_state_cached.cache_info(),
        "coupled_pore_equilibrium": cp.equilibrium_cache_info(),
        "wet_cell": ww.evaluate_cell_cache_info(),
        "binary_gas": bg.state_cache_info(),
        "hexane_state_tp": hx.state_tp_cache_info(),
        "water_state_tp": wa.state_tp_cache_info(),
    }


def assemble_backward_euler(
    before: CutTransportState,
    candidate: CutTransportUnknowns,
    dt_s: float,
    surface_boundary: ct.PoreBoundary,
    *,
    target_config: CutTransportConfig | None = None,
    enforce_reduced_film_thresholds: bool = True,
) -> CutTransportAssembly:
    """Assemble one square BE/ALE/RH residual or reject with exact rollback."""

    try:
        return _assemble_backward_euler(
            before,
            candidate,
            dt_s,
            surface_boundary,
            target_config,
            enforce_reduced_film_thresholds,
        )
    except CutTransportStepError:
        raise
    except Exception as exc:
        event_restart = isinstance(exc, CutTransportTopologyError) and ("event restart" in str(exc))
        raise CutTransportStepError(
            f"cut-cell ALE/RH residual rejected with exact rollback: {exc}",
            before,
            event_restart_required=event_restart,
        ) from exc


def _assemble_backward_euler(
    before: CutTransportState,
    candidate: CutTransportUnknowns,
    dt_s: float,
    surface_boundary: ct.PoreBoundary,
    target_config: CutTransportConfig | None,
    enforce_reduced_film_thresholds: bool,
) -> CutTransportAssembly:
    if not isinstance(before, CutTransportState):
        raise TypeError("before must be a CutTransportState")
    if not isinstance(candidate, CutTransportUnknowns):
        raise TypeError("candidate must be CutTransportUnknowns")
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("time step must be positive and finite")
    if not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("surface boundary must be Dirichlet or bounded reduced film")
    config = before.config if target_config is None else target_config
    validate_pressure_only_target_config(before.config, config)

    old_layout = before.layout
    rank_vector = candidate.equation_rank_vector()
    if len(rank_vector) != old_layout.unknown_count:
        raise ValueError("candidate does not match the current cut-chart rank")
    if not all(math.isfinite(value) for value in rank_vector) or not all(
        math.isfinite(value)
        for value in candidate.effective_wet_retained_water_capacity_duals_over_rt
    ):
        raise ValueError("candidate values must be finite")
    candidate_geometry = cg.partition_master_grid(
        before.geometry.master_grid,
        z=candidate.front_z,
    )
    try:
        new_layout = layout_for_geometry(candidate_geometry)
    except CutTransportTopologyError as exc:
        raise CutTransportTopologyError(
            f"front endpoint changed cut-chart rank; event restart required: {exc}"
        ) from exc
    if new_layout != old_layout:
        raise CutTransportTopologyError(
            "front reached or crossed a master face; event restart required before "
            "assembling the next fixed-rank chart"
        )
    swept = cg.swept_cut_geometry(before.geometry, candidate_geometry, dt_s)
    _validate_dry_boundary(surface_boundary, config)

    old_wet = _evaluate_wet_piece_states(
        old_layout,
        before.wet_temperatures_k,
        before.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        before.config,
        capacity_duals_over_rt=(before.effective_wet_retained_water_capacity_duals_over_rt),
    )
    old_dry = _evaluate_dry_piece_states(
        old_layout,
        before.dry_temperatures_k,
        before.dry_y_hexane,
        before.oil_fraction_labels,
        before.config,
    )
    new_wet = _evaluate_wet_piece_states(
        new_layout,
        candidate.wet_temperatures_k,
        candidate.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        config,
        capacity_duals_over_rt=(candidate.effective_wet_retained_water_capacity_duals_over_rt),
    )
    new_dry = _evaluate_dry_piece_states(
        new_layout,
        candidate.dry_temperatures_k,
        candidate.dry_y_hexane,
        before.oil_fraction_labels,
        config,
    )
    _require_length(
        "dry total Stefan fluxes",
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        new_layout.dry_face_count,
    )

    cut_index = new_layout.cut_cell_index
    interface = evaluate_interface_state(
        candidate.interface_temperature_k,
        config,
        before.historical_hexane_loadings[cut_index],
        before.oil_fraction_labels[cut_index],
    )
    wet_fluxes = _wet_face_fluxes(
        candidate_geometry,
        new_layout,
        new_wet,
        interface,
        config,
    )
    wet_retained_cap_certificate = wrc.certify_exact_graph(
        candidate.wet_retained_water_loadings,
        candidate.effective_wet_retained_water_capacity_duals_over_rt,
        config.wet.luikov,
    )
    dry_fluxes = _dry_face_fluxes(
        candidate_geometry,
        new_layout,
        new_dry,
        interface,
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        surface_boundary,
        before.oil_fraction_labels,
        config,
        enforce_reduced_film_thresholds=enforce_reduced_film_thresholds,
    )
    old_inventory = _piece_inventories(
        before.geometry,
        old_layout,
        old_wet,
        old_dry,
        before.historical_hexane_loadings,
        before.config,
    )
    new_inventory = _piece_inventories(
        candidate_geometry,
        new_layout,
        new_wet,
        new_dry,
        before.historical_hexane_loadings,
        config,
    )
    inventory_changes = _piece_inventory_changes(
        new_layout,
        old_inventory,
        new_inventory,
        swept,
    )
    residuals, wet_hexane_identity = _residual_blocks(
        before,
        candidate_geometry,
        new_layout,
        inventory_changes,
        interface,
        wet_fluxes,
        dry_fluxes,
        swept,
        dt_s,
        # O9a front-donor correction: the piece the front is consuming is the
        # OUTERMOST wet piece (the cut cell's wet part), and the donor is that
        # piece's own bulk retained-water concentration.
        #
        # THE TIME LEVEL IS THE MODULE'S OWN, NOT A NEW CHOICE.
        # ``_piece_inventory_changes`` books the matching geometric storage
        # term as ``Delta(V)*C_old``; the ALE sweep is that identical swept
        # volume seen from the flux side, so the donor is read at the same time
        # level - ``C_old`` of the consumed piece.  The two then cancel term by
        # term: the geometric mis-booking ``(-dV)*(W[j*] - donor)/V_new`` is
        # ``(-dV)*(W_old - W_old)/V_new``, the LITERAL zero, and the row keeps
        # exactly the ``V_new`` storage form the decomposition's discrete
        # identity is written in.  A front sweeping uniform material therefore
        # changes that piece's loading by exactly zero, and the row's
        # dependence on the front position - which lives in ``V_new`` - is
        # preserved rather than cancelled away.
        old_inventory.wet_water_concentration_mol_m3[-1],
        # O10a energy-donor completion (ruled 2026-08-21,
        # docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md).  The
        # Rankine-Hugoniot jump conditions must convect mass and enthalpy OF
        # THE SAME MATERIAL STATE - an internal-consistency requirement of the
        # balance laws - and that state is the bulk, undrained one the O9
        # traverse measured (Peclet 3.4e3).  So the enthalpy donor is the SAME
        # consumed piece read at the SAME time level as its mass donor one line
        # above: ``_piece_inventory_changes`` books the paired geometric energy
        # storage term as ``Delta(V)*rho_e_old`` from exactly this entry.
        old_inventory.wet_energy_density_j_m3[-1],
    )
    residuals, wet_cut_treatment_audit = _apply_vanishing_wet_cut_treatment(
        candidate_geometry,
        new_layout,
        new_wet,
        interface,
        wet_fluxes,
        residuals,
        config,
    )
    ledger = _ledger(
        inventory_changes,
        residuals,
        wet_hexane_identity,
        dry_fluxes[-1],
        before.geometry.master_grid.areas[-1],
        swept,
        dt_s,
        (*wet_fluxes, *dry_fluxes),
        before.config.dry.pressure_pa,
        config.dry.pressure_pa,
    )
    return CutTransportAssembly(
        before=before,
        target_config=config,
        candidate=candidate,
        candidate_geometry=candidate_geometry,
        swept_geometry=swept,
        layout=new_layout,
        interface=interface,
        wet_cells=new_wet,
        dry_cells=new_dry,
        wet_face_fluxes=wet_fluxes,
        dry_face_fluxes=dry_fluxes,
        wet_retained_cap_certificate=wet_retained_cap_certificate,
        vanishing_wet_cut_treatment_audit=wet_cut_treatment_audit,
        residuals=residuals,
        ledger=ledger,
        surface_boundary=surface_boundary,
        surface_film_audit=dry_fluxes[-1].surface_film_audit,
    )


def _evaluate_wet_piece_states(
    layout: CutTransportLayout,
    temperatures: Sequence[float],
    water: Sequence[float],
    historical_hexane: Sequence[float],
    oil_labels: Sequence[float],
    config: CutTransportConfig,
    *,
    capacity_duals_over_rt: Sequence[float] = (),
) -> tuple[ww.WetWaterCellState, ...]:
    _require_length("wet temperatures", temperatures, layout.wet_piece_count)
    _require_length("wet retained water", water, layout.wet_piece_count)
    capacity_duals = wrc.effective_capacity_duals(
        capacity_duals_over_rt,
        layout.wet_piece_count,
    )
    return tuple(
        ww.evaluate_cell(
            temperature,
            loading,
            historical_hexane[cell],
            oil_labels[cell],
            config.wet,
            retained_water_capacity_dual_over_rt=capacity_dual,
        )
        for cell, temperature, loading, capacity_dual in zip(
            layout.wet_cell_indices,
            temperatures,
            water,
            capacity_duals,
        )
    )


def _evaluate_dry_piece_states(
    layout: CutTransportLayout,
    temperatures: Sequence[float],
    compositions: Sequence[float],
    oil_labels: Sequence[float],
    config: CutTransportConfig,
) -> tuple[cp.EquilibriumPoreState, ...]:
    _require_length("dry temperatures", temperatures, layout.dry_piece_count)
    _require_length("dry compositions", compositions, layout.dry_piece_count)
    states: list[cp.EquilibriumPoreState] = []
    for cell, temperature, composition in zip(layout.dry_cell_indices, temperatures, compositions):
        pore = replace(config.dry.pore, w_o=oil_labels[cell])
        _validate_open_dry_band(temperature, composition, config, pore=pore)
        states.append(
            cp.evaluate_equilibrium(
                temperature,
                config.dry.pressure_pa,
                composition,
                pore,
            )
        )
    return tuple(states)


def _wet_face_fluxes(
    geometry: cg.CutGeometry,
    layout: CutTransportLayout,
    cells: Sequence[ww.WetWaterCellState],
    interface: CutInterfaceState,
    config: CutTransportConfig,
) -> tuple[ww.WetWaterFaceFlux, ...]:
    fluxes = [_zero_wet_flux("center symmetry")]
    centers = tuple(
        _wet_piece_center(geometry.cells[index], geometry.front.radius_m)
        for index in layout.wet_cell_indices
    )
    for left, right, left_center, right_center in zip(
        cells[:-1], cells[1:], centers[:-1], centers[1:]
    ):
        fluxes.append(
            _wet_flux(
                left.temperature_k,
                left.retained_water_loading,
                right.temperature_k,
                right.retained_water_loading,
                right_center - left_center,
                config,
                "wet cell",
                "wet cell",
                left_capacity_dual_over_rt=(left.retained_water_capacity_dual_over_rt),
                right_capacity_dual_over_rt=(right.retained_water_capacity_dual_over_rt),
            )
        )
    fluxes.append(
        _wet_flux(
            cells[-1].temperature_k,
            cells[-1].retained_water_loading,
            interface.temperature_k,
            interface.wet.retained_water_loading,
            geometry.front.radius_m - centers[-1],
            config,
            "wet cell",
            "interface",
            left_capacity_dual_over_rt=(cells[-1].retained_water_capacity_dual_over_rt),
            right_capacity_dual_over_rt=(interface.retained_water_trace_capacity_dual_over_rt),
        )
    )
    return tuple(fluxes)


def _wet_flux(
    left_temperature: float,
    left_water: float,
    right_temperature: float,
    right_water: float,
    distance_m: float,
    config: CutTransportConfig,
    left_name: str,
    right_name: str,
    *,
    left_capacity_dual_over_rt: float = 0.0,
    right_capacity_dual_over_rt: float = 0.0,
) -> ww.WetWaterFaceFlux:
    if not math.isfinite(distance_m) or distance_m <= 0.0:
        raise CutTransportTopologyError("wet face distance must be positive")
    face_temperature = 0.5 * (left_temperature + right_temperature)
    left_activity = sp.water_activity(left_water, config.wet.luikov)
    right_activity = sp.water_activity(right_water, config.wet.luikov)
    def log_fugacity(loading, activity, liquid):
        if isinstance(config.wet.luikov, sp.ContinuedPositiveLuikovParams):
            return sp.water_log_activity(loading, config.wet.luikov) + math.log(liquid.fugacity)
        return math.log(activity * liquid.fugacity)

    wrc.classify_exact_graph_point(
        left_water,
        left_capacity_dual_over_rt,
        config.wet.luikov,
    )
    wrc.classify_exact_graph_point(
        right_water,
        right_capacity_dual_over_rt,
        config.wet.luikov,
    )
    if (
        config.mass_force_mode
        is not nip.NonisothermalMassForceMode.COMPLETE_POTENTIAL_REFERENCE_GAUGE_COUNTERFACTUAL
    ):
        # Preserve the established branch exactly: both endpoints use the one
        # face-temperature liquid standard fugacity.
        liquid = wa.state_Tp(face_temperature, config.dry.pressure_pa, "liquid")
        left_potential = math.fsum(
            (
                log_fugacity(left_water, left_activity, liquid),
                left_capacity_dual_over_rt,
            )
        )
        right_potential = math.fsum(
            (
                log_fugacity(right_water, right_activity, liquid),
                right_capacity_dual_over_rt,
            )
        )
        gradient = (right_potential - left_potential) / distance_m
    else:
        if left_temperature == right_temperature:
            # Exact isothermal recovery without subtracting complete standard
            # potentials that cancel analytically.
            liquid = wa.state_Tp(face_temperature, config.dry.pressure_pa, "liquid")
            left_potential = math.fsum(
                (
                    log_fugacity(left_water, left_activity, liquid),
                    left_capacity_dual_over_rt,
                )
            )
            right_potential = math.fsum(
                (
                    log_fugacity(right_water, right_activity, liquid),
                    right_capacity_dual_over_rt,
                )
            )
        else:
            left_liquid = wa.state_Tp(
                left_temperature,
                config.dry.pressure_pa,
                "liquid",
            )
            right_liquid = wa.state_Tp(
                right_temperature,
                config.dry.pressure_pa,
                "liquid",
            )
            left_potential = nip.complete_retained_water_potential_over_rt(
                left_temperature,
                left_activity * left_liquid.fugacity,
                config.complete_potential_reference_gauge,
            )
            left_potential = math.fsum((left_potential, left_capacity_dual_over_rt))
            right_potential = nip.complete_retained_water_potential_over_rt(
                right_temperature,
                right_activity * right_liquid.fugacity,
                config.complete_potential_reference_gauge,
            )
            right_potential = math.fsum((right_potential, right_capacity_dual_over_rt))
        gradient = math.fsum((right_potential, -left_potential)) / distance_m
    chemical_gradient = gradient
    log_temperature_gradient = nip.log_temperature_gradient_m_inv(
        left_temperature,
        right_temperature,
        distance_m,
    )
    retained_thermal_gradient = 0.0
    if (
        config.mass_force_mode
        is nip.NonisothermalMassForceMode.REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS_COUNTERFACTUAL
    ):
        retained_thermal_gradient = nip.retained_water_thermal_force_gradient_m_inv(
            left_temperature,
            right_temperature,
            distance_m,
            config.retained_water_thermal_force_factor,
        )
        gradient = math.fsum((chemical_gradient, retained_thermal_gradient))
    water_flux = -config.wet.retained_water_mobility.value_mol_m_s * gradient
    fourier = config.wet.wet.conductivity * (left_temperature - right_temperature) / distance_m
    reciprocal_heat = 0.0
    if nip.reference_invariant_thermal_factor_mode(config.mass_force_mode):
        reciprocal_heat = nip.reciprocal_heat_of_transport_flux_w_m2(
            left_temperature_k=left_temperature,
            right_temperature_k=right_temperature,
            gas_counterflux_water_mol_m2_s=0.0,
            binary_thermal_diffusion_factor=config.binary_thermal_diffusion_factor,
            retained_water_flux_mol_m2_s=water_flux,
            retained_water_thermal_factor=(config.retained_water_thermal_force_factor),
        )
    reduced_heat = fourier if reciprocal_heat == 0.0 else math.fsum((fourier, reciprocal_heat))
    if water_flux > 0.0:
        donor_temperature, donor = left_temperature, left_name
    elif water_flux < 0.0:
        donor_temperature, donor = right_temperature, right_name
    else:
        donor_temperature, donor = face_temperature, "none (zero retained-water flux)"
    enthalpy = wa.state_Tp(donor_temperature, config.dry.pressure_pa, "liquid").h_mass * wa.M
    enthalpy_flux = water_flux * enthalpy
    total_energy = (
        reduced_heat + enthalpy_flux
        if reciprocal_heat == 0.0
        else math.fsum((reduced_heat, enthalpy_flux))
    )
    temperature_gradient = (right_temperature - left_temperature) / distance_m
    conduction_entropy = (
        config.wet.wet.conductivity
        * temperature_gradient
        * temperature_gradient
        / (left_temperature * right_temperature)
    )
    return ww.WetWaterFaceFlux(
        retained_water_flux_mol_m2_s=water_flux,
        conductive_heat_flux_w_m2=fourier,
        retained_water_energy_flux_w_m2=enthalpy_flux,
        total_energy_flux_w_m2=total_energy,
        dimensionless_force_gradient_m_inv=gradient,
        entropy_production_w_m3_k=(
            bg.R * config.wet.retained_water_mobility.value_mol_m_s * gradient * gradient
        ),
        enthalpy_donor=donor,
        used_explicit_nonisothermal_potential=(
            config.mass_force_mode
            is nip.NonisothermalMassForceMode.COMPLETE_POTENTIAL_REFERENCE_GAUGE_COUNTERFACTUAL
        ),
        retained_water_chemical_force_gradient_m_inv=chemical_gradient,
        retained_water_thermal_force_gradient_m_inv=retained_thermal_gradient,
        retained_water_thermal_force_factor=(config.retained_water_thermal_force_factor.alpha_ret),
        log_temperature_gradient_m_inv=log_temperature_gradient,
        fourier_heat_flux_w_m2=fourier,
        reciprocal_heat_of_transport_flux_w_m2=reciprocal_heat,
        reduced_heat_flux_w_m2=reduced_heat,
        conduction_entropy_production_w_m3_k=conduction_entropy,
        discrete_conjugate_temperature_k=(
            nip.discrete_conjugate_temperature_k(
                left_temperature,
                right_temperature,
            )
        ),
    )


def _zero_wet_flux(label: str) -> ww.WetWaterFaceFlux:
    return ww.WetWaterFaceFlux(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, label)


def _uses_nonzero_thermal_factors(config: CutTransportConfig) -> bool:
    """Return whether the selected datum-free counterfactual has nonzero heat cross-effects."""

    return nip.reference_invariant_thermal_factor_mode(config.mass_force_mode) and (
        config.binary_thermal_diffusion_factor.alpha_t != 0.0
        or config.retained_water_thermal_force_factor.alpha_ret != 0.0
    )


@lru_cache(maxsize=128)
def _fixed_composition_surface_temperature_bounds(
    pressure_pa: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
    conditioned_lower_k: float,
    conditioned_upper_k: float,
) -> tuple[float, float]:
    """Resolve the connected open gas-only interval for one surface trace."""

    lower = math.nextafter(conditioned_lower_k, conditioned_upper_k)
    upper = math.nextafter(conditioned_upper_k, conditioned_lower_k)
    if not lower < upper:
        raise CutTransportTopologyError(
            "pressure-conditioned dry domain has no open surface-temperature interval"
        )

    def feasible(temperature_k: float) -> bool:
        try:
            cp.encode_gas_only_y(
                temperature_k,
                pressure_pa,
                y_hexane,
                pore,
            )
        except ValueError, cp.CoupledPoreTopologyError:
            return False
        return True

    samples = tuple(lower + (upper - lower) * index / 64.0 for index in range(65))
    flags = tuple(feasible(temperature) for temperature in samples)
    feasible_indices = tuple(index for index, flag in enumerate(flags) if flag)
    if not feasible_indices:
        raise CutTransportTopologyError(
            "surface composition has no gas-only temperature in the conditioned domain"
        )
    first = feasible_indices[0]
    last = feasible_indices[-1]
    if any(not flags[index] for index in range(first, last + 1)):
        raise CutTransportTopologyError(
            "surface composition has a disconnected gas-only temperature domain"
        )

    def lower_feasible_endpoint(infeasible: float, feasible_value: float) -> float:
        lo = infeasible
        hi = feasible_value
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if mid == lo or mid == hi:
                break
            if feasible(mid):
                hi = mid
            else:
                lo = mid
        return hi

    def upper_feasible_endpoint(feasible_value: float, infeasible: float) -> float:
        lo = feasible_value
        hi = infeasible
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if mid == lo or mid == hi:
                break
            if feasible(mid):
                lo = mid
            else:
                hi = mid
        return lo

    fixed_lower = (
        lower if first == 0 else lower_feasible_endpoint(samples[first - 1], samples[first])
    )
    fixed_upper = (
        upper
        if last == len(samples) - 1
        else upper_feasible_endpoint(samples[last], samples[last + 1])
    )
    if not (fixed_lower < fixed_upper and feasible(fixed_lower) and feasible(fixed_upper)):
        raise CutTransportTopologyError(
            "surface gas-only temperature interval could not be resolved safely"
        )
    return fixed_lower, fixed_upper


def _non_fourier_surface_temperature_bounds(
    config: CutTransportConfig,
    surface_y_hexane: float,
    surface_pore: cp.CoupledPoreParams,
) -> tuple[float, float]:
    """Return strict caller-owned fixed-composition bounds for the surface root."""

    lower, upper = config.dry.conditioned_temperature_domain.solver_bounds_k
    return _fixed_composition_surface_temperature_bounds(
        config.dry.pressure_pa,
        surface_y_hexane,
        surface_pore,
        lower,
        upper,
    )


def _dry_face_fluxes(
    geometry: cg.CutGeometry,
    layout: CutTransportLayout,
    cells: Sequence[cp.EquilibriumPoreState],
    interface: CutInterfaceState,
    stefan_fluxes: Sequence[float],
    boundary: ct.PoreBoundary,
    oil_labels: Sequence[float],
    config: CutTransportConfig,
    *,
    newborn_temperature_gradient_k_m: float | None = None,
    newborn_composition_gradient_m_inv: float | None = None,
    enforce_reduced_film_thresholds: bool = True,
) -> tuple[DryFaceFlux, ...]:
    centers = tuple(
        _dry_piece_center(geometry.cells[index], geometry.front.radius_m)
        for index in layout.dry_cell_indices
    )
    fluxes: list[DryFaceFlux] = []
    first_cell_index = layout.dry_cell_indices[0]
    fluxes.append(
        _dry_flux(
            interface.temperature_k,
            interface.y_hexane,
            cells[0].temperature_k,
            cells[0].y_hexane,
            centers[0] - geometry.front.radius_m,
            stefan_fluxes[0],
            replace(config.dry.pore, w_o=oil_labels[first_cell_index]),
            replace(config.dry.pore, w_o=oil_labels[first_cell_index]),
            config,
            thermodynamic_force_temperature_k=interface.temperature_k,
            temperature_gradient_k_m=newborn_temperature_gradient_k_m,
            composition_gradient_m_inv=newborn_composition_gradient_m_inv,
        )
    )
    for face, (left, right, left_center, right_center) in enumerate(
        zip(cells[:-1], cells[1:], centers[:-1], centers[1:]),
        start=1,
    ):
        left_index = layout.dry_cell_indices[face - 1]
        right_index = layout.dry_cell_indices[face]
        fluxes.append(
            _dry_flux(
                left.temperature_k,
                left.y_hexane,
                right.temperature_k,
                right.y_hexane,
                right_center - left_center,
                stefan_fluxes[face],
                replace(config.dry.pore, w_o=oil_labels[left_index]),
                replace(config.dry.pore, w_o=oil_labels[right_index]),
                config,
            )
        )
    last_index = layout.dry_cell_indices[-1]
    outer_distance = geometry.master_grid.R - centers[-1]
    if isinstance(boundary, ct.DirichletPoreBoundary):
        # Preserve the original numerical-oracle surface path exactly.
        fluxes.append(
            _dry_flux(
                cells[-1].temperature_k,
                cells[-1].y_hexane,
                boundary.temperature_k,
                boundary.y_hexane,
                outer_distance,
                stefan_fluxes[-1],
                replace(config.dry.pore, w_o=oil_labels[last_index]),
                replace(config.dry.pore, w_o=oil_labels[last_index]),
                config,
            )
        )
    else:
        surface_reference = ct.prepare_reduced_film_surface_state(
            boundary,
            cell_temperature_k=cells[-1].temperature_k,
            pressure_pa=config.dry.pressure_pa,
            particle_radius_m=geometry.master_grid.R,
            outer_half_cell_distance_m=outer_distance,
            particle_thermal_conductivity_w_m_k=(config.dry.thermal_conductivity.value_w_m_k),
            binary_gas_interaction_k_wh=config.dry.pore.k_wh,
        )
        surface_pore = replace(config.dry.pore, w_o=oil_labels[last_index])

        if isinstance(boundary, csrf.CondensedSolventReducedFilmPoreBoundary):
            # BC-3: the ratified design's D4 seam, relocated by the BC-2
            # measurement into this sweep.  A sample strictly ABOVE the
            # hexane dew locus evaluates byte-identically to the legacy
            # branch below; a sample at/below the locus is in local
            # equilibrium with condensed solvent, so its composition becomes
            # the chart's own saturation-locus value y_dew(T_sample) and it
            # validates against that construction (the equilibrium gate's
            # inclusive admission) instead of the open gas-only band.  Typed
            # evidence for every condensed sample lands in the boundary's
            # caller-owned journal.  Never a clamp; no tolerance changes.
            def evaluate_surface_flux(
                surface_temperature_k: float,
                conductive_outward_w_m2: float,
            ) -> ct.ReducedFilmSurfaceFluxEvaluation:
                domain_lower_k, domain_upper_k = (
                    config.dry.conditioned_temperature_domain.solver_bounds_k
                )
                if not domain_lower_k < surface_temperature_k < domain_upper_k:
                    raise CutTransportTopologyError(
                        "dry piece temperature left the pressure-conditioned "
                        "open domain; reject, do not pin"
                    )
                sample_evidence = csrf.resolve_condensed_surface_sample(
                    boundary,
                    surface_temperature_k=surface_temperature_k,
                    cell_temperature_k=cells[-1].temperature_k,
                    pressure_pa=config.dry.pressure_pa,
                    pore=surface_pore,
                )
                if sample_evidence is None:
                    _validate_open_dry_band(
                        surface_temperature_k,
                        surface_reference.surface_y_hexane,
                        config,
                        pore=surface_pore,
                    )
                    surface_sample_y_hexane = surface_reference.surface_y_hexane
                else:
                    if not isinstance(
                        config.dry.primitive_band, ct.GasOnlyPrimitiveDomain
                    ):
                        raise CutTransportTopologyError(
                            "the condensed-solvent sweep construction requires "
                            "the exact gas-only chart; a rectangular primitive "
                            "band has no lawful dew-locus construction"
                        )
                    surface_sample_y_hexane = sample_evidence.constructed_y_hexane
                dry_face = _dry_flux(
                    cells[-1].temperature_k,
                    cells[-1].y_hexane,
                    surface_temperature_k,
                    surface_sample_y_hexane,
                    outer_distance,
                    stefan_fluxes[-1],
                    surface_pore,
                    surface_pore,
                    config,
                    conductive_heat_flux_override_w_m2=conductive_outward_w_m2,
                )
                return ct.ReducedFilmSurfaceFluxEvaluation(
                    component=dry_face.component,
                    payload=dry_face,
                    non_fourier_reduced_heat_flux_w_m2=(
                        dry_face.reciprocal_heat_of_transport_flux_w_m2
                    ),
                )

        else:

            def evaluate_surface_flux(
                surface_temperature_k: float,
                conductive_outward_w_m2: float,
            ) -> ct.ReducedFilmSurfaceFluxEvaluation:
                _validate_open_dry_band(
                    surface_temperature_k,
                    surface_reference.surface_y_hexane,
                    config,
                    pore=surface_pore,
                )
                dry_face = _dry_flux(
                    cells[-1].temperature_k,
                    cells[-1].y_hexane,
                    surface_temperature_k,
                    surface_reference.surface_y_hexane,
                    outer_distance,
                    stefan_fluxes[-1],
                    surface_pore,
                    surface_pore,
                    config,
                    conductive_heat_flux_override_w_m2=conductive_outward_w_m2,
                )
                return ct.ReducedFilmSurfaceFluxEvaluation(
                    component=dry_face.component,
                    payload=dry_face,
                    non_fourier_reduced_heat_flux_w_m2=(
                        dry_face.reciprocal_heat_of_transport_flux_w_m2
                    ),
                )

        surface = ct.solve_reduced_film_surface_closure(
            surface_reference,
            evaluate_surface_flux,
            admissible_surface_temperature_bounds_k=(
                _non_fourier_surface_temperature_bounds(
                    config,
                    surface_reference.surface_y_hexane,
                    surface_pore,
                )
                if _uses_nonzero_thermal_factors(config)
                else None
            ),
        )
        if isinstance(boundary, csrf.CondensedSolventReducedFilmPoreBoundary):
            # D4 admits the condensed construction "per trial": an
            # UN-thresholded trial assembly may converge its surface root on
            # the condensed branch and its surface equilibrium is then the
            # construction composition y_dew.  The ACCEPTED (threshold-
            # enforced) assembly requires the gas-only surface the D3 gate
            # certified, so a condensed root refuses typed BEFORE the
            # gas-only equilibrium evaluation below can obscure it.
            if enforce_reduced_film_thresholds:
                csrf.require_root_above_locus(
                    boundary,
                    surface_temperature_k=surface.surface_temperature_k,
                    cell_temperature_k=cells[-1].temperature_k,
                    pressure_pa=config.dry.pressure_pa,
                    pore=surface_pore,
                )
            surface_equilibrium_y_hexane = csrf.condensed_root_equilibrium_composition(
                boundary,
                surface_temperature_k=surface.surface_temperature_k,
                cell_temperature_k=cells[-1].temperature_k,
                pressure_pa=config.dry.pressure_pa,
                pore=surface_pore,
            )
        else:
            surface_equilibrium_y_hexane = surface.surface_y_hexane
        surface_flux = surface.surface_flux_evaluation.payload
        if not isinstance(surface_flux, DryFaceFlux):
            raise RuntimeError("cut surface closure lost its final dry-face payload")
        surface_thermo = cp.evaluate_equilibrium(
            surface.surface_temperature_k,
            config.dry.pressure_pa,
            surface_equilibrium_y_hexane,
            surface_pore,
        )
        surface_energy, audit = ct.reduced_film_surface_energy_and_audit(
            surface,
            water_gas_partial_enthalpy_j_mol=(surface_thermo.water_gas_partial_enthalpy_j_mol),
            hexane_gas_partial_enthalpy_j_mol=(surface_thermo.hexane_gas_partial_enthalpy_j_mol),
            enforce_fast_mass_threshold=enforce_reduced_film_thresholds,
        )
        if surface.fourier_heat_flux_outward_w_m2 != (surface_flux.fourier_heat_flux_w_m2):
            raise RuntimeError("reduced-film closure and particle face lost their Fourier identity")
        if surface.reciprocal_heat_of_transport_flux_outward_w_m2 != (
            surface_flux.reciprocal_heat_of_transport_flux_w_m2
        ):
            raise RuntimeError(
                "reduced-film closure and particle face lost their reciprocal-heat identity"
            )
        if surface_energy.conductive_heat_flux_w_m2 != (surface_flux.reduced_heat_flux_w_m2):
            raise RuntimeError(
                "reduced-film energy did not consume the locally closed reduced heat"
            )
        surface_flux = replace(surface_flux, energy=surface_energy)
        fluxes.append(replace(surface_flux, surface_film_audit=audit))
    return tuple(fluxes)


def _dry_flux(
    left_temperature: float,
    left_y_hexane: float,
    right_temperature: float,
    right_y_hexane: float,
    distance_m: float,
    total_stefan_flux: float,
    left_params: cp.CoupledPoreParams,
    right_params: cp.CoupledPoreParams,
    config: CutTransportConfig,
    *,
    thermodynamic_force_temperature_k: float | None = None,
    temperature_gradient_k_m: float | None = None,
    composition_gradient_m_inv: float | None = None,
    conductive_heat_flux_override_w_m2: float | None = None,
) -> DryFaceFlux:
    if not math.isfinite(distance_m) or distance_m <= 0.0:
        raise CutTransportTopologyError("dry face distance must be positive")
    if not math.isfinite(total_stefan_flux):
        raise ValueError("total Stefan flux must be finite")
    if conductive_heat_flux_override_w_m2 is not None and not math.isfinite(
        conductive_heat_flux_override_w_m2
    ):
        raise ValueError("conductive heat-flux override must be finite")
    uses_exact_trace = (
        temperature_gradient_k_m is not None or composition_gradient_m_inv is not None
    )
    if uses_exact_trace and (
        temperature_gradient_k_m is None or composition_gradient_m_inv is None
    ):
        raise ValueError("both newborn trace gradients must be supplied together")
    uses_actual_path = (
        thermodynamic_force_temperature_k is not None
        and config.moving_interface_composition_force_authority
        is MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
    )
    composition_force_path: acfp.ActualCompositionForcePath | None = None
    if thermodynamic_force_temperature_k is None:
        if uses_exact_trace:
            raise ValueError("newborn trace gradients require an explicit force temperature")
        face_temperature = 0.5 * (left_temperature + right_temperature)
        reconstruction = DryThermodynamicForceReconstruction.SYMMETRIC_ENDPOINT_TEMPERATURE
    else:
        face_temperature = float(thermodynamic_force_temperature_k)
        if not math.isfinite(face_temperature) or not (
            min(left_temperature, right_temperature)
            <= face_temperature
            <= max(left_temperature, right_temperature)
        ):
            raise CutTransportTopologyError(
                "explicit thermodynamic-force temperature must lie between the "
                "two endpoint temperatures"
            )
        if uses_actual_path:
            if face_temperature != left_temperature:
                raise CutTransportTopologyError(
                    "actual-path moving face must identify its physical interface endpoint"
                )
            reconstruction = DryThermodynamicForceReconstruction.MOVING_INTERFACE_ACTUAL_T_Y_PATH
        else:
            reconstruction = (
                DryThermodynamicForceReconstruction.MOVING_INTERFACE_EXACT_TRACE_SECANT
                if uses_exact_trace
                else DryThermodynamicForceReconstruction.MOVING_INTERFACE_T_GAMMA
            )
    if uses_exact_trace:
        assert temperature_gradient_k_m is not None
        assert composition_gradient_m_inv is not None
        if not all(
            math.isfinite(value)
            for value in (
                temperature_gradient_k_m,
                composition_gradient_m_inv,
            )
        ):
            raise ValueError("newborn trace gradients must be finite")
        if left_params != right_params:
            raise CutTransportTopologyError("newborn exact secant requires one material/pore label")
        if face_temperature != left_temperature:
            raise CutTransportTopologyError(
                "newborn exact secant force temperature must be T_Gamma"
            )
        if (
            math.fsum((left_temperature, temperature_gradient_k_m * distance_m))
            != right_temperature
        ):
            raise CutTransportTopologyError(
                "newborn temperature trace is inconsistent with its exact gradient"
            )
        if math.fsum((left_y_hexane, composition_gradient_m_inv * distance_m)) != right_y_hexane:
            raise CutTransportTopologyError(
                "newborn composition trace is inconsistent with its exact gradient"
            )
    if uses_actual_path:
        if left_params != right_params:
            raise CutTransportTopologyError(
                "actual-path moving face requires one unchanged material/pore label"
            )
        _validate_open_dry_band(
            right_temperature,
            right_y_hexane,
            config,
            pore=right_params,
        )
        composition_force_path = acfp.evaluate_actual_composition_force_path(
            left_temperature,
            left_y_hexane,
            right_temperature,
            right_y_hexane,
            config.dry.pressure_pa,
            distance_m,
            left_params,
            closed_hexane_derivative_authority=(
                acfp.ClosedHexaneDerivativeAuthority.NATIVE_LOCAL_ANALYTIC_ENVELOPE
            ),
        )
        left = composition_force_path.left_state
        right = composition_force_path.right_state
    else:
        left = cp.evaluate_equilibrium(
            face_temperature,
            config.dry.pressure_pa,
            left_y_hexane,
            left_params,
        )
        right = cp.evaluate_equilibrium(
            face_temperature,
            config.dry.pressure_pa,
            right_y_hexane,
            right_params,
        )
    face = tc.symmetric_binary_face_state(
        left.binary_gas.molar_density_mol_m3,
        left.y_water,
        left.y_hexane,
        right.binary_gas.molar_density_mol_m3,
        right.y_water,
        right.y_hexane,
    )
    mobility = tc.evaluate_binary_mobility_interval(
        config.dry.binary_diffusivity.interval,
        face,
    ).select(config.dry.binary_diffusivity.fraction)
    common_face_force_mode = config.mass_force_mode in (
        nip.NonisothermalMassForceMode.COMMON_FACE_TEMPERATURE_NEGLIGIBLE_SORET,
        nip.NonisothermalMassForceMode.REFERENCE_INVARIANT_BINARY_THERMAL_DIFFUSION_FACTOR_COUNTERFACTUAL,
        nip.NonisothermalMassForceMode.REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS_COUNTERFACTUAL,
    )
    if common_face_force_mode and uses_actual_path:
        assert composition_force_path is not None
        gas_force_gradient = composition_force_path.gas_force_gradient_m_inv
        retained_force_gradient = composition_force_path.retained_force_gradient_m_inv
        selected_binary_mobility = mobility.as_coupled_pore_selection()
        independent = cp.IndependentMaxwellStefanFlux(
            water_diffusive_flux_mol_m2_s=(
                -selected_binary_mobility.value_mol_m_s * gas_force_gradient
            ),
            hexane_diffusive_flux_mol_m2_s=(
                selected_binary_mobility.value_mol_m_s * gas_force_gradient
            ),
            dimensionless_force_gradient_m_inv=gas_force_gradient,
            mobility=selected_binary_mobility,
            entropy_production_w_m3_k=(
                bg.R
                * selected_binary_mobility.value_mol_m_s
                * gas_force_gradient
                * gas_force_gradient
            ),
            used_explicit_nonisothermal_potentials=False,
        )
        retained = cp.RetainedWaterFlux(
            flux_mol_m2_s=(
                -config.dry.retained_water_mobility.value_mol_m_s * retained_force_gradient
            ),
            dimensionless_force_gradient_m_inv=retained_force_gradient,
            mobility=config.dry.retained_water_mobility,
            entropy_production_w_m3_k=(
                bg.R
                * config.dry.retained_water_mobility.value_mol_m_s
                * retained_force_gradient
                * retained_force_gradient
            ),
            used_explicit_nonisothermal_potentials=False,
        )
    elif common_face_force_mode and uses_exact_trace:
        assert composition_gradient_m_inv is not None
        secants = cp.isothermal_potential_composition_secant_slopes(
            face_temperature,
            config.dry.pressure_pa,
            left_y_hexane,
            composition_gradient_m_inv * distance_m,
            left_params,
        )
        gas_force_gradient = (
            secants.gas_exchange_potential_secant_per_y_hexane * composition_gradient_m_inv
        )
        retained_force_gradient = (
            secants.retained_water_potential_secant_per_y_hexane * composition_gradient_m_inv
        )
        if (
            secants.retained_water_secant_branch is cp.RetainedWaterSecantBranch.RETAINED_CAP
            and retained_force_gradient != 0.0
        ):
            raise CutTransportTopologyError(
                "retained-cap exact trace produced a nonzero retained force"
            )
        selected_binary_mobility = mobility.as_coupled_pore_selection()
        independent = cp.IndependentMaxwellStefanFlux(
            water_diffusive_flux_mol_m2_s=(
                -selected_binary_mobility.value_mol_m_s * gas_force_gradient
            ),
            hexane_diffusive_flux_mol_m2_s=(
                selected_binary_mobility.value_mol_m_s * gas_force_gradient
            ),
            dimensionless_force_gradient_m_inv=gas_force_gradient,
            mobility=selected_binary_mobility,
            entropy_production_w_m3_k=(
                bg.R
                * selected_binary_mobility.value_mol_m_s
                * gas_force_gradient
                * gas_force_gradient
            ),
            used_explicit_nonisothermal_potentials=False,
        )
        retained = cp.RetainedWaterFlux(
            flux_mol_m2_s=(
                -config.dry.retained_water_mobility.value_mol_m_s * retained_force_gradient
            ),
            dimensionless_force_gradient_m_inv=retained_force_gradient,
            mobility=config.dry.retained_water_mobility,
            entropy_production_w_m3_k=(
                bg.R
                * config.dry.retained_water_mobility.value_mol_m_s
                * retained_force_gradient
                * retained_force_gradient
            ),
            used_explicit_nonisothermal_potentials=False,
        )
    elif (
        config.mass_force_mode
        is nip.NonisothermalMassForceMode.COMPLETE_POTENTIAL_REFERENCE_GAUGE_COUNTERFACTUAL
    ):
        actual_left = (
            left
            if uses_actual_path
            else cp.evaluate_equilibrium(
                left_temperature,
                config.dry.pressure_pa,
                left_y_hexane,
                left_params,
            )
        )
        actual_right = (
            right
            if uses_actual_path
            else cp.evaluate_equilibrium(
                right_temperature,
                config.dry.pressure_pa,
                right_y_hexane,
                right_params,
            )
        )
        force_difference = nip.complete_potential_force_difference(
            actual_left,
            actual_right,
            config.complete_potential_reference_gauge,
        )
        gas_force_gradient = force_difference.gas_exchange_right_minus_left / distance_m
        retained_force_gradient = force_difference.retained_water_right_minus_left / distance_m
        selected_binary_mobility = mobility.as_coupled_pore_selection()
        independent = cp.IndependentMaxwellStefanFlux(
            water_diffusive_flux_mol_m2_s=(
                -selected_binary_mobility.value_mol_m_s * gas_force_gradient
            ),
            hexane_diffusive_flux_mol_m2_s=(
                selected_binary_mobility.value_mol_m_s * gas_force_gradient
            ),
            dimensionless_force_gradient_m_inv=gas_force_gradient,
            mobility=selected_binary_mobility,
            entropy_production_w_m3_k=(
                bg.R
                * selected_binary_mobility.value_mol_m_s
                * gas_force_gradient
                * gas_force_gradient
            ),
            used_explicit_nonisothermal_potentials=True,
        )
        retained = cp.RetainedWaterFlux(
            flux_mol_m2_s=(
                -config.dry.retained_water_mobility.value_mol_m_s * retained_force_gradient
            ),
            dimensionless_force_gradient_m_inv=retained_force_gradient,
            mobility=config.dry.retained_water_mobility,
            entropy_production_w_m3_k=(
                bg.R
                * config.dry.retained_water_mobility.value_mol_m_s
                * retained_force_gradient
                * retained_force_gradient
            ),
            used_explicit_nonisothermal_potentials=True,
        )
        if uses_actual_path:
            reconstruction = DryThermodynamicForceReconstruction.MOVING_INTERFACE_COMPLETE_ENDPOINT_POTENTIALS_ACTUAL_STATES
        else:
            reconstruction = (
                DryThermodynamicForceReconstruction.MOVING_INTERFACE_COMPLETE_ENDPOINT_POTENTIALS
                if thermodynamic_force_temperature_k is not None
                else DryThermodynamicForceReconstruction.COMPLETE_ENDPOINT_POTENTIALS
            )
    else:
        independent = cp.independent_maxwell_stefan_flux(
            left,
            right,
            distance_m,
            mobility.as_coupled_pore_selection(),
        )
        retained = cp.retained_water_flux(
            left,
            right,
            distance_m,
            config.dry.retained_water_mobility,
        )
    gas_chemical_force_gradient = independent.dimensionless_force_gradient_m_inv
    retained_chemical_force_gradient = retained.dimensionless_force_gradient_m_inv
    log_temperature_gradient = nip.log_temperature_gradient_m_inv(
        left_temperature,
        right_temperature,
        distance_m,
    )
    binary_thermal_force_gradient = 0.0
    retained_thermal_force_gradient = 0.0
    if nip.reference_invariant_thermal_factor_mode(config.mass_force_mode):
        binary_thermal_force_gradient = nip.binary_thermal_force_gradient_m_inv(
            left_temperature,
            right_temperature,
            distance_m,
            config.binary_thermal_diffusion_factor,
        )
        effective_force_gradient = math.fsum(
            (gas_chemical_force_gradient, binary_thermal_force_gradient)
        )
        independent = cp.IndependentMaxwellStefanFlux(
            water_diffusive_flux_mol_m2_s=(
                -independent.mobility.value_mol_m_s * effective_force_gradient
            ),
            hexane_diffusive_flux_mol_m2_s=(
                independent.mobility.value_mol_m_s * effective_force_gradient
            ),
            dimensionless_force_gradient_m_inv=effective_force_gradient,
            mobility=independent.mobility,
            entropy_production_w_m3_k=(
                bg.R
                * independent.mobility.value_mol_m_s
                * effective_force_gradient
                * effective_force_gradient
            ),
            used_explicit_nonisothermal_potentials=False,
        )
        if (
            config.mass_force_mode
            is nip.NonisothermalMassForceMode.REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS_COUNTERFACTUAL
        ):
            retained_thermal_force_gradient = nip.retained_water_thermal_force_gradient_m_inv(
                left_temperature,
                right_temperature,
                distance_m,
                config.retained_water_thermal_force_factor,
            )
            effective_retained_force_gradient = math.fsum(
                (
                    retained_chemical_force_gradient,
                    retained_thermal_force_gradient,
                )
            )
            retained = cp.RetainedWaterFlux(
                flux_mol_m2_s=(
                    -retained.mobility.value_mol_m_s * effective_retained_force_gradient
                ),
                dimensionless_force_gradient_m_inv=(effective_retained_force_gradient),
                mobility=retained.mobility,
                entropy_production_w_m3_k=(
                    bg.R
                    * retained.mobility.value_mol_m_s
                    * effective_retained_force_gradient
                    * effective_retained_force_gradient
                ),
                used_explicit_nonisothermal_potentials=False,
            )
            if not uses_actual_path:
                reconstruction = (
                    DryThermodynamicForceReconstruction.REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS
                )
        elif not uses_actual_path:
            reconstruction = DryThermodynamicForceReconstruction.REFERENCE_INVARIANT_BINARY_THERMAL_DIFFUSION_FACTOR
    if total_stefan_flux > 0.0:
        advective_y = left.y_hexane
    elif total_stefan_flux < 0.0:
        advective_y = right.y_hexane
    else:
        advective_y = 0.5 * (left.y_hexane + right.y_hexane)
    component = cp.compose_component_fluxes(
        1.0 - advective_y,
        advective_y,
        total_stefan_flux,
        independent.water_diffusive_flux_mol_m2_s,
        retained.flux_mol_m2_s,
    )
    if conductive_heat_flux_override_w_m2 is not None:
        fourier_heat = conductive_heat_flux_override_w_m2
    elif uses_exact_trace:
        assert temperature_gradient_k_m is not None
        fourier_heat = -config.dry.thermal_conductivity.value_w_m_k * temperature_gradient_k_m
    else:
        fourier_heat = (
            config.dry.thermal_conductivity.value_w_m_k
            * (left_temperature - right_temperature)
            / distance_m
        )
    reciprocal_heat = 0.0
    if nip.reference_invariant_thermal_factor_mode(config.mass_force_mode):
        reciprocal_heat = nip.reciprocal_heat_of_transport_flux_w_m2(
            left_temperature_k=left_temperature,
            right_temperature_k=right_temperature,
            gas_counterflux_water_mol_m2_s=(independent.water_diffusive_flux_mol_m2_s),
            binary_thermal_diffusion_factor=(config.binary_thermal_diffusion_factor),
            retained_water_flux_mol_m2_s=retained.flux_mol_m2_s,
            retained_water_thermal_factor=(config.retained_water_thermal_force_factor),
        )
    reduced_heat = (
        fourier_heat if reciprocal_heat == 0.0 else math.fsum((fourier_heat, reciprocal_heat))
    )
    water_state = _upwind_state(component.gas_water_flux_mol_m2_s, left, right)
    hexane_state = _upwind_state(component.gas_hexane_flux_mol_m2_s, left, right)
    retained_state = _upwind_state(component.retained_water_flux_mol_m2_s, left, right)
    energy = cp.component_energy_flux(
        component,
        reduced_heat,
        water_state.water_gas_partial_enthalpy_j_mol,
        hexane_state.hexane_gas_partial_enthalpy_j_mol,
        retained_state.retained_water_enthalpy_j_mol,
    )
    gradient = (
        temperature_gradient_k_m
        if uses_exact_trace
        else (right_temperature - left_temperature) / distance_m
    )
    conduction_entropy = (
        config.dry.thermal_conductivity.value_w_m_k
        * gradient
        * gradient
        / (left_temperature * right_temperature)
    )
    return DryFaceFlux(
        component=component,
        energy=energy,
        independent_water_flux_mol_m2_s=(independent.water_diffusive_flux_mol_m2_s),
        retained_water_flux_mol_m2_s=retained.flux_mol_m2_s,
        conductive_heat_flux_w_m2=fourier_heat,
        gas_entropy_production_w_m3_k=independent.entropy_production_w_m3_k,
        retained_water_entropy_production_w_m3_k=(retained.entropy_production_w_m3_k),
        conduction_entropy_production_w_m3_k=conduction_entropy,
        thermodynamic_force_temperature_k=(None if uses_actual_path else face_temperature),
        thermodynamic_force_reconstruction=reconstruction,
        composition_force_path=composition_force_path,
        mass_force_mode=config.mass_force_mode,
        selected_binary_mobility_mol_m_s=(independent.mobility.value_mol_m_s),
        # A20 section 3.4 recomputability primitive; see the field comment.
        # This is the exact symmetric face state the selected binary mobility was
        # evaluated from at the top of this function.  It cannot be recovered from
        # ``independent.mobility``: every route adapts the binary selection through
        # ``as_coupled_pore_selection()``, which returns a generic
        # ``coupled_pore.MobilitySelection`` carrying only an interval and a
        # fraction, discarding the face.  That discard is precisely why a verifier
        # had to rebuild the state from serialized primitives and landed one ULP
        # away from the production value.
        selected_binary_mobility_face_state=face,
        gas_chemical_force_gradient_m_inv=gas_chemical_force_gradient,
        binary_thermal_force_gradient_m_inv=binary_thermal_force_gradient,
        binary_thermal_diffusion_factor=(config.binary_thermal_diffusion_factor.alpha_t),
        retained_water_chemical_force_gradient_m_inv=(retained_chemical_force_gradient),
        retained_water_thermal_force_gradient_m_inv=(retained_thermal_force_gradient),
        retained_water_thermal_force_factor=(config.retained_water_thermal_force_factor.alpha_ret),
        log_temperature_gradient_m_inv=log_temperature_gradient,
        fourier_heat_flux_w_m2=fourier_heat,
        reciprocal_heat_of_transport_flux_w_m2=reciprocal_heat,
        reduced_heat_flux_w_m2=reduced_heat,
        discrete_conjugate_temperature_k=(
            nip.discrete_conjugate_temperature_k(
                left_temperature,
                right_temperature,
            )
        ),
    )


def _piece_inventories(
    geometry: cg.CutGeometry,
    layout: CutTransportLayout,
    wet_cells: Sequence[ww.WetWaterCellState],
    dry_cells: Sequence[cp.EquilibriumPoreState],
    historical_hexane: Sequence[float],
    config: CutTransportConfig,
) -> _PieceInventories:
    wet_volumes = tuple(geometry.cells[index].wet_volume_m3 for index in layout.wet_cell_indices)
    dry_volumes = tuple(geometry.cells[index].dry_volume_m3 for index in layout.dry_cell_indices)
    wet_water_concentrations = tuple(cell.retained_water_concentration_mol_m3 for cell in wet_cells)
    wet_hexane_concentrations = tuple(
        config.wet.wet.rho_dm_p * historical_hexane[index] / hx.M
        for index in layout.wet_cell_indices
    )
    wet_energy_densities = tuple(cell.energy_density_j_m3 for cell in wet_cells)
    dry_water_concentrations = tuple(cell.total_water_concentration_mol_m3 for cell in dry_cells)
    dry_hexane_concentrations = tuple(cell.total_hexane_concentration_mol_m3 for cell in dry_cells)
    dry_energy_densities = tuple(cell.energy_density_j_m3 for cell in dry_cells)
    return _PieceInventories(
        wet_volumes_m3=wet_volumes,
        dry_volumes_m3=dry_volumes,
        wet_water_concentration_mol_m3=wet_water_concentrations,
        wet_hexane_concentration_mol_m3=wet_hexane_concentrations,
        wet_energy_density_j_m3=wet_energy_densities,
        dry_water_concentration_mol_m3=dry_water_concentrations,
        dry_hexane_concentration_mol_m3=dry_hexane_concentrations,
        dry_energy_density_j_m3=dry_energy_densities,
        wet_water_mol=tuple(
            volume * concentration
            for volume, concentration in zip(wet_volumes, wet_water_concentrations)
        ),
        wet_hexane_mol=tuple(
            volume * concentration
            for volume, concentration in zip(wet_volumes, wet_hexane_concentrations)
        ),
        wet_energy_j=tuple(
            volume * density for volume, density in zip(wet_volumes, wet_energy_densities)
        ),
        dry_water_mol=tuple(
            volume * concentration
            for volume, concentration in zip(dry_volumes, dry_water_concentrations)
        ),
        dry_hexane_mol=tuple(
            volume * concentration
            for volume, concentration in zip(dry_volumes, dry_hexane_concentrations)
        ),
        dry_energy_j=tuple(
            volume * density for volume, density in zip(dry_volumes, dry_energy_densities)
        ),
    )


def _piece_inventory_changes(
    layout: CutTransportLayout,
    old: _PieceInventories,
    new: _PieceInventories,
    swept: cg.SweptCutGeometry,
) -> _PieceInventoryChanges:
    """Use exact swept volumes without subtracting nearly equal inventories.

    For every piece, ``Delta(V*C)`` is evaluated as
    ``V_new*Delta(C) + Delta(V)*C_old``.  The signed volume change comes from
    the shared spherical sweep authority, so arbitrarily small within-cell
    motion does not lose the ALE storage term to cancellation.
    """

    wet_volume_changes = tuple(
        swept.cell_wet_volume_changes_m3[index] for index in layout.wet_cell_indices
    )
    dry_volume_changes = tuple(
        swept.cell_dry_volume_changes_m3[index] for index in layout.dry_cell_indices
    )

    def stable_changes(
        new_volumes: Sequence[float],
        volume_changes: Sequence[float],
        old_concentrations: Sequence[float],
        new_concentrations: Sequence[float],
    ) -> tuple[float, ...]:
        return tuple(
            math.fsum(
                (
                    new_volume * (new_concentration - old_concentration),
                    volume_change * old_concentration,
                )
            )
            for new_volume, volume_change, old_concentration, new_concentration in zip(
                new_volumes,
                volume_changes,
                old_concentrations,
                new_concentrations,
            )
        )

    return _PieceInventoryChanges(
        wet_water_mol=stable_changes(
            new.wet_volumes_m3,
            wet_volume_changes,
            old.wet_water_concentration_mol_m3,
            new.wet_water_concentration_mol_m3,
        ),
        wet_hexane_mol=stable_changes(
            new.wet_volumes_m3,
            wet_volume_changes,
            old.wet_hexane_concentration_mol_m3,
            new.wet_hexane_concentration_mol_m3,
        ),
        wet_energy_j=stable_changes(
            new.wet_volumes_m3,
            wet_volume_changes,
            old.wet_energy_density_j_m3,
            new.wet_energy_density_j_m3,
        ),
        dry_water_mol=stable_changes(
            new.dry_volumes_m3,
            dry_volume_changes,
            old.dry_water_concentration_mol_m3,
            new.dry_water_concentration_mol_m3,
        ),
        dry_hexane_mol=stable_changes(
            new.dry_volumes_m3,
            dry_volume_changes,
            old.dry_hexane_concentration_mol_m3,
            new.dry_hexane_concentration_mol_m3,
        ),
        dry_energy_j=stable_changes(
            new.dry_volumes_m3,
            dry_volume_changes,
            old.dry_energy_density_j_m3,
            new.dry_energy_density_j_m3,
        ),
    )


def _residual_blocks(
    before: CutTransportState,
    geometry: cg.CutGeometry,
    layout: CutTransportLayout,
    changes: _PieceInventoryChanges,
    interface: CutInterfaceState,
    wet_fluxes: Sequence[ww.WetWaterFaceFlux],
    dry_fluxes: Sequence[DryFaceFlux],
    swept: cg.SweptCutGeometry,
    dt_s: float,
    swept_wet_water_donor_concentration_mol_m3: float,
    swept_wet_energy_donor_density_j_m3: float,
) -> tuple[CutResidualBlocks, float]:
    grid = geometry.master_grid
    area_gamma = 4.0 * math.pi * geometry.front.radius_m**2
    q_gamma = swept.interface_swept_volume_rate_m3_s

    # O9a ALE FRONT-DONOR CORRECTION, ruled 2026-08-21 in
    # docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md (Class-B
    # edit to this frozen closure, authorized there together with the O9c
    # re-earn campaign).  The receding front consumes wet material at that
    # material's OWN BULK LOADING, so the ALE sweep is debited at the consumed
    # piece's bulk retained-water concentration and no longer at the
    # zero-volume equilibrium interface trace.
    #
    # REMOVED ASSUMPTION (this is a removal, not an addition): booking the
    # sweep at the trace implicitly asserted that retained water arriving at
    # the front has PRE-DRAINED to the dry-pore equilibrium loading before
    # conversion.  The traverse measured Peclet 3.4e3 (D/v = 2.6e-7 m against
    # an 8.85e-4 m particle), which forbids that migration; the shortfall was
    # re-averaged into the shrinking cut piece instead of being released.
    #
    # ROLE SEPARATION - only the MASS-BOOKKEEPING donor moves.  The interface
    # trace keeps its other role in full: it is still the local-equilibrium
    # potential the front-face diffusive flux is evaluated against
    # (``_wet_face_fluxes`` passes ``interface.wet.retained_water_loading`` and
    # ``interface.retained_water_trace_capacity_dual_over_rt`` as that face's
    # right state) and ``_retained_water_interface_trace_kkt`` and its
    # potential-continuity check are untouched.
    #
    # The released excess leaves through the UNCHANGED Rankine-Hugoniot
    # balance below: ``rh_water`` is the signed expansion of ``dry_h - wet_h``,
    # so its ``C_wet`` is this same donor by construction and the jump row
    # keeps exactly its existing structure - it simply now demands the larger
    # dry-side water flux that the richer arriving material implies.
    #
    # The donor is the consumed piece's bulk concentration at the time level
    # ``_piece_inventory_changes`` books the matching geometric storage term
    # at; see the call site in ``_assemble_backward_euler`` for that pairing.
    wet_water_donor_concentration = swept_wet_water_donor_concentration_mol_m3
    # O10a ENERGY-DONOR COMPLETION, ruled 2026-08-21 in
    # docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md.
    #
    # PHYSICS GROUND, as recorded in that ruling: two-phase moving-boundary
    # (Rankine-Hugoniot) jump conditions must convect mass and enthalpy OF THE
    # SAME MATERIAL STATE.  That is an internal-consistency requirement of the
    # balance laws themselves, not an empirical claim; WHICH material state it
    # is - the bulk, undrained one - is the O9-measured fact (Peclet 3.4e3
    # forbids pre-drainage).  With the mass donor already moved to the consumed
    # piece's bulk loading, leaving the energy donor on the zero-volume
    # equilibrium trace would have the same swept material carry water mass at
    # the bulk loading and water enthalpy at the trace loading.
    #
    # The donor is therefore the consumed piece's OWN bulk energy density, read
    # at the SAME time level as its mass donor - which is the level
    # ``_piece_inventory_changes`` books the paired geometric energy storage
    # term ``Delta(V)*rho_e_old`` at, so the energy mis-booking term cancels
    # term by term exactly as the mass one does.
    #
    # ROLE SEPARATION IS UNCHANGED.  The interface trace keeps its
    # potential/flux role in full on the energy side too: the front-face
    # conductive/advective energy flux ``wet_fluxes[-1].total_energy_flux_w_m2``
    # is still evaluated against the trace, and only the ALE convected density
    # moves.
    wet_energy_donor_density = swept_wet_energy_donor_density_j_m3
    wet_water_h = math.fsum(
        (
            area_gamma * wet_fluxes[-1].retained_water_flux_mol_m2_s,
            -q_gamma * wet_water_donor_concentration,
        )
    )
    wet_hexane_concentration = (
        before.config.wet.wet.rho_dm_p * interface.historical_hexane_loading / hx.M
    )
    wet_hexane_h = -q_gamma * wet_hexane_concentration
    wet_energy_h = math.fsum(
        (
            area_gamma * wet_fluxes[-1].total_energy_flux_w_m2,
            # O10a: the swept material's enthalpy travels at the same (bulk,
            # old-time) state as its mass; see the block above.
            -q_gamma * wet_energy_donor_density,
        )
    )
    dry_component = dry_fluxes[0].component
    dry_water_h = math.fsum(
        (
            area_gamma * dry_component.conserved_water_flux_mol_m2_s,
            -q_gamma * interface.dry.total_water_concentration_mol_m3,
        )
    )
    dry_hexane_h = math.fsum(
        (
            area_gamma * dry_component.conserved_hexane_flux_mol_m2_s,
            -q_gamma * interface.dry.total_hexane_concentration_mol_m3,
        )
    )
    dry_energy_h = math.fsum(
        (
            area_gamma * dry_fluxes[0].energy.total_energy_flux_w_m2,
            -q_gamma * interface.dry.energy_density_j_m3,
        )
    )

    wet_water_rates = [0.0]
    wet_energy_rates = [0.0]
    for face_position, cell_index in enumerate(layout.wet_cell_indices[1:], start=1):
        area = grid.areas[cell_index]
        wet_water_rates.append(area * wet_fluxes[face_position].retained_water_flux_mol_m2_s)
        wet_energy_rates.append(area * wet_fluxes[face_position].total_energy_flux_w_m2)
    wet_water_rates.append(wet_water_h)
    wet_energy_rates.append(wet_energy_h)
    wet_water_residuals = tuple(
        math.fsum(
            (
                change / dt_s,
                wet_water_rates[index + 1],
                -wet_water_rates[index],
            )
        )
        for index, change in enumerate(changes.wet_water_mol)
    )
    wet_energy_residuals = tuple(
        math.fsum(
            (
                change / dt_s,
                wet_energy_rates[index + 1],
                -wet_energy_rates[index],
            )
        )
        for index, change in enumerate(changes.wet_energy_j)
    )

    dry_water_rates = [dry_water_h]
    dry_hexane_rates = [dry_hexane_h]
    dry_energy_rates = [dry_energy_h]
    for face_position in range(1, layout.dry_face_count):
        if face_position == layout.dry_piece_count:
            area = grid.areas[-1]
        else:
            area = grid.areas[layout.dry_cell_indices[face_position]]
        flux = dry_fluxes[face_position]
        dry_water_rates.append(area * flux.component.conserved_water_flux_mol_m2_s)
        dry_hexane_rates.append(area * flux.component.conserved_hexane_flux_mol_m2_s)
        dry_energy_rates.append(area * flux.energy.total_energy_flux_w_m2)

    def regional_residuals(
        inventory_changes: Sequence[float],
        rates: Sequence[float],
    ) -> tuple[float, ...]:
        return tuple(
            math.fsum(
                (
                    change / dt_s,
                    rates[index + 1],
                    -rates[index],
                )
            )
            for index, change in enumerate(inventory_changes)
        )

    dry_water_residuals = regional_residuals(changes.dry_water_mol, dry_water_rates)
    dry_hexane_residuals = regional_residuals(changes.dry_hexane_mol, dry_hexane_rates)
    dry_energy_residuals = regional_residuals(changes.dry_energy_j, dry_energy_rates)
    wet_hexane_identity = math.fsum(
        (
            math.fsum(changes.wet_hexane_mol) / dt_s,
            wet_hexane_h,
        )
    )
    # Form the jump rows directly from their signed face/sweep terms.  The
    # mathematically equivalent ``dry_h - wet_h`` route first adds two large
    # sweep contributions on each side and then subtracts them, needlessly
    # losing digits precisely near a converged RH root.  This signed expansion
    # changes no equation or flux; it only preserves the cancellation already
    # present in the common ALE frame.
    rh_water = math.fsum(
        (
            area_gamma * dry_component.conserved_water_flux_mol_m2_s,
            -area_gamma * wet_fluxes[-1].retained_water_flux_mol_m2_s,
            -q_gamma * interface.dry.total_water_concentration_mol_m3,
            # O9a: the same corrected wet-side donor as ``wet_water_h`` above,
            # so this row stays exactly ``dry_h - wet_h`` and the global water
            # telescoping identity is untouched.
            q_gamma * wet_water_donor_concentration,
        )
    )
    rh_hexane = math.fsum(
        (
            area_gamma * dry_component.conserved_hexane_flux_mol_m2_s,
            -q_gamma * interface.dry.total_hexane_concentration_mol_m3,
            q_gamma * wet_hexane_concentration,
        )
    )
    rh_energy = math.fsum(
        (
            area_gamma * dry_fluxes[0].energy.total_energy_flux_w_m2,
            -area_gamma * wet_fluxes[-1].total_energy_flux_w_m2,
            -q_gamma * interface.dry.energy_density_j_m3,
            # O10a: the same corrected wet-side energy donor as
            # ``wet_energy_h`` above, so this row stays exactly
            # ``dry_h - wet_h`` and the global energy telescoping identity is
            # untouched.  Because this donor's water content is now the SAME
            # bulk concentration ``rh_water`` convects, the water caloric datum
            # is again an exact gauge freedom of this row - the acceptance
            # instrument the O10a ruling names.
            q_gamma * wet_energy_donor_density,
        )
    )
    return (
        CutResidualBlocks(
            wet_water_mol_s=wet_water_residuals,
            wet_energy_w=wet_energy_residuals,
            dry_water_mol_s=dry_water_residuals,
            dry_hexane_mol_s=dry_hexane_residuals,
            dry_energy_w=dry_energy_residuals,
            rh_water_mol_s=rh_water,
            rh_hexane_mol_s=rh_hexane,
            rh_energy_w=rh_energy,
        ),
        wet_hexane_identity,
    )


def _apply_vanishing_wet_cut_treatment(
    geometry: cg.CutGeometry,
    layout: CutTransportLayout,
    wet_cells: Sequence[ww.WetWaterCellState],
    interface: CutInterfaceState,
    wet_fluxes: Sequence[ww.WetWaterFaceFlux],
    original: CutResidualBlocks,
    config: CutTransportConfig,
) -> tuple[CutResidualBlocks, VanishingWetCutTreatmentAudit]:
    """Apply the fixed conservative Amendment-12 final-wet-pair scheme.

    This is the declared finite-volume operator, not a fallback selected after
    nonlinear failure.  It is evaluated on every residual call.  Only the last
    two wet retained-water and common-energy rows can change, and their pair
    sums remain equal to the original rows to roundoff.
    """

    cut_index = layout.cut_cell_index
    partition = geometry.cells[cut_index]
    full_volume = geometry.master_grid.volumes[cut_index]
    fraction = partition.wet_volume_m3 / full_volume
    authority = config.vanishing_wet_cut_threshold_authority
    geometric_weight = vanishing_wet_cut_weight(fraction, authority)
    full_neighbor_available = layout.wet_piece_count >= 2
    original_cut_water = original.wet_water_mol_s[-1]
    original_cut_energy = original.wet_energy_w[-1]

    if authority.target_wet_volume_fraction is None:
        audit = VanishingWetCutTreatmentAudit(
            threshold_authority=authority,
            target_wet_volume_fraction=None,
            treatment_enabled=False,
            wet_volume_fraction=fraction,
            geometric_weight=0.0,
            applied_weight=0.0,
            full_wet_neighbor_available=full_neighbor_available,
            active=False,
            activation_count=0,
            cut_master_cell_index=cut_index,
            adjacent_wet_piece_position=None,
            adjacent_master_cell_index=None,
            adjacent_to_interface_distance_m=None,
            interface_area_m2=None,
            original_adjacent_water_residual_mol_s=None,
            original_cut_water_residual_mol_s=original_cut_water,
            transformed_adjacent_water_residual_mol_s=None,
            transformed_cut_water_residual_mol_s=original_cut_water,
            original_adjacent_common_energy_residual_w=None,
            original_cut_common_energy_residual_w=original_cut_energy,
            transformed_adjacent_common_energy_residual_w=None,
            transformed_cut_common_energy_residual_w=original_cut_energy,
            original_cut_interface_water_flux_mol_m2_s=None,
            reconstructed_adjacent_interface_water_flux_mol_m2_s=None,
            original_cut_interface_common_energy_flux_w_m2=None,
            reconstructed_adjacent_interface_common_energy_flux_w_m2=None,
            water_flux_difference_area_weighted_mol_s=None,
            common_energy_flux_difference_area_weighted_w=None,
            original_water_pair_sum_mol_s=None,
            transformed_water_pair_sum_mol_s=None,
            water_pair_sum_difference_mol_s=0.0,
            original_common_energy_pair_sum_w=None,
            transformed_common_energy_pair_sum_w=None,
            common_energy_pair_sum_difference_w=0.0,
        )
        return original, audit

    if geometric_weight == 0.0:
        audit = VanishingWetCutTreatmentAudit(
            threshold_authority=authority,
            target_wet_volume_fraction=authority.target_wet_volume_fraction,
            treatment_enabled=True,
            wet_volume_fraction=fraction,
            geometric_weight=0.0,
            applied_weight=0.0,
            full_wet_neighbor_available=full_neighbor_available,
            active=False,
            activation_count=0,
            cut_master_cell_index=cut_index,
            adjacent_wet_piece_position=None,
            adjacent_master_cell_index=None,
            adjacent_to_interface_distance_m=None,
            interface_area_m2=None,
            original_adjacent_water_residual_mol_s=None,
            original_cut_water_residual_mol_s=original_cut_water,
            transformed_adjacent_water_residual_mol_s=None,
            transformed_cut_water_residual_mol_s=original_cut_water,
            original_adjacent_common_energy_residual_w=None,
            original_cut_common_energy_residual_w=original_cut_energy,
            transformed_adjacent_common_energy_residual_w=None,
            transformed_cut_common_energy_residual_w=original_cut_energy,
            original_cut_interface_water_flux_mol_m2_s=None,
            reconstructed_adjacent_interface_water_flux_mol_m2_s=None,
            original_cut_interface_common_energy_flux_w_m2=None,
            reconstructed_adjacent_interface_common_energy_flux_w_m2=None,
            water_flux_difference_area_weighted_mol_s=None,
            common_energy_flux_difference_area_weighted_w=None,
            original_water_pair_sum_mol_s=None,
            transformed_water_pair_sum_mol_s=None,
            water_pair_sum_difference_mol_s=0.0,
            original_common_energy_pair_sum_w=None,
            transformed_common_energy_pair_sum_w=None,
            common_energy_pair_sum_difference_w=0.0,
        )
        return original, audit

    applied_weight = geometric_weight if full_neighbor_available else 0.0

    if not full_neighbor_available:
        audit = VanishingWetCutTreatmentAudit(
            threshold_authority=authority,
            target_wet_volume_fraction=authority.target_wet_volume_fraction,
            treatment_enabled=True,
            wet_volume_fraction=fraction,
            geometric_weight=geometric_weight,
            applied_weight=0.0,
            full_wet_neighbor_available=False,
            active=False,
            activation_count=0,
            cut_master_cell_index=cut_index,
            adjacent_wet_piece_position=None,
            adjacent_master_cell_index=None,
            adjacent_to_interface_distance_m=None,
            interface_area_m2=None,
            original_adjacent_water_residual_mol_s=None,
            original_cut_water_residual_mol_s=original_cut_water,
            transformed_adjacent_water_residual_mol_s=None,
            transformed_cut_water_residual_mol_s=original_cut_water,
            original_adjacent_common_energy_residual_w=None,
            original_cut_common_energy_residual_w=original_cut_energy,
            transformed_adjacent_common_energy_residual_w=None,
            transformed_cut_common_energy_residual_w=original_cut_energy,
            original_cut_interface_water_flux_mol_m2_s=None,
            reconstructed_adjacent_interface_water_flux_mol_m2_s=None,
            original_cut_interface_common_energy_flux_w_m2=None,
            reconstructed_adjacent_interface_common_energy_flux_w_m2=None,
            water_flux_difference_area_weighted_mol_s=None,
            common_energy_flux_difference_area_weighted_w=None,
            original_water_pair_sum_mol_s=None,
            transformed_water_pair_sum_mol_s=None,
            water_pair_sum_difference_mol_s=0.0,
            original_common_energy_pair_sum_w=None,
            transformed_common_energy_pair_sum_w=None,
            common_energy_pair_sum_difference_w=0.0,
        )
        return original, audit

    adjacent_position = layout.wet_piece_count - 2
    adjacent_cell_index = layout.wet_cell_indices[adjacent_position]
    adjacent_center = _wet_piece_center(
        geometry.cells[adjacent_cell_index],
        geometry.front.radius_m,
    )
    distance = geometry.front.radius_m - adjacent_center
    adjacent = wet_cells[adjacent_position]
    reconstructed = _wet_flux(
        adjacent.temperature_k,
        adjacent.retained_water_loading,
        interface.temperature_k,
        interface.wet.retained_water_loading,
        distance,
        config,
        "adjacent full wet cell",
        "interface trace",
        left_capacity_dual_over_rt=(adjacent.retained_water_capacity_dual_over_rt),
        right_capacity_dual_over_rt=(interface.retained_water_trace_capacity_dual_over_rt),
    )
    area = 4.0 * math.pi * geometry.front.radius_m**2
    actual = wet_fluxes[-1]
    water_difference = area * math.fsum(
        (
            actual.retained_water_flux_mol_m2_s,
            -reconstructed.retained_water_flux_mol_m2_s,
        )
    )
    energy_difference = area * math.fsum(
        (
            actual.total_energy_flux_w_m2,
            -reconstructed.total_energy_flux_w_m2,
        )
    )

    wet_water = list(original.wet_water_mol_s)
    wet_energy = list(original.wet_energy_w)
    original_adjacent_water = wet_water[-2]
    original_adjacent_energy = wet_energy[-2]
    original_water_pair = math.fsum((original_adjacent_water, original_cut_water))
    original_energy_pair = math.fsum((original_adjacent_energy, original_cut_energy))
    wet_water[-2] = math.fsum(
        (
            original_adjacent_water,
            applied_weight * (original_cut_water - water_difference),
        )
    )
    wet_water[-1] = math.fsum(
        (
            (1.0 - applied_weight) * original_cut_water,
            applied_weight * water_difference,
        )
    )
    wet_energy[-2] = math.fsum(
        (
            original_adjacent_energy,
            applied_weight * (original_cut_energy - energy_difference),
        )
    )
    wet_energy[-1] = math.fsum(
        (
            (1.0 - applied_weight) * original_cut_energy,
            applied_weight * energy_difference,
        )
    )
    transformed_water_pair = math.fsum((wet_water[-2], wet_water[-1]))
    transformed_energy_pair = math.fsum((wet_energy[-2], wet_energy[-1]))
    transformed = replace(
        original,
        wet_water_mol_s=tuple(wet_water),
        wet_energy_w=tuple(wet_energy),
    )
    audit = VanishingWetCutTreatmentAudit(
        threshold_authority=authority,
        target_wet_volume_fraction=authority.target_wet_volume_fraction,
        treatment_enabled=True,
        wet_volume_fraction=fraction,
        geometric_weight=geometric_weight,
        applied_weight=applied_weight,
        full_wet_neighbor_available=True,
        active=applied_weight > 0.0,
        activation_count=int(applied_weight > 0.0),
        cut_master_cell_index=cut_index,
        adjacent_wet_piece_position=adjacent_position,
        adjacent_master_cell_index=adjacent_cell_index,
        adjacent_to_interface_distance_m=distance,
        interface_area_m2=area,
        original_adjacent_water_residual_mol_s=original_adjacent_water,
        original_cut_water_residual_mol_s=original_cut_water,
        transformed_adjacent_water_residual_mol_s=wet_water[-2],
        transformed_cut_water_residual_mol_s=wet_water[-1],
        original_adjacent_common_energy_residual_w=original_adjacent_energy,
        original_cut_common_energy_residual_w=original_cut_energy,
        transformed_adjacent_common_energy_residual_w=wet_energy[-2],
        transformed_cut_common_energy_residual_w=wet_energy[-1],
        original_cut_interface_water_flux_mol_m2_s=(actual.retained_water_flux_mol_m2_s),
        reconstructed_adjacent_interface_water_flux_mol_m2_s=(
            reconstructed.retained_water_flux_mol_m2_s
        ),
        original_cut_interface_common_energy_flux_w_m2=(actual.total_energy_flux_w_m2),
        reconstructed_adjacent_interface_common_energy_flux_w_m2=(
            reconstructed.total_energy_flux_w_m2
        ),
        water_flux_difference_area_weighted_mol_s=water_difference,
        common_energy_flux_difference_area_weighted_w=energy_difference,
        original_water_pair_sum_mol_s=original_water_pair,
        transformed_water_pair_sum_mol_s=transformed_water_pair,
        water_pair_sum_difference_mol_s=math.fsum((transformed_water_pair, -original_water_pair)),
        original_common_energy_pair_sum_w=original_energy_pair,
        transformed_common_energy_pair_sum_w=transformed_energy_pair,
        common_energy_pair_sum_difference_w=math.fsum(
            (transformed_energy_pair, -original_energy_pair)
        ),
    )
    return transformed, audit


def _ledger(
    changes: _PieceInventoryChanges,
    residuals: CutResidualBlocks,
    wet_hexane_identity: float,
    surface: DryFaceFlux,
    surface_area_m2: float,
    swept: cg.SweptCutGeometry,
    dt_s: float,
    all_fluxes: Sequence[object],
    pressure_before_pa: float,
    pressure_after_pa: float,
) -> CutTransportLedger:
    water_blocks = math.fsum(
        (
            *residuals.wet_water_mol_s,
            *residuals.dry_water_mol_s,
            residuals.rh_water_mol_s,
        )
    )
    hexane_blocks = math.fsum(
        (
            wet_hexane_identity,
            *residuals.dry_hexane_mol_s,
            residuals.rh_hexane_mol_s,
        )
    )
    energy_blocks = math.fsum(
        (
            *residuals.wet_energy_w,
            *residuals.dry_energy_w,
            residuals.rh_energy_w,
        )
    )
    water_direct = math.fsum(
        (
            math.fsum((*changes.wet_water_mol, *changes.dry_water_mol)) / dt_s,
            surface_area_m2 * surface.component.conserved_water_flux_mol_m2_s,
        )
    )
    hexane_direct = math.fsum(
        (
            math.fsum((*changes.wet_hexane_mol, *changes.dry_hexane_mol)) / dt_s,
            surface_area_m2 * surface.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    energy_direct = math.fsum(
        (
            math.fsum((*changes.wet_energy_j, *changes.dry_energy_j)) / dt_s,
            surface_area_m2 * surface.energy.total_energy_flux_w_m2,
        )
    )
    entropy_values: list[float] = []
    for flux in all_fluxes:
        if isinstance(flux, ww.WetWaterFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
        elif isinstance(flux, DryFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
    return CutTransportLedger(
        pressure_before_pa=pressure_before_pa,
        pressure_after_pa=pressure_after_pa,
        finite_pressure_transition=pressure_after_pa != pressure_before_pa,
        wet_hexane_sweep_identity_mol_s=wet_hexane_identity,
        water_global_from_blocks_mol_s=water_blocks,
        water_global_direct_mol_s=water_direct,
        water_telescoping_error_mol_s=water_blocks - water_direct,
        hexane_global_from_blocks_mol_s=hexane_blocks,
        hexane_global_direct_mol_s=hexane_direct,
        hexane_telescoping_error_mol_s=hexane_blocks - hexane_direct,
        energy_global_from_blocks_w=energy_blocks,
        energy_global_direct_w=energy_direct,
        energy_telescoping_error_w=energy_blocks - energy_direct,
        geometry_cell_gcl_residual_m3=swept.max_cell_gcl_residual_m3,
        minimum_entropy_production_w_m3_k=min(entropy_values, default=0.0),
    )


@dataclass(frozen=True)
class _AuditedInterfaceCompositionRoot:
    y_hexane: float
    log_fugacity_residual: float
    dry: cp.EquilibriumPoreState
    audit: InterfaceCompositionRootAudit


_UPPER_LOG_LIQUID_SIDE = "negative_log_fugacity_residual_liquid_side"
_UPPER_PRODUCT_SUPERSATURATED = "product_activity_supersaturated_on_log_tolerance_plateau"
_UPPER_RETURNED_CERTIFIED_ENDPOINT = "returned_fully_certified_upper_endpoint"


def _solve_interface_composition(
    temperature_k: float,
    config: CutTransportConfig,
    oil_fraction_label: float,
) -> tuple[float, float]:
    """Compatibility view of the audited representable interface root."""

    root = _solve_interface_composition_with_audit(
        temperature_k,
        config,
        oil_fraction_label,
    )
    return root.y_hexane, root.log_fugacity_residual


def _solve_interface_composition_with_audit(
    temperature_k: float,
    config: CutTransportConfig,
    oil_fraction_label: float,
) -> _AuditedInterfaceCompositionRoot:
    """Return a representable non-supersaturated PHY-039 interface root.

    The log-fugacity residual and its caller-owned tolerance remain the
    equation and acceptance authority.  Binary64 can nevertheless round that
    residual to zero for several adjacent compositions while the product-form
    activity arithmetic used by the coupled-pore state still gives
    ``a_h > 1`` by a few ulps.  Resolve only that discrete plateau inside the
    bisection: such a point remains the upper (liquid-side) endpoint and is
    never returned or repaired after the solve.
    """

    lower, upper = config.interface_composition.y_hexane_bounds
    tolerance = config.interface_composition.log_fugacity_tolerance
    pore = replace(config.dry.pore, w_o=oil_fraction_label)

    def residual(y_hexane: float) -> float:
        return front.log_fugacity_residual(
            temperature_k,
            config.dry.pressure_pa,
            1.0 - y_hexane,
            y_hexane,
            k_wh=config.dry.pore.k_wh,
        )[0]

    def certify_representable_trace(
        y_hexane: float,
    ) -> tuple[cp.EquilibriumPoreState | None, float]:
        """Certify the exact candidate or identify only hexane supersaturation."""

        try:
            dry = cp.evaluate_equilibrium(
                temperature_k,
                config.dry.pressure_pa,
                y_hexane,
                pore,
            )
        except cp.HexaneSupersaturationTopologyError:
            # This typed outcome is raised only after the authoritative pore
            # evaluator has checked parameters, primitives, mixture properties,
            # and water topology.  It is the sole failure that Amendment 18
            # retains as the semantic upper side of a near-root plateau.
            _, hexane_activity = cp._gas_component_activities(  # noqa: SLF001
                temperature_k,
                config.dry.pressure_pa,
                y_hexane,
                pore,
            )
            if not math.isfinite(hexane_activity) or hexane_activity <= 1.0:
                raise RuntimeError(
                    "typed hexane-supersaturation outcome disagrees with product activity"
                )
            return None, hexane_activity
        return dry, dry.hexane_activity

    def make_root(
        y_hexane: float,
        root_residual: float,
        dry: cp.EquilibriumPoreState,
        product_hexane_activity: float,
        *,
        final_lower: float,
        final_upper: float,
        iteration_count: int,
        rounding_plateau: bool,
        no_distinct_midpoint: bool,
        upper_classification: str,
        upper_residual: float,
        upper_product_activity: float | None,
    ) -> _AuditedInterfaceCompositionRoot:
        if (
            dry.temperature_k != temperature_k
            or dry.pressure_pa != config.dry.pressure_pa
            or dry.y_hexane != y_hexane
            or dry.hexane_activity != product_hexane_activity
        ):
            raise RuntimeError("interface root audit lost its exact coupled-pore certificate")
        audit = InterfaceCompositionRootAudit(
            caller_lower_y_hexane=lower,
            caller_upper_y_hexane=upper,
            final_lower_y_hexane=final_lower,
            final_upper_y_hexane=final_upper,
            returned_y_hexane=y_hexane,
            returned_log_fugacity_residual=root_residual,
            returned_product_hexane_activity=product_hexane_activity,
            actual_oil_fraction_label=oil_fraction_label,
            log_fugacity_tolerance=tolerance,
            bisection_iteration_count=iteration_count,
            rounding_plateau_encountered=rounding_plateau,
            no_distinct_midpoint_encountered=no_distinct_midpoint,
            final_upper_classification=upper_classification,
            final_upper_log_fugacity_residual=upper_residual,
            final_upper_product_hexane_activity=upper_product_activity,
        )
        return _AuditedInterfaceCompositionRoot(
            y_hexane,
            root_residual,
            dry,
            audit,
        )

    r_lower = residual(lower)
    r_upper = residual(upper)
    lower_dry: cp.EquilibriumPoreState | None = None
    lower_activity: float | None = None
    if 0.0 <= r_lower <= tolerance:
        lower_dry, lower_activity = certify_representable_trace(lower)
        if lower_dry is None:
            raise CutTransportTopologyError(
                "interface bracket lower endpoint is not a representable non-supersaturated root"
            )

    rounding_plateau_encountered = False
    upper_product_activity: float | None = None
    upper_dry: cp.EquilibriumPoreState | None = None
    if 0.0 <= r_upper <= tolerance:
        upper_dry, upper_product_activity = certify_representable_trace(upper)
        rounding_plateau_encountered = upper_dry is None

    if r_lower < 0.0:
        if r_upper <= 0.0:
            raise CutTransportTopologyError(
                "interface n-hexane fugacity root is not bracketed by the caller's "
                "declared interface interval"
            )
        raise CutTransportTopologyError(
            "interface fugacity bracket has the wrong branch orientation"
        )
    if upper_dry is not None and upper_product_activity is not None:
        return make_root(
            upper,
            r_upper,
            upper_dry,
            upper_product_activity,
            final_lower=lower,
            final_upper=upper,
            iteration_count=0,
            rounding_plateau=False,
            no_distinct_midpoint=False,
            upper_classification=_UPPER_RETURNED_CERTIFIED_ENDPOINT,
            upper_residual=r_upper,
            upper_product_activity=upper_product_activity,
        )
    if lower_dry is not None and lower_activity is not None and not (rounding_plateau_encountered):
        if r_upper >= 0.0:
            raise CutTransportTopologyError(
                "interface n-hexane fugacity root is not bracketed by the caller's "
                "declared interface interval"
            )
        return make_root(
            lower,
            r_lower,
            lower_dry,
            lower_activity,
            final_lower=lower,
            final_upper=upper,
            iteration_count=0,
            rounding_plateau=False,
            no_distinct_midpoint=False,
            upper_classification=_UPPER_LOG_LIQUID_SIDE,
            upper_residual=r_upper,
            upper_product_activity=None,
        )
    product_supersaturated_semantic_upper = (
        rounding_plateau_encountered and 0.0 <= r_upper <= tolerance
    )
    if r_lower * r_upper > 0.0 and not product_supersaturated_semantic_upper:
        raise CutTransportTopologyError(
            "interface n-hexane fugacity root is not bracketed by the caller's "
            "declared interface interval"
        )
    # Residual decreases with y_h on the admitted binary branch.  Preserve the
    # gas-side endpoint (f_g <= f_l) so roundoff cannot manufacture a_h>1.
    if r_upper > 0.0 and not product_supersaturated_semantic_upper:
        raise CutTransportTopologyError(
            "interface fugacity bracket has the wrong branch orientation"
        )
    lo, hi = lower, upper
    r_lo = r_lower
    r_hi = r_upper
    upper_classification = (
        _UPPER_PRODUCT_SUPERSATURATED if rounding_plateau_encountered else _UPPER_LOG_LIQUID_SIDE
    )
    representable_lo_in_tolerance = lower_dry is not None
    representable_lo_activity = lower_activity
    representable_lo_dry = lower_dry
    no_distinct_midpoint_encountered = False
    bisection_iteration_count = 0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            no_distinct_midpoint_encountered = True
            break
        bisection_iteration_count += 1
        r_mid = residual(mid)
        if r_mid > tolerance:
            lo, r_lo = mid, r_mid
            representable_lo_in_tolerance = False
            representable_lo_activity = None
            representable_lo_dry = None
        elif r_mid < 0.0:
            hi, r_hi = mid, r_mid
            upper_product_activity = None
            upper_classification = _UPPER_LOG_LIQUID_SIDE
        else:
            mid_dry, mid_activity = certify_representable_trace(mid)
            if mid_dry is not None:
                if not rounding_plateau_encountered:
                    return make_root(
                        mid,
                        r_mid,
                        mid_dry,
                        mid_activity,
                        final_lower=mid,
                        final_upper=hi,
                        iteration_count=bisection_iteration_count,
                        rounding_plateau=False,
                        no_distinct_midpoint=False,
                        upper_classification=upper_classification,
                        upper_residual=r_hi,
                        upper_product_activity=upper_product_activity,
                    )
                lo, r_lo = mid, r_mid
                representable_lo_in_tolerance = True
                representable_lo_activity = mid_activity
                representable_lo_dry = mid_dry
                continue
            # The unchanged log equation regards this float as gas-side, but
            # the authoritative product-form activity rounds just above one.
            # Keep it inside the live upper bracket; do not retreat or clip.
            hi, r_hi = mid, r_mid
            upper_product_activity = mid_activity
            upper_classification = _UPPER_PRODUCT_SUPERSATURATED
            rounding_plateau_encountered = True
    if (
        rounding_plateau_encountered
        and no_distinct_midpoint_encountered
        and math.nextafter(lo, math.inf) == hi
        and representable_lo_in_tolerance
        and representable_lo_activity is not None
        and representable_lo_dry is not None
        and 0.0 <= r_lo <= tolerance
    ):
        return make_root(
            lo,
            r_lo,
            representable_lo_dry,
            representable_lo_activity,
            final_lower=lo,
            final_upper=hi,
            iteration_count=bisection_iteration_count,
            rounding_plateau=True,
            no_distinct_midpoint=no_distinct_midpoint_encountered,
            upper_classification=upper_classification,
            upper_residual=r_hi,
            upper_product_activity=upper_product_activity,
        )
    if r_lo > tolerance:
        raise CutTransportTopologyError(
            "interface fugacity root did not contract to its numerical tolerance"
        )
    raise CutTransportTopologyError(
        "interface fugacity root has no jointly certified representable "
        "non-supersaturated float at the unchanged tolerance"
    )


def _validate_open_dry_band(
    temperature_k: float,
    y_hexane: float,
    config: CutTransportConfig,
    *,
    pore: cp.CoupledPoreParams | None = None,
) -> None:
    if not math.isfinite(temperature_k) or not math.isfinite(y_hexane):
        raise ValueError("dry piece primitives must be finite")
    t_lo, t_hi = config.dry.conditioned_temperature_domain.solver_bounds_k
    if not t_lo < temperature_k < t_hi:
        raise CutTransportTopologyError(
            "dry piece temperature left the pressure-conditioned open domain; reject, do not pin"
        )
    if isinstance(config.dry.primitive_band, ct.GasOnlyPrimitiveDomain):
        local_pore = config.dry.pore if pore is None else pore
        try:
            cp.encode_gas_only_y(
                temperature_k,
                config.dry.pressure_pa,
                y_hexane,
                local_pore,
            )
        except cp.CoupledPoreTopologyError as exc:
            # RR3: carry the underlying boundary identity through the wrap so
            # the surfaced refusal names WHICH gate fired (the erosion that
            # obscured the PART-01 blocker); the leading phrase is stable for
            # existing matchers.
            raise CutTransportTopologyError(
                "dry piece composition left the exact gas-only topology; "
                f"reject, do not pin or clip [{exc}]"
            ) from exc
        return
    y_lo, y_hi = config.dry.primitive_band.y_hexane_bounds
    if not y_lo < y_hexane < y_hi:
        raise CutTransportTopologyError(
            "dry piece primitives left the caller's open gas-only band; reject, do not clip"
        )


def _validate_dry_boundary(boundary: ct.PoreBoundary, config: CutTransportConfig) -> None:
    _validate_open_dry_band(
        boundary.temperature_k,
        boundary.y_hexane,
        config,
        pore=config.dry.pore,
    )
    cp.evaluate_equilibrium(
        boundary.temperature_k,
        config.dry.pressure_pa,
        boundary.y_hexane,
        config.dry.pore,
    )


def _wet_piece_center(cell: cg.MaterialCellPartition, front_radius_m: float) -> float:
    outer = min(cell.outer_radius_m, front_radius_m)
    if outer <= cell.inner_radius_m:
        raise CutTransportTopologyError("wet piece lost positive radial width")
    return 0.5 * (cell.inner_radius_m + outer)


def _dry_piece_center(cell: cg.MaterialCellPartition, front_radius_m: float) -> float:
    inner = max(cell.inner_radius_m, front_radius_m)
    if cell.outer_radius_m <= inner:
        raise CutTransportTopologyError("dry piece lost positive radial width")
    return 0.5 * (inner + cell.outer_radius_m)


def _upwind_state(
    flux: float,
    left: cp.EquilibriumPoreState,
    right: cp.EquilibriumPoreState,
) -> cp.EquilibriumPoreState:
    if flux > 0.0:
        return left
    if flux < 0.0:
        return right
    return left


def _require_length(name: str, values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} has length {len(values)}, expected {expected}")


__all__ = [
    "CutInterfaceState",
    "CutResidualBlocks",
    "CutTransportAssembly",
    "CutTransportConfig",
    "CutTransportLayout",
    "CutTransportLedger",
    "CutTransportState",
    "CutTransportStepError",
    "CutTransportTopologyError",
    "CutTransportUnknowns",
    "DryFaceFlux",
    "DryThermodynamicForceReconstruction",
    "InterfaceCompositionBracket",
    "InterfaceCompositionRootAudit",
    "MovingInterfaceCompositionForceAuthority",
    "assemble_backward_euler",
    "candidate_from_state",
    "clear_exact_input_caches",
    "evaluate_interface_state",
    "exact_input_cache_info",
    "initialize_cut_state",
    "layout_for_geometry",
    "pressure_only_target_config",
    "validate_pressure_only_target_config",
]
