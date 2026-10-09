r"""G1G-01 conservative water/n-hexane pore-transport foundation.

This module implements the fixed-domain building blocks selected by G1G-01:

* total water is an extensive stationary-material-shell inventory;
* dry-pore water is repartitioned locally between modified-Luikov retained
  water and the one shared binary pore gas;
* the evidence-cap active set anchors retained-water potential at its cap
  fugacity without clamping the actual pore-gas water activity;
* binary Maxwell--Stefan transport has one independent diffusive flux;
* the total molar Stefan flux is the algebraic multiplier supplied by summed
  component continuity at prescribed bed-node pressure; and
* energy is carried once as

  ``q + Nw_g*h_wg + Nh_g*h_hg + Nw_ret*h_wret``.

The module intentionally does not claim the moving n-hexane-front integration
is complete.  Its stationary-shell, component-flux, total-flux, energy-flux,
and birth-inversion APIs are the conservative boundary for that follow-on
work.  Every result therefore reports ``physically_qualifying=False`` until
the ALE/front, external active-set, mesh/time, and uncertainty gates pass.

No transport coefficient is given a default point value.  A calculation must
select a declared point from an explicit positive interval.  No activity,
composition, inventory, or phase state is clipped.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, field, replace
from enum import Enum
from functools import lru_cache
from typing import Sequence

import numpy as np

from dtdc_simulator.core2.exact_cache import float_bits_key
from dtdc_simulator.core2.particle import dry_thermo, wet_core
from dtdc_simulator.core2.particle.grid import SphericalGrid
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa


WATER_DIRECT_T_MIN_K = 323.15
WATER_DIRECT_T_MAX_K = 343.15
WATER_DIRECT_ACTIVITY_MIN = 0.109
WATER_DIRECT_ACTIVITY_MAX = 0.799
REFERENCE_FUGACITY_PA = 1.0
# Pure coordinate scaling: it keeps the cut solver's finite +/-30 numerical
# guard well separated from floating-point phase boundaries without changing
# the represented open thermodynamic domain.
_GAS_ONLY_BARRIER_COORDINATE_SCALE = 2.0


class CoupledPoreTopologyError(ValueError):
    """The requested state needs a phase/authority branch not owned here."""


class HexaneSupersaturationTopologyError(CoupledPoreTopologyError):
    """The exact dry-pore primitive lies above n-hexane saturation."""


class BirthRepartitionError(RuntimeError):
    """A conservative shell birth inversion found no admissible solution."""


class CoupledBirthStepRequiredError(BirthRepartitionError):
    """A zero-transfer birth is structurally insufficient at fixed pressure."""

    classification = "NEEDS_COUPLED_RH_STEP"


class CoincidentEndpointBranch(str, Enum):
    """A1 unilateral branches at ``X_f=0`` and ``s=R``.

    These values distinguish a converged physical branch from numerical
    recovery actions.  In particular, ``ROLLBACK_RETRY`` is never a contact
    state.
    """

    RECEDING_BIRTH = "receding_birth"
    ZERO_FLUX_CONTACT = "zero_flux_contact"
    ROLLBACK_RETRY = "rollback_retry"
    INWARD_ENDPOINT_FLAG = "inward_endpoint_flag"
    MODEL_TOPOLOGY_STOP = "model_topology_stop"


class EndpointFailureKind(str, Enum):
    """Non-physical outcomes of one coincident-endpoint solve attempt."""

    NONLINEAR_NONCONVERGENCE = "nonlinear_nonconvergence"
    ITERATION_EXHAUSTION = "iteration_exhaustion"
    SINGULAR_JACOBIAN = "singular_jacobian"
    OVERLARGE_EVENT_STEP = "overlarge_event_step"
    INADMISSIBLE_TRIAL_STATE = "inadmissible_trial_state"
    COUPLED_RESIDUAL_NOT_CONVERGED = "coupled_residual_not_converged"
    MISSING_COUPLED_BIRTH = "missing_coupled_birth"
    INCOMPATIBLE_BIRTH_LEDGER = "incompatible_birth_ledger"
    CONVERGED_MODEL_TOPOLOGY = "converged_model_topology"


class ExistingDryTransportBranch(str, Enum):
    """A1 audit branches after a dry shell already exists."""

    OUTWARD_OR_ZERO = "outward_or_zero"
    CONSERVATIVE_INWARD_DRY_TRANSPORT = "conservative_inward_dry_transport"
    REVERSE_CORE_STOP = "reverse_core_stop"


class RetainedWaterActiveSet(str, Enum):
    """Named production branches of the frozen retained-water authority."""

    LUIKOV_SMOOTH = "luikov_smooth"
    RETAINED_CAP = "retained_cap"


class RetainedWaterSecantBranch(str, Enum):
    """Exact active-set path followed by one finite isothermal face secant."""

    LUIKOV_SMOOTH = "luikov_smooth"
    RETAINED_CAP = "retained_cap"
    CROSSES_RETAINED_CAP = "crosses_retained_cap"


class GasOnlyCompositionConstraint(str, Enum):
    """Physical/authority boundary selecting one gas-only composition edge."""

    HEXANE_FREE_BINARY_ENDPOINT = "hexane_free_binary_endpoint"
    WATER_SATURATION = "water_saturation"
    HEXANE_SATURATION = "hexane_saturation"
    LUIKOV_LOWER = "luikov_lower"


@dataclass(frozen=True)
class GasOnlyCompositionInterval:
    """Exact open dry-pore composition interval at one fixed ``(T, P)``.

    The reported endpoints are physical phase/authority boundaries.  Solver
    coordinates must remain strictly between them.  In particular, the
    modified-Luikov evidence cap is not an endpoint: it remains the named
    retained-water active set while the actual pore-gas activity is unchanged.
    """

    lower_y_hexane: float
    upper_y_hexane: float
    lower_constraint: GasOnlyCompositionConstraint
    upper_constraint: GasOnlyCompositionConstraint
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in (self.lower_y_hexane, self.upper_y_hexane)):
            raise ValueError("gas-only composition endpoints must be finite")
        if not 0.0 <= self.lower_y_hexane < self.upper_y_hexane <= 1.0:
            raise ValueError("gas-only composition interval must be a nonempty subset of [0,1]")
        if not isinstance(self.lower_constraint, GasOnlyCompositionConstraint) or not isinstance(
            self.upper_constraint, GasOnlyCompositionConstraint
        ):
            raise TypeError("gas-only composition endpoints need named constraints")

    @property
    def y_hexane_bounds(self) -> tuple[float, float]:
        return self.lower_y_hexane, self.upper_y_hexane


@dataclass(frozen=True)
class MobilityInterval:
    """An explicit positive uncertainty interval, never a ground-truth point.

    The mobility unit for the dimensionless-force laws in this module is
    mol/(m s).  Numerical endpoints are deliberately required from the caller;
    G1G-01 approves the interval structure but no soybean-DT endpoint values.
    """

    lower_mol_m_s: float
    upper_mol_m_s: float
    label: str

    def __post_init__(self) -> None:
        values = (self.lower_mol_m_s, self.upper_mol_m_s)
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise ValueError("mobility interval endpoints must be positive and finite")
        if self.lower_mol_m_s > self.upper_mol_m_s:
            raise ValueError("mobility interval endpoints must be ordered")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("mobility interval needs a non-empty provenance label")

    def select(self, fraction: float) -> "MobilitySelection":
        """Select a declared interval coordinate without clipping it."""
        return MobilitySelection(self, fraction)


@dataclass(frozen=True)
class MobilitySelection:
    """One explicitly identified endpoint/interior uncertainty case."""

    interval: MobilityInterval
    fraction: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.fraction) or not 0.0 <= self.fraction <= 1.0:
            raise ValueError("mobility interval fraction must lie in [0, 1]")

    @property
    def value_mol_m_s(self) -> float:
        return self.interval.lower_mol_m_s + self.fraction * (
            self.interval.upper_mol_m_s - self.interval.lower_mol_m_s
        )

    @property
    def physically_qualifying(self) -> bool:
        return False


@dataclass(frozen=True)
class CoupledPoreParams:
    """Material/EOS parameters for the locally equilibrated dry-pore branch."""

    eps_g: float = 0.141
    rho_dm_p: float = 1159.65
    cp_dry_meal: float = sp.CP_NATIVE_DRY_MEAL
    cp_oil: float = sp.CP_OIL
    T_ref_solid: float = 298.15
    w_o: float = sp.W_O_REF
    gab: sp.GabParams = sp.GabParams()
    oil: sp.OilIsotherm = sp.OilIsotherm()
    luikov: sp.LuikovParams = sp.LuikovParams()
    k_wh: float = bg.K_WH_CENTRAL
    temperature_bounds_k: tuple[float, float] = (323.15, 433.0)
    # G1G-01 nominal is zero.  Values in [0,1] are the mandatory sensitivity
    # from zero to the pure-liquid retained-water specific volume.
    retained_water_liquid_volume_fraction: float = 0.0


@dataclass(frozen=True)
class EquilibriumPoreState:
    """One complete locally equilibrated state on the selected water active set."""

    temperature_k: float
    pressure_pa: float
    gas_accessible_fraction: float
    y_water: float
    y_hexane: float
    binary_gas: bg.BinaryVirialState
    water_activity: float
    hexane_activity: float
    retained_water_active_set: RetainedWaterActiveSet
    retained_water_loading: float
    retained_hexane_loading: float
    retained_water_fugacity_pa: float
    pore_water_concentration_kg_m3: float
    pore_hexane_concentration_kg_m3: float
    retained_water_concentration_kg_m3: float
    retained_hexane_concentration_kg_m3: float
    total_water_concentration_kg_m3: float
    total_hexane_concentration_kg_m3: float
    total_water_concentration_mol_m3: float
    total_hexane_concentration_mol_m3: float
    dry_meal_energy_density_j_m3: float
    retained_water_energy_density_j_m3: float
    retained_water_pv_sensitivity_j_m3: float
    retained_hexane_energy_density_j_m3: float
    pore_gas_energy_density_j_m3: float
    energy_density_j_m3: float
    water_gas_partial_enthalpy_j_mol: float
    hexane_gas_partial_enthalpy_j_mol: float
    retained_water_enthalpy_j_mol: float
    gas_exchange_potential_isothermal: float
    retained_water_potential_isothermal: float
    outside_project_pressure: bool
    below_cross_direct_evidence: bool
    outside_gab_direct_temperature_evidence: bool
    outside_water_direct_temperature_evidence: bool
    outside_water_direct_activity_evidence: bool
    uses_nominal_zero_retained_partial_volume: bool
    physically_qualifying: bool = False

    @property
    def has_source_domain_extrapolation(self) -> bool:
        return (
            self.outside_project_pressure
            or self.below_cross_direct_evidence
            or self.outside_gab_direct_temperature_evidence
            or self.outside_water_direct_temperature_evidence
            or self.outside_water_direct_activity_evidence
        )


@dataclass(frozen=True)
class IsothermalPotentialCompositionDerivatives:
    r"""Exact fixed-``(T,P)`` force slopes with respect to ``y_h``.

    The gas-exchange force is ``ln(f_w/f_h)``.  The retained-water force is
    ``ln(f_w,ret/f_ref)`` on the named retained-water active set.  Its slope
    is therefore the gas-water fugacity slope on the smooth Luikov branch and
    exactly zero on the retained-cap branch.  The latter is a real loss of the
    composition-to-retained-force map, not a small derivative to regularize.
    """

    gas_exchange_potential_derivative_per_y_hexane: float
    retained_water_potential_derivative_per_y_hexane: float
    retained_water_active_set: RetainedWaterActiveSet
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not all(
            math.isfinite(value)
            for value in (
                self.gas_exchange_potential_derivative_per_y_hexane,
                self.retained_water_potential_derivative_per_y_hexane,
            )
        ):
            raise ValueError("isothermal potential derivatives must be finite")
        if self.gas_exchange_potential_derivative_per_y_hexane >= 0.0:
            raise ValueError("binary gas-exchange force must decrease with y_hexane")
        if not isinstance(self.retained_water_active_set, RetainedWaterActiveSet):
            raise TypeError("potential derivatives require a retained-water active set")
        if self.retained_water_active_set is RetainedWaterActiveSet.LUIKOV_SMOOTH:
            if self.retained_water_potential_derivative_per_y_hexane >= 0.0:
                raise ValueError("smooth Luikov retained-water force must decrease with y_hexane")
        elif self.retained_water_potential_derivative_per_y_hexane != 0.0:
            raise ValueError("retained-cap force derivative must be exactly zero")

    @property
    def retained_force_map_is_nonsingular(self) -> bool:
        """Whether composition locally spans the retained-water force."""

        return self.retained_water_potential_derivative_per_y_hexane != 0.0


@dataclass(frozen=True)
class IsothermalPotentialCompositionSecantSlopes:
    """Cancellation-free fixed-``(T,P)`` force secants per ``delta y_h``.

    ``retained_water_active_set`` records the left endpoint and
    ``right_retained_water_active_set`` records the right endpoint.  A finite
    secant may cross the continuous retained-cap/Luikov kink; that is distinct
    from differentiating at the kink and requires no smoothed active set.
    """

    gas_exchange_potential_secant_per_y_hexane: float
    retained_water_potential_secant_per_y_hexane: float
    retained_water_active_set: RetainedWaterActiveSet
    right_retained_water_active_set: RetainedWaterActiveSet | None = None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.right_retained_water_active_set is None:
            object.__setattr__(
                self,
                "right_retained_water_active_set",
                self.retained_water_active_set,
            )
        slopes = (
            self.gas_exchange_potential_secant_per_y_hexane,
            self.retained_water_potential_secant_per_y_hexane,
        )
        if not all(math.isfinite(value) for value in slopes):
            raise ValueError("isothermal potential secant slopes must be finite")
        if self.gas_exchange_potential_secant_per_y_hexane >= 0.0:
            raise ValueError("binary gas-exchange secant must decrease with y_hexane")
        active_sets = (
            self.retained_water_active_set,
            self.right_retained_water_active_set,
        )
        if not all(isinstance(value, RetainedWaterActiveSet) for value in active_sets):
            raise TypeError("retained-water secant endpoints need named active sets")
        if all(value is RetainedWaterActiveSet.RETAINED_CAP for value in active_sets):
            if self.retained_water_potential_secant_per_y_hexane != 0.0:
                raise ValueError("retained-cap secant must be exactly zero")
        elif self.retained_water_potential_secant_per_y_hexane >= 0.0:
            raise ValueError(
                "smooth or cap-crossing retained-water secant must decrease with y_hexane"
            )

    @property
    def retained_water_secant_branch(self) -> RetainedWaterSecantBranch:
        """Return the exact constant, smooth, or kink-crossing face branch."""

        right = self.right_retained_water_active_set
        if right is self.retained_water_active_set:
            if right is RetainedWaterActiveSet.RETAINED_CAP:
                return RetainedWaterSecantBranch.RETAINED_CAP
            return RetainedWaterSecantBranch.LUIKOV_SMOOTH
        return RetainedWaterSecantBranch.CROSSES_RETAINED_CAP

    @property
    def crosses_retained_cap(self) -> bool:
        return self.retained_water_secant_branch is RetainedWaterSecantBranch.CROSSES_RETAINED_CAP


@dataclass(frozen=True)
class IndependentMaxwellStefanFlux:
    """The one rank-one binary pore-gas diffusive flux pair."""

    water_diffusive_flux_mol_m2_s: float
    hexane_diffusive_flux_mol_m2_s: float
    dimensionless_force_gradient_m_inv: float
    mobility: MobilitySelection
    entropy_production_w_m3_k: float
    used_explicit_nonisothermal_potentials: bool
    physically_qualifying: bool = False


@dataclass(frozen=True)
class RetainedWaterFlux:
    """Finite retained-matrix water flux, separate from pore-gas MS transport."""

    flux_mol_m2_s: float
    dimensionless_force_gradient_m_inv: float
    mobility: MobilitySelection
    entropy_production_w_m3_k: float
    used_explicit_nonisothermal_potentials: bool
    physically_qualifying: bool = False


@dataclass(frozen=True)
class ComponentMolarFluxes:
    """Gas-frame and conserved component fluxes at one radial face."""

    total_stefan_flux_mol_m2_s: float
    independent_water_flux_mol_m2_s: float
    retained_water_flux_mol_m2_s: float
    gas_water_flux_mol_m2_s: float
    gas_hexane_flux_mol_m2_s: float
    conserved_water_flux_mol_m2_s: float
    conserved_hexane_flux_mol_m2_s: float
    physically_qualifying: bool = False

    @property
    def conserved_total_flux_mol_m2_s(self) -> float:
        return self.conserved_water_flux_mol_m2_s + self.conserved_hexane_flux_mol_m2_s


@dataclass(frozen=True)
class ComponentEnergyFlux:
    """One common-datum component-plus-conduction energy flux."""

    conductive_heat_flux_w_m2: float
    gas_water_enthalpy_flux_w_m2: float
    gas_hexane_enthalpy_flux_w_m2: float
    retained_water_enthalpy_flux_w_m2: float
    total_energy_flux_w_m2: float
    physically_qualifying: bool = False


@lru_cache(maxsize=8192, typed=True)
def binary_gas_component_activities(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    params: CoupledPoreParams = CoupledPoreParams(),
) -> tuple[float, float]:
    """Return water/hexane gas activities without solving solid storage.

    This is the frozen pressure-virial gas/liquid phase diagnostic only.  It
    deliberately does not invert the modified-Luikov retained-water map or
    require the composition to lie inside the dry-pore storage chart.  A
    caller that needs an admissible dry-pore state must separately apply
    :func:`gas_only_composition_interval` or solve the full equilibrium state.
    """

    _validate_gas_only_temperature_pressure(temperature_k, pressure_pa, params)
    return _gas_component_activities(
        temperature_k,
        pressure_pa,
        y_hexane,
        params,
    )


@lru_cache(maxsize=8192, typed=True)
def gas_only_composition_interval(
    temperature_k: float,
    pressure_pa: float,
    params: CoupledPoreParams = CoupledPoreParams(),
) -> GasOnlyCompositionInterval:
    r"""Return the complete open dry-pore ``y_h`` interval at fixed ``T,P``.

    The interval is the intersection

    ``0 < y_h < 1``, ``a_ref <= a_w <= 1``, and ``a_h <= 1``.

    Here ``a_ref`` is the frozen modified-Luikov lower authority.  The Luikov
    cap activity is deliberately absent because it selects an existing storage
    active set; it is neither a phase boundary nor an activity clamp.  Empty
    intersections require another phase/topology and fail closed.
    """

    _validate_gas_only_temperature_pressure(temperature_k, pressure_pa, params)
    water_at_zero, hexane_at_zero = _gas_component_activities(
        temperature_k, pressure_pa, 0.0, params
    )
    water_at_one, hexane_at_one = _gas_component_activities(temperature_k, pressure_pa, 1.0, params)
    if not water_at_zero > water_at_one or not hexane_at_zero < hexane_at_one:
        raise CoupledPoreTopologyError(
            "frozen binary activities lost their required composition monotonicity"
        )

    activity_at_ref = sp.water_activity(params.luikov.W_ref, params.luikov)
    if water_at_zero < activity_at_ref:
        raise CoupledPoreTopologyError(
            "temperature/pressure point has no gas-only composition inside the "
            "qualified Luikov lower authority"
        )

    if water_at_zero <= 1.0:
        lower = 0.0
        lower_constraint = GasOnlyCompositionConstraint.HEXANE_FREE_BINARY_ENDPOINT
    else:
        lower = _bisect_activity_composition(
            temperature_k,
            pressure_pa,
            params,
            component="water",
            target=1.0,
        )
        lower_constraint = GasOnlyCompositionConstraint.WATER_SATURATION

    if hexane_at_one <= 1.0:
        hexane_upper = 1.0
    else:
        hexane_upper = _bisect_activity_composition(
            temperature_k,
            pressure_pa,
            params,
            component="hexane",
            target=1.0,
        )
    luikov_upper = _bisect_activity_composition(
        temperature_k,
        pressure_pa,
        params,
        component="water",
        target=activity_at_ref,
    )
    if hexane_upper < luikov_upper:
        upper = hexane_upper
        upper_constraint = GasOnlyCompositionConstraint.HEXANE_SATURATION
    else:
        upper = luikov_upper
        upper_constraint = GasOnlyCompositionConstraint.LUIKOV_LOWER

    if not lower < upper:
        raise CoupledPoreTopologyError(
            "temperature/pressure point has no gas-only water/hexane composition interval"
        )
    return GasOnlyCompositionInterval(lower, upper, lower_constraint, upper_constraint)


def smooth_luikov_composition_lower_bound(
    temperature_k: float,
    pressure_pa: float,
    params: CoupledPoreParams = CoupledPoreParams(),
) -> float:
    """Return the exact retained-cap/smooth composition seam at fixed ``(T,P)``."""

    _validate_gas_only_temperature_pressure(temperature_k, pressure_pa, params)
    return _bisect_activity_composition(
        temperature_k,
        pressure_pa,
        params,
        component="water",
        target=sp.water_activity(params.luikov.W_cap, params.luikov),
    )


def encode_gas_only_y(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    params: CoupledPoreParams = CoupledPoreParams(),
) -> float:
    r"""Map one strictly admissible gas-only composition to ``(-inf, inf)``.

    A smooth log barrier uses every physical/authority gap rather than a
    piecewise minimum of the two possible upper boundaries.  The coordinate
    therefore remains smooth where hexane saturation and the Luikov lower
    authority exchange which one limits the interval.
    """

    _validate_gas_only_temperature_pressure(temperature_k, pressure_pa, params)
    if not math.isfinite(y_hexane) or not 0.0 < y_hexane < 1.0:
        raise CoupledPoreTopologyError(
            "gas-only chart encoding requires 0 < y_hexane < 1; do not clip"
        )
    water_activity, hexane_activity = _gas_component_activities(
        temperature_k, pressure_pa, y_hexane, params
    )
    return _gas_only_barrier_coordinate(
        y_hexane,
        water_activity,
        hexane_activity,
        sp.water_activity(params.luikov.W_ref, params.luikov),
    )


@lru_cache(maxsize=8192, typed=True)
def decode_gas_only_y(
    temperature_k: float,
    pressure_pa: float,
    coordinate: float,
    params: CoupledPoreParams = CoupledPoreParams(),
) -> float:
    """Invert :func:`encode_gas_only_y` without clipping a trial state."""

    if not math.isfinite(coordinate):
        raise ValueError("gas-only composition coordinate must be finite")
    interval = gas_only_composition_interval(temperature_k, pressure_pa, params)
    lower, upper = interval.y_hexane_bounds
    activity_at_ref = sp.water_activity(params.luikov.W_ref, params.luikov)
    last_composition = math.nan
    best_composition = math.nan
    best_absolute_residual = math.inf
    saw_negative_residual = False
    saw_positive_residual = False
    composition = lower + 0.5 * (upper - lower)
    for _ in range(128):
        if not lower < composition < upper or composition == last_composition:
            break
        last_composition = composition
        (
            water_activity,
            hexane_activity,
            water_activity_derivative,
            hexane_activity_derivative,
        ) = _gas_component_activities_with_derivatives(
            temperature_k, pressure_pa, composition, params
        )
        try:
            encoded = _gas_only_barrier_coordinate(
                composition,
                water_activity,
                hexane_activity,
                activity_at_ref,
            )
        except CoupledPoreTopologyError:
            # Root-rounding can put an endpoint a few ulps outside its physical
            # side.  An interior trial identifies which exact inequality was
            # crossed and contracts the bracket; no inadmissible value is ever
            # returned to the caller.
            if water_activity >= 1.0:
                lower = composition
                composition = lower + 0.5 * (upper - lower)
                continue
            if hexane_activity >= 1.0 or water_activity <= activity_at_ref:
                upper = composition
                composition = lower + 0.5 * (upper - lower)
                continue
            raise
        last_residual = encoded - coordinate
        if abs(last_residual) < best_absolute_residual:
            best_composition = composition
            best_absolute_residual = abs(last_residual)
        # The chart is consumed inside nonlinear residuals certified at
        # 1e-10 after physical row scaling.  A loose coordinate-space stop
        # can otherwise become the residual floor when the barrier map is
        # steep.  Stop at a small ulp multiple; the adjacent-float branch
        # below remains the exact fallback when no closer composition exists.
        coordinate_tolerance = 32.0 * math.ulp(max(1.0, abs(coordinate)))
        if abs(last_residual) <= coordinate_tolerance:
            return composition
        if last_residual < 0.0:
            saw_negative_residual = True
            lower = composition
        else:
            saw_positive_residual = True
            upper = composition
        barrier_derivative = _gas_only_barrier_coordinate_derivative(
            composition,
            water_activity,
            hexane_activity,
            water_activity_derivative,
            hexane_activity_derivative,
            activity_at_ref,
        )
        newton_composition = composition - last_residual / barrier_derivative
        if lower < newton_composition < upper:
            next_composition = newton_composition
        else:
            next_composition = lower + 0.5 * (upper - lower)
        if next_composition == composition:
            next_composition = lower + 0.5 * (upper - lower)
        if not lower < next_composition < upper:
            break
        composition = next_composition
    if saw_negative_residual and saw_positive_residual and math.isfinite(best_composition):
        # The target lies between two adjacent representable interior states.
        # Returning the nearer one is ordinary correctly-directed rounding,
        # not projection onto a phase/composition boundary.
        return best_composition
    raise CoupledPoreTopologyError(
        "gas-only composition coordinate has no converged representable interior value; do not clip"
    )


def evaluate_equilibrium(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    params: CoupledPoreParams = CoupledPoreParams(),
    *,
    gas_accessible_fraction: float = 1.0,
) -> EquilibriumPoreState:
    """Evaluate local equilibrium, memoized for bit-exact repeated inputs."""

    return _evaluate_equilibrium_cached(
        temperature_k,
        pressure_pa,
        y_hexane,
        params,
        gas_accessible_fraction,
        float_bits_key(
            temperature_k,
            pressure_pa,
            y_hexane,
            gas_accessible_fraction,
        ),
        hx.caloric_datum_signature(),
        wa.caloric_datum_signature(),
    )


def isothermal_potential_composition_derivatives(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    params: CoupledPoreParams = CoupledPoreParams(),
) -> IsothermalPotentialCompositionDerivatives:
    r"""Return analytic ``d/dy_h`` slopes of both frozen dry-pore forces.

    Temperature and pressure are held fixed.  Binary-fugacity slopes come
    directly from the frozen pressure-second-virial activity derivatives.
    The production retained-water active set is selected by
    :func:`evaluate_equilibrium`; no finite difference crosses or smooths its
    boundary.
    """

    state = evaluate_equilibrium(temperature_k, pressure_pa, y_hexane, params)
    (
        water_activity,
        hexane_activity,
        water_activity_derivative,
        hexane_activity_derivative,
    ) = _gas_component_activities_with_derivatives(
        temperature_k,
        pressure_pa,
        y_hexane,
        params,
    )
    gas_derivative = math.fsum(
        (
            water_activity_derivative / water_activity,
            -hexane_activity_derivative / hexane_activity,
        )
    )
    if state.retained_water_active_set is RetainedWaterActiveSet.LUIKOV_SMOOTH:
        retained_derivative = water_activity_derivative / water_activity
    elif state.retained_water_active_set is RetainedWaterActiveSet.RETAINED_CAP:
        retained_derivative = 0.0
    else:  # pragma: no cover - exhaustive guard against a future active set
        raise CoupledPoreTopologyError("unrecognized retained-water active set in force derivative")
    try:
        return IsothermalPotentialCompositionDerivatives(
            gas_exchange_potential_derivative_per_y_hexane=gas_derivative,
            retained_water_potential_derivative_per_y_hexane=retained_derivative,
            retained_water_active_set=state.retained_water_active_set,
        )
    except (TypeError, ValueError) as exc:
        raise CoupledPoreTopologyError(
            "frozen isothermal potential derivative lost its required branch structure"
        ) from exc


def isothermal_potential_composition_secant_slopes(
    temperature_k: float,
    pressure_pa: float,
    left_y_hexane: float,
    delta_y_hexane: float,
    params: CoupledPoreParams = CoupledPoreParams(),
) -> IsothermalPotentialCompositionSecantSlopes:
    r"""Return exact frozen-virial force differences factored by ``delta y_h``.

    ``log1p`` keeps both ideal-mixture terms resolved as the endpoints
    coalesce.  The pressure-second-virial contribution is integrated
    analytically, so no nearly equal endpoint gas potentials are subtracted.

    The retained potential is the continuous active-set potential frozen by
    PHY-031: constant at the evidence cap and equal to the local gas-water
    fugacity on the smooth Luikov branch.  Cap-to-cap secants are therefore
    exactly zero.  A finite cap crossing has a unique piecewise secant, formed
    from the smooth endpoint and the exact cap activity without smoothing the
    kink.  Only a zero-length derivative at the exact seam is nonunique and
    fails closed through :func:`isothermal_potential_composition_derivatives`.
    """

    if not all(
        math.isfinite(value)
        for value in (
            temperature_k,
            pressure_pa,
            left_y_hexane,
            delta_y_hexane,
        )
    ):
        raise CoupledPoreTopologyError("isothermal secant inputs must be finite")
    right_y_hexane = math.fsum((left_y_hexane, delta_y_hexane))
    if not 0.0 < left_y_hexane < 1.0 or not 0.0 < right_y_hexane < 1.0:
        raise CoupledPoreTopologyError("isothermal secant endpoints require 0 < y_hexane < 1")
    left = evaluate_equilibrium(
        temperature_k,
        pressure_pa,
        left_y_hexane,
        params,
    )
    right = evaluate_equilibrium(
        temperature_k,
        pressure_pa,
        right_y_hexane,
        params,
    )
    if delta_y_hexane == 0.0:
        seam = smooth_luikov_composition_lower_bound(
            temperature_k,
            pressure_pa,
            params,
        )
        if left_y_hexane == seam:
            raise CoupledPoreTopologyError(
                "exact retained-cap/Luikov seam has two one-sided force derivatives; "
                "select an unambiguous active-set side rather than smoothing the kink"
            )
        derivative = isothermal_potential_composition_derivatives(
            temperature_k,
            pressure_pa,
            left_y_hexane,
            params,
        )
        return IsothermalPotentialCompositionSecantSlopes(
            derivative.gas_exchange_potential_derivative_per_y_hexane,
            derivative.retained_water_potential_derivative_per_y_hexane,
            derivative.retained_water_active_set,
            derivative.retained_water_active_set,
        )
    if right_y_hexane == left_y_hexane:
        raise CoupledPoreTopologyError(
            "newborn composition increment is below floating-point representability"
        )
    mixture = bg.state(
        temperature_k,
        pressure_pa,
        1.0 - left_y_hexane,
        left_y_hexane,
        k_wh=params.k_wh,
    )
    delta_virial = mixture.B_ww - 2.0 * mixture.B_wh + mixture.B_hh
    virial_factor = 2.0 * pressure_pa * delta_virial / (bg.R * temperature_k)
    water_log_slope = math.log1p(-delta_y_hexane / (1.0 - left_y_hexane)) / delta_y_hexane
    hexane_log_slope = math.log1p(delta_y_hexane / left_y_hexane) / delta_y_hexane
    gas_slope = math.fsum((water_log_slope, -hexane_log_slope, -virial_factor))
    left_active_set = left.retained_water_active_set
    right_active_set = right.retained_water_active_set
    if (
        left_active_set is RetainedWaterActiveSet.RETAINED_CAP
        and right_active_set is RetainedWaterActiveSet.RETAINED_CAP
    ):
        retained_slope = 0.0
    elif (
        left_active_set is RetainedWaterActiveSet.LUIKOV_SMOOTH
        and right_active_set is RetainedWaterActiveSet.LUIKOV_SMOOTH
    ):
        retained_slope = math.fsum(
            (
                water_log_slope,
                -virial_factor * math.fsum((left_y_hexane, 0.5 * delta_y_hexane)),
            )
        )
    else:
        activity_at_cap = sp.water_activity(params.luikov.W_cap, params.luikov)
        if left_active_set is RetainedWaterActiveSet.RETAINED_CAP:
            smooth_activity = right.water_activity
            retained_difference = math.log1p(
                math.fsum((smooth_activity, -activity_at_cap)) / activity_at_cap
            )
        else:
            smooth_activity = left.water_activity
            retained_difference = -math.log1p(
                math.fsum((smooth_activity, -activity_at_cap)) / activity_at_cap
            )
        retained_slope = retained_difference / delta_y_hexane
    try:
        return IsothermalPotentialCompositionSecantSlopes(
            gas_slope,
            retained_slope,
            left_active_set,
            right_active_set,
        )
    except (TypeError, ValueError) as exc:
        raise CoupledPoreTopologyError(
            "frozen isothermal potential secant lost its required branch structure"
        ) from exc


@lru_cache(maxsize=8192, typed=True)
def _evaluate_equilibrium_cached(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    params: CoupledPoreParams,
    gas_accessible_fraction: float,
    _scalar_key: tuple[bytes, ...],
    _hexane_datum_key: tuple[bytes, ...],
    _water_datum_key: tuple[bytes, ...],
) -> EquilibriumPoreState:
    """Evaluate the G1G-01 local dry-pore storage and common-datum energy.

    Below the retained cap, Luikov equality is used.  At and above the cap
    activity, retained loading is exactly ``W_cap`` while the actual gas water
    activity is retained unchanged.  That is the named capacity active set,
    not an activity clamp.
    """

    _validate_params(params)
    _validate_primitives(
        temperature_k,
        pressure_pa,
        y_hexane,
        gas_accessible_fraction,
        params,
    )
    y_water = 1.0 - y_hexane
    mixture = bg.state(
        temperature_k,
        pressure_pa,
        y_water,
        y_hexane,
        k_wh=params.k_wh,
    )
    liquid_water = wa.state_Tp(temperature_k, pressure_pa, "liquid")
    liquid_hexane = hx.state_Tp(temperature_k, pressure_pa, "liquid")
    water_activity = mixture.fugacity_water_pa / liquid_water.fugacity
    hexane_activity = mixture.fugacity_hexane_pa / liquid_hexane.fugacity
    if not math.isfinite(water_activity) or not 0.0 < water_activity <= 1.0:
        raise CoupledPoreTopologyError(
            f"dry-pore water activity must lie in (0,1], got {water_activity}; "
            "activate the external/free-water topology, do not clamp"
        )

    # Amendment 18 permits only *hexane-only* supersaturation to identify the
    # semantic upper side of a representable interface-root plateau.  A state
    # that has already left the qualified Luikov domain remains a fatal
    # coupled-pore failure even when its hexane activity is also above one.
    activity_at_ref = sp.water_activity(params.luikov.W_ref, params.luikov)
    activity_at_cap = sp.water_activity(params.luikov.W_cap, params.luikov)
    lower_tolerance = 2.0e-13 * max(activity_at_ref, 1.0)
    if water_activity < activity_at_ref - lower_tolerance:
        raise CoupledPoreTopologyError(
            f"water activity {water_activity} lies below the qualified Luikov "
            f"lower-bound activity {activity_at_ref}; no continuation is authorized"
        )
    if not math.isfinite(hexane_activity) or not 0.0 < hexane_activity <= 1.0:
        raise HexaneSupersaturationTopologyError(
            f"dry-pore hexane activity must lie in (0,1], got {hexane_activity}; "
            "activate the mobile-liquid topology, do not clamp"
        )
    if water_activity < activity_at_cap:
        try:
            retained_water = sp.water_retained(water_activity, params.luikov)
        except ValueError as exc:
            raise CoupledPoreTopologyError(
                "water state left the smooth modified-Luikov branch"
            ) from exc
        water_active_set = RetainedWaterActiveSet.LUIKOV_SMOOTH
        retained_water_fugacity = mixture.fugacity_water_pa
    else:
        retained_water = params.luikov.W_cap
        water_active_set = RetainedWaterActiveSet.RETAINED_CAP
        retained_water_fugacity = activity_at_cap * liquid_water.fugacity

    try:
        retained_hexane = sp.retained_hexane(
            hexane_activity,
            temperature_k,
            params.w_o,
            params.gab,
            params.oil,
        )
    except ValueError as exc:
        raise CoupledPoreTopologyError("hexane state left the frozen retained dry branch") from exc
    if not math.isfinite(retained_hexane) or retained_hexane < 0.0:
        raise CoupledPoreTopologyError("retained hexane is negative or non-finite")

    effective_eps = params.eps_g * gas_accessible_fraction
    pore_water = effective_eps * mixture.water_density_kg_m3
    pore_hexane = effective_eps * mixture.hexane_density_kg_m3
    retained_water_concentration = params.rho_dm_p * retained_water
    retained_hexane_concentration = params.rho_dm_p * retained_hexane
    total_water = retained_water_concentration + pore_water
    total_hexane = retained_hexane_concentration + pore_hexane

    dry_meal_energy = params.rho_dm_p * sp.composite_dry_meal_sensible_energy(
        temperature_k,
        params.T_ref_solid,
        params.w_o,
        params.cp_dry_meal,
        params.cp_oil,
    )
    retained_water_zero_volume = retained_water_concentration * liquid_water.h_mass
    retained_water_pv = retained_water_concentration * pressure_pa / liquid_water.rho_mass
    retained_water_energy = retained_water_zero_volume - (
        params.retained_water_liquid_volume_fraction * retained_water_pv
    )
    retained_hexane_energy = params.rho_dm_p * sp.retained_internal_energy(
        retained_hexane,
        temperature_k,
        params.w_o,
        params.gab,
        params.oil,
        pressure_pa=pressure_pa,
        activity=hexane_activity,
    )
    ideal_gas_enthalpy = y_water * bg.water_ideal_enthalpy_molar(
        temperature_k
    ) + y_hexane * hx.h_ideal(temperature_k)
    gas_internal_energy = (
        ideal_gas_enthalpy
        + mixture.residual_enthalpy_molar
        - pressure_pa / mixture.molar_density_mol_m3
    )
    pore_gas_energy = effective_eps * mixture.molar_density_mol_m3 * gas_internal_energy
    energy = dry_meal_energy + retained_water_energy + retained_hexane_energy + pore_gas_energy

    h_w_g = (
        bg.water_ideal_enthalpy_molar(temperature_k) + mixture.partial_residual_enthalpy_water_molar
    )
    h_h_g = hx.h_ideal(temperature_k) + mixture.partial_residual_enthalpy_hexane_molar
    h_w_ret = liquid_water.h_mass * wa.M
    values = (
        water_activity,
        hexane_activity,
        retained_water,
        retained_hexane,
        retained_water_fugacity,
        pore_water,
        pore_hexane,
        total_water,
        total_hexane,
        dry_meal_energy,
        retained_water_energy,
        retained_hexane_energy,
        pore_gas_energy,
        energy,
        h_w_g,
        h_h_g,
        h_w_ret,
    )
    if not all(math.isfinite(value) for value in values):
        raise CoupledPoreTopologyError("coupled-pore equilibrium produced a non-finite state")

    return EquilibriumPoreState(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        gas_accessible_fraction=gas_accessible_fraction,
        y_water=y_water,
        y_hexane=y_hexane,
        binary_gas=mixture,
        water_activity=water_activity,
        hexane_activity=hexane_activity,
        retained_water_active_set=water_active_set,
        retained_water_loading=retained_water,
        retained_hexane_loading=retained_hexane,
        retained_water_fugacity_pa=retained_water_fugacity,
        pore_water_concentration_kg_m3=pore_water,
        pore_hexane_concentration_kg_m3=pore_hexane,
        retained_water_concentration_kg_m3=retained_water_concentration,
        retained_hexane_concentration_kg_m3=retained_hexane_concentration,
        total_water_concentration_kg_m3=total_water,
        total_hexane_concentration_kg_m3=total_hexane,
        total_water_concentration_mol_m3=total_water / wa.M,
        total_hexane_concentration_mol_m3=total_hexane / hx.M,
        dry_meal_energy_density_j_m3=dry_meal_energy,
        retained_water_energy_density_j_m3=retained_water_energy,
        retained_water_pv_sensitivity_j_m3=retained_water_pv,
        retained_hexane_energy_density_j_m3=retained_hexane_energy,
        pore_gas_energy_density_j_m3=pore_gas_energy,
        energy_density_j_m3=energy,
        water_gas_partial_enthalpy_j_mol=h_w_g,
        hexane_gas_partial_enthalpy_j_mol=h_h_g,
        retained_water_enthalpy_j_mol=h_w_ret,
        gas_exchange_potential_isothermal=math.log(
            mixture.fugacity_water_pa / mixture.fugacity_hexane_pa
        ),
        retained_water_potential_isothermal=math.log(
            retained_water_fugacity / REFERENCE_FUGACITY_PA
        ),
        outside_project_pressure=mixture.outside_project_pressure,
        below_cross_direct_evidence=mixture.below_cross_direct_evidence,
        outside_gab_direct_temperature_evidence=not (
            dry_thermo.GAB_DIRECT_T_MIN_K <= temperature_k <= dry_thermo.GAB_DIRECT_T_MAX_K
        ),
        outside_water_direct_temperature_evidence=not (
            WATER_DIRECT_T_MIN_K <= temperature_k <= WATER_DIRECT_T_MAX_K
        ),
        outside_water_direct_activity_evidence=not (
            WATER_DIRECT_ACTIVITY_MIN <= water_activity <= WATER_DIRECT_ACTIVITY_MAX
        ),
        uses_nominal_zero_retained_partial_volume=(
            params.retained_water_liquid_volume_fraction == 0.0
        ),
    )


def clear_equilibrium_cache() -> None:
    """Drop every cached local-equilibrium state and its configuration keys."""

    _evaluate_equilibrium_cached.cache_clear()


def equilibrium_cache_info():
    """Return bounded-cache hit/miss/size diagnostics."""

    return _evaluate_equilibrium_cached.cache_info()


def ideal_binary_mobility_mol_m_s(
    molar_density_mol_m3: float,
    diffusivity_m2_s: float,
    y_water: float,
    y_hexane: float,
) -> float:
    """Return ``c D y_w y_h`` for the ideal binary Maxwell--Stefan oracle."""
    values = (molar_density_mol_m3, diffusivity_m2_s, y_water, y_hexane)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("ideal binary mobility inputs must be finite")
    if molar_density_mol_m3 <= 0.0 or diffusivity_m2_s <= 0.0:
        raise ValueError("ideal binary density and diffusivity must be positive")
    if y_water <= 0.0 or y_hexane <= 0.0 or abs(y_water + y_hexane - 1.0) > 1.0e-12:
        raise ValueError("ideal binary mobility requires positive mole fractions summing to one")
    return molar_density_mol_m3 * diffusivity_m2_s * y_water * y_hexane


def independent_maxwell_stefan_flux(
    left: EquilibriumPoreState,
    right: EquilibriumPoreState,
    distance_m: float,
    mobility: MobilitySelection,
    *,
    left_chemical_potential_difference_over_rt: float | None = None,
    right_chemical_potential_difference_over_rt: float | None = None,
) -> IndependentMaxwellStefanFlux:
    r"""Return the one independent binary diffusive flux ``J_w=-J_h``.

    At an isothermal face, ``ln(f_w/f_h)`` is the exactly equivalent fugacity
    force.  A non-isothermal caller must provide both complete
    ``(mu_w-mu_h)/(R T)`` values; this foundation refuses to invent a Soret or
    standard-chemical-potential term.
    """
    _validate_distance(distance_m)
    explicit = (
        left_chemical_potential_difference_over_rt is not None
        or right_chemical_potential_difference_over_rt is not None
    )
    if explicit:
        if (
            left_chemical_potential_difference_over_rt is None
            or right_chemical_potential_difference_over_rt is None
        ):
            raise ValueError("both non-isothermal chemical potentials are required")
        potentials = (
            left_chemical_potential_difference_over_rt,
            right_chemical_potential_difference_over_rt,
        )
        if not all(math.isfinite(value) for value in potentials):
            raise ValueError("chemical-potential forces must be finite")
        phi_left, phi_right = potentials
    else:
        if left.temperature_k != right.temperature_k:
            raise CoupledPoreTopologyError(
                "non-isothermal Maxwell-Stefan flux needs explicit complete "
                "chemical potentials; a log-fugacity-only Soret law is not authorized"
            )
        phi_left = left.gas_exchange_potential_isothermal
        phi_right = right.gas_exchange_potential_isothermal
    gradient = (phi_right - phi_left) / distance_m
    water_flux = -mobility.value_mol_m_s * gradient
    entropy_production = bg.R * mobility.value_mol_m_s * gradient * gradient
    return IndependentMaxwellStefanFlux(
        water_diffusive_flux_mol_m2_s=water_flux,
        hexane_diffusive_flux_mol_m2_s=-water_flux,
        dimensionless_force_gradient_m_inv=gradient,
        mobility=mobility,
        entropy_production_w_m3_k=entropy_production,
        used_explicit_nonisothermal_potentials=explicit,
    )


def retained_water_flux(
    left: EquilibriumPoreState,
    right: EquilibriumPoreState,
    distance_m: float,
    mobility: MobilitySelection,
    *,
    left_retained_chemical_potential_over_rt: float | None = None,
    right_retained_chemical_potential_over_rt: float | None = None,
) -> RetainedWaterFlux:
    r"""Return ``N_w,ret=-L_w,ret grad(mu_w,ret/(R T))``.

    As for the pore-gas force, log retained fugacity is sufficient only on an
    isothermal face.  The retained-cap potential uses the cap fugacity while
    leaving the actual gas activity untouched.  Thus a cap-to-cap face has
    exactly zero retained flux even when its independent binary-gas flux is
    nonzero; this is the PHY-031 capacity active set, not an activity clamp.
    """
    _validate_distance(distance_m)
    explicit = (
        left_retained_chemical_potential_over_rt is not None
        or right_retained_chemical_potential_over_rt is not None
    )
    if explicit:
        if (
            left_retained_chemical_potential_over_rt is None
            or right_retained_chemical_potential_over_rt is None
        ):
            raise ValueError("both retained-water chemical potentials are required")
        potentials = (
            left_retained_chemical_potential_over_rt,
            right_retained_chemical_potential_over_rt,
        )
        if not all(math.isfinite(value) for value in potentials):
            raise ValueError("retained-water chemical-potential forces must be finite")
        phi_left, phi_right = potentials
    else:
        if left.temperature_k != right.temperature_k:
            raise CoupledPoreTopologyError(
                "non-isothermal retained-water flux needs explicit complete "
                "chemical potentials; no Soret law is authorized"
            )
        phi_left = left.retained_water_potential_isothermal
        phi_right = right.retained_water_potential_isothermal
    gradient = (phi_right - phi_left) / distance_m
    flux = -mobility.value_mol_m_s * gradient
    entropy_production = bg.R * mobility.value_mol_m_s * gradient * gradient
    return RetainedWaterFlux(
        flux_mol_m2_s=flux,
        dimensionless_force_gradient_m_inv=gradient,
        mobility=mobility,
        entropy_production_w_m3_k=entropy_production,
        used_explicit_nonisothermal_potentials=explicit,
    )


def compose_component_fluxes(
    y_water: float,
    y_hexane: float,
    total_stefan_flux_mol_m2_s: float,
    independent_water_flux_mol_m2_s: float,
    retained_water_flux_mol_m2_s: float,
) -> ComponentMolarFluxes:
    """Compose gas-frame and conserved fluxes without adding another rate law."""
    values = (
        y_water,
        y_hexane,
        total_stefan_flux_mol_m2_s,
        independent_water_flux_mol_m2_s,
        retained_water_flux_mol_m2_s,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("component-flux inputs must be finite")
    if y_water < 0.0 or y_hexane < 0.0 or abs(y_water + y_hexane - 1.0) > 1.0e-12:
        raise ValueError("face mole fractions must be non-negative and sum to one")
    # FLUX CONVENTION (annotation only; no law, gate, tolerance or value changes here).
    # ``total_stefan_flux_mol_m2_s`` is the TOTAL pore-gas mixture molar flux
    # ``N_t = N_water + N_hexane`` at the face -- in the birth chart it is the unknown
    # ``dry_total_stefan_fluxes_mol_m2_s[i]``, with ``i = 0`` the receding front and
    # ``i = 1`` the outer particle surface -- positive OUTWARD along ``+r`` in the
    # solid-matrix-fixed (laboratory) frame.  It is NOT the hexane rate: adding the two
    # lines below and using ``y_water + y_hexane = 1`` gives
    # ``gas_water + gas_hexane = total_stefan_flux`` identically, so the unknown is the
    # SUM of the two species gas fluxes (advective drag ``y_k * N_t`` plus the
    # equal-and-opposite binary diffusion ``+/- J_w``).  The single independent
    # diffusive quantity is the rank-one binary pair
    # ``independent_water_flux_mol_m2_s = J_w = -J_h``, in the molar-average (mixture)
    # frame and oriented by ``independent_maxwell_stefan_flux`` above
    # (``water_flux = -mobility * (phi_right - phi_left) / distance_m`` with ``right``
    # the outer state, so positive is transport toward increasing radius).  A NEGATIVE
    # total at an accepted birth root means net INWARD water uptake by the newborn dry
    # pore on the ``RETAINED_CAP`` retained-water active set, not an inward hexane flux;
    # measured and explained in
    # docs/GT_PS2_Q_B1_4A_FLUX_SIGN_DIAGNOSIS_2026-09-04.md (item R-Q4A-1).
    gas_water = y_water * total_stefan_flux_mol_m2_s + independent_water_flux_mol_m2_s
    gas_hexane = y_hexane * total_stefan_flux_mol_m2_s - independent_water_flux_mol_m2_s
    return ComponentMolarFluxes(
        total_stefan_flux_mol_m2_s=total_stefan_flux_mol_m2_s,
        independent_water_flux_mol_m2_s=independent_water_flux_mol_m2_s,
        retained_water_flux_mol_m2_s=retained_water_flux_mol_m2_s,
        gas_water_flux_mol_m2_s=gas_water,
        gas_hexane_flux_mol_m2_s=gas_hexane,
        conserved_water_flux_mol_m2_s=gas_water + retained_water_flux_mol_m2_s,
        conserved_hexane_flux_mol_m2_s=gas_hexane,
    )


def component_energy_flux(
    component_fluxes: ComponentMolarFluxes,
    conductive_heat_flux_w_m2: float,
    water_gas_partial_enthalpy_j_mol: float,
    hexane_gas_partial_enthalpy_j_mol: float,
    retained_water_enthalpy_j_mol: float,
) -> ComponentEnergyFlux:
    """Assemble the G1G-01 energy flux with every component carried once."""
    values = (
        conductive_heat_flux_w_m2,
        water_gas_partial_enthalpy_j_mol,
        hexane_gas_partial_enthalpy_j_mol,
        retained_water_enthalpy_j_mol,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("component energy-flux inputs must be finite")
    water_gas = component_fluxes.gas_water_flux_mol_m2_s * water_gas_partial_enthalpy_j_mol
    hexane_gas = component_fluxes.gas_hexane_flux_mol_m2_s * hexane_gas_partial_enthalpy_j_mol
    retained_water = component_fluxes.retained_water_flux_mol_m2_s * retained_water_enthalpy_j_mol
    total = math.fsum((conductive_heat_flux_w_m2, water_gas, hexane_gas, retained_water))
    return ComponentEnergyFlux(
        conductive_heat_flux_w_m2=conductive_heat_flux_w_m2,
        gas_water_enthalpy_flux_w_m2=water_gas,
        gas_hexane_enthalpy_flux_w_m2=hexane_gas,
        retained_water_enthalpy_flux_w_m2=retained_water,
        total_energy_flux_w_m2=total,
    )


@dataclass(frozen=True)
class TotalStefanFluxRecovery:
    """Algebraic fixed-pressure total flux recovered inside a DAE residual."""

    combined_component_fluxes_mol_m2_s: tuple[float, ...]
    total_stefan_fluxes_mol_m2_s: tuple[float, ...]
    retained_water_fluxes_mol_m2_s: tuple[float, ...]
    combined_component_rates_mol_s: tuple[float, ...]
    summed_storage_changes_mol_s: tuple[float, ...]
    telescope_residual_mol_s: float
    normalized_telescope_residual: float
    physically_qualifying: bool = False


@dataclass(frozen=True)
class ComponentBalanceResiduals:
    """Per-shell component residuals after composing the recovered total flux."""

    water_residuals_mol_s: tuple[float, ...]
    hexane_residuals_mol_s: tuple[float, ...]
    summed_residuals_mol_s: tuple[float, ...]
    max_normalized_water_residual: float
    max_normalized_hexane_residual: float
    max_normalized_summed_residual: float
    physically_qualifying: bool = False


def recover_total_stefan_flux(
    grid: SphericalGrid,
    old_water_inventory_mol: Sequence[float],
    new_water_inventory_mol: Sequence[float],
    old_hexane_inventory_mol: Sequence[float],
    new_hexane_inventory_mol: Sequence[float],
    retained_water_face_fluxes_mol_m2_s: Sequence[float],
    dt_s: float,
    *,
    inner_combined_component_rate_mol_s: float = 0.0,
) -> TotalStefanFluxRecovery:
    r"""Recover ``N_t`` cumulatively from summed discrete continuity.

    For each stationary spherical cell,

    ``Delta(n_w+n_h)/dt + A_R G_R - A_L G_L = 0``,

    where ``G=N_t+N_w,ret``.  The recurrence is evaluated from the known inner
    integration constant.  ``N_t=G-N_w,ret`` is therefore a DAE/Lagrange face
    flux of the current candidate inventories, not a post-hoc balance repair.

    A partial-front/ALE integrator must supply its exact swept-storage changes
    in the inventory arguments and the conservative inner/front integration
    constant.  This function does not invent that front condition.
    """
    _validate_dt(dt_s)
    sequences = (
        old_water_inventory_mol,
        new_water_inventory_mol,
        old_hexane_inventory_mol,
        new_hexane_inventory_mol,
    )
    for values in sequences:
        _validate_cell_inventory(values, grid.n)
    _validate_face_values(
        retained_water_face_fluxes_mol_m2_s,
        grid.n + 1,
        "retained-water face fluxes",
    )
    if not math.isfinite(inner_combined_component_rate_mol_s):
        raise ValueError("inner combined component rate must be finite")
    if grid.areas[0] == 0.0:
        if inner_combined_component_rate_mol_s != 0.0:
            raise ValueError("a spherical center has exactly zero inner component rate")
        if retained_water_face_fluxes_mol_m2_s[0] != 0.0:
            raise ValueError("a spherical center has exactly zero retained-water flux")

    storage_changes = tuple(
        math.fsum((new_w, -old_w, new_h, -old_h)) / dt_s
        for old_w, new_w, old_h, new_h in zip(
            old_water_inventory_mol,
            new_water_inventory_mol,
            old_hexane_inventory_mol,
            new_hexane_inventory_mol,
        )
    )
    rates = [inner_combined_component_rate_mol_s]
    for storage_change in storage_changes:
        rates.append(math.fsum((rates[-1], -storage_change)))

    combined_fluxes: list[float] = []
    stefan_fluxes: list[float] = []
    for face, (area, rate, retained_flux) in enumerate(
        zip(grid.areas, rates, retained_water_face_fluxes_mol_m2_s)
    ):
        if area == 0.0:
            combined_flux = 0.0
            stefan_flux = 0.0
        else:
            combined_flux = rate / area
            stefan_flux = combined_flux - retained_flux
        if not math.isfinite(stefan_flux):
            raise RuntimeError(f"recovered total Stefan flux is non-finite at face {face}")
        combined_fluxes.append(combined_flux)
        stefan_fluxes.append(stefan_flux)

    telescope_residual = math.fsum((rates[-1], -rates[0], *storage_changes))
    scale = max(
        abs(rates[-1]),
        abs(rates[0]),
        *(abs(value) for value in storage_changes),
        1.0e-300,
    )
    return TotalStefanFluxRecovery(
        combined_component_fluxes_mol_m2_s=tuple(combined_fluxes),
        total_stefan_fluxes_mol_m2_s=tuple(stefan_fluxes),
        retained_water_fluxes_mol_m2_s=tuple(retained_water_face_fluxes_mol_m2_s),
        combined_component_rates_mol_s=tuple(rates),
        summed_storage_changes_mol_s=storage_changes,
        telescope_residual_mol_s=telescope_residual,
        normalized_telescope_residual=abs(telescope_residual) / scale,
    )


def component_balance_residuals(
    grid: SphericalGrid,
    old_water_inventory_mol: Sequence[float],
    new_water_inventory_mol: Sequence[float],
    old_hexane_inventory_mol: Sequence[float],
    new_hexane_inventory_mol: Sequence[float],
    face_fluxes: Sequence[ComponentMolarFluxes],
    dt_s: float,
) -> ComponentBalanceResiduals:
    """Return independent water/hexane FV residuals on a stationary grid."""
    _validate_dt(dt_s)
    sequences = (
        old_water_inventory_mol,
        new_water_inventory_mol,
        old_hexane_inventory_mol,
        new_hexane_inventory_mol,
    )
    for values in sequences:
        _validate_cell_inventory(values, grid.n)
    if len(face_fluxes) != grid.n + 1:
        raise ValueError("component face fluxes must align with all grid faces")

    water_rates = tuple(
        area * flux.conserved_water_flux_mol_m2_s for area, flux in zip(grid.areas, face_fluxes)
    )
    hexane_rates = tuple(
        area * flux.conserved_hexane_flux_mol_m2_s for area, flux in zip(grid.areas, face_fluxes)
    )
    water_residuals: list[float] = []
    hexane_residuals: list[float] = []
    water_scales: list[float] = []
    hexane_scales: list[float] = []
    for cell in range(grid.n):
        water_change = (
            math.fsum((new_water_inventory_mol[cell], -old_water_inventory_mol[cell])) / dt_s
        )
        hexane_change = (
            math.fsum((new_hexane_inventory_mol[cell], -old_hexane_inventory_mol[cell])) / dt_s
        )
        water_residual = math.fsum((water_change, water_rates[cell + 1], -water_rates[cell]))
        hexane_residual = math.fsum((hexane_change, hexane_rates[cell + 1], -hexane_rates[cell]))
        water_residuals.append(water_residual)
        hexane_residuals.append(hexane_residual)
        water_scales.append(
            max(
                abs(water_change),
                abs(water_rates[cell]),
                abs(water_rates[cell + 1]),
                1.0e-300,
            )
        )
        hexane_scales.append(
            max(
                abs(hexane_change),
                abs(hexane_rates[cell]),
                abs(hexane_rates[cell + 1]),
                1.0e-300,
            )
        )
    summed = tuple(
        math.fsum((water, hexane)) for water, hexane in zip(water_residuals, hexane_residuals)
    )
    summed_scales = tuple(
        max(water_scale, hexane_scale)
        for water_scale, hexane_scale in zip(water_scales, hexane_scales)
    )
    return ComponentBalanceResiduals(
        water_residuals_mol_s=tuple(water_residuals),
        hexane_residuals_mol_s=tuple(hexane_residuals),
        summed_residuals_mol_s=summed,
        max_normalized_water_residual=max(
            abs(value) / scale for value, scale in zip(water_residuals, water_scales)
        ),
        max_normalized_hexane_residual=max(
            abs(value) / scale for value, scale in zip(hexane_residuals, hexane_scales)
        ),
        max_normalized_summed_residual=max(
            abs(value) / scale for value, scale in zip(summed, summed_scales)
        ),
    )


def conservative_component_update(
    grid: SphericalGrid,
    old_water_inventory_mol: Sequence[float],
    old_hexane_inventory_mol: Sequence[float],
    face_fluxes: Sequence[ComponentMolarFluxes],
    dt_s: float,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Apply one exact flux-form inventory update for manufactured checks."""
    _validate_dt(dt_s)
    _validate_cell_inventory(old_water_inventory_mol, grid.n)
    _validate_cell_inventory(old_hexane_inventory_mol, grid.n)
    if len(face_fluxes) != grid.n + 1:
        raise ValueError("component face fluxes must align with all grid faces")
    water_rates = tuple(
        area * flux.conserved_water_flux_mol_m2_s for area, flux in zip(grid.areas, face_fluxes)
    )
    hexane_rates = tuple(
        area * flux.conserved_hexane_flux_mol_m2_s for area, flux in zip(grid.areas, face_fluxes)
    )
    new_water = tuple(
        old - dt_s * (water_rates[cell + 1] - water_rates[cell])
        for cell, old in enumerate(old_water_inventory_mol)
    )
    new_hexane = tuple(
        old - dt_s * (hexane_rates[cell + 1] - hexane_rates[cell])
        for cell, old in enumerate(old_hexane_inventory_mol)
    )
    _validate_cell_inventory(new_water, grid.n)
    _validate_cell_inventory(new_hexane, grid.n)
    return new_water, new_hexane


@dataclass(frozen=True)
class StationaryMaterialInventory:
    """Extensive component/energy states on fixed material shells.

    Water owns this stationary topology.  N-hexane is retained here as the
    coupling ledger required at birth; its subsequent wet/dry evolution is
    owned by the moving-front solver.
    """

    grid: SphericalGrid
    water_masses_kg: tuple[float, ...]
    hexane_masses_kg: tuple[float, ...]
    energies_j: tuple[float, ...]
    physically_qualifying: bool = False

    def __post_init__(self) -> None:
        if len(self.water_masses_kg) != self.grid.n:
            raise ValueError("stationary water masses must align with material shells")
        if len(self.hexane_masses_kg) != self.grid.n:
            raise ValueError("stationary hexane masses must align with material shells")
        if len(self.energies_j) != self.grid.n:
            raise ValueError("stationary energies must align with material shells")
        if not all(math.isfinite(value) and value >= 0.0 for value in self.water_masses_kg):
            raise ValueError("stationary shell water masses must be finite and non-negative")
        if not all(math.isfinite(value) and value >= 0.0 for value in self.hexane_masses_kg):
            raise ValueError("stationary shell hexane masses must be finite and non-negative")
        if not all(math.isfinite(value) for value in self.energies_j):
            raise ValueError("stationary shell energies must be finite")

    @property
    def total_water_mass_kg(self) -> float:
        return math.fsum(self.water_masses_kg)

    @property
    def total_energy_j(self) -> float:
        return math.fsum(self.energies_j)

    @property
    def total_hexane_mass_kg(self) -> float:
        return math.fsum(self.hexane_masses_kg)


def inventory_from_activation(
    activation: object,
    wet_params: wet_core.WetCoreParams,
) -> StationaryMaterialInventory:
    """Lift an exact radial activation payload onto stationary material shells.

    The function intentionally accepts the activation protocol rather than
    importing ``sphere`` and creating a module cycle.  Every required field is
    validated and the reconstructed extensive totals must match the activation
    payload before the inventory is returned.
    """
    required = (
        "grid",
        "temperatures_k",
        "oil_fraction_labels",
        "retained_water_loadings",
        "wet_inventory",
        "retained_water_mass_kg",
        "internal_hexane_mass_kg",
        "total_internal_energy_j",
    )
    if any(not hasattr(activation, name) for name in required):
        raise TypeError("activation payload lacks stationary-shell material fields")
    material_grid = activation.grid
    if not isinstance(material_grid, SphericalGrid):
        raise TypeError("activation grid is not a spherical material grid")
    temperatures = tuple(activation.temperatures_k)
    oil_labels = tuple(activation.oil_fraction_labels)
    water_loadings = tuple(activation.retained_water_loadings)
    hexane_loadings = tuple(activation.wet_inventory.loadings)
    if any(
        len(values) != material_grid.n
        for values in (temperatures, oil_labels, water_loadings, hexane_loadings)
    ):
        raise ValueError("activation material fields do not align with its grid")
    if any(label != oil_labels[0] for label in oil_labels):
        raise ValueError("activation oil labels are not one stationary material basis")
    if any(label != water_loadings[0] for label in water_loadings):
        raise ValueError("activation water labels are not one stationary material basis")
    if not math.isclose(oil_labels[0], wet_params.w_o, rel_tol=2.0e-14):
        raise ValueError("activation and stationary inventory oil bases differ")
    if not math.isclose(water_loadings[0], wet_params.X_water, rel_tol=2.0e-14, abs_tol=1.0e-16):
        raise ValueError("activation and stationary inventory water bases differ")
    material_wet_params = replace(
        wet_params,
        w_o=oil_labels[0],
        X_water=water_loadings[0],
    )
    if material_grid.R <= 0.0:
        raise ValueError("activation material radius must be positive")
    water_masses = tuple(
        volume * material_wet_params.rho_dm_p * loading
        for volume, loading in zip(material_grid.volumes, water_loadings)
    )
    hexane_masses = tuple(
        volume * material_wet_params.rho_dm_p * loading
        for volume, loading in zip(material_grid.volumes, hexane_loadings)
    )
    energies = tuple(
        volume * wet_core.energy_density(temperature, hexane_loading, material_wet_params)
        for volume, temperature, hexane_loading in zip(
            material_grid.volumes,
            temperatures,
            hexane_loadings,
        )
    )
    water_residual = math.fsum((*water_masses, -float(activation.retained_water_mass_kg)))
    energy_residual = math.fsum((*energies, -float(activation.total_internal_energy_j)))
    hexane_residual = math.fsum((*hexane_masses, -float(activation.internal_hexane_mass_kg)))
    water_scale = max(
        abs(float(activation.retained_water_mass_kg)),
        abs(math.fsum(water_masses)),
        1.0e-300,
    )
    energy_scale = max(
        abs(float(activation.total_internal_energy_j)),
        abs(math.fsum(energies)),
        1.0e-300,
    )
    hexane_scale = max(
        abs(float(activation.internal_hexane_mass_kg)),
        abs(math.fsum(hexane_masses)),
        1.0e-300,
    )
    if abs(water_residual) / water_scale > 2.0e-12:
        raise ValueError("activation water payload fails stationary-shell reconstruction")
    if abs(energy_residual) / energy_scale > 2.0e-12:
        raise ValueError("activation energy payload fails stationary-shell reconstruction")
    if abs(hexane_residual) / hexane_scale > 2.0e-12:
        raise ValueError("activation hexane payload fails stationary-shell reconstruction")
    return StationaryMaterialInventory(
        material_grid,
        water_masses,
        hexane_masses,
        energies,
    )


@dataclass(frozen=True)
class BirthCouplingIncrements:
    """Signed conservative increments from the coupled front/DAE step.

    Positive values enter the shell.  ``energy_increment_j`` is the single
    common-datum integral of conduction, both gas-component enthalpy fluxes,
    retained-water enthalpy flux, and exact swept/ALE terms.  These increments
    are not fitted background sources.
    """

    water_mass_increment_kg: float
    hexane_mass_increment_kg: float
    energy_increment_j: float
    provenance: str

    def __post_init__(self) -> None:
        values = (
            self.water_mass_increment_kg,
            self.hexane_mass_increment_kg,
            self.energy_increment_j,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("birth coupling increments must be finite")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("birth coupling increments need explicit provenance")

    @property
    def is_zero_transfer(self) -> bool:
        return (
            self.water_mass_increment_kg == 0.0
            and self.hexane_mass_increment_kg == 0.0
            and self.energy_increment_j == 0.0
        )


@dataclass(frozen=True)
class BirthRepartitionLedger:
    """Independent extensive audit of one zero-pore-volume shell birth."""

    shell_index: int
    dry_pore_volume_before_m3: float
    dry_pore_volume_after_m3: float
    water_mass_before_kg: float
    hexane_mass_before_kg: float
    retained_water_mass_after_kg: float
    pore_water_mass_after_kg: float
    water_mass_after_kg: float
    hexane_mass_after_kg: float
    energy_before_j: float
    energy_after_j: float
    coupled_water_increment_kg: float
    coupled_hexane_increment_kg: float
    coupled_energy_increment_j: float
    artificial_water_source_kg: float
    artificial_energy_source_j: float
    water_residual_kg: float
    hexane_residual_kg: float
    energy_residual_j: float
    normalized_water_residual: float
    normalized_hexane_residual: float
    normalized_energy_residual: float
    max_normalized_residual: float
    physically_qualifying: bool = False


@dataclass(frozen=True)
class BirthRepartitionResult:
    """One dry-pore equilibrium state recovered from existing shell totals."""

    stationary_inventory: StationaryMaterialInventory
    coupling: BirthCouplingIncrements
    shell_index: int
    state: EquilibriumPoreState
    resulting_hexane_mass_kg: float
    iterations: int
    max_scaled_inversion_residual: float
    ledger: BirthRepartitionLedger
    physically_qualifying: bool = False


@dataclass(frozen=True)
class CoincidentEndpointResiduals:
    """Normalized residual evidence from the complete endpoint problem.

    A contact decision needs the same component, energy, equilibrium, and
    external-active-set evidence as a receding decision.  A small signed
    surface flux is deliberately not included here and is never converted to
    zero by a tolerance.
    """

    water: float
    hexane: float
    common_datum_energy: float
    local_equilibrium: float
    external_active_set: float
    property_bounds_admissible: bool
    phase_topology_admissible: bool
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.water,
            self.hexane,
            self.common_datum_energy,
            self.local_equilibrium,
            self.external_active_set,
        )
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            raise ValueError("endpoint residuals must be finite and non-negative")
        if type(self.property_bounds_admissible) is not bool:
            raise TypeError("property-bounds admissibility must be boolean")
        if type(self.phase_topology_admissible) is not bool:
            raise TypeError("phase-topology admissibility must be boolean")

    @property
    def max_normalized_residual(self) -> float:
        return max(
            self.water,
            self.hexane,
            self.common_datum_energy,
            self.local_equilibrium,
            self.external_active_set,
        )


@dataclass(frozen=True)
class CoincidentBirthALELedger:
    """Exact extensive n-hexane ledger for one finite newborn dry shell.

    All transfers are integrated over the accepted event substep and use an
    outward-positive laboratory-frame convention.  The dry-domain balance is

    ``M_after - M_before = interface_transfer + swept_dry_storage - surface_out``.

    The time-integrated Rankine--Hugoniot identity is
    ``interface_transfer = swept_wet_inventory - swept_dry_storage``.

    This deliberately keeps finite-event surface loss separate from the
    instantaneous coincident-interface flux used by the onset RH condition.
    """

    dry_hexane_inventory_before_mol: float
    dry_hexane_inventory_after_mol: float
    interface_hexane_transfer_to_dry_mol: float
    swept_wet_hexane_inventory_mol: float
    swept_dry_hexane_storage_mol: float
    surface_hexane_out_mol: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.dry_hexane_inventory_before_mol,
            self.dry_hexane_inventory_after_mol,
            self.swept_wet_hexane_inventory_mol,
            self.swept_dry_hexane_storage_mol,
        )
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            raise ValueError("birth ALE inventories and swept storage must be non-negative")
        signed_values = (
            self.interface_hexane_transfer_to_dry_mol,
            self.surface_hexane_out_mol,
        )
        if not all(math.isfinite(value) for value in signed_values):
            raise ValueError("birth ALE transfers must be finite")

    @property
    def dry_balance_residual_mol(self) -> float:
        return math.fsum(
            (
                self.dry_hexane_inventory_after_mol,
                -self.dry_hexane_inventory_before_mol,
                -self.interface_hexane_transfer_to_dry_mol,
                -self.swept_dry_hexane_storage_mol,
                self.surface_hexane_out_mol,
            )
        )

    @property
    def normalized_dry_balance_residual(self) -> float:
        scale = max(
            abs(self.dry_hexane_inventory_before_mol),
            abs(self.dry_hexane_inventory_after_mol),
            abs(self.interface_hexane_transfer_to_dry_mol),
            abs(self.swept_wet_hexane_inventory_mol),
            abs(self.swept_dry_hexane_storage_mol),
            abs(self.surface_hexane_out_mol),
            1.0e-300,
        )
        return abs(self.dry_balance_residual_mol) / scale

    @property
    def extensive_rh_residual_mol(self) -> float:
        return math.fsum(
            (
                self.interface_hexane_transfer_to_dry_mol,
                -self.swept_wet_hexane_inventory_mol,
                self.swept_dry_hexane_storage_mol,
            )
        )

    @property
    def normalized_extensive_rh_residual(self) -> float:
        scale = max(
            abs(self.interface_hexane_transfer_to_dry_mol),
            abs(self.swept_wet_hexane_inventory_mol),
            abs(self.swept_dry_hexane_storage_mol),
            1.0e-300,
        )
        return abs(self.extensive_rh_residual_mol) / scale


@dataclass(frozen=True)
class CoincidentEndpointCandidate:
    """One converged candidate from the simultaneous surface/front solve.

    Onset force and onset flux are independently recorded and
    outward-positive.  Access and step-mean active area are likewise explicit
    rather than inferred from a zero flux.  The finite-step mean surface flux
    is separate from the onset/interface flux because the newborn dry shell
    can accumulate.  ``mean_recession_speed_m_s`` is the exact integrated
    thickness change divided by ``dt_s`` and may not be negative on the frozen
    primary-drainage topology.  The dry thickness is actual accepted-step
    geometry, not an epsilon initializer.  ``accepted_stationary_inventory``
    is the candidate's conservative post-step water/hexane/energy state; it
    is exposed only if the commit barrier accepts the branch.
    """

    dt_s: float
    onset_outward_hexane_thermodynamic_force: float
    hexane_access_fraction: float
    active_surface_area_m2: float
    onset_surface_hexane_flux_mol_m2_s: float
    onset_interface_hexane_flux_mol_m2_s: float
    onset_recession_speed_m_s: float
    mean_surface_hexane_flux_mol_m2_s: float
    mean_recession_speed_m_s: float
    dry_shell_thickness_after_m: float
    wet_hexane_concentration_mol_m3: float
    dry_hexane_concentration_mol_m3: float
    residuals: CoincidentEndpointResiduals
    accepted_stationary_inventory: StationaryMaterialInventory
    birth: BirthRepartitionResult | None
    birth_ale_ledger: CoincidentBirthALELedger | None
    provenance: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.dt_s,
            self.onset_outward_hexane_thermodynamic_force,
            self.hexane_access_fraction,
            self.active_surface_area_m2,
            self.onset_surface_hexane_flux_mol_m2_s,
            self.onset_interface_hexane_flux_mol_m2_s,
            self.onset_recession_speed_m_s,
            self.mean_surface_hexane_flux_mol_m2_s,
            self.mean_recession_speed_m_s,
            self.dry_shell_thickness_after_m,
            self.wet_hexane_concentration_mol_m3,
            self.dry_hexane_concentration_mol_m3,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("coincident-endpoint candidate values must be finite")
        if self.dt_s <= 0.0:
            raise ValueError("coincident-endpoint time step must be positive")
        if not 0.0 <= self.hexane_access_fraction <= 1.0:
            raise ValueError("hexane access fraction must lie in [0,1]")
        if self.active_surface_area_m2 < 0.0:
            raise ValueError("active endpoint surface area may not be negative")
        if self.dry_shell_thickness_after_m < 0.0:
            raise ValueError("dry-shell thickness may not be negative")
        if (
            min(
                self.wet_hexane_concentration_mol_m3,
                self.dry_hexane_concentration_mol_m3,
            )
            < 0.0
        ):
            raise ValueError("endpoint hexane concentrations may not be negative")
        if not isinstance(self.residuals, CoincidentEndpointResiduals):
            raise TypeError("endpoint candidate needs a residual ledger")
        if not isinstance(self.accepted_stationary_inventory, StationaryMaterialInventory):
            raise TypeError("endpoint candidate needs its accepted material inventory")
        if self.birth is not None and not isinstance(self.birth, BirthRepartitionResult):
            raise TypeError("endpoint birth payload has the wrong type")
        if self.birth_ale_ledger is not None and not isinstance(
            self.birth_ale_ledger, CoincidentBirthALELedger
        ):
            raise TypeError("endpoint birth ALE payload has the wrong type")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("endpoint candidate needs explicit provenance")


@dataclass(frozen=True)
class CoincidentEndpointFailedAttempt:
    """A failed nonlinear attempt that owns no accepted state mutation."""

    kind: EndpointFailureKind
    attempted_dt_s: float
    positive_evaporation_drive: bool
    detail: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.kind, EndpointFailureKind):
            raise TypeError("endpoint failure kind must be explicit")
        if not math.isfinite(self.attempted_dt_s) or self.attempted_dt_s <= 0.0:
            raise ValueError("failed endpoint attempt needs a positive time step")
        if type(self.positive_evaporation_drive) is not bool:
            raise TypeError("positive-drive failure evidence must be boolean")
        if not isinstance(self.detail, str) or not self.detail.strip():
            raise ValueError("failed endpoint attempt needs a diagnostic detail")


@dataclass(frozen=True)
class CoincidentEndpointResolution:
    """Fail-closed A1 branch resolution with exact rollback semantics."""

    branch: CoincidentEndpointBranch
    stationary_inventory: StationaryMaterialInventory
    candidate: CoincidentEndpointCandidate | None
    failed_attempt: CoincidentEndpointFailedAttempt | None
    failure_kind: EndpointFailureKind | None
    recession_speed_m_s: float
    dry_shell_thickness_before_m: float
    dry_shell_thickness_after_m: float
    rh_normalized_residual: float
    kinematic_normalized_residual: float
    detail: str
    geometry_normalized_residual: float = 0.0
    ale_normalized_residual: float = 0.0
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def commit_allowed(self) -> bool:
        return self.branch in (
            CoincidentEndpointBranch.RECEDING_BIRTH,
            CoincidentEndpointBranch.ZERO_FLUX_CONTACT,
        )

    @property
    def rollback_required(self) -> bool:
        return self.branch is CoincidentEndpointBranch.ROLLBACK_RETRY

    @property
    def stop_or_flag_required(self) -> bool:
        return self.branch in (
            CoincidentEndpointBranch.INWARD_ENDPOINT_FLAG,
            CoincidentEndpointBranch.MODEL_TOPOLOGY_STOP,
        )

    @property
    def birth(self) -> BirthRepartitionResult | None:
        if self.branch is not CoincidentEndpointBranch.RECEDING_BIRTH or self.candidate is None:
            return None
        return self.candidate.birth


@dataclass(frozen=True)
class ExistingDryTransportAudit:
    """Distinguish conservative inward dry transport from reverse-core demand."""

    branch: ExistingDryTransportBranch
    dry_shell_thickness_m: float
    surface_hexane_flux_mol_m2_s: float
    recession_speed_m_s: float
    water_only_condensation_requested: bool
    stop_or_flag_required: bool
    physically_qualifying: bool = field(default=False, init=False)


def resolve_coincident_endpoint(
    inventory: StationaryMaterialInventory,
    attempt: CoincidentEndpointCandidate | CoincidentEndpointFailedAttempt,
    *,
    acceptance_tolerance: float,
) -> CoincidentEndpointResolution:
    """Resolve the A1 endpoint without a flux deadband or epsilon dry shell.

    The caller owns the simultaneous external-film, water, energy, and
    ALE/Rankine--Hugoniot solve.  This function is its bounded
    production-facing commit audit, not that solver.  Only a fully checked
    receding or exact-zero-flux contact candidate can commit.  With the
    current whole-material-shell repartition primitive, receding commit is
    deliberately limited to complete exposure of the outer stationary shell:
    gas-accessible fraction must be one and its exact shell volume must equal
    the volume swept by the accepted thickness.  Partial-cell evidence is
    rejected until a mixed wet/dry stationary-cell update exists.  A failed
    attempt returns the original immutable inventory by identity and is
    classified for rollback/retry, never as contact.
    """
    if not isinstance(inventory, StationaryMaterialInventory):
        raise TypeError("coincident endpoint needs a stationary material inventory")
    if (
        not math.isfinite(acceptance_tolerance)
        or acceptance_tolerance <= 0.0
        or acceptance_tolerance > 1.0e-10
    ):
        raise ValueError("endpoint acceptance tolerance must lie in (0, 1e-10]")
    if isinstance(attempt, CoincidentEndpointFailedAttempt):
        return _endpoint_resolution(
            CoincidentEndpointBranch.ROLLBACK_RETRY,
            inventory,
            None,
            attempt,
            attempt.kind,
            0.0,
            0.0,
            0.0,
            "failed nonlinear attempt rolled back; refine or continue and retry",
        )
    if not isinstance(attempt, CoincidentEndpointCandidate):
        raise TypeError("endpoint attempt has the wrong type")
    if attempt.accepted_stationary_inventory.grid != inventory.grid:
        raise ValueError("endpoint candidate changed the stationary material grid")

    residuals = attempt.residuals
    if residuals.max_normalized_residual > acceptance_tolerance:
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.COUPLED_RESIDUAL_NOT_CONVERGED,
            "component, energy, equilibrium, or active-set residual did not converge",
        )
    if not (residuals.property_bounds_admissible and residuals.phase_topology_admissible):
        return _endpoint_resolution(
            CoincidentEndpointBranch.MODEL_TOPOLOGY_STOP,
            inventory,
            attempt,
            None,
            EndpointFailureKind.CONVERGED_MODEL_TOPOLOGY,
            0.0,
            0.0,
            0.0,
            "converged endpoint candidate is outside a property or phase topology",
        )

    onset_surface_flux = attempt.onset_surface_hexane_flux_mol_m2_s
    onset_interface_flux = attempt.onset_interface_hexane_flux_mol_m2_s
    onset_speed = attempt.onset_recession_speed_m_s
    mean_surface_flux = attempt.mean_surface_hexane_flux_mol_m2_s
    mean_speed = attempt.mean_recession_speed_m_s
    thickness = attempt.dry_shell_thickness_after_m
    if onset_surface_flux == 0.0:
        if (
            onset_interface_flux != 0.0
            or onset_speed != 0.0
            or mean_surface_flux != 0.0
            or mean_speed != 0.0
            or thickness != 0.0
            or attempt.birth is not None
            or attempt.birth_ale_ledger is not None
        ):
            return _retry_endpoint(
                inventory,
                attempt,
                EndpointFailureKind.COUPLED_RESIDUAL_NOT_CONVERGED,
                "zero-flux contact candidate moved or created the dry domain",
            )
        if attempt.accepted_stationary_inventory.hexane_masses_kg != (inventory.hexane_masses_kg):
            return _retry_endpoint(
                inventory,
                attempt,
                EndpointFailureKind.INCOMPATIBLE_BIRTH_LEDGER,
                "zero-flux contact changed a shell-conserved wet hexane inventory",
            )
        return _endpoint_resolution(
            CoincidentEndpointBranch.ZERO_FLUX_CONTACT,
            inventory,
            attempt,
            None,
            None,
            0.0,
            0.0,
            0.0,
            "exact zero-flux contact accepted; wet heat and water equations remain active",
        )

    if onset_surface_flux < 0.0:
        if (
            onset_interface_flux != onset_surface_flux
            or onset_speed != 0.0
            or mean_surface_flux >= 0.0
            or mean_speed != 0.0
            or thickness != 0.0
            or attempt.birth is not None
            or attempt.birth_ale_ledger is not None
        ):
            return _endpoint_resolution(
                CoincidentEndpointBranch.MODEL_TOPOLOGY_STOP,
                inventory,
                attempt,
                None,
                EndpointFailureKind.CONVERGED_MODEL_TOPOLOGY,
                0.0,
                0.0,
                0.0,
                "inward endpoint candidate attempted front motion or dry-domain birth",
            )
        return _endpoint_resolution(
            CoincidentEndpointBranch.INWARD_ENDPOINT_FLAG,
            inventory,
            attempt,
            None,
            None,
            0.0,
            0.0,
            0.0,
            "inward endpoint tendency is outside the qualified monotonic branch",
        )

    if (
        attempt.onset_outward_hexane_thermodynamic_force <= 0.0
        or attempt.hexane_access_fraction <= 0.0
        or attempt.active_surface_area_m2 <= 0.0
    ):
        return _endpoint_resolution(
            CoincidentEndpointBranch.MODEL_TOPOLOGY_STOP,
            inventory,
            attempt,
            None,
            EndpointFailureKind.CONVERGED_MODEL_TOPOLOGY,
            0.0,
            0.0,
            0.0,
            "positive flux lacks positive force, access, or active area evidence",
        )

    birth = attempt.birth
    ale = attempt.birth_ale_ledger
    if birth is None or ale is None:
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.MISSING_COUPLED_BIRTH,
            "positive outward flux requires simultaneous birth and ALE ledgers",
        )
    if birth.stationary_inventory is not inventory:
        raise ValueError("endpoint birth was not solved from this inventory object")
    if birth.coupling.is_zero_transfer:
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.INCOMPATIBLE_BIRTH_LEDGER,
            "positive birth may not reuse the isolated zero-transfer map",
        )
    if (
        birth.ledger.dry_pore_volume_before_m3 != 0.0
        or birth.ledger.dry_pore_volume_after_m3 <= 0.0
        or birth.ledger.artificial_water_source_kg != 0.0
        or birth.ledger.artificial_energy_source_j != 0.0
        or birth.ledger.max_normalized_residual > acceptance_tolerance
    ):
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.INCOMPATIBLE_BIRTH_LEDGER,
            "birth payload failed zero-old-volume or conservative-ledger checks",
        )
    if birth.shell_index != inventory.grid.n - 1:
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.INCOMPATIBLE_BIRTH_LEDGER,
            "first dry material shell is not the stationary surface shell",
        )
    if birth.state.gas_accessible_fraction != 1.0:
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.INCOMPATIBLE_BIRTH_LEDGER,
            "partial-cell whole-shell repartition is not A1 swept-volume evidence",
        )
    expected_dry = birth.state.total_hexane_concentration_mol_m3
    dry_scale = max(abs(expected_dry), abs(attempt.dry_hexane_concentration_mol_m3), 1.0)
    dry_mismatch = abs(attempt.dry_hexane_concentration_mol_m3 - expected_dry) / dry_scale
    concentration_jump = math.fsum(
        (
            attempt.wet_hexane_concentration_mol_m3,
            -attempt.dry_hexane_concentration_mol_m3,
        )
    )
    if concentration_jump <= 0.0 or dry_mismatch > acceptance_tolerance:
        return _endpoint_resolution(
            CoincidentEndpointBranch.MODEL_TOPOLOGY_STOP,
            inventory,
            attempt,
            None,
            EndpointFailureKind.CONVERGED_MODEL_TOPOLOGY,
            0.0,
            dry_mismatch,
            0.0,
            "receding candidate lacks a positive, birth-consistent storage jump",
        )
    if onset_interface_flux <= 0.0 or onset_speed <= 0.0 or mean_speed <= 0.0 or thickness <= 0.0:
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.COUPLED_RESIDUAL_NOT_CONVERGED,
            "positive outward flux requires positive recession and dry thickness",
        )
    if thickness >= inventory.grid.R:
        return _retry_endpoint(
            inventory,
            attempt,
            EndpointFailureKind.OVERLARGE_EVENT_STEP,
            "event substep crossed the entire particle and must be refined",
        )

    radius = inventory.grid.R
    swept_bulk_volume = 4.0 / 3.0 * math.pi * (radius**3 - (radius - thickness) ** 3)
    exposed_material_volume = inventory.grid.volumes[birth.shell_index]
    geometry_scale = max(
        abs(swept_bulk_volume),
        abs(exposed_material_volume),
        1.0e-300,
    )
    geometry_residual = (
        abs(math.fsum((swept_bulk_volume, -exposed_material_volume))) / geometry_scale
    )
    accepted_inventory = attempt.accepted_stationary_inventory
    ledger = birth.ledger
    shell = birth.shell_index
    unchanged_shells = all(
        (
            accepted_inventory.water_masses_kg[index] == inventory.water_masses_kg[index]
            and accepted_inventory.hexane_masses_kg[index] == inventory.hexane_masses_kg[index]
            and accepted_inventory.energies_j[index] == inventory.energies_j[index]
        )
        for index in range(inventory.grid.n)
        if index != shell
    )
    accepted_values = (
        accepted_inventory.water_masses_kg[shell],
        accepted_inventory.hexane_masses_kg[shell],
        accepted_inventory.energies_j[shell],
    )
    ledger_values = (
        ledger.water_mass_after_kg,
        ledger.hexane_mass_after_kg,
        ledger.energy_after_j,
    )
    inventory_scales = tuple(
        max(abs(accepted), abs(expected), 1.0e-300)
        for accepted, expected in zip(accepted_values, ledger_values)
    )
    inventory_residual = max(
        abs(math.fsum((accepted, -expected))) / scale
        for accepted, expected, scale in zip(
            accepted_values,
            ledger_values,
            inventory_scales,
        )
    )
    dry_after_expected = ledger.hexane_mass_after_kg / hx.M
    dry_inventory_scale = max(
        abs(ale.dry_hexane_inventory_after_mol),
        abs(dry_after_expected),
        1.0e-300,
    )
    dry_inventory_residual = (
        abs(math.fsum((ale.dry_hexane_inventory_after_mol, -dry_after_expected)))
        / dry_inventory_scale
    )
    swept_wet_expected = attempt.wet_hexane_concentration_mol_m3 * swept_bulk_volume
    swept_dry_expected = attempt.dry_hexane_concentration_mol_m3 * swept_bulk_volume
    swept_scale = max(
        abs(ale.swept_wet_hexane_inventory_mol),
        abs(ale.swept_dry_hexane_storage_mol),
        abs(swept_wet_expected),
        abs(swept_dry_expected),
        1.0e-300,
    )
    swept_residual = max(
        abs(math.fsum((ale.swept_wet_hexane_inventory_mol, -swept_wet_expected))) / swept_scale,
        abs(math.fsum((ale.swept_dry_hexane_storage_mol, -swept_dry_expected))) / swept_scale,
    )
    surface_out_from_mean = mean_surface_flux * attempt.active_surface_area_m2 * attempt.dt_s
    surface_scale = max(
        abs(ale.surface_hexane_out_mol),
        abs(surface_out_from_mean),
        1.0e-300,
    )
    mean_surface_residual = (
        abs(math.fsum((ale.surface_hexane_out_mol, -surface_out_from_mean))) / surface_scale
    )
    whole_before_mol = inventory.total_hexane_mass_kg / hx.M
    whole_after_mol = accepted_inventory.total_hexane_mass_kg / hx.M
    whole_scale = max(
        abs(whole_before_mol),
        abs(whole_after_mol),
        abs(ale.surface_hexane_out_mol),
        1.0e-300,
    )
    whole_residual = (
        abs(math.fsum((whole_after_mol, -whole_before_mol, ale.surface_hexane_out_mol)))
        / whole_scale
    )
    coupling_surface_residual = abs(
        math.fsum(
            (
                birth.coupling.hexane_mass_increment_kg / hx.M,
                ale.surface_hexane_out_mol,
            )
        )
    ) / max(
        abs(birth.coupling.hexane_mass_increment_kg / hx.M),
        abs(ale.surface_hexane_out_mol),
        1.0e-300,
    )
    ale_residual = max(
        ale.normalized_dry_balance_residual,
        ale.normalized_extensive_rh_residual,
        dry_inventory_residual,
        swept_residual,
        mean_surface_residual,
        whole_residual,
        coupling_surface_residual,
    )
    if (
        geometry_residual > acceptance_tolerance
        or inventory_residual > acceptance_tolerance
        or ale.dry_hexane_inventory_before_mol != 0.0
        or ale_residual > acceptance_tolerance
        or not unchanged_shells
    ):
        return _endpoint_resolution(
            CoincidentEndpointBranch.ROLLBACK_RETRY,
            inventory,
            attempt,
            None,
            EndpointFailureKind.INCOMPATIBLE_BIRTH_LEDGER,
            0.0,
            0.0,
            0.0,
            "birth geometry, ALE ledger, or accepted inventory is inconsistent",
            geometry_residual=geometry_residual,
            ale_residual=ale_residual,
        )

    onset_rh_flux = onset_speed * concentration_jump
    rh_scale = max(
        abs(onset_surface_flux),
        abs(onset_interface_flux),
        abs(onset_rh_flux),
        1.0e-300,
    )
    rh_residual = max(
        abs(math.fsum((onset_surface_flux, -onset_interface_flux))) / rh_scale,
        abs(math.fsum((onset_interface_flux, -onset_rh_flux))) / rh_scale,
    )
    kinematic_scale = max(
        abs(thickness),
        abs(attempt.dt_s * mean_speed),
        inventory.grid.R,
        1.0e-300,
    )
    kinematic_residual = abs(math.fsum((thickness, -attempt.dt_s * mean_speed))) / kinematic_scale
    if max(rh_residual, kinematic_residual) > acceptance_tolerance:
        return _endpoint_resolution(
            CoincidentEndpointBranch.ROLLBACK_RETRY,
            inventory,
            attempt,
            None,
            EndpointFailureKind.COUPLED_RESIDUAL_NOT_CONVERGED,
            rh_residual,
            kinematic_residual,
            0.0,
            "Rankine--Hugoniot or endpoint kinematic residual did not converge",
        )
    return _endpoint_resolution(
        CoincidentEndpointBranch.RECEDING_BIRTH,
        inventory,
        attempt,
        None,
        None,
        rh_residual,
        kinematic_residual,
        mean_speed,
        "positive outward flux accepted on the simultaneous zero-old-volume birth branch",
        geometry_residual=geometry_residual,
        ale_residual=ale_residual,
    )


def audit_existing_dry_transport(
    *,
    dry_shell_thickness_m: float,
    surface_hexane_flux_mol_m2_s: float,
    recession_speed_m_s: float,
    hexane_condensation_or_rewetting_requested: bool,
    water_only_condensation_requested: bool,
    attached_inventory_recreation_requested: bool,
    surface_reset_requested: bool,
    replacement_core_requested: bool,
) -> ExistingDryTransportAudit:
    """Apply A1/PHY-014 after recession without banning inward gas transport."""
    values = (
        dry_shell_thickness_m,
        surface_hexane_flux_mol_m2_s,
        recession_speed_m_s,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("existing-dry transport audit values must be finite")
    if dry_shell_thickness_m <= 0.0:
        raise ValueError("existing-dry audit requires a positive dry-shell thickness")
    switches = (
        hexane_condensation_or_rewetting_requested,
        water_only_condensation_requested,
        attached_inventory_recreation_requested,
        surface_reset_requested,
        replacement_core_requested,
    )
    if any(type(value) is not bool for value in switches):
        raise TypeError("existing-dry topology requests must be boolean")
    reverse_requested = any(
        (
            recession_speed_m_s < 0.0,
            hexane_condensation_or_rewetting_requested,
            attached_inventory_recreation_requested,
            surface_reset_requested,
            replacement_core_requested,
        )
    )
    if reverse_requested:
        branch = ExistingDryTransportBranch.REVERSE_CORE_STOP
    elif surface_hexane_flux_mol_m2_s < 0.0:
        branch = ExistingDryTransportBranch.CONSERVATIVE_INWARD_DRY_TRANSPORT
    else:
        branch = ExistingDryTransportBranch.OUTWARD_OR_ZERO
    return ExistingDryTransportAudit(
        branch=branch,
        dry_shell_thickness_m=dry_shell_thickness_m,
        surface_hexane_flux_mol_m2_s=surface_hexane_flux_mol_m2_s,
        recession_speed_m_s=recession_speed_m_s,
        water_only_condensation_requested=water_only_condensation_requested,
        stop_or_flag_required=(branch is ExistingDryTransportBranch.REVERSE_CORE_STOP),
    )


def _retry_endpoint(
    inventory: StationaryMaterialInventory,
    candidate: CoincidentEndpointCandidate,
    kind: EndpointFailureKind,
    detail: str,
) -> CoincidentEndpointResolution:
    return _endpoint_resolution(
        CoincidentEndpointBranch.ROLLBACK_RETRY,
        inventory,
        candidate,
        None,
        kind,
        0.0,
        0.0,
        0.0,
        detail,
    )


def _endpoint_resolution(
    branch: CoincidentEndpointBranch,
    inventory: StationaryMaterialInventory,
    candidate: CoincidentEndpointCandidate | None,
    failed_attempt: CoincidentEndpointFailedAttempt | None,
    failure_kind: EndpointFailureKind | None,
    rh_residual: float,
    kinematic_residual: float,
    speed: float,
    detail: str,
    *,
    geometry_residual: float = 0.0,
    ale_residual: float = 0.0,
) -> CoincidentEndpointResolution:
    thickness_after = 0.0 if candidate is None else candidate.dry_shell_thickness_after_m
    if branch is not CoincidentEndpointBranch.RECEDING_BIRTH:
        thickness_after = 0.0
        speed = 0.0
    resolved_inventory = inventory
    if (
        branch
        in (
            CoincidentEndpointBranch.RECEDING_BIRTH,
            CoincidentEndpointBranch.ZERO_FLUX_CONTACT,
        )
        and candidate is not None
    ):
        resolved_inventory = candidate.accepted_stationary_inventory
    return CoincidentEndpointResolution(
        branch=branch,
        stationary_inventory=resolved_inventory,
        candidate=candidate,
        failed_attempt=failed_attempt,
        failure_kind=failure_kind,
        recession_speed_m_s=speed,
        dry_shell_thickness_before_m=0.0,
        dry_shell_thickness_after_m=thickness_after,
        rh_normalized_residual=rh_residual,
        kinematic_normalized_residual=kinematic_residual,
        detail=detail,
        geometry_normalized_residual=geometry_residual,
        ale_normalized_residual=ale_residual,
    )


def repartition_birth_shell(
    inventory: StationaryMaterialInventory,
    shell_index: int,
    pressure_pa: float,
    gas_accessible_fraction: float,
    coupling: BirthCouplingIncrements,
    initial_temperature_k: float,
    initial_y_hexane: float,
    params: CoupledPoreParams = CoupledPoreParams(),
    *,
    tolerance: float = 2.0e-11,
    max_iterations: int = 28,
) -> BirthRepartitionResult:
    """Birth dry pore volume from the shell's existing water and energy.

    The old dry pore volume is exactly zero.  The two fixed-pressure local
    primitives ``(T,y_h)`` are recovered from coupled post-transfer water and
    energy.  The resulting n-hexane inventory must also match the explicit
    front/DAE increment.  This preserves the fixed-pressure rank count instead
    of pretending three conserved totals can be inverted from two primitives
    or that an actual activation birth is a zero-transfer algebraic map.
    """
    if not isinstance(coupling, BirthCouplingIncrements):
        raise TypeError("birth requires explicit conservative coupling increments")

    def inversion_failure(message: str) -> BirthRepartitionError:
        if coupling.is_zero_transfer:
            return CoupledBirthStepRequiredError(
                "NEEDS_COUPLED_RH_STEP: zero-transfer activation birth is not a "
                f"closed fixed-pressure map; {message}"
            )
        return BirthRepartitionError(message)

    try:
        index = operator.index(shell_index)
    except TypeError as exc:
        raise ValueError("birth shell index must be an integer") from exc
    if isinstance(shell_index, bool):
        raise ValueError("birth shell index must be an integer")
    if not 0 <= index < inventory.grid.n:
        raise ValueError("birth shell index is outside the stationary grid")
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("birth inversion tolerance must be positive and finite")
    if max_iterations < 1:
        raise ValueError("birth inversion iteration budget must be positive")
    if not math.isfinite(gas_accessible_fraction) or not 0.0 < gas_accessible_fraction <= 1.0:
        raise ValueError("birth requires a positive gas-accessible fraction in (0,1]")

    volume = inventory.grid.volumes[index]
    target_water_mass = math.fsum(
        (inventory.water_masses_kg[index], coupling.water_mass_increment_kg)
    )
    target_hexane_mass = math.fsum(
        (inventory.hexane_masses_kg[index], coupling.hexane_mass_increment_kg)
    )
    target_energy_j = math.fsum((inventory.energies_j[index], coupling.energy_increment_j))
    if target_water_mass <= 0.0:
        raise BirthRepartitionError("coupled birth target water must remain positive")
    if target_hexane_mass < 0.0:
        raise BirthRepartitionError("coupled birth target hexane became negative")
    target_water = target_water_mass / volume
    target_energy = target_energy_j / volume
    scales = np.array(
        (
            max(abs(target_water), params.rho_dm_p * params.luikov.W_ref, 1.0),
            max(
                abs(target_energy),
                params.rho_dm_p
                * sp.composite_dry_meal_heat_capacity(params.w_o, params.cp_dry_meal, params.cp_oil)
                * max(initial_temperature_k, 1.0),
                1.0e6,
            ),
        ),
        dtype=float,
    )
    transformed = np.array(
        (
            _temperature_to_raw(initial_temperature_k, params),
            _logit(initial_y_hexane),
        ),
        dtype=float,
    )
    last_norm = math.inf
    accepted_state: EquilibriumPoreState | None = None

    def evaluate_raw(values: np.ndarray) -> EquilibriumPoreState:
        temperature, y_hexane = _decode_birth_coordinates(values, params)
        return evaluate_equilibrium(
            temperature,
            pressure_pa,
            y_hexane,
            params,
            gas_accessible_fraction=gas_accessible_fraction,
        )

    def residual_for(state: EquilibriumPoreState) -> np.ndarray:
        return (
            np.array(
                (
                    state.total_water_concentration_kg_m3 - target_water,
                    state.energy_density_j_m3 - target_energy,
                ),
                dtype=float,
            )
            / scales
        )

    for iteration in range(max_iterations + 1):
        try:
            state = evaluate_raw(transformed)
        except (ValueError, RuntimeError) as exc:
            raise inversion_failure(
                "initial or accepted birth iterate is outside the approved topology"
            ) from exc
        residual = residual_for(state)
        last_norm = float(np.max(np.abs(residual)))
        if last_norm <= tolerance:
            accepted_state = state
            break
        if iteration == max_iterations:
            break

        jacobian = np.empty((2, 2), dtype=float)
        for column in range(2):
            derivative_found = False
            step = 2.0e-5
            for _ in range(8):
                for direction in (1.0, -1.0):
                    trial_coordinates = transformed.copy()
                    trial_coordinates[column] += direction * step
                    try:
                        trial_state = evaluate_raw(trial_coordinates)
                    except ValueError, RuntimeError:
                        continue
                    jacobian[:, column] = (residual_for(trial_state) - residual) / (
                        direction * step
                    )
                    derivative_found = True
                    break
                if derivative_found:
                    break
                step *= 0.25
            if not derivative_found:
                raise inversion_failure("birth Jacobian has no admissible local perturbation")
        if not np.all(np.isfinite(jacobian)):
            raise inversion_failure("birth Jacobian is non-finite")
        try:
            update = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError as exc:
            raise inversion_failure("birth Jacobian is singular") from exc
        if not np.all(np.isfinite(update)):
            raise inversion_failure("birth Newton update is non-finite")
        update_norm = float(np.linalg.norm(update))
        if update_norm > 4.0:
            update *= 4.0 / update_norm

        alpha = 1.0
        installed = False
        while alpha >= 2.0**-18:
            trial_coordinates = transformed + alpha * update
            try:
                trial_state = evaluate_raw(trial_coordinates)
            except ValueError, RuntimeError:
                alpha *= 0.5
                continue
            trial_norm = float(np.max(np.abs(residual_for(trial_state))))
            if trial_norm < last_norm:
                transformed = trial_coordinates
                installed = True
                break
            alpha *= 0.5
        if not installed:
            raise inversion_failure(
                "birth Newton line search found no admissible residual decrease"
            )
    else:  # pragma: no cover - loop always exits through range exhaustion
        iteration = max_iterations

    if accepted_state is None:
        raise inversion_failure(
            "birth repartition did not converge within the iteration budget "
            f"(max scaled residual={last_norm:.6e})"
        )

    retained_water_mass = volume * accepted_state.retained_water_concentration_kg_m3
    pore_water_mass = volume * accepted_state.pore_water_concentration_kg_m3
    water_after = math.fsum((retained_water_mass, pore_water_mass))
    hexane_after = volume * accepted_state.total_hexane_concentration_kg_m3
    energy_after = volume * accepted_state.energy_density_j_m3
    water_before = inventory.water_masses_kg[index]
    hexane_before = inventory.hexane_masses_kg[index]
    energy_before = inventory.energies_j[index]
    water_residual = math.fsum((water_after, -target_water_mass))
    hexane_residual = math.fsum((hexane_after, -target_hexane_mass))
    energy_residual = math.fsum((energy_after, -target_energy_j))
    water_scale = max(
        abs(target_water_mass),
        abs(retained_water_mass),
        abs(pore_water_mass),
        abs(coupling.water_mass_increment_kg),
        1.0e-300,
    )
    hexane_scale = max(
        abs(target_hexane_mass),
        abs(hexane_after),
        abs(coupling.hexane_mass_increment_kg),
        1.0e-300,
    )
    energy_scale = max(
        abs(target_energy_j),
        abs(energy_after),
        abs(coupling.energy_increment_j),
        volume
        * params.rho_dm_p
        * sp.composite_dry_meal_heat_capacity(params.w_o, params.cp_dry_meal, params.cp_oil)
        * max(accepted_state.temperature_k, 1.0),
        1.0e-300,
    )
    normalized_water = abs(water_residual) / water_scale
    normalized_hexane = abs(hexane_residual) / hexane_scale
    normalized_energy = abs(energy_residual) / energy_scale
    ledger = BirthRepartitionLedger(
        shell_index=index,
        dry_pore_volume_before_m3=0.0,
        dry_pore_volume_after_m3=(volume * params.eps_g * gas_accessible_fraction),
        water_mass_before_kg=water_before,
        hexane_mass_before_kg=hexane_before,
        retained_water_mass_after_kg=retained_water_mass,
        pore_water_mass_after_kg=pore_water_mass,
        water_mass_after_kg=water_after,
        hexane_mass_after_kg=hexane_after,
        energy_before_j=energy_before,
        energy_after_j=energy_after,
        coupled_water_increment_kg=coupling.water_mass_increment_kg,
        coupled_hexane_increment_kg=coupling.hexane_mass_increment_kg,
        coupled_energy_increment_j=coupling.energy_increment_j,
        artificial_water_source_kg=0.0,
        artificial_energy_source_j=0.0,
        water_residual_kg=water_residual,
        hexane_residual_kg=hexane_residual,
        energy_residual_j=energy_residual,
        normalized_water_residual=normalized_water,
        normalized_hexane_residual=normalized_hexane,
        normalized_energy_residual=normalized_energy,
        max_normalized_residual=max(normalized_water, normalized_hexane, normalized_energy),
    )
    if ledger.max_normalized_residual > 1.0e-10:
        raise inversion_failure(
            "birth repartition failed the independent 1e-10 component/energy ledger; "
            "the supplied front/DAE increments are incompatible with the fixed-pressure state"
        )
    return BirthRepartitionResult(
        stationary_inventory=inventory,
        coupling=coupling,
        shell_index=index,
        state=accepted_state,
        resulting_hexane_mass_kg=hexane_after,
        iterations=iteration,
        max_scaled_inversion_residual=last_norm,
        ledger=ledger,
    )


def _validate_params(params: CoupledPoreParams) -> None:
    values = (
        params.eps_g,
        params.rho_dm_p,
        params.cp_dry_meal,
        params.cp_oil,
        params.T_ref_solid,
        params.w_o,
        params.k_wh,
        params.retained_water_liquid_volume_fraction,
        *params.temperature_bounds_k,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("coupled-pore parameters must be finite")
    if not 0.0 < params.eps_g < 1.0:
        raise ValueError("pore-gas fraction must lie strictly inside (0,1)")
    if min(params.rho_dm_p, params.cp_dry_meal, params.cp_oil) <= 0.0:
        raise ValueError("material density and heat capacities must be positive")
    if not sp.W_O_MIN <= params.w_o <= sp.W_O_MAX:
        raise ValueError("residual-oil label is outside the frozen envelope")
    if not bg.K_WH_MIN <= params.k_wh <= bg.K_WH_MAX:
        raise ValueError("k_wh is outside the frozen qualification interval")
    if not 0.0 <= params.retained_water_liquid_volume_fraction <= 1.0:
        raise ValueError("retained-water partial-volume sensitivity fraction must lie in [0,1]")
    lower, upper = params.temperature_bounds_k
    if not 0.0 < lower < upper < hx.TC:
        raise ValueError("coupled-pore temperature bounds must be ordered and subcritical")
    luikov_values = (
        params.luikov.A1,
        params.luikov.A2,
        params.luikov.W_ref,
        params.luikov.W_cap,
    )
    continued = isinstance(params.luikov, sp.ContinuedPositiveLuikovParams)
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (luikov_values[:2] + luikov_values[3:] if continued else luikov_values)
    ):
        raise ValueError("modified-Luikov parameters must be positive and finite")
    if params.luikov.W_ref >= params.luikov.W_cap:
        raise ValueError("modified-Luikov retained-water bounds must be ordered")
    # Exercise the exact frozen inverse at both active-set boundaries.
    sp.water_activity(params.luikov.W_ref, params.luikov)
    sp.water_activity(params.luikov.W_cap, params.luikov)


def _validate_gas_only_temperature_pressure(
    temperature_k: float,
    pressure_pa: float,
    params: CoupledPoreParams,
) -> None:
    if not isinstance(params, CoupledPoreParams):
        raise TypeError("gas-only composition utilities require CoupledPoreParams")
    _validate_params(params)
    if not math.isfinite(temperature_k) or not math.isfinite(pressure_pa):
        raise CoupledPoreTopologyError("gas-only temperature and pressure must be finite")
    lower, upper = params.temperature_bounds_k
    if not lower <= temperature_k <= upper:
        raise CoupledPoreTopologyError(
            f"temperature {temperature_k} K is outside [{lower},{upper}] K"
        )
    if pressure_pa <= 0.0:
        raise CoupledPoreTopologyError("gas-only pressure must be positive")


def _gas_component_activities(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    params: CoupledPoreParams,
) -> tuple[float, float]:
    if not math.isfinite(y_hexane) or not 0.0 <= y_hexane <= 1.0:
        raise CoupledPoreTopologyError("activity evaluation requires 0 <= y_hexane <= 1")
    mixture = bg.state(
        temperature_k,
        pressure_pa,
        1.0 - y_hexane,
        y_hexane,
        k_wh=params.k_wh,
    )
    liquid_water = wa.state_Tp(temperature_k, pressure_pa, "liquid")
    liquid_hexane = hx.state_Tp(temperature_k, pressure_pa, "liquid")
    water_activity = mixture.fugacity_water_pa / liquid_water.fugacity
    hexane_activity = mixture.fugacity_hexane_pa / liquid_hexane.fugacity
    if not all(
        math.isfinite(value) and value >= 0.0 for value in (water_activity, hexane_activity)
    ):
        raise CoupledPoreTopologyError("gas-only component activities must be nonnegative")
    return water_activity, hexane_activity


def _gas_component_activities_with_derivatives(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    params: CoupledPoreParams,
) -> tuple[float, float, float, float]:
    r"""Return activities and their exact frozen-virial ``d/dy_h`` values."""

    if not math.isfinite(y_hexane) or not 0.0 <= y_hexane <= 1.0:
        raise CoupledPoreTopologyError("activity evaluation requires 0 <= y_hexane <= 1")
    mixture = bg.state(
        temperature_k,
        pressure_pa,
        1.0 - y_hexane,
        y_hexane,
        k_wh=params.k_wh,
    )
    liquid_water = wa.state_Tp(temperature_k, pressure_pa, "liquid")
    liquid_hexane = hx.state_Tp(temperature_k, pressure_pa, "liquid")
    water_activity = mixture.fugacity_water_pa / liquid_water.fugacity
    hexane_activity = mixture.fugacity_hexane_pa / liquid_hexane.fugacity

    # For the pressure-second-virial binary closure, writing
    # DeltaB = B_ww - 2 B_wh + B_hh gives
    #
    # d ln(phi_w)/dy_h = -2 P y_h DeltaB/(R T),
    # d ln(phi_h)/dy_h =  2 P y_w DeltaB/(R T).
    #
    # The product forms below remain finite at y_h=0 and y_h=1, unlike a
    # logarithmic derivative of the vanishing component activity.
    delta_virial = mixture.B_ww - 2.0 * mixture.B_wh + mixture.B_hh
    virial_factor = 2.0 * pressure_pa * delta_virial / (bg.R * temperature_k)
    water_fugacity_factor = math.exp(mixture.ln_phi_water) * pressure_pa / liquid_water.fugacity
    hexane_fugacity_factor = math.exp(mixture.ln_phi_hexane) * pressure_pa / liquid_hexane.fugacity
    water_activity_derivative = water_fugacity_factor * (
        -1.0 - (1.0 - y_hexane) * virial_factor * y_hexane
    )
    hexane_activity_derivative = hexane_fugacity_factor * (
        1.0 + y_hexane * virial_factor * (1.0 - y_hexane)
    )
    values = (
        water_activity,
        hexane_activity,
        water_activity_derivative,
        hexane_activity_derivative,
    )
    if not all(math.isfinite(value) for value in values):
        raise CoupledPoreTopologyError(
            "gas-only component activities and derivatives must be finite"
        )
    if (
        water_activity < 0.0
        or hexane_activity < 0.0
        or water_activity_derivative >= 0.0
        or hexane_activity_derivative <= 0.0
    ):
        raise CoupledPoreTopologyError(
            "frozen binary activities lost their required composition monotonicity"
        )
    return values


def _bisect_activity_composition(
    temperature_k: float,
    pressure_pa: float,
    params: CoupledPoreParams,
    *,
    component: str,
    target: float,
) -> float:
    if component not in {"water", "hexane"}:
        raise ValueError("activity root component must be water or hexane")
    component_index = 0 if component == "water" else 1

    def residual_and_derivative(y_hexane: float) -> tuple[float, float]:
        values = _gas_component_activities_with_derivatives(
            temperature_k, pressure_pa, y_hexane, params
        )
        return values[component_index] - target, values[component_index + 2]

    lower = 0.0
    upper = 1.0
    lower_residual, _ = residual_and_derivative(lower)
    upper_residual, _ = residual_and_derivative(upper)
    if lower_residual == 0.0:
        return lower
    if upper_residual == 0.0:
        return upper
    if lower_residual * upper_residual > 0.0:
        raise CoupledPoreTopologyError(
            f"{component} activity boundary is not bracketed at this temperature/pressure"
        )
    best_composition = lower
    best_absolute_residual = abs(lower_residual)
    if abs(upper_residual) < best_absolute_residual:
        best_composition = upper
        best_absolute_residual = abs(upper_residual)
    composition = lower - lower_residual * (upper - lower) / (upper_residual - lower_residual)
    if not lower < composition < upper:
        composition = lower + 0.5 * (upper - lower)
    for _ in range(64):
        if not lower < composition < upper:
            break
        composition_residual, composition_derivative = residual_and_derivative(composition)
        if abs(composition_residual) < best_absolute_residual:
            best_composition = composition
            best_absolute_residual = abs(composition_residual)
        if abs(composition_residual) <= 16.0 * math.ulp(max(1.0, abs(target))):
            return composition
        if lower_residual * composition_residual <= 0.0:
            upper = composition
            upper_residual = composition_residual
        else:
            lower = composition
            lower_residual = composition_residual
        newton_composition = composition - composition_residual / composition_derivative
        if lower < newton_composition < upper:
            next_composition = newton_composition
        else:
            next_composition = lower + 0.5 * (upper - lower)
        if next_composition == composition:
            next_composition = lower + 0.5 * (upper - lower)
        if next_composition == lower or next_composition == upper:
            break
        composition = next_composition
    root = best_composition
    if not math.isfinite(root) or not 0.0 <= root <= 1.0:
        raise CoupledPoreTopologyError("activity boundary inversion failed")
    return root


def _gas_only_barrier_coordinate(
    y_hexane: float,
    water_activity: float,
    hexane_activity: float,
    activity_at_ref: float,
) -> float:
    gaps = (
        y_hexane,
        1.0 - y_hexane,
        1.0 - water_activity,
        1.0 - hexane_activity,
        water_activity - activity_at_ref,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in gaps):
        # RR3 (owner-accepted 2026-08-21): name the violated boundary so the
        # Amendment-18 semantics survive the wrapping on the live path.  The
        # exception CLASS stays the generic one on purpose - retyping the
        # barrier's hexane branch could newly satisfy the plateau catch in
        # cut_transport (a certified-path behavior change reserved for the
        # RR1 packet); only the diagnostic is sharpened here.
        if math.isfinite(hexane_activity) and 1.0 - hexane_activity <= 0.0:
            raise CoupledPoreTopologyError(
                "gas-only chart state violates a binary, saturation, or "
                "Luikov-lower boundary; do not clip "
                f"[hexane activity {hexane_activity!r} at or above saturation: "
                "a supersaturated state is outside this open chart; the "
                "equilibrium gate's typed semantics apply - activate the "
                "mobile-liquid topology, do not clamp]"
            )
        raise CoupledPoreTopologyError(
            "gas-only chart state violates a binary, saturation, or Luikov-lower "
            "boundary; do not clip "
            f"[gaps: y_h={gaps[0]!r}, 1-y_h={gaps[1]!r}, 1-a_w={gaps[2]!r}, "
            f"1-a_h={gaps[3]!r}, a_w-a_ref={gaps[4]!r}]"
        )
    coordinate = _GAS_ONLY_BARRIER_COORDINATE_SCALE * math.fsum(
        (
            math.log(gaps[0]),
            -math.log(gaps[1]),
            math.log(gaps[2]),
            -math.log(gaps[3]),
            -math.log(gaps[4]),
        )
    )
    if not math.isfinite(coordinate):
        raise CoupledPoreTopologyError("gas-only barrier coordinate is non-finite")
    return coordinate


def _gas_only_barrier_coordinate_derivative(
    y_hexane: float,
    water_activity: float,
    hexane_activity: float,
    water_activity_derivative: float,
    hexane_activity_derivative: float,
    activity_at_ref: float,
) -> float:
    """Return the analytic derivative of the exact barrier coordinate."""

    derivative = _GAS_ONLY_BARRIER_COORDINATE_SCALE * math.fsum(
        (
            1.0 / y_hexane,
            1.0 / (1.0 - y_hexane),
            -water_activity_derivative / (1.0 - water_activity),
            hexane_activity_derivative / (1.0 - hexane_activity),
            -water_activity_derivative / (water_activity - activity_at_ref),
        )
    )
    if not math.isfinite(derivative) or derivative <= 0.0:
        raise CoupledPoreTopologyError(
            "gas-only barrier lost its required strict composition monotonicity"
        )
    return derivative


def _validate_primitives(
    temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    gas_accessible_fraction: float,
    params: CoupledPoreParams,
) -> None:
    values = (
        temperature_k,
        pressure_pa,
        y_hexane,
        gas_accessible_fraction,
    )
    if not all(math.isfinite(value) for value in values):
        raise CoupledPoreTopologyError("pore-state primitives must be finite")
    lower, upper = params.temperature_bounds_k
    if not lower <= temperature_k <= upper:
        raise CoupledPoreTopologyError(
            f"temperature {temperature_k} K is outside [{lower},{upper}] K"
        )
    if pressure_pa <= 0.0:
        raise CoupledPoreTopologyError("pore pressure must be positive")
    if not 0.0 < y_hexane < 1.0:
        raise CoupledPoreTopologyError("coupled binary pore state requires 0 < y_hexane < 1")
    if not 0.0 <= gas_accessible_fraction <= 1.0:
        raise CoupledPoreTopologyError("gas-accessible pore fraction must lie in [0,1]")


def _validate_distance(distance_m: float) -> None:
    if not math.isfinite(distance_m) or distance_m <= 0.0:
        raise ValueError("face-state distance must be positive and finite")


def _validate_dt(dt_s: float) -> None:
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("component-balance time step must be positive and finite")


def _validate_cell_inventory(values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError("component inventories must align with grid cells")
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        raise ValueError("component inventories must be finite and non-negative")


def _validate_face_values(values: Sequence[float], expected: int, name: str) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} must align with grid faces")
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{name} must be finite")


def _logit(value: float) -> float:
    if not math.isfinite(value) or not 0.0 < value < 1.0:
        raise ValueError("birth initial fraction must lie strictly inside (0,1)")
    return math.log(value / (1.0 - value))


def _logistic(value: float) -> float:
    if value >= 0.0:
        exponential = math.exp(-value)
        return 1.0 / (1.0 + exponential)
    exponential = math.exp(value)
    return exponential / (1.0 + exponential)


def _temperature_to_raw(temperature_k: float, params: CoupledPoreParams) -> float:
    lower, upper = params.temperature_bounds_k
    if not math.isfinite(temperature_k) or not lower < temperature_k < upper:
        raise ValueError("birth initial temperature must lie strictly inside the declared bounds")
    return _logit((temperature_k - lower) / (upper - lower))


def _decode_birth_coordinates(
    values: Sequence[float], params: CoupledPoreParams
) -> tuple[float, float]:
    lower, upper = params.temperature_bounds_k
    temperature = lower + (upper - lower) * _logistic(float(values[0]))
    y_hexane = _logistic(float(values[1]))
    return temperature, y_hexane


__all__ = [
    "BirthCouplingIncrements",
    "BirthRepartitionError",
    "BirthRepartitionLedger",
    "BirthRepartitionResult",
    "CoincidentBirthALELedger",
    "CoincidentEndpointBranch",
    "CoincidentEndpointCandidate",
    "CoincidentEndpointFailedAttempt",
    "CoincidentEndpointResiduals",
    "CoincidentEndpointResolution",
    "ComponentBalanceResiduals",
    "ComponentEnergyFlux",
    "ComponentMolarFluxes",
    "CoupledPoreParams",
    "CoupledPoreTopologyError",
    "HexaneSupersaturationTopologyError",
    "CoupledBirthStepRequiredError",
    "EndpointFailureKind",
    "EquilibriumPoreState",
    "ExistingDryTransportAudit",
    "ExistingDryTransportBranch",
    "GasOnlyCompositionConstraint",
    "GasOnlyCompositionInterval",
    "IndependentMaxwellStefanFlux",
    "IsothermalPotentialCompositionDerivatives",
    "IsothermalPotentialCompositionSecantSlopes",
    "MobilityInterval",
    "MobilitySelection",
    "RetainedWaterActiveSet",
    "RetainedWaterFlux",
    "RetainedWaterSecantBranch",
    "StationaryMaterialInventory",
    "TotalStefanFluxRecovery",
    "audit_existing_dry_transport",
    "clear_equilibrium_cache",
    "component_balance_residuals",
    "component_energy_flux",
    "compose_component_fluxes",
    "conservative_component_update",
    "evaluate_equilibrium",
    "equilibrium_cache_info",
    "decode_gas_only_y",
    "encode_gas_only_y",
    "gas_only_composition_interval",
    "ideal_binary_mobility_mol_m_s",
    "independent_maxwell_stefan_flux",
    "isothermal_potential_composition_derivatives",
    "isothermal_potential_composition_secant_slopes",
    "inventory_from_activation",
    "recover_total_stefan_flux",
    "repartition_birth_shell",
    "resolve_coincident_endpoint",
    "retained_water_flux",
    "smooth_luikov_composition_lower_bound",
]
