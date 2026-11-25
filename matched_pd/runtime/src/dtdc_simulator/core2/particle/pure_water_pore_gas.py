r"""The water vertex (``y_hexane = 0``) of the binary dry-pore gas: pure steam in a dry shell.

Why this module exists
----------------------
The certified dry shell carries a **binary** pore gas, and every entry into that topology refuses a
composition on the boundary of the open simplex (``coupled_pore._validate_primitives``: "coupled
binary pore state requires 0 < y_hexane < 1"; ``transport_coefficients``: "an exact-purity face has
zero mobility and a singular logit force").  A dry particle stripped of n-hexane in hexane-free
superheated steam drives its cells towards that vertex without limit; the measured refusal is
``docs/GT_PS2_UNIVERSAL_KERNEL_STRIPPED_PARTICLE_IN_SUPERHEATED_STEAM_DIAGNOSIS_2026-11-01.md`` and
the owner ruled (2026-11-01) a declared n-hexane-exhausted regime.  This module is the water analogue
of :mod:`dtdc_simulator.core2.particle.single_component_pore_gas` (the pure n-hexane vertex): a
**separate, typed, non-qualifying** authority.  It edits nothing: every frozen closure it uses
(``props.binary_gas.state`` at ``(y_w, y_h) = (1, 0)``, the IAPWS water states, the modified-Luikov
retained-water isotherm and its cap, the GAB/oil retained n-hexane law and its caloric, the dry-meal
sensible energy, ``coupled_pore.retained_water_flux``, ``compose_component_fluxes`` and
``component_energy_flux``) is imported and called unchanged.

What the vertex is
------------------
*   **The cell state.**  With one component at fixed total pressure the pore gas has no composition
    degree of freedom: it is steam at the pore pressure in every cell.  The cell state is the binary
    state's own expressions evaluated with ``y_water = 1`` and ``y_hexane = 0`` exactly: the virial gas
    at the vertex (which the frozen closure admits; ``B_mix = B_ww``, the n-hexane fugacity exactly 0),
    the retained water on the unchanged Luikov/cap active set at ``a_w = f_w / f_w,liq``, the retained
    n-hexane loading and its energy from the frozen GAB/oil law at activity 0 (exactly 0 and -0.0), the
    common-datum energy density with every n-hexane term zero.  Nothing is extrapolated.
*   **No mutual diffusion.**  The independent binary flux ``J_w = -J_h`` is identically zero by the
    absence of the second component, not as a limit of the binary formula (whose exchange force
    ``ln(f_w/f_h)`` diverges at the vertex while its mobility ``c D y_w y_h`` vanishes).
*   **No composition force on the retained water.**  The kernel's dry faces use the common-face-
    temperature, negligible-Soret mode (``alpha_T = alpha_ret = 0``): an interior face evaluates both
    endpoint states at the face temperature and an outer face takes only the composition part of the
    force along its straight (T, y) path.  Two vertex endpoints have the same composition, so both
    forces are exactly zero, and the frozen ``retained_water_flux`` returns zero from the two identical
    potentials.  Water released by sorption leaves as the total (Stefan) molar flux of steam; heat moves
    by Fourier conduction and by the enthalpy that steam carries.  Any other mass-force mode is refused.

Scope and honesty
-----------------
``physically_qualifying`` is ``False`` on every object here.  The vertex is declared, not measured;
its agreement with the binary authority at small ``y_hexane`` is the regime's stated error (the
harness build record ``docs/GT_PS2_UNIVERSAL_KERNEL_HEXANE_EXHAUSTED_REGIME_BUILD_2026-11-01.md``).  A
caller holding any n-hexane must use the binary authority.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import nonisothermal_potential as nip
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa

RECORD = "docs/GT_PS2_UNIVERSAL_KERNEL_HEXANE_EXHAUSTED_REGIME_BUILD_2026-11-01.md"

#: The only composition this module admits, exactly.
ADMITTED_WATER_MOLE_FRACTION = 1.0
ADMITTED_HEXANE_MOLE_FRACTION = 0.0

#: The kernel's declared dry-face mass-force mode, the only one the vertex reduction is stated for.
VERTEX_MASS_FORCE_MODE = nip.NonisothermalMassForceMode.COMMON_FACE_TEMPERATURE_NEGLIGIBLE_SORET


class PureWaterVertexError(cp.CoupledPoreTopologyError):
    """A state, face or configuration the pure-water vertex does not carry (refused, never clipped)."""


@dataclass(frozen=True)
class PureWaterPoreState:
    """One locally equilibrated dry-pore state at the water vertex, ``y_hexane = 0`` exactly.

    The field names are those of :class:`coupled_pore.EquilibriumPoreState` for every quantity the
    vertex carries, so the frozen flux closures that read them (``retained_water_flux``) are called
    unchanged.  It is not that class: no binary validator is asked to accept the vertex.
    """

    temperature_k: float
    pressure_pa: float
    gas_accessible_fraction: float
    binary_gas: bg.BinaryVirialState
    water_activity: float
    retained_water_active_set: cp.RetainedWaterActiveSet
    retained_water_loading: float
    retained_water_fugacity_pa: float
    pore_water_concentration_kg_m3: float
    retained_water_concentration_kg_m3: float
    total_water_concentration_kg_m3: float
    total_water_concentration_mol_m3: float
    dry_meal_energy_density_j_m3: float
    retained_water_energy_density_j_m3: float
    retained_water_pv_sensitivity_j_m3: float
    retained_hexane_energy_density_j_m3: float
    pore_gas_energy_density_j_m3: float
    energy_density_j_m3: float
    water_gas_partial_enthalpy_j_mol: float
    hexane_gas_partial_enthalpy_j_mol: float
    retained_water_enthalpy_j_mol: float
    retained_water_potential_isothermal: float
    y_water: float = ADMITTED_WATER_MOLE_FRACTION
    y_hexane: float = ADMITTED_HEXANE_MOLE_FRACTION
    hexane_activity: float = 0.0
    retained_hexane_loading: float = 0.0
    pore_hexane_concentration_kg_m3: float = 0.0
    retained_hexane_concentration_kg_m3: float = 0.0
    total_hexane_concentration_kg_m3: float = 0.0
    total_hexane_concentration_mol_m3: float = 0.0
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.y_water != ADMITTED_WATER_MOLE_FRACTION or self.y_hexane != 0.0:
            raise PureWaterVertexError(
                "the pure-water vertex holds y_water = 1, y_hexane = 0 exactly"
            )
        hexane = (
            self.hexane_activity,
            self.retained_hexane_loading,
            self.pore_hexane_concentration_kg_m3,
            self.retained_hexane_concentration_kg_m3,
            self.total_hexane_concentration_kg_m3,
            self.total_hexane_concentration_mol_m3,
        )
        if any(value != 0.0 for value in hexane):
            raise PureWaterVertexError("the pure-water vertex carries no n-hexane")

    @property
    def has_pore_gas_composition_degree_of_freedom(self) -> bool:
        """False: at fixed total pressure a one-component gas has none."""

        return False


def evaluate_pure_water_state(
    temperature_k: float,
    pressure_pa: float,
    params: cp.CoupledPoreParams = cp.CoupledPoreParams(),
    *,
    gas_accessible_fraction: float = 1.0,
) -> PureWaterPoreState:
    """The binary local-equilibrium expressions (``coupled_pore._evaluate_equilibrium_cached``) at
    ``y_water = 1``, ``y_hexane = 0`` exactly, every frozen closure called unchanged.

    The binary validator's temperature, pressure and accessibility checks and its water-activity
    topology checks are applied with the binary messages; the composition check is replaced by the
    vertex itself (no composition is an input).
    """

    cp._validate_params(params)
    values = (temperature_k, pressure_pa, gas_accessible_fraction)
    if not all(math.isfinite(value) for value in values):
        raise PureWaterVertexError("pore-state primitives must be finite")
    lower, upper = params.temperature_bounds_k
    if not lower <= temperature_k <= upper:
        raise PureWaterVertexError(f"temperature {temperature_k} K is outside [{lower},{upper}] K")
    if pressure_pa <= 0.0:
        raise PureWaterVertexError("pore pressure must be positive")
    if not 0.0 <= gas_accessible_fraction <= 1.0:
        raise PureWaterVertexError("gas-accessible pore fraction must lie in [0,1]")
    y_water, y_hexane = ADMITTED_WATER_MOLE_FRACTION, ADMITTED_HEXANE_MOLE_FRACTION
    mixture = bg.state(temperature_k, pressure_pa, y_water, y_hexane, k_wh=params.k_wh)
    liquid_water = wa.state_Tp(temperature_k, pressure_pa, "liquid")
    water_activity = mixture.fugacity_water_pa / liquid_water.fugacity
    if not math.isfinite(water_activity) or not 0.0 < water_activity <= 1.0:
        raise PureWaterVertexError(
            f"dry-pore water activity must lie in (0,1], got {water_activity}; "
            "activate the external/free-water topology, do not clamp"
        )
    activity_at_ref = sp.water_activity(params.luikov.W_ref, params.luikov)
    activity_at_cap = sp.water_activity(params.luikov.W_cap, params.luikov)
    # The binary validator's own lower-authority tolerance expression (coupled_pore.py, the
    # Luikov lower-bound guard of _evaluate_equilibrium_cached), unchanged.
    lower_tolerance = 2.0e-13 * max(activity_at_ref, 1.0)
    if water_activity < activity_at_ref - lower_tolerance:
        raise PureWaterVertexError(
            f"water activity {water_activity} lies below the qualified Luikov "
            f"lower-bound activity {activity_at_ref}; no continuation is authorized"
        )
    if mixture.fugacity_hexane_pa != 0.0 or mixture.hexane_density_kg_m3 != 0.0:
        raise PureWaterVertexError("the frozen virial state at the water vertex carries n-hexane")
    if water_activity < activity_at_cap:
        try:
            retained_water = sp.water_retained(water_activity, params.luikov)
        except ValueError as exc:
            raise PureWaterVertexError(
                "water state left the smooth modified-Luikov branch"
            ) from exc
        water_active_set = cp.RetainedWaterActiveSet.LUIKOV_SMOOTH
        retained_water_fugacity = mixture.fugacity_water_pa
    else:
        retained_water = params.luikov.W_cap
        water_active_set = cp.RetainedWaterActiveSet.RETAINED_CAP
        retained_water_fugacity = activity_at_cap * liquid_water.fugacity
    # The frozen retained n-hexane law at activity 0: the vertex requires its zero.
    retained_hexane = sp.retained_hexane(0.0, temperature_k, params.w_o, params.gab, params.oil)
    if retained_hexane != 0.0:
        raise PureWaterVertexError("the frozen retained n-hexane law is not zero at activity 0")

    effective_eps = params.eps_g * gas_accessible_fraction
    pore_water = effective_eps * mixture.water_density_kg_m3
    retained_water_concentration = params.rho_dm_p * retained_water
    total_water = retained_water_concentration + pore_water

    dry_meal_energy = params.rho_dm_p * sp.composite_dry_meal_sensible_energy(
        temperature_k,
        params.T_ref_solid,
        params.w_o,
        params.cp_dry_meal,
        params.cp_oil,
        oil_law=params.oil_heat_capacity_law,
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
        activity=0.0,
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
    finite = (
        water_activity,
        retained_water,
        retained_water_fugacity,
        pore_water,
        total_water,
        dry_meal_energy,
        retained_water_energy,
        retained_hexane_energy,
        pore_gas_energy,
        energy,
        h_w_g,
        h_h_g,
        h_w_ret,
    )
    if not all(math.isfinite(value) for value in finite):
        raise PureWaterVertexError("pure-water vertex equilibrium produced a non-finite state")
    return PureWaterPoreState(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        gas_accessible_fraction=gas_accessible_fraction,
        binary_gas=mixture,
        water_activity=water_activity,
        retained_water_active_set=water_active_set,
        retained_water_loading=retained_water,
        retained_water_fugacity_pa=retained_water_fugacity,
        pore_water_concentration_kg_m3=pore_water,
        retained_water_concentration_kg_m3=retained_water_concentration,
        total_water_concentration_kg_m3=total_water,
        total_water_concentration_mol_m3=total_water / wa.M,
        dry_meal_energy_density_j_m3=dry_meal_energy,
        retained_water_energy_density_j_m3=retained_water_energy,
        retained_water_pv_sensitivity_j_m3=retained_water_pv,
        retained_hexane_energy_density_j_m3=retained_hexane_energy,
        pore_gas_energy_density_j_m3=pore_gas_energy,
        energy_density_j_m3=energy,
        water_gas_partial_enthalpy_j_mol=h_w_g,
        hexane_gas_partial_enthalpy_j_mol=h_h_g,
        retained_water_enthalpy_j_mol=h_w_ret,
        retained_water_potential_isothermal=math.log(
            retained_water_fugacity / cp.REFERENCE_FUGACITY_PA
        ),
    )


@dataclass(frozen=True)
class PureWaterFaceFlux:
    """The fluxes at one dry face between two vertex states (the dry face's component and energy
    objects, the frozen ``ComponentMolarFluxes`` and ``ComponentEnergyFlux``)."""

    component: cp.ComponentMolarFluxes
    energy: cp.ComponentEnergyFlux
    independent_water_flux_mol_m2_s: float
    retained_water_flux_mol_m2_s: float
    conductive_heat_flux_w_m2: float
    gas_entropy_production_w_m3_k: float
    retained_water_entropy_production_w_m3_k: float
    conduction_entropy_production_w_m3_k: float
    force_reconstruction: str
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def total_entropy_production_w_m3_k(self) -> float:
        return (
            self.gas_entropy_production_w_m3_k
            + self.retained_water_entropy_production_w_m3_k
            + self.conduction_entropy_production_w_m3_k
        )


def _require_vertex_mode(config) -> None:
    if config.mass_force_mode is not VERTEX_MASS_FORCE_MODE:
        raise PureWaterVertexError(
            f"the pure-water vertex is stated for {VERTEX_MASS_FORCE_MODE.value} only, "
            f"not {config.mass_force_mode.value}"
        )
    if (
        config.binary_thermal_diffusion_factor.alpha_t != 0.0
        or config.retained_water_thermal_force_factor.alpha_ret != 0.0
    ):
        raise PureWaterVertexError("the pure-water vertex is stated for zero thermal force factors")


def _upwind(flux: float, left: PureWaterPoreState, right: PureWaterPoreState) -> PureWaterPoreState:
    # cut_transport._upwind_state, unchanged in rule: outward flux takes the left state, zero the left.
    if flux > 0.0:
        return left
    if flux < 0.0:
        return right
    return left


def _face(left, right, left_t, right_t, distance_m, total_stefan_flux, params, config, label):
    retained = cp.retained_water_flux(left, right, distance_m, config.dry.retained_water_mobility)
    component = cp.compose_component_fluxes(
        ADMITTED_WATER_MOLE_FRACTION,
        ADMITTED_HEXANE_MOLE_FRACTION,
        total_stefan_flux,
        0.0,
        retained.flux_mol_m2_s,
    )
    fourier = config.dry.thermal_conductivity.value_w_m_k * (left_t - right_t) / distance_m
    water_state = _upwind(component.gas_water_flux_mol_m2_s, left, right)
    hexane_state = _upwind(component.gas_hexane_flux_mol_m2_s, left, right)
    retained_state = _upwind(component.retained_water_flux_mol_m2_s, left, right)
    energy = cp.component_energy_flux(
        component,
        fourier,
        water_state.water_gas_partial_enthalpy_j_mol,
        hexane_state.hexane_gas_partial_enthalpy_j_mol,
        retained_state.retained_water_enthalpy_j_mol,
    )
    gradient = (right_t - left_t) / distance_m
    conduction_entropy = (
        config.dry.thermal_conductivity.value_w_m_k * gradient * gradient / (left_t * right_t)
    )
    return PureWaterFaceFlux(
        component=component,
        energy=energy,
        independent_water_flux_mol_m2_s=0.0,
        retained_water_flux_mol_m2_s=retained.flux_mol_m2_s,
        conductive_heat_flux_w_m2=fourier,
        gas_entropy_production_w_m3_k=0.0,
        retained_water_entropy_production_w_m3_k=retained.entropy_production_w_m3_k,
        conduction_entropy_production_w_m3_k=conduction_entropy,
        force_reconstruction=label,
    )


def interior_face_flux(
    left_temperature_k: float,
    right_temperature_k: float,
    distance_m: float,
    total_stefan_flux_mol_m2_s: float,
    params: cp.CoupledPoreParams,
    config,
) -> PureWaterFaceFlux:
    """An interior dry face (``cut_transport._dry_flux`` without a force temperature): both endpoint
    states at the face temperature ``(T_L + T_R)/2``, the independent flux zero by absence, the
    retained flux from the two identical vertex potentials (exactly zero), steam carried by the total
    flux, Fourier heat between the endpoint temperatures."""

    if not math.isfinite(distance_m) or distance_m <= 0.0:
        raise PureWaterVertexError("dry face distance must be positive")
    if not math.isfinite(total_stefan_flux_mol_m2_s):
        raise ValueError("total Stefan flux must be finite")
    _require_vertex_mode(config)
    face_temperature = 0.5 * (left_temperature_k + right_temperature_k)
    state = evaluate_pure_water_state(face_temperature, config.dry.pressure_pa, params)
    return _face(
        state,
        state,
        left_temperature_k,
        right_temperature_k,
        distance_m,
        total_stefan_flux_mol_m2_s,
        params,
        config,
        "symmetric_endpoint_temperature_vertex",
    )


def surface_face_flux(
    cell_temperature_k: float,
    surface_temperature_k: float,
    distance_m: float,
    total_stefan_flux_mol_m2_s: float,
    params: cp.CoupledPoreParams,
    config,
):
    """The outer half-cell face (the kernel's actual linear (T, y) path with the force temperature at
    the cell): actual endpoint states at the cell and surface temperatures, the composition part of
    both forces exactly zero (equal endpoint compositions), the independent and retained fluxes zero,
    steam carried by the total flux with the upwind endpoint's enthalpy, Fourier heat.  Returns the
    face flux and the two endpoint states."""

    if not math.isfinite(distance_m) or distance_m <= 0.0:
        raise PureWaterVertexError("dry face distance must be positive")
    if not math.isfinite(total_stefan_flux_mol_m2_s):
        raise ValueError("total Stefan flux must be finite")
    _require_vertex_mode(config)
    left = evaluate_pure_water_state(cell_temperature_k, config.dry.pressure_pa, params)
    right = evaluate_pure_water_state(surface_temperature_k, config.dry.pressure_pa, params)
    # The actual path's composition-only forces: zero along a constant-composition path.  The
    # retained flux is assembled as the kernel's actual-path branch assembles it (mobility times
    # the force gradient), with the gradient zero.
    retained_gradient = 0.0
    mobility = config.dry.retained_water_mobility.value_mol_m_s
    retained_flux = -mobility * retained_gradient
    component = cp.compose_component_fluxes(
        ADMITTED_WATER_MOLE_FRACTION,
        ADMITTED_HEXANE_MOLE_FRACTION,
        total_stefan_flux_mol_m2_s,
        0.0,
        retained_flux,
    )
    fourier = (
        config.dry.thermal_conductivity.value_w_m_k
        * (cell_temperature_k - surface_temperature_k)
        / distance_m
    )
    water_state = _upwind(component.gas_water_flux_mol_m2_s, left, right)
    hexane_state = _upwind(component.gas_hexane_flux_mol_m2_s, left, right)
    retained_state = _upwind(component.retained_water_flux_mol_m2_s, left, right)
    energy = cp.component_energy_flux(
        component,
        fourier,
        water_state.water_gas_partial_enthalpy_j_mol,
        hexane_state.hexane_gas_partial_enthalpy_j_mol,
        retained_state.retained_water_enthalpy_j_mol,
    )
    gradient = (surface_temperature_k - cell_temperature_k) / distance_m
    conduction_entropy = (
        config.dry.thermal_conductivity.value_w_m_k
        * gradient
        * gradient
        / (cell_temperature_k * surface_temperature_k)
    )
    face = PureWaterFaceFlux(
        component=component,
        energy=energy,
        independent_water_flux_mol_m2_s=0.0,
        retained_water_flux_mol_m2_s=retained_flux,
        conductive_heat_flux_w_m2=fourier,
        gas_entropy_production_w_m3_k=0.0,
        retained_water_entropy_production_w_m3_k=bg.R
        * mobility
        * retained_gradient
        * retained_gradient,
        conduction_entropy_production_w_m3_k=conduction_entropy,
        force_reconstruction="actual_linear_t_y_path_vertex",
    )
    return face, left, right


__all__ = [
    "ADMITTED_HEXANE_MOLE_FRACTION",
    "ADMITTED_WATER_MOLE_FRACTION",
    "PureWaterFaceFlux",
    "PureWaterPoreState",
    "PureWaterVertexError",
    "RECORD",
    "VERTEX_MASS_FORCE_MODE",
    "evaluate_pure_water_state",
    "interior_face_flux",
    "surface_face_flux",
]
