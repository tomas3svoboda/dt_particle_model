r"""Complete-potential counterfactual for the P1EF mass-force ablation.

The frozen G1G-01 pore laws use dimensionless thermodynamic forces

``grad((mu_w-mu_h)/(R*T))`` and ``grad(mu_w,ret/(R*T))``.

At one temperature their standard-state terms cancel between face endpoints,
so the existing log-fugacity implementation is exact.  Across a temperature
gradient they do not cancel.  This module constructs those missing terms from
the same IAPWS-95, Span--Wagner, and PHY-053 authorities used by storage and
energy.  It adds no fitted heat of transport.

There is an unavoidable model-form qualification.  A diagonal non-isothermal
law ``J=-L*grad(mu/(R*T))`` is not invariant when an arbitrary constant is
added to a component enthalpy datum unless the conjugate heat flux / Onsager
cross coefficient is transformed with it.  No such Soret coefficient is
identified for soybean DT meal.  The counterfactual therefore fixes an
explicit reference gauge: each ideal-gas component has zero normalized
enthalpy and entropy at ``(T_ref,p_ref)``.  Raw EOS datum changes cancel from
this construction, but changing the declared reference gauge remains a
mandatory sensitivity, not a new physical parameter estimate.

The production-comparison mode is selected immutably by the cut-transport
configuration.  The existing common-face-temperature negligible-Soret branch
does not call this module and remains numerically unchanged.  A separate,
reference-invariant engineering ablation uses explicit thermal-force factors
and their reciprocal heat-of-transport contribution.  Its bounds are frozen
before trajectory evaluation: ``alpha_wh=+-0.5`` and ``alpha_ret=+-25`` are
representative qualification stresses, while ``alpha_ret=+-40`` is an
anomalous convergence/conservation stress only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Protocol

from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import water as wa


class NonisothermalMassForceMode(str, Enum):
    """Immutable model-form selection for particle mass forces."""

    COMMON_FACE_TEMPERATURE_NEGLIGIBLE_SORET = (
        "common_face_temperature_negligible_soret"
    )
    COMPLETE_POTENTIAL_REFERENCE_GAUGE_COUNTERFACTUAL = (
        "complete_potential_reference_gauge_counterfactual"
    )
    REFERENCE_INVARIANT_BINARY_THERMAL_DIFFUSION_FACTOR_COUNTERFACTUAL = (
        "reference_invariant_binary_thermal_diffusion_factor_counterfactual"
    )
    REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS_COUNTERFACTUAL = (
        "reference_invariant_thermal_force_factors_counterfactual"
    )


@dataclass(frozen=True)
class BinaryThermalDiffusionFactorSelection:
    r"""One explicit dimensionless ``alpha_T`` counterfactual selection.

    The associated binary-water flux increment is

    ``J_w^T = -L_wh * alpha_T * grad(ln(T))``.

    In the ideal Maxwell--Stefan limit ``L_wh=c*D*y_w*y_h``, this is exactly
    the classical ``-c*D*alpha_T*y_w*y_h*grad(ln(T))`` form.  It is component-
    reference invariant and conserves ``J_w^T=-J_h^T``.  This generic value
    object accepts manufactured selections; the P1EF engineering builder
    separately freezes its owner-authorized grid to ``{-0.5,0,+0.5}``.
    Every nonzero selection remains a non-identifying model-form stress until
    direct or validated kinetic-theory evidence for water-vapour/n-hexane
    exists.
    """

    alpha_t: float
    label: str
    direct_binary_evidence: bool = False
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.alpha_t):
            raise ValueError("binary thermal-diffusion factor must be finite")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("binary thermal-diffusion factor needs provenance")
        if not isinstance(self.direct_binary_evidence, bool):
            raise TypeError("direct-binary-evidence flag must be boolean")


ZERO_BINARY_THERMAL_DIFFUSION_FACTOR = BinaryThermalDiffusionFactorSelection(
    alpha_t=0.0,
    label=(
        "P1EF production negligible-Soret point alpha_T=0; no direct "
        "water-vapour/n-hexane thermal-diffusion identification"
    ),
    direct_binary_evidence=False,
)


@dataclass(frozen=True)
class RetainedWaterThermalForceFactorSelection:
    r"""One explicit retained-water thermal-force-factor selection.

    The selected retained-water flux is

    ``J_ret=-L_ret*(grad(phi_ret)+alpha_ret*grad(ln(T)))``.

    ``alpha_ret`` can be read as a dimensionless heat-of-transport stress,
    ``Q*_ret/(R*T)``, only inside this phenomenological ablation.  The animal-
    feed activation energy and isosteric heat used to size the brackets do not
    identify ``Q*_ret`` or its sign for desolventizing soybean meal.  Every
    nonzero value therefore remains non-production and non-identifying.
    """

    alpha_ret: float
    label: str
    directly_identified_heat_of_transport: bool = False
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.alpha_ret):
            raise ValueError("retained-water thermal-force factor must be finite")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("retained-water thermal-force factor needs provenance")
        if not isinstance(self.directly_identified_heat_of_transport, bool):
            raise TypeError("heat-of-transport identification flag must be boolean")


ZERO_RETAINED_WATER_THERMAL_FORCE_FACTOR = (
    RetainedWaterThermalForceFactorSelection(
        alpha_ret=0.0,
        label=(
            "P1EF production common-face-temperature point alpha_ret=0; no "
            "retained-water heat of transport is identified"
        ),
    )
)


# These grids are deliberately immutable and declared before any P1EF
# trajectory result.  The representative and anomalous retained-water sets
# must never be pooled into one acceptance interval.
P1EF_GAS_REPRESENTATIVE_ALPHA_T_VALUES = (-0.5, 0.0, 0.5)
P1EF_RETAINED_REPRESENTATIVE_ALPHA_VALUES = (-25.0, 0.0, 25.0)
P1EF_RETAINED_ANOMALOUS_STRESS_ALPHA_VALUES = (-40.0, 40.0)

TOUFFET_DIFFUSION_ACTIVATION_ENERGY_KJ_MOL = 57.17
TOUFFET_DIFFUSION_ACTIVATION_ENERGY_UNCERTAINTY_KJ_MOL = 5.02
P1EF_RETAINED_BRACKET_REFERENCE_TEMPERATURE_K = 323.15
TOUFFET_UPPER_ACTIVATION_ENERGY_OVER_RT = (
    (
        TOUFFET_DIFFUSION_ACTIVATION_ENERGY_KJ_MOL
        + TOUFFET_DIFFUSION_ACTIVATION_ENERGY_UNCERTAINTY_KJ_MOL
    )
    * 1000.0
    / (bg.R * P1EF_RETAINED_BRACKET_REFERENCE_TEMPERATURE_K)
)
TOUFFET_MAX_TOTAL_ISOSTERIC_HEAT_KJ_KG_WATER = 5500.0
TOUFFET_MAX_TOTAL_ISOSTERIC_HEAT_OVER_RT = (
    TOUFFET_MAX_TOTAL_ISOSTERIC_HEAT_KJ_KG_WATER
    * 1000.0
    * wa.M
    / (bg.R * P1EF_RETAINED_BRACKET_REFERENCE_TEMPERATURE_K)
)
HESS_ANALOG_MASS_CONTRAST = (hx.M - wa.M) / (hx.M + wa.M)
HESS_ANALOG_ALPHA_T_ESTIMATE = 0.38 * HESS_ANALOG_MASS_CONTRAST


def reference_invariant_thermal_factor_mode(
    mode: NonisothermalMassForceMode,
) -> bool:
    """Return whether ``mode`` selects explicit datum-free thermal factors."""

    if not isinstance(mode, NonisothermalMassForceMode):
        raise TypeError("thermal-factor mode query requires an immutable mode")
    return mode in (
        NonisothermalMassForceMode.REFERENCE_INVARIANT_BINARY_THERMAL_DIFFUSION_FACTOR_COUNTERFACTUAL,
        NonisothermalMassForceMode.REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS_COUNTERFACTUAL,
    )


@dataclass(frozen=True)
class CompletePotentialReferenceGauge:
    """Declared component gauge used only by the complete-potential ablation.

    Both component ideal-gas enthalpies and entropies are re-zeroed at this
    state.  This removes dependence on the arbitrary numerical zeros embedded
    in the two pure-fluid EOS implementations.  It does *not* identify the
    missing heat-of-transport coefficient; reference-temperature sensitivity
    must therefore remain visible in the ablation report.
    """

    reference_temperature_k: float
    reference_pressure_pa: float
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not (
            math.isfinite(self.reference_temperature_k)
            and self.reference_temperature_k > 0.0
        ):
            raise ValueError("complete-potential reference temperature must be positive")
        if not (
            math.isfinite(self.reference_pressure_pa)
            and self.reference_pressure_pa > 0.0
        ):
            raise ValueError("complete-potential reference pressure must be positive")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("complete-potential reference gauge needs a provenance label")


DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE = CompletePotentialReferenceGauge(
    reference_temperature_k=hx.T_REF,
    reference_pressure_pa=hx.P_REF,
    label=(
        "P1EF counterfactual gauge: each component ideal-gas h and s is zero "
        "at the Span-Wagner reporting reference (298.15 K, 101325 Pa); this "
        "is not an identified soybean Soret heat of transport"
    ),
)


class _PorePotentialState(Protocol):
    temperature_k: float
    retained_water_fugacity_pa: float
    gas_exchange_potential_isothermal: float
    retained_water_potential_isothermal: float
    binary_gas: bg.BinaryVirialState


@dataclass(frozen=True)
class CompletePotentialState:
    """Gauge-fixed component potentials at one dry-pore equilibrium state."""

    water_gas_chemical_potential_over_rt: float
    hexane_gas_chemical_potential_over_rt: float
    gas_exchange_chemical_potential_over_rt: float
    retained_water_chemical_potential_over_rt: float
    normalized_water_ideal_enthalpy_j_mol: float
    normalized_hexane_ideal_enthalpy_j_mol: float
    reference_gauge: CompletePotentialReferenceGauge
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.water_gas_chemical_potential_over_rt,
            self.hexane_gas_chemical_potential_over_rt,
            self.gas_exchange_chemical_potential_over_rt,
            self.retained_water_chemical_potential_over_rt,
            self.normalized_water_ideal_enthalpy_j_mol,
            self.normalized_hexane_ideal_enthalpy_j_mol,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("complete chemical-potential state must be finite")


@dataclass(frozen=True)
class CompletePotentialForceDifference:
    """Cancellation-aware dimensionless force differences across one face."""

    gas_exchange_right_minus_left: float
    retained_water_right_minus_left: float
    used_exact_isothermal_recovery: bool
    reference_gauge: CompletePotentialReferenceGauge
    physically_qualifying: bool = field(default=False, init=False)


def water_ideal_entropy_molar(temperature_k: float, pressure_pa: float) -> float:
    """IAPWS-95 ideal-water entropy, J/(mol K), on its native datum.

    The public IAPWS state contains ideal plus residual Helmholtz terms.  At
    ``rho=p/(R*T)`` the ideal part corresponds exactly to the requested ideal
    pressure.  Subtracting the analytic residual entropy therefore exposes the
    same ideal authority without introducing another correlation or standard
    state.
    """

    _validate_temperature_pressure(temperature_k, pressure_pa)
    density_kg_m3 = pressure_pa / (wa.R * temperature_k)
    state = wa.state(temperature_k, density_kg_m3)
    residual = wa.residual(density_kg_m3 / wa.RHOC, wa.TC / temperature_k)
    residual_entropy_j_kg_k = wa.R * (
        residual.tau * residual.phir_t - residual.phir
    )
    return (state.s_mass - residual_entropy_j_kg_k) * wa.M


def normalized_water_ideal_standard_state(
    temperature_k: float,
    gauge: CompletePotentialReferenceGauge = (
        DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE
    ),
) -> tuple[float, float]:
    """Return ``(mu_w^0/(RT), h_w^0-h_w^0_ref)`` in the fixed gauge."""

    _validate_gauge_and_temperature(gauge, temperature_k)
    h = bg.water_ideal_enthalpy_molar(temperature_k)
    h_ref = bg.water_ideal_enthalpy_molar(gauge.reference_temperature_k)
    s = water_ideal_entropy_molar(temperature_k, gauge.reference_pressure_pa)
    s_ref = water_ideal_entropy_molar(
        gauge.reference_temperature_k,
        gauge.reference_pressure_pa,
    )
    normalized_h = h - h_ref
    normalized_s = s - s_ref
    potential = normalized_h / (bg.R * temperature_k) - normalized_s / bg.R
    return potential, normalized_h


def normalized_hexane_ideal_standard_state(
    temperature_k: float,
    gauge: CompletePotentialReferenceGauge = (
        DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE
    ),
) -> tuple[float, float]:
    """Return ``(mu_h^0/(RT), h_h^0-h_h^0_ref)`` in the fixed gauge."""

    _validate_gauge_and_temperature(gauge, temperature_k)
    h = hx.h_ideal(temperature_k)
    h_ref = hx.h_ideal(gauge.reference_temperature_k)
    s = hx.s_ideal(temperature_k, gauge.reference_pressure_pa)
    s_ref = hx.s_ideal(
        gauge.reference_temperature_k,
        gauge.reference_pressure_pa,
    )
    normalized_h = h - h_ref
    normalized_s = s - s_ref
    potential = normalized_h / (bg.R * temperature_k) - normalized_s / bg.R
    return potential, normalized_h


def complete_pore_potential_state(
    state: _PorePotentialState,
    gauge: CompletePotentialReferenceGauge = (
        DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE
    ),
) -> CompletePotentialState:
    """Construct complete gas-exchange and retained-water potentials.

    ``ln(f_i/p_ref)`` uses the exact PHY-053 fugacity.  The retained potential
    uses the frozen active-set fugacity already selected by the pore state:
    actual gas water fugacity on the smooth Luikov branch and cap activity
    times IAPWS liquid fugacity on the retained-cap branch.
    """

    water_standard, water_h = normalized_water_ideal_standard_state(
        state.temperature_k,
        gauge,
    )
    hexane_standard, hexane_h = normalized_hexane_ideal_standard_state(
        state.temperature_k,
        gauge,
    )
    mixture = state.binary_gas
    water = water_standard + math.log(
        mixture.fugacity_water_pa / gauge.reference_pressure_pa
    )
    hexane = hexane_standard + math.log(
        mixture.fugacity_hexane_pa / gauge.reference_pressure_pa
    )
    retained = water_standard + math.log(
        state.retained_water_fugacity_pa / gauge.reference_pressure_pa
    )
    return CompletePotentialState(
        water_gas_chemical_potential_over_rt=water,
        hexane_gas_chemical_potential_over_rt=hexane,
        gas_exchange_chemical_potential_over_rt=math.fsum((water, -hexane)),
        retained_water_chemical_potential_over_rt=retained,
        normalized_water_ideal_enthalpy_j_mol=water_h,
        normalized_hexane_ideal_enthalpy_j_mol=hexane_h,
        reference_gauge=gauge,
    )


def complete_retained_water_potential_over_rt(
    temperature_k: float,
    retained_water_fugacity_pa: float,
    gauge: CompletePotentialReferenceGauge = (
        DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE
    ),
) -> float:
    """Complete retained-water potential for wet or dry material."""

    if not math.isfinite(retained_water_fugacity_pa) or retained_water_fugacity_pa <= 0.0:
        raise ValueError("retained-water fugacity must be positive and finite")
    standard, _ = normalized_water_ideal_standard_state(temperature_k, gauge)
    return standard + math.log(
        retained_water_fugacity_pa / gauge.reference_pressure_pa
    )


def complete_potential_force_difference(
    left: _PorePotentialState,
    right: _PorePotentialState,
    gauge: CompletePotentialReferenceGauge = (
        DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE
    ),
) -> CompletePotentialForceDifference:
    """Return complete right-minus-left forces with exact isothermal recovery.

    When endpoint temperatures are bit equal, the standard-state terms cancel
    identically.  Reusing the existing stored isothermal differences makes the
    recovery bit exact and avoids subtracting two large standard potentials.
    """

    if left.temperature_k == right.temperature_k:
        return CompletePotentialForceDifference(
            gas_exchange_right_minus_left=(
                right.gas_exchange_potential_isothermal
                - left.gas_exchange_potential_isothermal
            ),
            retained_water_right_minus_left=(
                right.retained_water_potential_isothermal
                - left.retained_water_potential_isothermal
            ),
            used_exact_isothermal_recovery=True,
            reference_gauge=gauge,
        )
    left_complete = complete_pore_potential_state(left, gauge)
    right_complete = complete_pore_potential_state(right, gauge)
    return CompletePotentialForceDifference(
        gas_exchange_right_minus_left=math.fsum(
            (
                right_complete.gas_exchange_chemical_potential_over_rt,
                -left_complete.gas_exchange_chemical_potential_over_rt,
            )
        ),
        retained_water_right_minus_left=math.fsum(
            (
                right_complete.retained_water_chemical_potential_over_rt,
                -left_complete.retained_water_chemical_potential_over_rt,
            )
        ),
        used_exact_isothermal_recovery=False,
        reference_gauge=gauge,
    )


def gauge_shifted_potential_over_rt(
    potential_over_rt: float,
    temperature_k: float,
    *,
    enthalpy_shift_j_mol: float,
    entropy_shift_j_mol_k: float,
) -> float:
    """Apply ``mu' = mu + C_h - T*C_s`` for covariance audits."""

    values = (
        potential_over_rt,
        temperature_k,
        enthalpy_shift_j_mol,
        entropy_shift_j_mol_k,
    )
    if not all(math.isfinite(value) for value in values) or temperature_k <= 0.0:
        raise ValueError("gauge-shift inputs must be finite and T must be positive")
    return math.fsum(
        (
            potential_over_rt,
            enthalpy_shift_j_mol / (bg.R * temperature_k),
            -entropy_shift_j_mol_k / bg.R,
        )
    )


def binary_thermal_force_gradient_m_inv(
    left_temperature_k: float,
    right_temperature_k: float,
    distance_m: float,
    selection: BinaryThermalDiffusionFactorSelection,
) -> float:
    r"""Return ``alpha_T*grad(ln(T))`` without a component reference datum."""

    if not isinstance(selection, BinaryThermalDiffusionFactorSelection):
        raise TypeError("binary thermal force requires an explicit alpha_T selection")
    return selection.alpha_t * log_temperature_gradient_m_inv(
        left_temperature_k,
        right_temperature_k,
        distance_m,
    )


def retained_water_thermal_force_gradient_m_inv(
    left_temperature_k: float,
    right_temperature_k: float,
    distance_m: float,
    selection: RetainedWaterThermalForceFactorSelection,
) -> float:
    r"""Return ``alpha_ret*grad(ln(T))`` without a component datum."""

    if not isinstance(selection, RetainedWaterThermalForceFactorSelection):
        raise TypeError(
            "retained-water thermal force requires an explicit alpha_ret selection"
        )
    return selection.alpha_ret * log_temperature_gradient_m_inv(
        left_temperature_k,
        right_temperature_k,
        distance_m,
    )


def log_temperature_gradient_m_inv(
    left_temperature_k: float,
    right_temperature_k: float,
    distance_m: float,
) -> float:
    """Return the exact endpoint secant ``grad(ln(T))``."""

    values = (left_temperature_k, right_temperature_k, distance_m)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("log-temperature-gradient inputs must be finite")
    if left_temperature_k <= 0.0 or right_temperature_k <= 0.0 or distance_m <= 0.0:
        raise ValueError("log-temperature gradient requires positive T and distance")
    return math.log1p(
        (right_temperature_k - left_temperature_k) / left_temperature_k
    ) / distance_m


def discrete_conjugate_temperature_k(
    left_temperature_k: float,
    right_temperature_k: float,
) -> float:
    r"""Return the face temperature conjugate to endpoint ``grad(1/T)``.

    For unequal endpoints this is

    ``T_c=T_l*T_r*ln(T_r/T_l)/(T_r-T_l)``.

    Consequently ``T_c*grad(1/T)=-grad(ln(T))`` exactly up to roundoff.
    The equality branch returns the common endpoint exactly.
    """

    if not all(
        math.isfinite(value) for value in (left_temperature_k, right_temperature_k)
    ):
        raise ValueError("conjugate-temperature endpoints must be finite")
    if left_temperature_k <= 0.0 or right_temperature_k <= 0.0:
        raise ValueError("conjugate-temperature endpoints must be positive")
    if left_temperature_k == right_temperature_k:
        return left_temperature_k
    difference = right_temperature_k - left_temperature_k
    return (
        left_temperature_k
        * right_temperature_k
        * math.log1p(difference / left_temperature_k)
        / difference
    )


def reciprocal_heat_of_transport_flux_w_m2(
    *,
    left_temperature_k: float,
    right_temperature_k: float,
    gas_counterflux_water_mol_m2_s: float,
    binary_thermal_diffusion_factor: BinaryThermalDiffusionFactorSelection,
    retained_water_flux_mol_m2_s: float,
    retained_water_thermal_factor: RetainedWaterThermalForceFactorSelection,
) -> float:
    r"""Return the reciprocal reduced heat flux for the selected mass laws.

    With ``g_T=grad(ln(T))`` and
    ``J_i=-L_i(g_i+alpha_i*g_T)``, the reduced heat flux is represented as

    ``q=q_F + R*T_c*sum(alpha_i*J_i)``.

    The finite-volume face keeps its existing endpoint Fourier law
    ``q_F=-k*(T_r-T_l)/d`` (whose continuum limit is ``-k*T*g_T``).

    ``T_c`` is the exact discrete conjugate temperature above.  Therefore the
    endpoint entropy-production bilinear form reduces to

    ``R*sum(L_i*(g_i+alpha_i*g_T)^2) + q_F*grad(1/T)``.

    This zero-mass-flux conductivity representation is positive for every
    finite ``alpha_i`` when ``k`` and the mobilities are positive; no clipping
    or postulated Schur bound is required.  Numerical/property-domain failures
    at stress endpoints must still reject the state unchanged.
    """

    if not isinstance(
        binary_thermal_diffusion_factor,
        BinaryThermalDiffusionFactorSelection,
    ):
        raise TypeError("reciprocal heat flux requires an explicit binary alpha_T")
    if not isinstance(
        retained_water_thermal_factor,
        RetainedWaterThermalForceFactorSelection,
    ):
        raise TypeError("reciprocal heat flux requires an explicit retained alpha")
    fluxes = (gas_counterflux_water_mol_m2_s, retained_water_flux_mol_m2_s)
    if not all(math.isfinite(value) for value in fluxes):
        raise ValueError("reciprocal heat-of-transport mass fluxes must be finite")
    conjugate_temperature = discrete_conjugate_temperature_k(
        left_temperature_k,
        right_temperature_k,
    )
    weighted_flux = math.fsum(
        (
            binary_thermal_diffusion_factor.alpha_t
            * gas_counterflux_water_mol_m2_s,
            retained_water_thermal_factor.alpha_ret
            * retained_water_flux_mol_m2_s,
        )
    )
    if weighted_flux == 0.0:
        return 0.0
    return bg.R * conjugate_temperature * weighted_flux


def gauge_covariant_entropy_production_w_m3_k(
    *,
    left_temperature_k: float,
    right_temperature_k: float,
    distance_m: float,
    total_energy_flux_w_m2: float,
    gas_water_flux_mol_m2_s: float,
    gas_hexane_flux_mol_m2_s: float,
    retained_water_flux_mol_m2_s: float,
    left_water_gas_potential_over_rt: float,
    right_water_gas_potential_over_rt: float,
    left_hexane_gas_potential_over_rt: float,
    right_hexane_gas_potential_over_rt: float,
    left_retained_water_potential_over_rt: float,
    right_retained_water_potential_over_rt: float,
) -> float:
    r"""Return the discrete entropy-production bilinear form.

    This is ``J_E grad(1/T) - sum_i N_i grad(mu_i/T)``.  It is algebraically
    invariant when component chemical potentials and the common-datum energy
    flux are transformed together.  Its covariance is distinct from (and
    does not identify) a positive-definite full heat/mass Onsager matrix.
    """

    values = (
        left_temperature_k,
        right_temperature_k,
        distance_m,
        total_energy_flux_w_m2,
        gas_water_flux_mol_m2_s,
        gas_hexane_flux_mol_m2_s,
        retained_water_flux_mol_m2_s,
        left_water_gas_potential_over_rt,
        right_water_gas_potential_over_rt,
        left_hexane_gas_potential_over_rt,
        right_hexane_gas_potential_over_rt,
        left_retained_water_potential_over_rt,
        right_retained_water_potential_over_rt,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("entropy-production audit inputs must be finite")
    if left_temperature_k <= 0.0 or right_temperature_k <= 0.0 or distance_m <= 0.0:
        raise ValueError("entropy-production audit requires positive T and distance")
    inverse_temperature_gradient = (
        1.0 / right_temperature_k - 1.0 / left_temperature_k
    ) / distance_m
    water_force = (
        right_water_gas_potential_over_rt - left_water_gas_potential_over_rt
    ) / distance_m
    hexane_force = (
        right_hexane_gas_potential_over_rt - left_hexane_gas_potential_over_rt
    ) / distance_m
    retained_force = (
        right_retained_water_potential_over_rt
        - left_retained_water_potential_over_rt
    ) / distance_m
    return math.fsum(
        (
            total_energy_flux_w_m2 * inverse_temperature_gradient,
            -bg.R * gas_water_flux_mol_m2_s * water_force,
            -bg.R * gas_hexane_flux_mol_m2_s * hexane_force,
            -bg.R * retained_water_flux_mol_m2_s * retained_force,
        )
    )


def _validate_gauge_and_temperature(
    gauge: CompletePotentialReferenceGauge,
    temperature_k: float,
) -> None:
    if not isinstance(gauge, CompletePotentialReferenceGauge):
        raise TypeError("complete-potential force requires an explicit reference gauge")
    if not math.isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("complete-potential temperature must be positive and finite")


def _validate_temperature_pressure(temperature_k: float, pressure_pa: float) -> None:
    if not all(math.isfinite(value) for value in (temperature_k, pressure_pa)):
        raise ValueError("ideal-water standard-state inputs must be finite")
    if temperature_k <= 0.0 or pressure_pa <= 0.0:
        raise ValueError("ideal-water standard state requires positive T and p")


__all__ = [
    "BinaryThermalDiffusionFactorSelection",
    "CompletePotentialForceDifference",
    "CompletePotentialReferenceGauge",
    "CompletePotentialState",
    "DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE",
    "HESS_ANALOG_ALPHA_T_ESTIMATE",
    "HESS_ANALOG_MASS_CONTRAST",
    "NonisothermalMassForceMode",
    "P1EF_GAS_REPRESENTATIVE_ALPHA_T_VALUES",
    "P1EF_RETAINED_ANOMALOUS_STRESS_ALPHA_VALUES",
    "P1EF_RETAINED_BRACKET_REFERENCE_TEMPERATURE_K",
    "P1EF_RETAINED_REPRESENTATIVE_ALPHA_VALUES",
    "RetainedWaterThermalForceFactorSelection",
    "TOUFFET_DIFFUSION_ACTIVATION_ENERGY_KJ_MOL",
    "TOUFFET_DIFFUSION_ACTIVATION_ENERGY_UNCERTAINTY_KJ_MOL",
    "TOUFFET_MAX_TOTAL_ISOSTERIC_HEAT_KJ_KG_WATER",
    "TOUFFET_MAX_TOTAL_ISOSTERIC_HEAT_OVER_RT",
    "TOUFFET_UPPER_ACTIVATION_ENERGY_OVER_RT",
    "ZERO_BINARY_THERMAL_DIFFUSION_FACTOR",
    "ZERO_RETAINED_WATER_THERMAL_FORCE_FACTOR",
    "binary_thermal_force_gradient_m_inv",
    "complete_pore_potential_state",
    "complete_potential_force_difference",
    "complete_retained_water_potential_over_rt",
    "discrete_conjugate_temperature_k",
    "gauge_covariant_entropy_production_w_m3_k",
    "gauge_shifted_potential_over_rt",
    "log_temperature_gradient_m_inv",
    "normalized_hexane_ideal_standard_state",
    "normalized_water_ideal_standard_state",
    "reciprocal_heat_of_transport_flux_w_m2",
    "reference_invariant_thermal_factor_mode",
    "retained_water_thermal_force_gradient_m_inv",
    "water_ideal_entropy_molar",
]
