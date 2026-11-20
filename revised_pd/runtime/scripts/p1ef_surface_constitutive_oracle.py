"""Independent P1EF surface Maxwell--Stefan constitutive evidence.

This module is an evidence oracle, not a transport implementation.  It consumes
only a self-contained mapping of serialized primitive values, reevaluates the
frozen PHY-053 binary-gas authority, and independently reconstructs the surface
binary mobility, force, n-hexane counterflux, and entropy production.

The bounded Phase-1 Particle Engineering Foundation (P1EF) remains
nonqualifying and non-predictive.  Passing this oracle changes no balance,
coefficient, state, event, solver tolerance, or normative F2/F3 status.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
import math
from typing import Any

from dtdc_simulator.core2.props import binary_gas as bg


AMENDMENT_ID = "GT-PS-2-P1E-15"
RAW_SAMPLE_SCHEMA_VERSION = 1
RAW_SAMPLE_KIND = "p1ef_surface_constitutive_primitives_v1"
ORACLE_OUTPUT_SCHEMA_VERSION = 1
COMPARISON_ULPS = 8

COMMON_FACE_NEGLIGIBLE_SORET = "common_face_temperature_negligible_soret"
REFERENCE_INVARIANT_BINARY_THERMAL_FACTOR = (
    "reference_invariant_binary_thermal_diffusion_factor_counterfactual"
)
REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS = (
    "reference_invariant_thermal_force_factors_counterfactual"
)
SUPPORTED_MASS_FORCE_MODES = (
    COMMON_FACE_NEGLIGIBLE_SORET,
    REFERENCE_INVARIANT_BINARY_THERMAL_FACTOR,
    REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS,
)

TOP_LEVEL_FIELDS = (
    "schema_version",
    "sample_kind",
    "outer_dry_state",
    "resolved_surface_state",
    "center_to_surface_distance_m",
    "pressure_pa",
    "binary_gas_interaction_k_wh",
    "binary_diffusivity",
    "mass_force",
    "production_values",
)
OUTER_DRY_STATE_FIELDS = (
    "temperature_k",
    "y_hexane",
    "center_radius_m",
)
RESOLVED_SURFACE_STATE_FIELDS = (
    "temperature_k",
    "y_hexane",
    "radius_m",
)
BINARY_DIFFUSIVITY_FIELDS = (
    "lower_m2_s",
    "upper_m2_s",
    "fraction",
    "selected_m2_s",
)
MASS_FORCE_FIELDS = (
    "mode",
    "binary_thermal_diffusion_factor",
)
PRODUCTION_VALUE_FIELDS = (
    "thermodynamic_force_temperature_k",
    "selected_binary_mobility_mol_m_s",
    "gas_chemical_force_gradient_m_inv",
    "binary_thermal_force_gradient_m_inv",
    "effective_binary_force_gradient_m_inv",
    "independent_hexane_flux_mol_m2_s",
    "binary_gas_entropy_production_w_m3_k",
)
HISTORY_SAMPLE_ROLES = (
    "baseline_end",
    "pulse_first",
    "pulse_end",
    "return_first",
)

_FORBIDDEN_SURROGATE_KEYS = frozenset(
    {
        "nh",
        "nt",
        "conservedhexaneflux",
        "conservedhexanefluxmolm2s",
        "surfaceconservedhexanefluxmolm2s",
        "totalstefanflux",
        "totalstefanfluxmolm2s",
        "surfacetotalstefanfluxmolm2s",
        "drytotalstefanfluxesmolm2s",
    }
)


class SurfaceConstitutiveOracleError(ValueError):
    """A primitive sample cannot support the frozen A15 evidence claim."""


@dataclass(frozen=True)
class SurfaceConstitutiveReconstruction:
    """Primitive-only reconstruction of one surface constitutive state."""

    common_face_temperature_k: float
    outer_water_fugacity_pa: float
    outer_hexane_fugacity_pa: float
    surface_water_fugacity_pa: float
    surface_hexane_fugacity_pa: float
    outer_molar_density_mol_m3: float
    surface_molar_density_mol_m3: float
    face_molar_density_mol_m3: float
    binary_chain_rule_mean: float
    selected_binary_diffusivity_m2_s: float
    selected_binary_mobility_mol_m_s: float
    gas_chemical_force_gradient_m_inv: float
    binary_thermal_force_gradient_m_inv: float
    effective_binary_force_gradient_m_inv: float
    independent_hexane_flux_mol_m2_s: float
    binary_gas_entropy_production_w_m3_k: float


@dataclass(frozen=True)
class _ParsedSample:
    outer_temperature_k: float
    outer_y_hexane: float
    outer_center_radius_m: float
    surface_temperature_k: float
    surface_y_hexane: float
    surface_radius_m: float
    distance_m: float
    pressure_pa: float
    k_wh: float
    diffusivity_lower_m2_s: float
    diffusivity_upper_m2_s: float
    diffusivity_fraction: float
    serialized_selected_diffusivity_m2_s: float
    mass_force_mode: str
    binary_thermal_diffusion_factor: float
    production_values: Mapping[str, Any]


def raw_sample_schema() -> dict[str, Any]:
    """Return the exact schema a backend must serialize for every accepted root."""

    return {
        "schema_version": RAW_SAMPLE_SCHEMA_VERSION,
        "sample_kind": RAW_SAMPLE_KIND,
        "required_fields": {
            "top_level": list(TOP_LEVEL_FIELDS),
            "outer_dry_state": list(OUTER_DRY_STATE_FIELDS),
            "resolved_surface_state": list(RESOLVED_SURFACE_STATE_FIELDS),
            "binary_diffusivity": list(BINARY_DIFFUSIVITY_FIELDS),
            "mass_force": list(MASS_FORCE_FIELDS),
            "production_values": list(PRODUCTION_VALUE_FIELDS),
        },
        "supported_mass_force_modes": list(SUPPORTED_MASS_FORCE_MODES),
        "comparison_ulps": COMPARISON_ULPS,
        "forbidden_flux_substitutes": [
            "N_h or any conserved n-hexane component flux",
            "N_t or any total Stefan flux",
        ],
        "claim_boundary": {
            "engineering_foundation_only": True,
            "physically_qualifying": False,
            "plant_predictive": False,
            "normative_F2_solvent_particle": False,
            "normative_F3_water_interface": False,
        },
    }


def _normalized_key(value: str) -> str:
    return "".join(character for character in value.lower() if character.isalnum())


def _reject_flux_surrogates(value: object, path: str = "sample") -> None:
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if not isinstance(key, str):
                raise SurfaceConstitutiveOracleError(f"{path} contains a non-string key")
            if _normalized_key(key) in _FORBIDDEN_SURROGATE_KEYS:
                raise SurfaceConstitutiveOracleError(
                    f"{path}.{key} is a forbidden N_h/N_t flux surrogate; "
                    "A15 requires independent J_h"
                )
            _reject_flux_surrogates(nested, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, nested in enumerate(value):
            _reject_flux_surrogates(nested, f"{path}[{index}]")


def _strict_mapping(
    value: object,
    expected_fields: tuple[str, ...],
    label: str,
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SurfaceConstitutiveOracleError(f"{label} must be an object")
    actual = set(value)
    expected = set(expected_fields)
    if actual != expected:
        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected, key=str)
        raise SurfaceConstitutiveOracleError(
            f"{label} fields are incomplete or unexpected; "
            f"missing={missing}, unexpected={unexpected}"
        )
    if not all(isinstance(key, str) for key in value):
        raise SurfaceConstitutiveOracleError(f"{label} keys must be strings")
    return value


def _finite_number(value: object, label: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
    ):
        raise SurfaceConstitutiveOracleError(f"{label} must be a finite number")
    return float(value)


def _positive_number(value: object, label: str) -> float:
    result = _finite_number(value, label)
    if result <= 0.0:
        raise SurfaceConstitutiveOracleError(f"{label} must be strictly positive")
    return result


def _interior_composition(value: object, label: str) -> float:
    result = _finite_number(value, label)
    if not 0.0 < result < 1.0:
        raise SurfaceConstitutiveOracleError(
            f"{label} must lie strictly inside the binary composition interval"
        )
    return result


def _parse_sample(value: object) -> _ParsedSample:
    _reject_flux_surrogates(value)
    sample = _strict_mapping(value, TOP_LEVEL_FIELDS, "sample")
    schema_version = sample["schema_version"]
    if (
        isinstance(schema_version, bool)
        or not isinstance(schema_version, int)
        or schema_version != RAW_SAMPLE_SCHEMA_VERSION
    ):
        raise SurfaceConstitutiveOracleError(
            f"sample schema must equal {RAW_SAMPLE_SCHEMA_VERSION}"
        )
    if sample["sample_kind"] != RAW_SAMPLE_KIND:
        raise SurfaceConstitutiveOracleError("sample kind is not the frozen A15 primitive kind")

    outer = _strict_mapping(
        sample["outer_dry_state"],
        OUTER_DRY_STATE_FIELDS,
        "outer_dry_state",
    )
    surface = _strict_mapping(
        sample["resolved_surface_state"],
        RESOLVED_SURFACE_STATE_FIELDS,
        "resolved_surface_state",
    )
    diffusivity = _strict_mapping(
        sample["binary_diffusivity"],
        BINARY_DIFFUSIVITY_FIELDS,
        "binary_diffusivity",
    )
    mass_force = _strict_mapping(
        sample["mass_force"],
        MASS_FORCE_FIELDS,
        "mass_force",
    )
    production = _strict_mapping(
        sample["production_values"],
        PRODUCTION_VALUE_FIELDS,
        "production_values",
    )

    outer_temperature = _positive_number(
        outer["temperature_k"],
        "outer_dry_state.temperature_k",
    )
    surface_temperature = _positive_number(
        surface["temperature_k"],
        "resolved_surface_state.temperature_k",
    )
    outer_y_hexane = _interior_composition(
        outer["y_hexane"],
        "outer_dry_state.y_hexane",
    )
    surface_y_hexane = _interior_composition(
        surface["y_hexane"],
        "resolved_surface_state.y_hexane",
    )
    outer_center = _finite_number(
        outer["center_radius_m"],
        "outer_dry_state.center_radius_m",
    )
    surface_radius = _positive_number(
        surface["radius_m"],
        "resolved_surface_state.radius_m",
    )
    distance = _positive_number(
        sample["center_to_surface_distance_m"],
        "center_to_surface_distance_m",
    )
    if outer_center < 0.0 or outer_center >= surface_radius:
        raise SurfaceConstitutiveOracleError(
            "outer dry center radius must lie inside the resolved particle surface"
        )

    lower = _positive_number(diffusivity["lower_m2_s"], "binary_diffusivity.lower_m2_s")
    upper = _positive_number(diffusivity["upper_m2_s"], "binary_diffusivity.upper_m2_s")
    if lower > upper:
        raise SurfaceConstitutiveOracleError("binary diffusivity endpoints are unordered")
    fraction = _finite_number(diffusivity["fraction"], "binary_diffusivity.fraction")
    if not 0.0 <= fraction <= 1.0:
        raise SurfaceConstitutiveOracleError("binary diffusivity fraction must lie in [0, 1]")
    selected = _positive_number(
        diffusivity["selected_m2_s"],
        "binary_diffusivity.selected_m2_s",
    )

    mode = mass_force["mode"]
    if not isinstance(mode, str) or mode not in SUPPORTED_MASS_FORCE_MODES:
        raise SurfaceConstitutiveOracleError(
            "mass-force mode is unsupported by the reference-invariant A15 oracle"
        )
    alpha_t = _finite_number(
        mass_force["binary_thermal_diffusion_factor"],
        "mass_force.binary_thermal_diffusion_factor",
    )
    if mode == COMMON_FACE_NEGLIGIBLE_SORET and alpha_t != 0.0:
        raise SurfaceConstitutiveOracleError(
            "negligible-Soret production mode requires an exact zero thermal factor"
        )

    for name in PRODUCTION_VALUE_FIELDS:
        _finite_number(production[name], f"production_values.{name}")
    if float(production["selected_binary_mobility_mol_m_s"]) < 0.0:
        raise SurfaceConstitutiveOracleError("production binary mobility is negative")
    if float(production["binary_gas_entropy_production_w_m3_k"]) < 0.0:
        raise SurfaceConstitutiveOracleError("production binary entropy is negative")

    return _ParsedSample(
        outer_temperature_k=outer_temperature,
        outer_y_hexane=outer_y_hexane,
        outer_center_radius_m=outer_center,
        surface_temperature_k=surface_temperature,
        surface_y_hexane=surface_y_hexane,
        surface_radius_m=surface_radius,
        distance_m=distance,
        pressure_pa=_positive_number(sample["pressure_pa"], "pressure_pa"),
        k_wh=_finite_number(
            sample["binary_gas_interaction_k_wh"],
            "binary_gas_interaction_k_wh",
        ),
        diffusivity_lower_m2_s=lower,
        diffusivity_upper_m2_s=upper,
        diffusivity_fraction=fraction,
        serialized_selected_diffusivity_m2_s=selected,
        mass_force_mode=mode,
        binary_thermal_diffusion_factor=alpha_t,
        production_values=production,
    )


def _binary_chain_rule_mean(left_y_water: float, right_y_water: float) -> float:
    """Independent entropy-variable mean; no production transport helper is used."""

    if left_y_water == right_y_water:
        return left_y_water * (1.0 - left_y_water)
    if (
        left_y_water == 0.0
        or left_y_water == 1.0
        or right_y_water == 0.0
        or right_y_water == 1.0
    ):
        return 0.0
    delta_y = right_y_water - left_y_water
    delta_logit = math.log1p(delta_y / left_y_water) - math.log1p(
        -delta_y / (1.0 - left_y_water)
    )
    result = delta_y / delta_logit
    if not math.isfinite(result) or not 0.0 < result <= 0.25:
        raise SurfaceConstitutiveOracleError(
            "independent binary chain-rule mean left its analytic interior bounds"
        )
    return result


def reconstruct_surface_constitutive_sample(
    value: object,
) -> SurfaceConstitutiveReconstruction:
    """Reconstruct PHY-053/L/G/J/sigma solely from one raw primitive sample."""

    sample = _parse_sample(value)
    face_temperature = 0.5 * (
        sample.outer_temperature_k + sample.surface_temperature_k
    )
    outer = bg.state(
        face_temperature,
        sample.pressure_pa,
        1.0 - sample.outer_y_hexane,
        sample.outer_y_hexane,
        k_wh=sample.k_wh,
    )
    surface = bg.state(
        face_temperature,
        sample.pressure_pa,
        1.0 - sample.surface_y_hexane,
        sample.surface_y_hexane,
        k_wh=sample.k_wh,
    )
    fugacities = (
        outer.fugacity_water_pa,
        outer.fugacity_hexane_pa,
        surface.fugacity_water_pa,
        surface.fugacity_hexane_pa,
    )
    if not all(math.isfinite(item) and item > 0.0 for item in fugacities):
        raise SurfaceConstitutiveOracleError("PHY-053 produced a nonpositive fugacity")

    left_y_water = 1.0 - sample.outer_y_hexane
    right_y_water = 1.0 - sample.surface_y_hexane
    chain_mean = _binary_chain_rule_mean(left_y_water, right_y_water)
    face_density = 0.5 * outer.molar_density_mol_m3 + 0.5 * (
        surface.molar_density_mol_m3
    )
    selected_diffusivity = sample.diffusivity_lower_m2_s + sample.diffusivity_fraction * (
        sample.diffusivity_upper_m2_s - sample.diffusivity_lower_m2_s
    )
    mobility = face_density * selected_diffusivity * chain_mean
    outer_force = math.log(outer.fugacity_water_pa / outer.fugacity_hexane_pa)
    surface_force = math.log(
        surface.fugacity_water_pa / surface.fugacity_hexane_pa
    )
    chemical_gradient = (surface_force - outer_force) / sample.distance_m
    thermal_gradient = 0.0
    if sample.mass_force_mode in (
        REFERENCE_INVARIANT_BINARY_THERMAL_FACTOR,
        REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS,
    ):
        log_temperature_gradient = math.log1p(
            (sample.surface_temperature_k - sample.outer_temperature_k)
            / sample.outer_temperature_k
        ) / sample.distance_m
        thermal_gradient = (
            sample.binary_thermal_diffusion_factor * log_temperature_gradient
        )
    effective_gradient = math.fsum((chemical_gradient, thermal_gradient))
    j_hexane = mobility * effective_gradient
    entropy = bg.R * mobility * effective_gradient * effective_gradient
    reconstructed = (
        face_temperature,
        face_density,
        chain_mean,
        selected_diffusivity,
        mobility,
        chemical_gradient,
        thermal_gradient,
        effective_gradient,
        j_hexane,
        entropy,
    )
    if not all(math.isfinite(item) for item in reconstructed):
        raise SurfaceConstitutiveOracleError(
            "surface constitutive reconstruction produced a nonfinite value"
        )
    if mobility <= 0.0 or entropy < 0.0:
        raise SurfaceConstitutiveOracleError(
            "surface constitutive reconstruction lost positive mobility/entropy"
        )

    return SurfaceConstitutiveReconstruction(
        common_face_temperature_k=face_temperature,
        outer_water_fugacity_pa=outer.fugacity_water_pa,
        outer_hexane_fugacity_pa=outer.fugacity_hexane_pa,
        surface_water_fugacity_pa=surface.fugacity_water_pa,
        surface_hexane_fugacity_pa=surface.fugacity_hexane_pa,
        outer_molar_density_mol_m3=outer.molar_density_mol_m3,
        surface_molar_density_mol_m3=surface.molar_density_mol_m3,
        face_molar_density_mol_m3=face_density,
        binary_chain_rule_mean=chain_mean,
        selected_binary_diffusivity_m2_s=selected_diffusivity,
        selected_binary_mobility_mol_m_s=mobility,
        gas_chemical_force_gradient_m_inv=chemical_gradient,
        binary_thermal_force_gradient_m_inv=thermal_gradient,
        effective_binary_force_gradient_m_inv=effective_gradient,
        independent_hexane_flux_mol_m2_s=j_hexane,
        binary_gas_entropy_production_w_m3_k=entropy,
    )


def _eight_ulp_comparison(serialized: float, reconstructed: float) -> dict[str, Any]:
    def binary64_ulp(value: float) -> float:
        """Return the IEEE-754 binary64 ULP without requiring Python 3.9."""

        magnitude = abs(value)
        if magnitude == 0.0 or magnitude < float.fromhex("0x1.0p-1022"):
            return float.fromhex("0x0.0000000000001p-1022")
        _fraction, exponent = math.frexp(magnitude)
        return math.ldexp(1.0, exponent - 53)

    bound = COMPARISON_ULPS * max(
        binary64_ulp(serialized),
        binary64_ulp(reconstructed),
    )
    try:
        residual = math.fsum((serialized, -reconstructed))
    except OverflowError:
        return {
            "serialized": serialized,
            "reconstructed": reconstructed,
            "residual": None,
            "absolute_residual": None,
            "eight_ulp_roundoff_bound": bound,
            "comparison_overflowed": True,
            "within_eight_ulps": False,
        }
    return {
        "serialized": serialized,
        "reconstructed": reconstructed,
        "residual": residual,
        "absolute_residual": abs(residual),
        "eight_ulp_roundoff_bound": bound,
        "comparison_overflowed": False,
        "within_eight_ulps": abs(residual) <= bound,
    }


def _failure_audit(errors: list[str]) -> dict[str, Any]:
    return {
        "schema_version": ORACLE_OUTPUT_SCHEMA_VERSION,
        "amendment_id": AMENDMENT_ID,
        "passed": False,
        "raw_sample_schema_version": RAW_SAMPLE_SCHEMA_VERSION,
        "comparison_ulps": COMPARISON_ULPS,
        "primitive_only_reconstruction": True,
        "N_h_or_N_t_substitution_allowed": False,
        "physics_or_numerical_equations_changed": False,
        "engineering_foundation_only": True,
        "physically_qualifying": False,
        "plant_predictive": False,
        "normative_F2_solvent_particle": False,
        "normative_F3_water_interface": False,
        "reconstruction": None,
        "comparisons": {},
        "errors": errors,
    }


def audit_surface_constitutive_sample(value: object) -> dict[str, Any]:
    """Return a strict, finite A15 assessment without mutating ``value``."""

    try:
        sample = _parse_sample(value)
        reconstruction = reconstruct_surface_constitutive_sample(value)
    except (SurfaceConstitutiveOracleError, TypeError, ValueError, OverflowError) as exc:
        return _failure_audit([str(exc)])

    production = sample.production_values
    observed_distance = sample.distance_m
    expected_distance = math.fsum(
        (sample.surface_radius_m, -sample.outer_center_radius_m)
    )
    comparisons = {
        "center_to_surface_distance_m": _eight_ulp_comparison(
            observed_distance,
            expected_distance,
        ),
        "selected_binary_diffusivity_m2_s": _eight_ulp_comparison(
            sample.serialized_selected_diffusivity_m2_s,
            reconstruction.selected_binary_diffusivity_m2_s,
        ),
        "thermodynamic_force_temperature_k": _eight_ulp_comparison(
            float(production["thermodynamic_force_temperature_k"]),
            reconstruction.common_face_temperature_k,
        ),
        "selected_binary_mobility_mol_m_s": _eight_ulp_comparison(
            float(production["selected_binary_mobility_mol_m_s"]),
            reconstruction.selected_binary_mobility_mol_m_s,
        ),
        "gas_chemical_force_gradient_m_inv": _eight_ulp_comparison(
            float(production["gas_chemical_force_gradient_m_inv"]),
            reconstruction.gas_chemical_force_gradient_m_inv,
        ),
        "binary_thermal_force_gradient_m_inv": _eight_ulp_comparison(
            float(production["binary_thermal_force_gradient_m_inv"]),
            reconstruction.binary_thermal_force_gradient_m_inv,
        ),
        "effective_binary_force_gradient_m_inv": _eight_ulp_comparison(
            float(production["effective_binary_force_gradient_m_inv"]),
            reconstruction.effective_binary_force_gradient_m_inv,
        ),
        "independent_hexane_flux_mol_m2_s": _eight_ulp_comparison(
            float(production["independent_hexane_flux_mol_m2_s"]),
            reconstruction.independent_hexane_flux_mol_m2_s,
        ),
        "binary_gas_entropy_production_w_m3_k": _eight_ulp_comparison(
            float(production["binary_gas_entropy_production_w_m3_k"]),
            reconstruction.binary_gas_entropy_production_w_m3_k,
        ),
    }
    failed = [name for name, item in comparisons.items() if not item["within_eight_ulps"]]
    errors = [f"{name} differs from primitive reconstruction by more than eight ULPs" for name in failed]
    return {
        **_failure_audit(errors),
        "passed": not errors,
        "reconstruction": asdict(reconstruction),
        "comparisons": comparisons,
    }


def require_surface_constitutive_sample(value: object) -> dict[str, Any]:
    """Return the audit or raise when the primitive sample fails closed."""

    audit = audit_surface_constitutive_sample(value)
    if audit["passed"] is not True:
        raise SurfaceConstitutiveOracleError("; ".join(audit["errors"]))
    return audit


def audit_history_aware_flux_response(value: object) -> dict[str, Any]:
    """Gate the A15 accepted-endpoint J_h response without prescribing return sign.

    The caller selects the four chronological accepted roots.  Full trajectory
    chronology remains a backend/release concern; this function independently
    audits each selected root and derives every response from reconstructed J_h.
    """

    base = {
        "schema_version": ORACLE_OUTPUT_SCHEMA_VERSION,
        "amendment_id": AMENDMENT_ID,
        "passed": False,
        "sample_roles": list(HISTORY_SAMPLE_ROLES),
        "baseline_end_outward": False,
        "pulse_first_inward": False,
        "pulse_strictly_decreased_J_h": False,
        "return_strictly_increased_J_h": False,
        "absolute_return_sign_prescribed": False,
        "pulse_minus_baseline_J_h_mol_m2_s": None,
        "return_minus_pulse_J_h_mol_m2_s": None,
        "reconstructed_J_h_mol_m2_s": None,
        "sample_audits": None,
        "N_h_or_N_t_substitution_allowed": False,
        "engineering_foundation_only": True,
        "physically_qualifying": False,
        "plant_predictive": False,
        "normative_F2_solvent_particle": False,
        "normative_F3_water_interface": False,
        "errors": [],
    }
    try:
        _reject_flux_surrogates(value, "history_samples")
        samples = _strict_mapping(value, HISTORY_SAMPLE_ROLES, "history_samples")
    except SurfaceConstitutiveOracleError as exc:
        return {**base, "errors": [str(exc)]}

    audits = {role: audit_surface_constitutive_sample(samples[role]) for role in HISTORY_SAMPLE_ROLES}
    failed_roles = [role for role, audit in audits.items() if audit["passed"] is not True]
    if failed_roles:
        return {
            **base,
            "sample_audits": audits,
            "errors": [f"history sample {role} failed its primitive oracle" for role in failed_roles],
        }

    fluxes = {
        role: float(audits[role]["reconstruction"]["independent_hexane_flux_mol_m2_s"])
        for role in HISTORY_SAMPLE_ROLES
    }
    pulse_change = math.fsum((fluxes["pulse_first"], -fluxes["baseline_end"]))
    return_change = math.fsum((fluxes["return_first"], -fluxes["pulse_end"]))
    baseline_outward = fluxes["baseline_end"] > 0.0
    pulse_inward = fluxes["pulse_first"] < 0.0
    pulse_decreased = pulse_change < 0.0
    return_increased = return_change > 0.0
    passed = baseline_outward and pulse_inward and pulse_decreased and return_increased
    errors: list[str] = []
    if not baseline_outward:
        errors.append("baseline-end reconstructed J_h is not strictly outward")
    if not pulse_inward:
        errors.append("pulse-first reconstructed J_h is not strictly inward")
    if not pulse_decreased:
        errors.append("pulse-first reconstructed J_h did not strictly decrease")
    if not return_increased:
        errors.append("return-first reconstructed J_h did not strictly increase from pulse-end")
    return {
        **base,
        "passed": passed,
        "baseline_end_outward": baseline_outward,
        "pulse_first_inward": pulse_inward,
        "pulse_strictly_decreased_J_h": pulse_decreased,
        "return_strictly_increased_J_h": return_increased,
        "pulse_minus_baseline_J_h_mol_m2_s": pulse_change,
        "return_minus_pulse_J_h_mol_m2_s": return_change,
        "reconstructed_J_h_mol_m2_s": fluxes,
        "sample_audits": audits,
        "errors": errors,
    }


__all__ = [
    "AMENDMENT_ID",
    "BINARY_DIFFUSIVITY_FIELDS",
    "COMMON_FACE_NEGLIGIBLE_SORET",
    "COMPARISON_ULPS",
    "HISTORY_SAMPLE_ROLES",
    "MASS_FORCE_FIELDS",
    "ORACLE_OUTPUT_SCHEMA_VERSION",
    "OUTER_DRY_STATE_FIELDS",
    "PRODUCTION_VALUE_FIELDS",
    "RAW_SAMPLE_KIND",
    "RAW_SAMPLE_SCHEMA_VERSION",
    "REFERENCE_INVARIANT_BINARY_THERMAL_FACTOR",
    "REFERENCE_INVARIANT_THERMAL_FORCE_FACTORS",
    "RESOLVED_SURFACE_STATE_FIELDS",
    "SUPPORTED_MASS_FORCE_MODES",
    "SurfaceConstitutiveOracleError",
    "SurfaceConstitutiveReconstruction",
    "TOP_LEVEL_FIELDS",
    "audit_history_aware_flux_response",
    "audit_surface_constitutive_sample",
    "raw_sample_schema",
    "reconstruct_surface_constitutive_sample",
    "require_surface_constitutive_sample",
]
