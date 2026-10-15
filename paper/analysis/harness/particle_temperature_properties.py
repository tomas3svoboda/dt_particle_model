"""Explicit paper-only temperature-property variant; native source stays read-only.

Soybean oil: Garcia Rojas et al. (2013), Table 7, DOI
10.1080/10942912.2011.604889, cp=1024+3*T J/(kg K).
Binary gas: Fuller et al. (1966), Eq. 4, DOI 10.1021/ie50677a007,
normalized T**1.75/P shape with a declared effective-pore reference value.
The reference magnitude is not a converted Cardarelli coefficient.

Only four native functions receive exact, checked in-memory substitutions:
wet-core capacity, dry-thermo capacity, and the two binary face intervals.
All storage paths reuse the native reference-anchored oil weights and a common
integrated sensible-energy helper. Existing numerical normalization scales,
retention laws, water mobilities, proof budgets and acceptance gates stay fixed.
"""

from contextlib import contextmanager
from dataclasses import replace
import inspect
import math
import textwrap

from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.particle import transport_coefficients as tc

REFERENCE_T_K = 378.15
REFERENCE_P_PA = 101325.0
OIL_SOURCE_RANGE_K = (299.15, 433.15)


def oil_heat_capacity(temperature_k):
    if not math.isfinite(temperature_k) or temperature_k <= 0:
        raise ValueError("oil heat capacity requires a positive finite Kelvin temperature")
    return 1024.0 + 3.0 * temperature_k


def composite_capacity(temperature_k, w_o, cp_native=sp.CP_NATIVE_DRY_MEAL, cp_oil=sp.CP_OIL):
    """Native composite construction with the published oil cp(T)."""
    return sp.composite_dry_meal_heat_capacity(w_o, cp_native, oil_heat_capacity(temperature_k))


def composite_energy(T, T_ref, w_o, cp_native=sp.CP_NATIVE_DRY_MEAL, cp_oil=sp.CP_OIL):
    """Exact integral, including the native signed off-anchor oil correction."""
    oil_heat_capacity(T)
    oil_heat_capacity(T_ref)
    alpha, beta = sp._oil_weights(w_o)
    return (T - T_ref) * (alpha * cp_native + beta * (1024.0 + 1.5 * (T + T_ref)))


def effective_diffusivity(T, P, reference_m2_s):
    if not all(math.isfinite(x) and x > 0 for x in (T, P, reference_m2_s)):
        raise ValueError("diffusivity requires positive finite T, P and reference magnitude")
    return reference_m2_s * (T / REFERENCE_T_K) ** 1.75 * REFERENCE_P_PA / P


def face_interval(interval, T, P):
    """Evaluate an already pore-effective interval once at native common face T."""
    scale = effective_diffusivity(T, P, 1.0)
    return tc.EffectiveBinaryDiffusivityInterval(
        interval.lower_m2_s * scale,
        interval.upper_m2_s * scale,
        interval.provenance
        + "; Fuller T^1.75/P shape, fixed effective hindrance; no extra porosity",
    )


def declaration(reference_m2_s, oil_law=True, temperature_law=True):
    return {
        "variant": "literature-informed unvalidated effective-pore reference",
        "D_b_reference_m2_s": reference_m2_s,
        "reference_temperature_k": REFERENCE_T_K,
        "reference_pressure_pa": REFERENCE_P_PA,
        "binary_temperature_law": "Fuller T^1.75/P" if temperature_law else "constant",
        "oil_heat_capacity_law": (
            "Garcia Rojas 2013 Table 7: 1024+3*T J/(kg K)" if oil_law else "constant 2000 J/(kg K)"
        ),
        "oil_source_temperature_range_k": OIL_SOURCE_RANGE_K,
        "oil_energy": "analytic integral with native reference-anchored signed weights",
        "effective_hindrance": "assumed independent of temperature, loading and moisture",
        "source_coefficient_convention_conversion_completed": False,
        "retained_water_mobility": "held at original shared selection; not retargeted to outcome",
        "native_acceptance_or_numerical_normalization_changed": False,
        "physically_qualifying": False,
    }


@contextmanager
def installed(reference_m2_s, *, oil_law=True, temperature_law=True):
    """Install the variant before initialization; restore every changed object."""
    from dtdc_simulator.core2.particle import wet_core as wc, dry_thermo as dt
    from dtdc_simulator.core2.particle import coupled_transport as ct, cut_transport as cut
    import positive_water_continuation as continued

    changes = []
    missing = object()

    def assign(module, name, value):
        changes.append((module, name, getattr(module, name, missing)))
        setattr(module, name, value)

    def substitute(module, name, old, new):
        original = getattr(module, name)
        source = textwrap.dedent(inspect.getsource(original))
        if source.count(old) != 1:
            raise ValueError(
                f"native function changed: expected one substitution in {module.__name__}.{name}"
            )
        assign(module, "_paper_cp", composite_capacity)
        assign(module, "_paper_face_interval", face_interval)
        changes.append((module, name, original))
        exec(
            compile(source.replace(old, new), inspect.getsourcefile(original), "exec"),
            module.__dict__,
        )

    try:
        if oil_law:
            assign(sp, "composite_dry_meal_sensible_energy", composite_energy)
            substitute(
                wc,
                "specific_heat_capacity",
                "sp.composite_dry_meal_heat_capacity(p.w_o, p.cp_dry_meal, p.cp_oil)",
                "_paper_cp(T, p.w_o, p.cp_dry_meal, p.cp_oil)",
            )
            substitute(
                dt,
                "evaluate",
                "sp.composite_dry_meal_heat_capacity(p.w_o, p.cp_dry_meal, p.cp_oil)",
                "_paper_cp(T, p.w_o, p.cp_dry_meal, p.cp_oil)",
            )
        if temperature_law:
            substitute(
                ct,
                "_face_constitutive_data",
                "config.binary_diffusivity.interval,",
                "_paper_face_interval(config.binary_diffusivity.interval, face_temperature, config.pressure_pa),",
            )
            substitute(
                cut,
                "_dry_flux",
                "config.dry.binary_diffusivity.interval,",
                "_paper_face_interval(config.dry.binary_diffusivity.interval, face_temperature, config.dry.pressure_pa),",
            )
        original_configuration = continued.continued_configuration

        def configuration(config):
            new = original_configuration(config)
            reference = tc.EffectiveBinaryDiffusivityInterval(
                reference_m2_s,
                reference_m2_s,
                "Declared effective binary pore reference at 378.15 K, 101325 Pa; Cardarelli-scale engineering assumption, not measured/mapped binary data",
            ).select(0.0)
            return replace(new, dry=replace(new.dry, binary_diffusivity=reference))

        assign(continued, "continued_configuration", configuration)
        wc.specific_energy.exact_memo_cache.clear()
        cut.clear_exact_input_caches()
        yield declaration(reference_m2_s, oil_law, temperature_law)
    finally:
        for module, name, original in reversed(changes):
            if original is missing:
                delattr(module, name)
            else:
                setattr(module, name, original)
        wc.specific_energy.exact_memo_cache.clear()
        cut.clear_exact_input_caches()
