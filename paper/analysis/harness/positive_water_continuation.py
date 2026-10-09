"""Shared, owner-approved continued positive-water configuration for both baths."""

from dataclasses import dataclass, fields, replace
import math

from scipy.optimize import brentq
from dtdc_simulator.core2.props import sorption as sp, binary_gas as bg, hexane as hx
from dtdc_simulator.core2.particle import wet_core, cut_integrator as ci

APPROVED_NATIVE_CHANGES = {
    "props/sorption.py",
    "particle/wet_core.py",
    "particle/wet_water.py",
    "particle/coupled_pore.py",
    "particle/cut_transport.py",
}


@dataclass(frozen=True)
class RegionalRootStudyControls(ci.CutSolverControls):
    """Owner RULE 2: local root resolution with strict native ledgers.

    This declaration does not change state, constitutive laws, phase gates,
    or the 1e-10 component/energy acceptance bound. The native validation
    probe checks all other inherited fields without advancing a particle.
    """

    def __post_init__(self):
        if not 0 < self.nonlinear_residual_tolerance <= 1e-8:
            raise ValueError("regional root study resolution must be <=1e-8")
        if not 10 <= self.maximum_logit_magnitude <= 70:
            raise ValueError("regional root study chart must lie in [10,70]")
        kwargs = {f.name: getattr(self, f.name) for f in fields(ci.CutSolverControls) if f.init}
        kwargs["nonlinear_residual_tolerance"] = min(self.nonlinear_residual_tolerance, 1e-10)
        kwargs["maximum_logit_magnitude"] = min(self.maximum_logit_magnitude, 35)
        ci.CutSolverControls(**kwargs)


@dataclass(frozen=True)
class VanishingCoreStudyControls(ci.CutSolverControls):
    """Declared local resolution for negligible remaining wet volume only.

    Owner RULE 2; inherited raw per-step/cumulative conservation <=1e-10
    remains binding. The runner additionally requires incoming z<=1e-7.
    No physics, material state, or accepted native codec identity is reset.
    """

    def __post_init__(self):
        if not 0 < self.nonlinear_residual_tolerance <= 2e-8:
            raise ValueError("vanishing-core numerical resolution must be <=2e-8")
        if not 10 <= self.maximum_logit_magnitude <= 70:
            raise ValueError("vanishing-core encoded chart must lie in [10,70]")
        kwargs = {f.name: getattr(self, f.name) for f in fields(ci.CutSolverControls) if f.init}
        kwargs["nonlinear_residual_tolerance"] = min(self.nonlinear_residual_tolerance, 1e-10)
        kwargs["maximum_logit_magnitude"] = min(self.maximum_logit_magnitude, 35)
        ci.CutSolverControls(**kwargs)


def continued_configuration(config):
    """Explicitly select one identical water law on both sides of the front."""
    old = config.wet.luikov
    law = sp.ContinuedPositiveLuikovParams(A1=old.A1, A2=old.A2, W_cap=old.W_cap)
    old_wet = config.wet.wet
    wet = wet_core.ContinuedWaterWetCoreParams(
        **{f.name: getattr(old_wet, f.name) for f in fields(wet_core.WetCoreParams)}, luikov=law
    )
    return replace(
        config,
        wet=replace(config.wet, wet=wet, luikov=law),
        dry=replace(config.dry, pore=replace(config.dry.pore, luikov=law)),
        interface_composition=replace(
            config.interface_composition,
            y_hexane_bounds=(
                config.interface_composition.y_hexane_bounds[0],
                math.nextafter(1.0, 0.0),
            ),
        ),
    )


def continued_controls(controls, config):
    """Positive-water chart; upper interface endpoint from the native pure-H EOS."""
    p = config.dry.pressure_pa
    endpoint = brentq(
        lambda t: math.log(
            bg.state(t, p, 0.0, 1.0, k_wh=config.dry.pore.k_wh).fugacity_hexane_pa
            / hx.state_Tp(t, p, "liquid").fugacity
        ),
        330.0,
        350.0,
        xtol=1e-11,
    )
    return replace(
        controls,
        wet_water_bounds=(0.0, controls.wet_water_bounds[1]),
        interface_temperature_bounds_k=(
            controls.interface_temperature_bounds_k[0],
            math.nextafter(endpoint, -math.inf),
        ),
    )


def adopt_checkpoint(checkpoint):
    """Requalify an incoming physical state under the approved new parameter type.

    This is explicit model adoption, not exact-identity substitution or a
    continuation claimed to use the parent's immutable executable. All native
    reference inventories and cumulative balances are carried without reset.
    """
    before = checkpoint["before"]
    old = ci.inventory_snapshot(before.transport)
    config = continued_configuration(checkpoint["config"])
    transport = replace(before.transport, config=config)
    new = ci.inventory_snapshot(transport)
    if old != new:
        raise ValueError("continued-law adoption changed native particle inventory/energy")
    controls = continued_controls(before.controls, config)
    after = replace(before, transport=transport, controls=controls)
    checkpoint.update(before=after, config=config, controls=controls)
    return {
        "authorization": "owner approval 2026-10-06",
        "same_physical_primitives": True,
        "inventory_and_energy_bit_identical": True,
        "native_reference_and_cumulative_ledgers_carried": True,
        "same_immutable_executable_claimed": False,
        "continued_law": type(config.wet.luikov).__name__,
        "qualified_lower_anchor": sp.LuikovParams().W_ref,
        "open_positive_water_lower_endpoint": 0.0,
        "physically_qualifying": False,
    }
